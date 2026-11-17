# -*- coding: utf-8 -*-
"""分子进化分析：
1. 成对距离矩阵（p-distance / Jukes-Cantor 校正）
2. 邻接法 (Neighbor-Joining) 建树
3. k-mer 距离 + PCoA (经典多维标度) 降维
4. k-means 聚类
5. Newick 导出
"""

import math
import random
import re
from collections import Counter

# ---------------------------------------------------------------- 距离矩阵

def _pairwise(p: str, q: str) -> tuple:
    """返回 (有效比较位点数, 差异数)"""
    valid = diff = 0
    for a, b in zip(p, q):
        if a in "ACGTU" and b in "ACGTU":
            valid += 1
            if a != b:
                diff += 1
    return valid, diff


def distance_matrix(seqs: list, model: str = "jc69") -> dict:
    """计算成对距离矩阵。

    model: 'p' = p-distance; 'jc69' = Jukes-Cantor 校正
    返回 {names, matrix, model, note}
    """
    n = len(seqs)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            valid, diff = _pairwise(seqs[i], seqs[j])
            if valid == 0:
                d = 0.0
            else:
                p = diff / valid
                if model == "jc69":
                    if p >= 0.75:
                        d = 10.0  # 饱和，给一个大值
                    else:
                        d = round(-0.75 * math.log(1.0 - 4.0 * p / 3.0), 6)
                else:
                    d = round(p, 6)
            mat[i][j] = mat[j][i] = d
    return {"matrix": mat, "model": model}


# ---------------------------------------------------------------- 邻接法

class _Node:
    __slots__ = ("name", "children", "length", "idx")

    def __init__(self, name=None, idx=-1):
        self.name = name
        self.children = []
        self.length = 0.0   # 到父节点的枝长
        self.idx = idx


def neighbor_joining(names: list, dist: list) -> dict:
    """邻接法建树（Saitou & Nei 1987）。

    约定：D 恒为 m×m 方阵，行/列下标 = active 列表中的位置；
    active[k] 才是 nodes 列表中的节点下标。
    返回 {tree: 嵌套节点 JSON, newick: str}
    """
    n = len(names)
    if n == 0:
        return {"error": "没有序列可建树。"}
    if n == 1:
        node = {"name": names[0], "length": 0.0, "children": [], "leaf": True}
        return {"tree": node, "newick": f"{names[0]};"}

    D = [row[:] for row in dist]              # 深拷贝，避免污染输入
    nodes = [_Node(names[i], i) for i in range(n)]
    active = list(range(n))                   # 活动节点在 nodes 中的下标

    while len(active) > 2:
        m = len(active)
        # 1) 净分歧率 r_i
        r = [sum(D[i][j] for j in range(m)) / (m - 2) for i in range(m)]
        # 2) Q 矩阵最小的节点对（位置 i, j）
        best_i, best_j, best_q = 0, 1, float("inf")
        for i in range(m):
            for j in range(i + 1, m):
                q = (m - 2) * D[i][j] - r[i] - r[j]
                if q < best_q:
                    best_q, best_i, best_j = q, i, j
        duv = D[best_i][best_j]
        # 3) 枝长（负值截断为 0）
        lu = 0.5 * duv + 0.5 * (r[best_i] - r[best_j])
        lv = duv - lu
        lu, lv = max(lu, 0.0), max(lv, 0.0)

        # 4) 新建父节点，挂上 u、v
        parent = _Node()
        nodes[active[best_i]].length = lu
        nodes[active[best_j]].length = lv
        parent.children = [nodes[active[best_i]], nodes[active[best_j]]]
        nodes.append(parent)
        parent_idx = len(nodes) - 1

        # 5) 压缩矩阵：保留位置（除 best_i/best_j），末位为新节点
        keep = [k for k in range(m) if k != best_i and k != best_j]
        new_m = m - 1
        new_mat = [[0.0] * new_m for _ in range(new_m)]
        for a_i, ka in enumerate(keep):
            for b_i in range(a_i + 1, len(keep)):
                v = D[ka][keep[b_i]]
                new_mat[a_i][b_i] = new_mat[b_i][a_i] = v
        last = new_m - 1
        for a_i, ka in enumerate(keep):
            dk = 0.5 * (D[ka][best_i] + D[ka][best_j] - duv)
            dk = max(dk, 0.0)
            new_mat[a_i][last] = new_mat[last][a_i] = dk

        active = [active[k] for k in keep] + [parent_idx]
        D = new_mat

    # 剩余两个节点：等分连接到根
    a_idx, b_idx = active[0], active[1]
    root = _Node()
    half = max(D[0][1] / 2.0, 0.0)
    nodes[a_idx].length = half
    nodes[b_idx].length = half
    root.children = [nodes[a_idx], nodes[b_idx]]

    def to_json(nd: _Node) -> dict:
        if not nd.children:
            return {"name": nd.name, "length": round(nd.length, 5), "leaf": True}
        return {
            "name": nd.name or "",
            "length": round(nd.length, 5),
            "leaf": False,
            "children": [to_json(c) for c in nd.children],
        }

    tree = to_json(root)
    newick = _to_newick(root)
    return {"tree": tree, "newick": newick}


def _to_newick(nd: _Node) -> str:
    if not nd.children:
        label = nd.name or "Leaf"
        return f"{label}:{nd.length:.5f}"
    inner = ",".join(_to_newick(c) for c in nd.children)
    label = nd.name or ""
    return f"({inner}){label}:{nd.length:.5f}"


# ---------------------------------------------------------------- PCoA

def jacobi_eigen(a: list, tol: float = 1e-10, max_iter: int = 100):
    """Jacobi 旋转法求实对称矩阵全部特征值/特征向量。

    返回 (eigvals, eigvecs) — eigvecs[k] 为第 k 个特征向量（列向量存为行）
    """
    n = len(a)
    m = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _ in range(max_iter):
        # 找最大非对角元
        off, p, q = 0.0, 0, 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(m[i][j]) > off:
                    off, p, q = abs(m[i][j]), i, j
        if off < tol:
            break
        app, aqq, apq = m[p][p], m[q][q], m[p][q]
        if abs(apq) < 1e-300:
            continue
        theta = (aqq - app) / (2.0 * apq)
        t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
        c = 1.0 / math.sqrt(t * t + 1.0)
        s = t * c
        # 旋转
        for k in range(n):
            if k != p and k != q:
                mkp, mkq = m[k][p], m[k][q]
                m[k][p] = m[p][k] = c * mkp - s * mkq
                m[k][q] = m[q][k] = s * mkp + c * mkq
        m[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq
        m[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq
        m[p][q] = m[q][p] = 0.0
        for k in range(n):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq

    eigvals = [m[i][i] for i in range(n)]
    pairs = sorted(range(n), key=lambda i: -eigvals[i])
    eigvals_sorted = [eigvals[i] for i in pairs]
    eigvecs_sorted = [[v[r][i] for r in range(n)] for i in pairs]
    return eigvals_sorted, eigvecs_sorted


def pcoa(dist: list) -> dict:
    """经典多维标度 (Principal Coordinates Analysis)。

    B = -1/2 · J D² J，取前两大特征向量作为二维坐标。
    返回 {coords:[[x,y]...], explained:[pc1,pc2], all_eigen:[...]}
    """
    n = len(dist)
    if n == 0:
        return {"coords": [], "explained": []}
    if n == 1:
        return {"coords": [[0.0, 0.0]], "explained": [0.0, 0.0]}

    # D²
    d2 = [[dist[i][j] ** 2 for j in range(n)] for i in range(n)]
    # 行/列均值与总均值
    row_mean = [sum(row) / n for row in d2]
    grand = sum(row_mean) / n
    # B = -0.5 * (d2 - row - col + grand)
    B = [[-0.5 * (d2[i][j] - row_mean[i] - row_mean[j] + grand)
          for j in range(n)] for i in range(n)]

    vals, vecs = jacobi_eigen(B)
    total = sum(v for v in vals if v > 0) or 1.0
    coords = []
    for i in range(n):
        x = vecs[0][i] * math.sqrt(max(vals[0], 0.0)) if len(vals) > 0 else 0.0
        y = vecs[1][i] * math.sqrt(max(vals[1], 0.0)) if len(vals) > 1 else 0.0
        coords.append([round(x, 5), round(y, 5)])
    explained = [
        round(max(vals[0], 0.0) * 100.0 / total, 2) if len(vals) > 0 else 0.0,
        round(max(vals[1], 0.0) * 100.0 / total, 2) if len(vals) > 1 else 0.0,
    ]
    return {
        "coords": coords,
        "explained": explained,
        "eigenvalues": [round(v, 5) for v in vals[: min(8, len(vals))]],
    }


# ---------------------------------------------------------------- k-means

def kmeans(points: list, k: int = 3, seed: int = 42, max_iter: int = 100) -> dict:
    """二维点集 k-means++ 初始化 + Lloyd 迭代；固定随机种子保证结果可复现。"""
    n = len(points)
    if n == 0:
        return {"clusters": [], "centers": []}
    k = max(1, min(k, n))
    rng = random.Random(seed)
    pts = [(p[0], p[1]) for p in points]

    # k-means++ 初始化
    centers = [pts[rng.randrange(n)]]
    while len(centers) < k:
        d2min = []
        for p in pts:
            best = min((p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 for c in centers)
            d2min.append(best)
        total = sum(d2min)
        if total <= 0:
            centers.append(pts[rng.randrange(n)])
            continue
        target = rng.random() * total
        acc = 0.0
        for i, d in enumerate(d2min):
            acc += d
            if acc >= target:
                centers.append(pts[i])
                break

    labels = [0] * n
    for _ in range(max_iter):
        changed = False
        # 分配
        for i, p in enumerate(pts):
            bi, bd = 0, float("inf")
            for ci, c in enumerate(centers):
                d = (p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2
                if d < bd:
                    bd, bi = d, ci
            if labels[i] != bi:
                labels[i] = bi
                changed = True
        # 更新中心
        for ci in range(k):
            members = [pts[i] for i in range(n) if labels[i] == ci]
            if members:
                cx = sum(p[0] for p in members) / len(members)
                cy = sum(p[1] for p in members) / len(members)
                centers[ci] = (cx, cy)
        if not changed:
            break

    return {
        "labels": labels,
        "centers": [[round(c[0], 5), round(c[1], 5)] for c in centers],
        "k": k,
        "seed": seed,
    }


# ---------------------------------------------------------------- 汇总入口

def phylo_analysis(names: list, seqs: list, model: str = "jc69") -> dict:
    """距离矩阵 + NJ 建树"""
    err = ""
    if len(names) < 3:
        err = "系统发育建树至少需要 3 条序列（建议 4 条以上）。"
    if len(names) != len(seqs):
        return {"error": "序列名称与数量不一致。"}
    dm = distance_matrix(seqs, model)
    out = {"names": names, "model": model, "dist_matrix": dm["matrix"]}
    if err:
        out["error"] = err
        return out
    nj = neighbor_joining(names, dm["matrix"])
    out.update(nj)
    return out


def cluster_analysis(names: list, seqs: list, k: int = 3, model: str = "p") -> dict:
    """k-mer Jaccard 距离 + PCoA + k-means 聚类（无监督机器学习视角）"""
    if len(names) != len(seqs):
        return {"error": "序列名称与数量不一致。"}
    if len(names) < 3:
        return {"error": "聚类分析至少需要 3 条序列。"}
    km = kmer_jaccard_matrix(seqs, k=3)
    pc = pcoa(km)
    kmres = kmeans(pc["coords"], k=k)
    return {
        "names": names,
        "dist_matrix": km,
        "coords": pc["coords"],
        "explained": pc["explained"],
        "eigenvalues": pc.get("eigenvalues", []),
        "labels": kmres["labels"],
        "centers": kmres["centers"],
        "k": kmres["k"],
        "method": "3-mer 加权 Jaccard → PCoA → k-means",
    }


def kmer_jaccard_matrix(seqs: list, k: int = 3) -> list:
    """k-mer 计数向量的加权 Jaccard（Ruzicka）距离。

    d = 1 - Σ min(f_a, f_b) / Σ max(f_a, f_b)

    相比集合 Jaccard，计数版能反映 k-mer 频率（组成）差异：
    短序列上集合 Jaccard 会因 k-mer 空间饱和而失效。
    """
    counts = []
    for s in seqs:
        s = re.sub(r"[^ACGTU]", "", s.upper())
        if len(s) < k:
            counts.append(Counter())
        else:
            counts.append(Counter(s[i:i + k] for i in range(len(s) - k + 1)))
    n = len(counts)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            keys = set(counts[i]) | set(counts[j])
            inter = sum(min(counts[i][t], counts[j][t]) for t in keys)
            union = sum(max(counts[i][t], counts[j][t]) for t in keys)
            d = 1.0 - inter / union if union else 0.0
            mat[i][j] = mat[j][i] = round(d, 6)
    return mat
