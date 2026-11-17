# -*- coding: utf-8 -*-
"""双序列比对：Needleman-Wunsch 全局比对 / Smith-Waterman 局部比对（线性 gap 罚分）。

打分模型:
    s(a,b) = +match    (a == b)
             -mismatch (a != b)
    gap    = -gap_penalty (每个空位)

算法: 动态规划打分矩阵 + 回溯, 时间/空间 O(m*n)。
"""

from array import array


def _score_fn(match: float, mismatch: float):
    def s(a: str, b: str) -> float:
        return match if a == b else -mismatch
    return s


def needleman_wunsch(a: str, b: str, match=1.0, mismatch=1.0, gap=1.0) -> dict:
    """全局比对：两条序列首尾都要对齐"""
    return _align(a, b, match, mismatch, gap, mode="global")


def smith_waterman(a: str, b: str, match=1.0, mismatch=1.0, gap=1.0) -> dict:
    """局部比对：只对齐得分最高的相似片段"""
    return _align(a, b, match, mismatch, gap, mode="local")


def _align(a: str, b: str, match: float, mismatch: float, gap: float, mode: str) -> dict:
    a = a.upper().replace(" ", "").replace("\n", "")
    b = b.upper().replace(" ", "").replace("\n", "")
    m, n = len(a), len(b)
    if m == 0 or n == 0:
        return {"error": "两条序列均不能为空。", "mode": mode}
    if m * n > 25_000_000:
        return {"error": f"序列过长（{m}×{n} 超过 2.5×10⁷ 个矩阵单元），请截短后再比对。", "mode": mode}

    sc = _score_fn(match, mismatch)
    neg_gap = -gap
    local = mode == "local"

    # DP 分数矩阵与回溯方向矩阵（0=对角 1=来自上方 2=来自左方）
    H = [array("d", [0.0] * (n + 1)) for _ in range(m + 1)]
    ptr = [bytearray(n + 1) for _ in range(m + 1)]

    # 边界初始化：全局比对 H[0][j] = -j*gap, H[i][0] = -i*gap；局部比对全 0
    if not local:
        for j in range(1, n + 1):
            H[0][j] = neg_gap * j
            ptr[0][j] = 2
        for i in range(1, m + 1):
            H[i][0] = neg_gap * i
            ptr[i][0] = 1

    best, bi, bj = 0.0, 0, 0
    for i in range(1, m + 1):
        row_cur, row_up = H[i], H[i - 1]
        ptr_cur = ptr[i]
        ai = a[i - 1]
        for j in range(1, n + 1):
            d = row_up[j - 1] + sc(ai, b[j - 1])
            u = row_up[j] + neg_gap
            l = row_cur[j - 1] + neg_gap
            # 平局时优先对角 → 上 → 左，保证回溯结果稳定
            if d >= u and d >= l:
                row_cur[j], ptr_cur[j] = d, 0
            elif u >= l:
                row_cur[j], ptr_cur[j] = u, 1
            else:
                row_cur[j], ptr_cur[j] = l, 2
            if local:
                if row_cur[j] <= 0.0:
                    row_cur[j], ptr_cur[j] = 0.0, 0
                elif row_cur[j] > best:
                    best, bi, bj = row_cur[j], i, j
            else:
                pass  # 全局分数取 H[m][n]

    # 回溯
    if local:
        if best <= 0.0:
            return _empty_result(mode, a, b)
        i, j = bi, bj
    else:
        i, j = m, n

    chars_a, chars_b = [], []
    while True:
        if local and H[i][j] == 0.0:
            break
        p = ptr[i][j]
        if p == 0:
            chars_a.append(a[i - 1]); chars_b.append(b[j - 1]); i -= 1; j -= 1
        elif p == 1:
            chars_a.append(a[i - 1]); chars_b.append("-"); i -= 1
        else:
            chars_a.append("-"); chars_b.append(b[j - 1]); j -= 1
        if i <= 0 and j <= 0:
            break

    # 全局比对补齐前缀（单侧先到 0 的情况）
    if not local:
        while i > 0:
            chars_a.append(a[i - 1]); chars_b.append("-"); i -= 1
        while j > 0:
            chars_a.append("-"); chars_b.append(b[j - 1]); j -= 1

    chars_a.reverse(); chars_b.reverse()
    aa, bb = "".join(chars_a), "".join(chars_b)

    matches = mismatches = gaps = 0
    for x, y in zip(aa, bb):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matches += 1
        else:
            mismatches += 1
    aligned_len = len(aa)
    ident = round(matches * 100.0 / aligned_len, 2) if aligned_len else 0.0

    if local:
        a_chars = sum(1 for c in aa if c != "-")
        b_chars = sum(1 for c in bb if c != "-")
        start_a, end_a = bi - a_chars + 1, bi
        start_b, end_b = bj - b_chars + 1, bj
        score = best
    else:
        start_a, end_a, start_b, end_b = 1, m, 1, n
        score = H[m][n]

    return {
        "mode": mode,
        "score": round(score, 3),
        "aligned_a": aa,
        "aligned_b": bb,
        "identity": ident,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "aligned_len": aligned_len,
        "coverage_a": round((end_a - start_a + 1) * 100.0 / m, 2),
        "coverage_b": round((end_b - start_b + 1) * 100.0 / n, 2),
        "start_a": start_a, "end_a": end_a,
        "start_b": start_b, "end_b": end_b,
        "len_a": m, "len_b": n,
        "matrix_preview": _preview(a, b, H, mode),
    }


def _empty_result(mode: str, a: str, b: str) -> dict:
    return {
        "mode": mode, "score": 0.0, "aligned_a": "", "aligned_b": "",
        "identity": 0.0, "matches": 0, "mismatches": 0, "gaps": 0,
        "aligned_len": 0, "coverage_a": 0.0, "coverage_b": 0.0,
        "start_a": 0, "end_a": 0, "start_b": 0, "end_b": 0,
        "len_a": len(a), "len_b": len(b),
        "matrix_preview": _preview(a, b, [[0.0]], mode),
        "note": "未找到正得分的局部比对区域。",
    }


def _preview(a: str, b: str, H, mode: str, max_dim: int = 28) -> dict:
    """DP 矩阵左上角预览（供前端热图展示）"""
    m = min(len(a), max_dim)
    n = min(len(b), max_dim)
    grid, lo, hi = [], float("inf"), float("-inf")
    for i in range(m + 1):
        if i < len(H):
            row_src = H[i]
        else:
            row_src = None
        row = []
        for j in range(n + 1):
            v = row_src[j] if row_src is not None and j < len(row_src) else 0.0
            v = 0.0 if v in (float("-inf"),) else v
            row.append(round(v, 2))
            lo, hi = min(lo, v), max(hi, v)
        grid.append(row)
    if lo == float("inf"):
        lo, hi = 0.0, 0.0
    return {"a": a[:m], "b": b[:n], "grid": grid, "min": lo, "max": hi, "mode": mode}


def align(a: str, b: str, mode: str = "global", match=1.0, mismatch=1.0, gap=1.0) -> dict:
    """统一入口：mode = global | local"""
    if mode == "local":
        return smith_waterman(a, b, match, mismatch, gap)
    return needleman_wunsch(a, b, match, mismatch, gap)
