# -*- coding: utf-8 -*-
"""限制性内切酶酶切图谱：位点搜索、切割位置计算、片段长度统计。

酶切位点表达式: 以 ^ 标记切割位置（如 G^AATTC），
支持 IUPAC 简并碱基（N/R/Y/S/W/K/M/B/D/H/V）。
"""

# 常用限制性内切酶（^ 为切割位点）
ENZYMES = [
    {"name": "EcoRI",    "site": "G^AATTC",   "palindromic": True},
    {"name": "BamHI",    "site": "G^GATCC",   "palindromic": True},
    {"name": "HindIII",  "site": "A^AGCTT",   "palindromic": True},
    {"name": "PstI",     "site": "CTGCA^G",   "palindromic": True},
    {"name": "NotI",     "site": "GC^GGCCGC", "palindromic": True},
    {"name": "XbaI",     "site": "T^CTAGA",   "palindromic": True},
    {"name": "SpeI",     "site": "A^ACTAGT",  "palindromic": True},
    {"name": "SmaI",     "site": "CCC^GGG",   "palindromic": True},
    {"name": "KpnI",     "site": "GGTAC^C",   "palindromic": True},
    {"name": "SacI",     "site": "GAGCT^C",   "palindromic": True},
    {"name": "EcoRV",    "site": "GAT^ATC",   "palindromic": True},
    {"name": "HaeIII",   "site": "GG^CC",     "palindromic": True},
    {"name": "AluI",     "site": "AG^CT",     "palindromic": True},
    {"name": "HinfI",    "site": "G^ANTC",    "palindromic": True},
    {"name": "MboI",     "site": "^GATC",     "palindromic": True},
    {"name": "Sau3AI",   "site": "^GATC",     "palindromic": True},
    {"name": "DpnI",     "site": "GA^TC",     "palindromic": True},
    {"name": "NcoI",     "site": "C^CATGG",   "palindromic": True},
    {"name": "NdeI",     "site": "CA^TATG",   "palindromic": True},
    {"name": "XhoI",     "site": "C^TCGAG",   "palindromic": True},
    {"name": "BglII",    "site": "A^GATCT",   "palindromic": True},
    {"name": "ClaI",     "site": "AT^CGAT",   "palindromic": True},
    {"name": "SalI",     "site": "G^TCGAC",   "palindromic": True},
    {"name": "ScaI",     "site": "AGT^ACT",   "palindromic": True},
    {"name": "PvuII",    "site": "CAG^CTG",   "palindromic": True},
    {"name": "TaqI",     "site": "T^CGA",     "palindromic": True},
]

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def parse_site(expr: str):
    """解析 'G^AATTC' → (pattern: str, cut_offset: int)
    cut_offset = 切割位点在 pattern 中的插入位置（0..len）"""
    expr = expr.strip().upper().replace(" ", "")
    if "^" not in expr:
        raise ValueError(f"酶切位点缺少 ^ 标记: {expr}")
    left, right = expr.split("^", 1)
    pattern = left + right
    return pattern, len(left)


def match_at(dna: str, pos: int, pattern: str) -> bool:
    """判断 dna[pos:] 是否匹配 pattern（支持简并碱基）"""
    if pos + len(pattern) > len(dna):
        return False
    for k, ch in enumerate(pattern):
        allowed = IUPAC.get(ch)
        if allowed is None or dna[pos + k] not in allowed:
            return False
    return True


def find_sites(dna: str, expr: str) -> list:
    """找出所有切割位置。

    pos = 切点左侧碱基数（0-based 切点索引），即在第 pos 个与第 pos+1 个
    碱基之间切开；片段长度 = 相邻切点之差。线性分子两端无效切点已剔除。
    返回 [{pos, site_start(1-based), site_end, matched}]
    """
    dna = dna.upper()
    pattern, cut_off = parse_site(expr)
    L = len(pattern)
    n = len(dna)
    out = []
    for i in range(n - L + 1):
        if not match_at(dna, i, pattern):
            continue
        pos = i + cut_off
        if pos <= 0 or pos >= n:
            continue  # 线性分子两端的切点不产生片段
        out.append({
            "pos": pos,
            "site_start": i + 1,
            "site_end": i + L,
            "matched": dna[i:i + L],
        })
    return out


def digest(seq: str, enzyme_names: list) -> dict:
    """对序列进行酶切分析。

    返回 {
      enzymes: [{name, site, count, positions, fragments}],
      fragment_hist: [{len, count}]   # 所有酶共同作用下的片段长度分布（并集切点）
      linear: 序列视为线性分子
    }
    """
    dna = seq.upper().replace("-", "")
    n = len(dna)
    table = {e["name"]: e for e in ENZYMES}
    result_enzymes = []
    all_cuts = set()

    for name in enzyme_names:
        e = table.get(name)
        if not e:
            continue
        sites = find_sites(dna, e["site"])
        positions = [s["pos"] for s in sites]
        # 单酶切：线性分子片段
        frags = _fragments(n, positions)
        result_enzymes.append({
            "name": name,
            "site": e["site"],
            "count": len(sites),
            "positions": positions,
            "sites": sites,
            "fragments": frags,
        })
        all_cuts.update(positions)

    # 多酶共同酶切的片段分布
    combined = _fragments(n, sorted(all_cuts))
    hist = {}
    for L in combined:
        hist[L] = hist.get(L, 0) + 1
    fragment_hist = sorted(
        ({"len": L, "count": c} for L, c in hist.items()),
        key=lambda x: -x["len"],
    )

    return {
        "length": n,
        "enzymes": result_enzymes,
        "total_cuts": len(all_cuts),
        "combined_fragments": combined,
        "fragment_hist": fragment_hist[:30],
        "linear": True,
    }


def _fragments(n: int, cuts: list) -> list:
    """线性序列在切点集合下的片段长度"""
    cuts = sorted(c for c in cuts if 0 < c < n)
    frags, prev = [], 0
    for c in cuts:
        frags.append(c - prev)
        prev = c
    frags.append(n - prev)
    return frags
