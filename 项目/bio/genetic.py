# -*- coding: utf-8 -*-
"""遗传密码：六框翻译、ORF 查找、蛋白质理化性质（分子量 / 等电点 / 消光系数）。"""

from .fasta import revcomp

# 标准遗传密码表（NCBI transl_table=1）
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODONS = ["ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"]

# 氨基酸平均残基质量（Da）
AA_MASS = {
    "G": 57.05, "A": 71.08, "S": 87.08, "P": 97.12, "V": 99.13,
    "T": 101.11, "C": 103.14, "L": 113.16, "I": 113.16, "N": 114.11,
    "D": 115.09, "Q": 128.13, "K": 128.17, "E": 129.12, "M": 131.19,
    "H": 137.14, "F": 147.17, "R": 156.19, "Y": 163.17, "W": 186.21,
    "U": 150.04, "X": 110.0,
}
WATER = 18.015

# 等电点计算用 pKa（Bjellqvist / EMBOSS 风格）
N_TERM_PKA = {"A": 7.59, "M": 6.93, "S": 6.93, "P": 8.36, "T": 6.82,
              "V": 7.44, "E": 7.70, "L": 6.70, "I": 7.52, "K": 6.69,
              "D": 6.77, "F": 6.77, "Q": 6.73, "Y": 6.77, "N": 6.65,
              "G": 7.75, "H": 6.77, "R": 6.73, "C": 6.73, "W": 6.73}
C_TERM_PKA = {"D": 2.55, "E": 3.55, "C": 3.55, "Y": 3.55, "H": 5.98,
              "K": 3.55, "R": 3.55, "S": 3.55, "T": 3.55, "W": 3.55}
SIDE_PKA = {
    "D": (3.9, "acid"), "E": (4.1, "acid"),
    "C": (8.5, "acid"), "Y": (10.1, "acid"),
    "H": (6.0, "base"), "K": (10.5, "base"),
    "R": (12.5, "base"),
}
# 组氨酸在 pH 6.0 附近的正电荷修正
IMIDAZOLIUM = (6.0, "base")

AA_CHARGE = {  # 各残基对总电荷的贡献（正碱性 / 负酸性）
    "K": 1, "R": 1, "H": 1,
    "D": -1, "E": -1, "C": -1, "Y": -1,
}


def translate(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    """按读框翻译（frame: 0/1/2；负链请调用方先 revcomp）。未知/简并密码子 → X"""
    s = seq.upper()[frame:]
    out = []
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "X")
        out.append(aa)
        if to_stop and aa == "*":
            break
    return "".join(out)


def six_frame(seq: str) -> list:
    """六框翻译结果 [{frame: +1..+3/-1..-3, strand, frame_idx, protein}]"""
    s = seq.upper()
    rc = revcomp(s)
    res = []
    for f in range(3):
        res.append({"frame": f"+{f + 1}", "strand": "+", "frame_idx": f,
                    "protein": translate(s, f)})
    for f in range(3):
        res.append({"frame": f"-{f + 1}", "strand": "-", "frame_idx": f,
                    "protein": translate(rc, f)})
    return res


def find_orfs(seq: str, min_aa: int = 30, start_codons=None, both_strands: bool = True,
              to_stop: bool = True) -> list:
    """查找开放阅读框。

    返回 [{frame, strand, start, end, length_nt, length_aa, protein, start_codon, stop_codon}]
    start/end 为正链 1-based 坐标。
    """
    s = seq.upper().replace("-", "")
    starts = set(start_codons or ["ATG"])
    out = []

    def scan(dna: str, strand: str):
        n = len(dna)
        for f in range(3):
            prot = translate(dna, f, to_stop=False)  # 预先整框翻译，便于定位起止
            i = 0
            while i < len(prot):
                start_codon = dna[f + i * 3:f + i * 3 + 3]
                if start_codon in starts and prot[i] != "*":
                    stop_at = -1
                    for j in range(i + 1, len(prot)):
                        if prot[j] == "*":
                            stop_at = j
                            break
                    if stop_at == -1:
                        if to_stop:
                            break  # 该读框后续均无终止子
                        stop_at = len(prot)  # 记为延伸到末端的不完整 ORF
                    aa_len = stop_at - i
                    if aa_len >= min_aa:
                        # 有终止子时 ORF 含终止密码子；不完整 ORF 延伸至读框末端
                        last_codon = stop_at if stop_at < len(prot) else len(prot) - 1
                        if strand == "+":
                            nt_start = f + i * 3 + 1
                            nt_end = f + (last_codon + 1) * 3
                        else:
                            # 负链：反向互补坐标 → 正链 1-based 坐标
                            rc_start = f + i * 3
                            rc_end_excl = f + (last_codon + 1) * 3
                            nt_start = n - rc_end_excl + 1
                            nt_end = n - rc_start
                        protein = prot[i:stop_at]
                        out.append({
                            "frame": (f"+{f + 1}" if strand == "+" else f"-{f + 1}"),
                            "strand": strand,
                            "frame_idx": f,
                            "start": nt_start,
                            "end": nt_end,
                            "length_nt": nt_end - nt_start + 1,
                            "length_aa": len(protein),
                            "start_codon": start_codon,
                            "stop_codon": "*" if stop_at < len(prot) else "-",
                            "protein": protein,
                        })
                    i = stop_at + 1 if stop_at >= 0 else i + 1
                else:
                    i += 1

    scan(s, "+")
    if both_strands:
        scan(revcomp(s), "-")

    out.sort(key=lambda x: -x["length_aa"])
    # 去重（同一区域不同读框重复计数较少，保留即可）
    return out


def molecular_weight(protein: str) -> float:
    """蛋白质平均分子量（Da），忽略 * 与 X 未知残基"""
    p = protein.upper().replace("*", "").replace("-", "")
    if not p:
        return 0.0
    total = WATER
    for aa in p:
        total += AA_MASS.get(aa, 110.0)
    return round(total, 2)


def isoelectric_point(protein: str, pH_low: float = 0.0, pH_high: float = 14.0,
                      tol: float = 1e-4) -> float:
    """二分法求等电点（净电荷 = 0 时的 pH）"""
    p = protein.upper().replace("*", "").replace("-", "")
    if not p:
        return 0.0

    def charge(pH: float) -> float:
        # N 端
        pka_n = N_TERM_PKA.get(p[0], 7.5)
        q = 1.0 / (1.0 + 10 ** (pH - pka_n))
        # C 端
        pka_c = C_TERM_PKA.get(p[-1], 3.55)
        q -= 1.0 / (1.0 + 10 ** (pka_c - pH))
        # 侧链
        for aa in p[1:]:
            if aa in SIDE_PKA:
                pka, kind = SIDE_PKA[aa]
                if kind == "acid":
                    q -= 1.0 / (1.0 + 10 ** (pka - pH))
                else:
                    q += 1.0 / (1.0 + 10 ** (pH - pka))
        return q

    lo, hi = pH_low, pH_high
    if charge(lo) < 0:
        return round(lo, 2)
    if charge(hi) > 0:
        return round(hi, 2)
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return round((lo + hi) / 2, 2)


def aromaticity(protein: str) -> float:
    """芳香族氨基酸（F+W+Y）占比"""
    p = [c for c in protein.upper() if c in "ACDEFGHIKLMNPQRSTVWY"]
    if not p:
        return 0.0
    return round(sum(1 for c in p if c in "FWY") * 100.0 / len(p), 2)


def gravy(protein: str) -> float:
    """GRAVY 疏水性指数（Kyte-Doolittle）"""
    kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
          "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
          "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
          "Y": -1.3, "V": 4.2}
    p = [c for c in protein.upper() if c in kd]
    if not p:
        return 0.0
    return round(sum(kd[c] for c in p) / len(p), 3)


def protein_properties(protein: str) -> dict:
    """蛋白质理化性质汇总"""
    p = protein.upper().replace("*", "").replace("-", "")
    return {
        "length_aa": len(p),
        "molecular_weight": molecular_weight(p),
        "isoelectric_point": isoelectric_point(p),
        "aromaticity": aromaticity(p),
        "gravy": gravy(p),
        "charge_ph7": round(_net_charge(p, 7.0), 3),
    }


def _net_charge(protein: str, pH: float) -> float:
    if not protein:
        return 0.0
    p = protein
    pka_n = N_TERM_PKA.get(p[0], 7.5)
    q = 1.0 / (1.0 + 10 ** (pH - pka_n))
    pka_c = C_TERM_PKA.get(p[-1], 3.55)
    q -= 1.0 / (1.0 + 10 ** (pka_c - pH))
    for aa in p[1:]:
        if aa in SIDE_PKA:
            pka, kind = SIDE_PKA[aa]
            if kind == "acid":
                q -= 1.0 / (1.0 + 10 ** (pka - pH))
            else:
                q += 1.0 / (1.0 + 10 ** (pH - pka))
    return q
