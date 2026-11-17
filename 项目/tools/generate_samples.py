# -*- coding: utf-8 -*-
"""生成内置示例数据（samples/ 目录）。

- 示例A：真实序列 —— 人 β-珠蛋白基因 (HBB) CDS 前 93 bp（经典教学片段）
- 示例B：模拟编码基因 —— 用于 ORF / 翻译 / 酶切演示
- 示例C：模拟病毒演化数据 —— 用于邻接法建树演示
- 示例D：模拟同源基因对 —— 用于双序列比对演示
- 示例E：模拟微生物群落（3 个类群）—— 用于 k-mer + PCoA + k-means 聚类演示

所有“模拟”序列均由固定随机种子生成，仅用于教学演示。
"""

import os
import random

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "samples")
os.makedirs(OUT, exist_ok=True)

rng = random.Random(20261113)


def rand_dna(n: int, gc: float = 0.5) -> str:
    p_gc = gc / 2
    p_at = (1 - gc) / 2
    return "".join(
        rng.choices(["G", "C", "A", "T"], weights=[p_gc, p_gc, p_at, p_at], k=n)
    )


def mutate(seq: str, rate: float) -> str:
    out = []
    for ch in seq:
        if rng.random() < rate:
            out.append(rng.choice([c for c in "ACGT" if c != ch]))
        else:
            out.append(ch)
    return "".join(out)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------- A: 真实 HBB 片段
HBB = (
    ">HBB_fragment 人β-珠蛋白基因(HBB) CDS 前93bp [真实序列, 教学片段]\n"
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTG\n"
    "GATGAAGTTGGTGGTGAGGCCCTGGGCAGG\n"
)

# ---------------------------------------------------------------- B: 模拟编码基因
def make_cds(n_aa: int) -> str:
    bases = {"F": "TTT", "L": "CTG", "I": "ATT", "M": "ATG", "V": "GTG",
             "S": "AGC", "P": "CCT", "T": "ACC", "A": "GCC", "Y": "TAT",
             "H": "CAT", "Q": "CAG", "N": "AAT", "K": "AAG", "D": "GAT",
             "E": "GAA", "C": "TGT", "W": "TGG", "R": "CGT", "G": "GGT"}
    keys = list(bases)
    codons = ["ATG"] + [bases[rng.choice(keys)] for _ in range(n_aa - 2)] + ["TAA"]
    return "".join(codons)


cds = make_cds(180)                      # 540 bp CDS
up = rand_dna(60, 0.42)
down = rand_dna(60, 0.58)
# 确保含常见酶切位点（直接植入）
up = up[:10] + "GAATTC" + up[16:25] + "GGATCC" + up[31:]
gene = up + cds + down
B = (f">SimGene_toxinX 模拟编码基因 toxinX (540bp CDS+侧翼) [教学用模拟序列]\n"
     + "\n".join(gene[i:i + 60] for i in range(0, len(gene), 60)))

# ---------------------------------------------------------------- C: 稡拟病毒演化
ancestor = rand_dna(420, 0.38)
rates = [0.02, 0.05, 0.08, 0.12, 0.16, 0.22]
labels = ["Virus_A", "Virus_B", "Virus_C", "Virus_D", "Virus_E", "Virus_F"]
# 两支系：A/B/C 一支，D/E/F 另一支（先各自派生中间祖先）
branch1 = mutate(ancestor, 0.06)
branch2 = mutate(ancestor, 0.07)
seqs_c = []
for i, lab in enumerate(labels):
    base = branch1 if i < 3 else branch2
    seqs_c.append(mutate(base, rates[i]))
C = "".join(f">{lab}_isolate{i+1} 模拟病毒分离株\n"
            + "\n".join(s[j:j + 60] for j in range(0, len(s), 60)) + "\n"
            for i, (lab, s) in enumerate(zip(labels, seqs_c)))

# ---------------------------------------------------------------- D: 同源基因对
anc = rand_dna(480, 0.52)
hom1 = mutate(anc, 0.04)
hom2 = mutate(anc, 0.13)
# 制造一小段插入/缺失，展示 gap
hom2 = hom2[:200] + "TTTTTT" + hom2[200:]
D = (">Homolog_gene_1 模拟同源基因1 [教学用模拟序列]\n"
     + "\n".join(hom1[i:i + 60] for i in range(0, len(hom1), 60)) + "\n"
     + ">Homolog_gene_2 模拟同源基因2 (含6bp插入) [教学用模拟序列]\n"
     + "\n".join(hom2[i:i + 60] for i in range(0, len(hom2), 60)))

# ---------------------------------------------------------------- E: 微生物群落 (3类群)
# (类群名, GC 含量, 突变率基数) —— 组内低突变、组间 GC 差异明显，保证聚类可分
groups = [("Firmicutes_like", 0.30, 0.045), ("Proteobacteria_like", 0.55, 0.055),
          ("GC_high_clade", 0.72, 0.05)]
entries = []
gi = 0
for gname, gc, rate in groups:
    core = rand_dna(360, gc)
    for m in range(3):
        gi += 1
        s = mutate(core, rate * (m + 1) * 0.5)
        entries.append((f"{gname}_sp{m+1}", s))
E = "".join(f">{name} 模拟宏基因组组装基因 (类群: {name.rsplit('_sp',1)[0]})\n"
            + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n"
            for name, s in entries)

files = {
    "示例A_人β珠蛋白HBB片段.fas": HBB,
    "示例B_模拟编码基因ORF.fas": B,
    "示例C_模拟病毒演化树.fas": C,
    "示例D_模拟同源基因对.fas": D,
    "示例E_模拟微生物群落聚类.fas": E,
}
for fn, content in files.items():
    with open(os.path.join(OUT, fn), "w", encoding="utf-8", newline="\n") as f:
        f.write(content.strip() + "\n")
    print(f"  ✓ {fn}  ({len(content)} bytes)")

print(f"\n示例数据已写入: {OUT}")
