# -*- coding: utf-8 -*-
"""后端算法冒烟测试：python tools/test_smoke.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bio import align as A
from bio import enzyme as E
from bio import fasta as F
from bio import genetic as G
from bio import phylo as P
from bio import stats as S

def main():
    # 1. FASTA 解析
    recs = F.parse_fasta(">a desc\nATGCATGC\n>b\nTTTTGGGG")
    assert len(recs) == 2 and recs[0]["seq"] == "ATGCATGC", recs
    print("1  fasta OK")

    # 2. 统计
    b = S.composition_bundle(recs)
    assert b["summary"]["count"] == 2
    assert abs(S.gc_content("GGCCAT") - 66.667) < 0.01
    print("2  stats OK, gc_curve points:", len(b["gc_curve"]))

    # 3. NW / SW 比对（用回溯结果反算分数，校验自洽）
    r = A.align("GATTACA", "GCATCG", "global", 1, 1, 1)
    assert r["aligned_a"] and len(r["aligned_a"]) == len(r["aligned_b"])
    re_score = (r["matches"] * 1 - r["mismatches"] * 1 - r["gaps"] * 1)
    assert abs(re_score - r["score"]) < 1e-6, (re_score, r["score"])
    print("3  NW score:", r["score"], "identity:", r["identity"])
    print("     ", r["aligned_a"], "\n     ", r["aligned_b"])

    # 暴力枚举所有比对方式，验证 NW 分数确为最优
    def brute(a, b, match=1, mism=1, gap=1):
        best = [float("-inf")]
        def rec(i, j, sc):
            if i == len(a) and j == len(b):
                best[0] = max(best[0], sc)
                return
            if i < len(a) and j < len(b):
                rec(i + 1, j + 1, sc + (match if a[i] == b[j] else -mism))
            if i < len(a):
                rec(i + 1, j, sc - gap)
            if j < len(b):
                rec(i, j + 1, sc - gap)
        rec(0, 0, 0)
        return best[0]
    bf = brute("GATTACA", "GCATCG")
    assert abs(bf - r["score"]) < 1e-9, (bf, r["score"])
    print("   brute-force optimum:", bf, "✓ DP 最优")

    r2 = A.align("ACGTACGTACGT", "TTTTACGTACGTAAAA", "local", 2, 1, 1)
    assert r2["score"] > 0
    re2 = r2["matches"] * 2 - r2["mismatches"] * 1 - r2["gaps"] * 1
    assert abs(re2 - r2["score"]) < 1e-6, (re2, r2["score"])
    print("   SW score:", r2["score"], "span:",
          r2["start_a"], "-", r2["end_a"])

    # 4. 翻译 / ORF / 蛋白质性质
    prot = G.translate(
        "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTG")
    assert prot.startswith("MVHLTPEEKSAVTALWGKVNV"), prot
    orf = G.find_orfs("ATGAAACCCGGGTTTTAG", min_aa=3, both_strands=False)
    assert orf and orf[0]["length_aa"] == 5, orf
    props = G.protein_properties("MKVLAAAG")
    assert props["molecular_weight"] > 500
    print("4  translate:", prot)
    print("   ORF:", orf[0]["start"], "-", orf[0]["end"],
          f"({orf[0]['length_aa']} aa), MW:", props["molecular_weight"],
          "pI:", props["isoelectric_point"])

    # 5. 酶切（切点 = 左侧碱基数；AGG|AATTC 在第 3 个碱基后切开）
    sites = E.find_sites("AGGAATTCCCT", "G^AATTC")
    assert sites and sites[0]["pos"] == 3, sites
    d = E.digest("AGGAATTCGGGATCCTTT", ["EcoRI", "BamHI"])
    assert d["total_cuts"] == 2, d
    assert d["combined_fragments"] == [3, 7, 8], d["combined_fragments"]
    print("5  enzyme OK, cuts:", d["total_cuts"],
          "fragments:", d["combined_fragments"])

    # 6. 距离矩阵 + NJ 建树
    names = ["s1", "s2", "s3", "s4"]
    seqs = ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGA",
            "ACGTACGTACGTACGTTTTT", "TTGTACGTACGTACGTACGT"]
    dm = P.distance_matrix(seqs, "jc69")
    nj = P.neighbor_joining(names, dm["matrix"])
    assert nj["newick"].count("(") >= 2
    print("6  newick:", nj["newick"])

    # 7. PCoA + k-means
    km = P.kmer_jaccard_matrix(seqs, 3)
    pc = P.pcoa(km)
    cl = P.kmeans(pc["coords"], k=2)
    print("7  pcoa:", pc["coords"], "explained:", pc["explained"],
          "labels:", cl["labels"])

    # 8. 示例文件全链路
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sdir = os.path.join(here, "samples")
    for fn in sorted(os.listdir(sdir)):
        txt = open(os.path.join(sdir, fn), encoding="utf-8").read()
        rs = F.parse_fasta(txt)
        print(f"   sample {fn}: {len(rs)} seqs, lens="
              f"{[len(x['seq']) for x in rs]}")

    # 9. 用真实示例数据跑完整分析
    txt = open(os.path.join(sdir, "示例C_模拟病毒演化树.fas"),
               encoding="utf-8").read()
    rs = F.parse_fasta(txt)
    out = P.phylo_analysis([x["name"] for x in rs],
                           [x["seq"] for x in rs], model="jc69")
    assert "tree" in out, out
    print("9  phylo on sample C OK")

    txt = open(os.path.join(sdir, "示例B_模拟编码基因ORF.fas"),
               encoding="utf-8").read()
    rs = F.parse_fasta(txt)
    orfs = G.find_orfs(rs[0]["seq"], min_aa=50, both_strands=False)
    assert orfs, "示例B 应能找到 ORF"
    print("   ORF on sample B:", orfs[0]["start"], "-",
          orfs[0]["end"], orfs[0]["length_aa"], "aa")

    print("\nALL BACKEND TESTS PASSED ✓")


if __name__ == "__main__":
    main()
