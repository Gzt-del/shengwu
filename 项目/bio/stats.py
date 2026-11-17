# -*- coding: utf-8 -*-
"""序列统计分析：碱基/氨基酸组成、GC 含量、滑窗 GC、GC-skew、密码子使用偏性等。"""

import re
from collections import Counter, OrderedDict

from .fasta import infer_seq_type

NUCS = ["A", "C", "G", "T", "U", "N"]
AA20 = list("ACDEFGHIKLMNPQRSTVWY")


def clean(seq: str) -> str:
    """去掉缺口与非字母字符，转大写"""
    return re.sub(r"[^A-Za-z*]", "", seq).upper()


def gc_content(seq: str) -> float:
    """GC 含量百分比（忽略 N 与缺口）"""
    s = clean(seq)
    g, c = s.count("G"), s.count("C")
    atgc = g + c + s.count("A") + s.count("T") + s.count("U")
    return round((g + c) * 100.0 / atgc, 3) if atgc else 0.0


def atgc_percent(seq: str) -> dict:
    s = clean(seq)
    total = sum(s.count(x) for x in "ACGTU") or 1
    return {
        "A": round(s.count("A") * 100.0 / total, 3),
        "T": round((s.count("T") + s.count("U")) * 100.0 / total, 3),
        "G": round(s.count("G") * 100.0 / total, 3),
        "C": round(s.count("C") * 100.0 / total, 3),
    }


def composition(seq: str, kind: str = "nucleotide") -> list:
    """组成统计 [{symbol, count, pct}]"""
    s = clean(seq)
    cnt = Counter(s)
    total = sum(cnt.values()) or 1
    if kind == "nucleotide":
        keys = [k for k in NUCS if k in cnt] + sorted(k for k in cnt if k not in NUCS)
    else:
        keys = [k for k in AA20 if k in cnt] + sorted(k for k in cnt if k not in AA20)
    return [
        {"symbol": k, "count": cnt[k], "pct": round(cnt[k] * 100.0 / total, 3)}
        for k in keys
    ]


def sliding_gc(seq: str, window: int = 100, step: int = 50) -> list:
    """滑窗 GC 曲线 [{pos, gc}]，pos 为窗口中心（1-based）；前缀和实现 O(n)"""
    s = clean(seq)
    n = len(s)
    if n == 0:
        return []
    window = max(1, min(window, n))
    step = max(1, step)
    pref = [0] * (n + 1)
    for i, ch in enumerate(s):
        pref[i + 1] = pref[i] + (1 if ch in "GC" else 0)
    out = []
    for start in range(0, n - window + 1, step):
        gc = pref[start + window] - pref[start]
        out.append({"pos": start + window // 2 + 1, "gc": round(gc * 100.0 / window, 2)})
    if not out:
        out.append({"pos": n, "gc": gc_content(s)})
    return out


def gc_skew(seq: str, window: int = 100, step: int = 50) -> list:
    """GC-skew = (G-C)/(G+C)，用于复制起点/基因岛分析"""
    s = clean(seq)
    n = len(s)
    if n == 0:
        return []
    window = max(1, min(window, n))
    step = max(1, step)
    pref_g, pref_c = [0] * (n + 1), [0] * (n + 1)
    for i, ch in enumerate(s):
        pref_g[i + 1] = pref_g[i] + (ch == "G")
        pref_c[i + 1] = pref_c[i] + (ch == "C")
    out = []
    start = 0
    while start < n:
        end = min(start + window, n)
        g = pref_g[end] - pref_g[start]
        c = pref_c[end] - pref_c[start]
        skew = (g - c) / (g + c) if (g + c) else 0.0
        out.append({"pos": (start + end) // 2 + 1, "skew": round(skew, 4)})
        if end == n:
            break
        start += step
    return out


def seq_summary(records: list) -> dict:
    """多序列总览统计"""
    n = len(records)
    lengths = [len(r["seq"]) for r in records]
    types = Counter(infer_seq_type(r["seq"]) for r in records)
    dna = [r for r in records if infer_seq_type(r["seq"]) == "nucleotide"]
    detail = []
    for r in records:
        seq = r["seq"]
        t = infer_seq_type(seq)
        detail.append({
            "name": r["name"],
            "desc": r.get("desc", ""),
            "type": t,
            "length": len(seq),
            "gc": gc_content(seq) if t == "nucleotide" else None,
            "n_count": seq.upper().count("N") if t == "nucleotide" else None,
            "gap_count": seq.count("-"),
        })
    return {
        "count": n,
        "total_length": sum(lengths),
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "avg_length": round(sum(lengths) / n, 1) if n else 0,
        "avg_gc": round(sum(gc_content(r["seq"]) for r in dna) / len(dna), 3) if dna else None,
        "types": dict(types),
        "detail": detail,
    }


def codon_usage(seq: str) -> list:
    """密码子使用频率表（读框 1，仅统计完整密码子）"""
    s = clean(seq)
    cnt = Counter(s[i:i + 3] for i in range(0, len(s) - 2, 3))
    total = sum(cnt.values()) or 1
    return [
        {"codon": c, "count": k, "pct": round(k * 100.0 / total, 3)}
        for c, k in sorted(cnt.items(), key=lambda kv: -kv[1])
    ]


def kmer_freq(seq: str, k: int = 3, top: int = 20) -> list:
    """k-mer 频率表（去除非 ACGT 字符）"""
    s = re.sub(r"[^ACGTU]", "", clean(seq))
    if not s:
        return []
    k = max(1, min(k, len(s)))
    cnt = Counter(s[i:i + k] for i in range(len(s) - k + 1))
    total = sum(cnt.values()) or 1
    rows = sorted(cnt.items(), key=lambda kv: -kv[1])[:top]
    return [{"kmer": km, "count": c, "pct": round(c * 100.0 / total, 3)} for km, c in rows]


def composition_bundle(records: list) -> dict:
    """一次请求返回：总览 + 每条序列组成 + 首条序列滑窗 GC / GC-skew"""
    if not records:
        return {}
    first = records[0]["seq"]
    t = infer_seq_type(first)
    win = 100 if len(first) >= 200 else max(10, len(first) // 5)
    step = max(1, win // 2)
    return {
        "summary": seq_summary(records),
        "composition": [
            {"name": r["name"], "kind": infer_seq_type(r["seq"]),
             "items": composition(r["seq"], infer_seq_type(r["seq"]))}
            for r in records
        ],
        "gc_curve": sliding_gc(first, win, step) if t == "nucleotide" else [],
        "gc_skew": gc_skew(first, win, step) if t == "nucleotide" else [],
        "focus": {"name": first and records[0]["name"], "window": win, "step": step},
    }
