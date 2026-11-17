# -*- coding: utf-8 -*-
"""FASTA / 纯序列文本解析模块。"""

import re

DNA_RE = re.compile(r"[ACGTNRYSWKMBDHV]", re.I)
AA_RE = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBXZOU*]", re.I)

# 常见蛋白质氨基酸（用于判定序列类型）
PROTEIN_ONLY = set("EFILPQXZ*")


def _clean_seq(raw: str) -> str:
    """去掉空白、数字与常见 FASTA 外部符号，保留字母与 * - ."""
    return re.sub(r"[\s0-9\-\.'\"`~]", "", raw).upper()


def parse_fasta(text: str) -> list:
    """解析 FASTA 或裸序列文本，返回 [{name, desc, seq}, ...]

    - 支持标准 FASTA（>header）
    - 也支持直接粘贴的裸序列（自动命名 Seq-1 / Seq-2 ...）
    """
    if text is None:
        return []
    text = text.replace("\ufeff", "").strip()
    if not text:
        return []

    records = []
    if ">" in text:
        cur_name, cur_desc, cur_buf = None, "", []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None or cur_buf:
                    records.append(_make_record(cur_name, cur_desc, cur_buf))
                header = line[1:].strip()
                parts = header.split(None, 1)
                cur_name = parts[0] if parts else "unnamed"
                cur_desc = parts[1] if len(parts) > 1 else ""
                cur_buf = []
            else:
                if cur_name is None:
                    cur_name, cur_desc = "Seq-1", ""
                cur_buf.append(line)
        records.append(_make_record(cur_name, cur_desc, cur_buf))
    else:
        # 裸序列：按空行分块，或整体视为一条
        blocks = [b for b in re.split(r"\n\s*\n", text) if b.strip()]
        if len(blocks) > 1:
            for i, b in enumerate(blocks, 1):
                records.append(_make_record(f"Seq-{i}", "", [b]))
        else:
            records.append(_make_record("Seq-1", "", [text]))

    records = [r for r in records if r["seq"]]
    return records


def _make_record(name, desc, buf) -> dict:
    seq = _clean_seq("\n".join(buf))
    return {"name": name or "unnamed", "desc": desc or "", "seq": seq}


def infer_seq_type(seq: str) -> str:
    """判定序列类型: nucleotide / protein / unknown"""
    s = seq.upper().replace("*", "").replace("-", "")
    if not s:
        return "unknown"
    letters = set(s)
    if letters <= set("ACGTUN"):
        return "nucleotide"
    if letters <= set("ACGTURYKMSWBDHVN"):
        # 只有核酸简并碱基 → 核酸
        return "nucleotide"
    if letters & PROTEIN_ONLY:
        return "protein"
    # 无法区分时按核酸简并处理
    dna_frac = sum(s.count(c) for c in "ACGTUN") / len(s)
    return "nucleotide" if dna_frac > 0.85 else "protein"


def revcomp(seq: str) -> str:
    """DNA 反向互补"""
    table = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")
    return seq.upper().translate(table)[::-1]


def validate(records: list, need_dna: bool = False) -> str:
    """校验解析结果，返回错误信息（空串表示通过）"""
    if not records:
        return "未检测到任何序列，请检查输入格式（FASTA 或裸序列均可）。"
    for r in records:
        if not r["seq"]:
            return f"序列 {r['name']} 内容为空。"
        if len(r["seq"]) > 200000:
            return f"序列 {r['name']} 长度过长（>200kb），请拆分后分析。"
        if need_dna:
            t = infer_seq_type(r["seq"])
            if t != "nucleotide":
                return f"序列 {r['name']} 不是核酸序列（判定为 {t}），该分析仅支持 DNA/RNA。"
    return ""
