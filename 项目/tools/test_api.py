# -*- coding: utf-8 -*-
"""API 端点集成测试：python tools/test_api.py [base_url]"""
import json
import os
import sys
import urllib.request

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8765"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=30) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def post(path, payload):
    req = urllib.request.Request(
        BASE + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.status, json.loads(r.read().decode("utf-8"))


def sample(name_part):
    for fn in os.listdir(os.path.join(HERE, "samples")):
        if name_part in fn:
            return open(os.path.join(HERE, "samples", fn), encoding="utf-8").read()
    raise FileNotFoundError(name_part)


def main():
    ok = 0

    # 静态资源
    for path, expect in [("/", "text/html"), ("/static/css/style.css", "text/css"),
                         ("/static/js/app.js", "javascript"),
                         ("/static/index.html", "text/html")]:
        status, ctype, body = get(path)
        assert status == 200 and expect in ctype and len(body) > 1000, (path, status, ctype, len(body))
        ok += 1
        print(f"  GET {path} -> {status} {ctype.split(';')[0]} ({len(body)} bytes)")

    # health / samples（中文文件名需正确 UTF-8）
    _, _, body = get("/api/health")
    assert json.loads(body)["ok"]
    ok += 1
    _, _, body = get("/api/samples")
    samples = json.loads(body.decode("utf-8"))["samples"]
    assert len(samples) == 5, samples
    assert any("示例A" in s["id"] for s in samples), samples[0]["id"]
    print(f"  GET /api/samples -> 5 items, 中文名 OK: {samples[0]['id']}")
    ok += 1

    # stats（真实 HBB 示例）
    st, d = post("/api/stats", {"text": sample("示例A"), "window": 30, "step": 15})
    assert d["summary"]["count"] == 1 and d["summary"]["avg_gc"] > 40, d.get("summary")
    assert d["gc_curve"] and d["composition"]
    print(f"  POST /api/stats -> GC={d['summary']['avg_gc']}%, "
          f"gc_curve={len(d['gc_curve'])} pts, {d['_elapsed_ms']} ms")
    ok += 1

    # align（同源基因对）
    def split_fasta(text):
        out, name, buf = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    out.append(">" + name + "\n" + "\n".join(buf))
                name, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            out.append(">" + name + "\n" + "\n".join(buf))
        return out

    blocks = split_fasta(sample("示例D"))
    assert len(blocks) == 2, len(blocks)
    st, d = post("/api/align", {"a": blocks[0], "b": blocks[1], "mode": "global",
                                "match": 2, "mismatch": 1, "gap": 2})
    assert d["identity"] > 70, d["identity"]
    assert len(d["aligned_a"]) == len(d["aligned_b"])
    print(f"  POST /api/align -> score={d['score']}, identity={d['identity']}%, "
          f"aln_len={d['aligned_len']}, {d['_elapsed_ms']} ms")
    ok += 1

    st, d = post("/api/align", {"a": "ACGTACGTACGT", "b": "TTTTACGTACGTAAAA",
                                "mode": "local", "match": 2, "mismatch": 1, "gap": 1})
    assert d["score"] == 18, d["score"]
    print(f"  POST /api/align local -> score={d['score']} (期望 18)")
    ok += 1

    # orf（模拟编码基因）
    st, d = post("/api/orf", {"text": sample("示例B"), "min_aa": 50,
                              "starts": ["ATG"], "both_strands": True})
    assert d["orf_count"] >= 1, d["orf_count"]
    top = d["orfs"][0]
    assert top["length_aa"] >= 100, top
    assert d["top_properties"]["isoelectric_point"] > 0
    print(f"  POST /api/orf -> {d['orf_count']} ORFs, top={top['length_aa']} aa "
          f"({top['start']}-{top['end']}), MW={d['top_properties']['molecular_weight']}")
    ok += 1

    # enzyme
    st, d = post("/api/enzyme", {"text": sample("示例B"),
                                 "enzymes": ["EcoRI", "BamHI", "HindIII", "PstI"]})
    dg = d["digest"]
    assert dg["total_cuts"] >= 2, dg
    assert dg["combined_fragments"] and sum(dg["combined_fragments"]) == dg["length"]
    print(f"  POST /api/enzyme -> {dg['total_cuts']} cuts, "
          f"{len(dg['combined_fragments'])} frags, sum OK")
    ok += 1

    # phylo（模拟病毒 6 条）
    st, d = post("/api/phylo", {"text": sample("示例C"), "model": "jc69"})
    assert "tree" in d and "newick" in d, d.get("error")
    assert len(d["names"]) == 6
    assert d["newick"].count(",") == 5
    print(f"  POST /api/phylo -> {d['newick'][:70]}...")
    ok += 1

    # cluster（微生物群落 9 条 → 3 簇；三个类群应各自聚成一簇）
    st, d = post("/api/cluster", {"text": sample("示例E"), "k": 3})
    assert len(d["labels"]) == 9 and len(set(d["labels"])) == 3, d["labels"]
    labels = d["labels"]
    groups_ok = (labels[0] == labels[1] == labels[2] and
                 labels[3] == labels[4] == labels[5] and
                 labels[6] == labels[7] == labels[8] and
                 len({labels[0], labels[3], labels[6]}) == 3)
    assert groups_ok, f"三个类群未各自聚簇: {labels}"
    assert d["explained"][0] > 0
    print(f"  POST /api/cluster -> labels={labels}, 类群各自成簇 ✓, "
          f"explained={d['explained']}, {d['_elapsed_ms']} ms")
    ok += 1

    # translate
    st, d = post("/api/translate", {"text": sample("示例A")})
    assert d["properties"]["molecular_weight"] > 0
    print(f"  POST /api/translate -> best frame {d['best_frame']}, "
          f"MW={d['properties']['molecular_weight']}")
    ok += 1

    # 错误处理：业务校验错误以 HTTP 200 + error 字段返回
    st, d = post("/api/align", {"a": "", "b": "ACGT"})
    assert "error" in d, d
    st, d = post("/api/phylo", {"text": ">a\nACGT\n>b\nACGT", "model": "jc69"})
    assert "error" in d, d
    print(f"  错误处理 -> 空序列/少于3条序列均返回 error 字段 ✓")
    ok += 1

    print(f"\nALL {ok} API TESTS PASSED ✓")


if __name__ == "__main__":
    main()
