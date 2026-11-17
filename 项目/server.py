# -*- coding: utf-8 -*-
"""
BioNova · 生物信息学智能分析平台
纯 Python 标准库 HTTP 服务器（零第三方依赖），双击 启动项目.bat 即可运行。

路由:
  GET  /            → 前端单页应用
  GET  /api/samples → 内置示例数据列表
  GET  /api/health  → 健康检查
  POST /api/<module> → 各分析模块 JSON 接口
"""

import json
import os
import sys
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from bio import align as A          # noqa: E402
from bio import enzyme as E         # noqa: E402
from bio import fasta as F          # noqa: E402
from bio import genetic as G        # noqa: E402
from bio import phylo as P          # noqa: E402
from bio import stats as S          # noqa: E402

# 控制台中文输出兼容 Windows GBK
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

PORT_RANGE = (8765, 8790)
STATIC_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
    ".woff2": "font/woff2",
    ".md": "text/markdown; charset=utf-8",
    ".fas": "text/plain; charset=utf-8",
    ".fasta": "text/plain; charset=utf-8",
    ".txt": "text/plain; charset=utf-8",
}


# ----------------------------------------------------------------- API 处理

def api_samples(payload):
    """列出内置示例"""
    d = os.path.join(ROOT, "samples")
    items = []
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            if fn.lower().endswith((".fas", ".fasta", ".fa", ".txt")):
                path = os.path.join(d, fn)
                with open(path, "r", encoding="utf-8") as f:
                    content = f.read()
                desc = ""
                for line in content.splitlines():
                    if line.startswith("#") or line.startswith(">"):
                        desc = line.lstrip("#> ").strip()
                        break
                items.append({"id": fn, "name": os.path.splitext(fn)[0], "desc": desc,
                              "content": content})
    return {"samples": items}


def api_stats(payload):
    text = payload.get("text", "")
    records = F.parse_fasta(text)
    err = F.validate(records)
    if err:
        return {"error": err}
    win = int(payload.get("window", 100) or 100)
    step = int(payload.get("step", 50) or 50)
    out = S.composition_bundle(records)
    # 用户可覆盖滑窗参数；序列过短时自动缩小窗口，保证曲线有足够数据点
    first = records[0]["seq"]
    if F.infer_seq_type(first) == "nucleotide" and win > 0:
        if win >= len(first):
            win = max(10, len(first) // 4)
            step = max(1, win // 2)
        out["gc_curve"] = S.sliding_gc(first, win, step)
        out["gc_skew"] = S.gc_skew(first, win, step)
        out["focus"] = {"name": records[0]["name"], "window": win, "step": step}
    return out


def api_align(payload):
    a = payload.get("a", "").strip()
    b = payload.get("b", "").strip()
    mode = payload.get("mode", "global")
    try:
        match = float(payload.get("match", 1))
        mismatch = float(payload.get("mismatch", 1))
        gap = float(payload.get("gap", 1))
    except (TypeError, ValueError):
        return {"error": "打分参数必须是数字。"}
    if not a or not b:
        return {"error": "请提供两条待比对序列。"}
    # 若传入 FASTA，取第一条
    if ">" in a:
        ra = F.parse_fasta(a)
        a = ra[0]["seq"] if ra else ""
    if ">" in b:
        rb = F.parse_fasta(b)
        b = rb[0]["seq"] if rb else ""
    return A.align(a, b, mode=mode, match=match, mismatch=mismatch, gap=gap)


def api_orf(payload):
    text = payload.get("text", "")
    records = F.parse_fasta(text)
    err = F.validate(records, need_dna=True)
    if err:
        return {"error": err}
    rec = records[0]
    seq = rec["seq"].upper()
    try:
        min_aa = int(payload.get("min_aa", 30))
    except (TypeError, ValueError):
        min_aa = 30
    min_aa = max(5, min(min_aa, 2000))
    starts = payload.get("starts") or ["ATG"]
    if isinstance(starts, str):
        starts = [s.strip().upper() for s in starts.split(",") if s.strip()]
    both = bool(payload.get("both_strands", True))

    orfs = G.find_orfs(seq, min_aa=min_aa, start_codons=starts, both_strands=both)
    frames = G.six_frame(seq)

    # 主 ORF（最长）附蛋白质理化性质
    top_props = None
    if orfs:
        top_props = G.protein_properties(orfs[0]["protein"])

    return {
        "name": rec["name"],
        "length": len(seq),
        "orf_count": len(orfs),
        "orfs": orfs[:40],
        "frames": frames,
        "top_properties": top_props,
        "codon_usage": S.codon_usage(seq)[:30],
    }


def api_enzyme(payload):
    text = payload.get("text", "")
    records = F.parse_fasta(text)
    err = F.validate(records, need_dna=True)
    if err:
        return {"error": err}
    seq = records[0]["seq"].upper()
    enzymes = payload.get("enzymes") or []
    if isinstance(enzymes, str):
        enzymes = [x.strip() for x in enzymes.split(",") if x.strip()]
    if not enzymes:
        return {"error": "请至少选择一种限制性内切酶。"}
    return {
        "name": records[0]["name"],
        "digest": E.digest(seq, enzymes),
        "catalog": E.ENZYMES,
    }


def api_enzyme_catalog(payload):
    return {"catalog": E.ENZYMES}


def api_phylo(payload):
    records = F.parse_fasta(payload.get("text", ""))
    err = F.validate(records, need_dna=True)
    if err:
        return {"error": err}
    if len(records) < 3:
        return {"error": "系统发育建树至少需要 3 条序列（建议 4 条以上，可多选示例数据）。"}
    model = payload.get("model", "jc69")
    if model not in ("jc69", "p"):
        model = "jc69"
    names = [r["name"] for r in records]
    seqs = [r["seq"].upper() for r in records]
    return P.phylo_analysis(names, seqs, model=model)


def api_cluster(payload):
    records = F.parse_fasta(payload.get("text", ""))
    err = F.validate(records, need_dna=True)
    if err:
        return {"error": err}
    if len(records) < 3:
        return {"error": "聚类分析至少需要 3 条序列。"}
    try:
        k = int(payload.get("k", 3))
    except (TypeError, ValueError):
        k = 3
    k = max(2, min(k, 8))
    names = [r["name"] for r in records]
    seqs = [r["seq"].upper() for r in records]
    return P.cluster_analysis(names, seqs, k=k)


def api_translate(payload):
    """整段序列六框翻译 + 蛋白质性质（序列类型可为 DNA）"""
    records = F.parse_fasta(payload.get("text", ""))
    err = F.validate(records, need_dna=True)
    if err:
        return {"error": err}
    rec = records[0]
    frames = G.six_frame(rec["seq"])
    # 选最长且不含大量 X 的框作为代表
    best = max(frames, key=lambda f: len(f["protein"].replace("X", "")))
    props = G.protein_properties(best["protein"])
    return {"name": rec["name"], "frames": frames, "best_frame": best["frame"],
            "properties": props}


ROUTES = {
    "/api/stats": api_stats,
    "/api/align": api_align,
    "/api/orf": api_orf,
    "/api/enzyme": api_enzyme,
    "/api/enzyme_catalog": api_enzyme_catalog,
    "/api/phylo": api_phylo,
    "/api/cluster": api_cluster,
    "/api/translate": api_translate,
}


# ----------------------------------------------------------------- HTTP

class Handler(BaseHTTPRequestHandler):
    server_version = "BioNova/1.0"
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):
        # 精简访问日志
        pass

    def _send(self, code: int, body: bytes, ctype: str, cache: bool = False):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        if cache:
            self.send_header("Cache-Control", "public, max-age=3600")
        else:
            self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        try:
            self.wfile.write(body)
        except BrokenPipeError:
            pass

    def _json(self, obj, code: int = 200):
        data = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self._send(code, data, "application/json; charset=utf-8")

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/samples":
            try:
                return self._json(api_samples({}))
            except Exception as ex:
                return self._json({"error": str(ex)}, 500)
        if path == "/api/health":
            return self._json({"ok": True, "app": "BioNova", "version": "1.0.0"})
        if path == "/api/enzyme_catalog":
            return self._json(api_enzyme_catalog({}))
        return self._static(path)

    def do_POST(self):
        path = urlparse(self.path).path
        fn = ROUTES.get(path)
        if not fn:
            return self._json({"error": f"未知接口: {path}"}, 404)
        try:
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b"{}"
            payload = json.loads(raw.decode("utf-8") or "{}")
        except Exception:
            return self._json({"error": "请求体必须是合法的 JSON。"}, 400)
        try:
            t0 = time.time()
            result = fn(payload)
            result = result if isinstance(result, dict) else {"data": result}
            result["_elapsed_ms"] = int((time.time() - t0) * 1000)
            return self._json(result)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            return self._json({"error": f"分析失败: {ex}"}, 500)

    def _static(self, path: str):
        if path in ("/", "/index.html"):
            rel = os.path.join("static", "index.html")
        else:
            rel = os.path.normpath(path.lstrip("/")).replace("\\", "/")
            if rel.startswith("..") or rel.startswith("/"):
                return self._send(403, b"forbidden", "text/plain")
            # 仅允许访问 static/ 与 samples/ 目录
            top = rel.split("/")[0]
            if top not in ("static", "samples"):
                rel = os.path.join("static", rel)
        fp = os.path.join(ROOT, rel)
        if not os.path.isfile(fp):
            return self._send(404, b"404 Not Found", "text/plain; charset=utf-8")
        ext = os.path.splitext(fp)[1].lower()
        ctype = STATIC_TYPES.get(ext, "application/octet-stream")
        with open(fp, "rb") as f:
            body = f.read()
        self._send(200, body, ctype, cache=ext not in (".html",))


def pick_port() -> int:
    import socket
    for p in range(PORT_RANGE[0], PORT_RANGE[1] + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    raise RuntimeError("端口 8765-8790 均被占用，请关闭占用进程后重试。")


BANNER = r"""
============================================================
   BioNova · 生物信息学智能分析平台  v1.0
   Sequence Analysis & Molecular Evolution Workbench
------------------------------------------------------------
   纯 Python 标准库实现 · 零第三方依赖 · 离线可用
   算法: Needleman-Wunsch / Smith-Waterman / 邻接法(NJ)
         六框翻译 / ORF 查找 / 酶切图谱 / PCoA + k-means
============================================================
"""


def main():
    print(BANNER)
    port = pick_port()
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}"
    print(f"  服务已启动: {url}")
    print(f"  按 Ctrl+C 停止服务\n")

    def open_browser():
        time.sleep(0.8)
        try:
            webbrowser.open(url)
        except Exception:
            pass

    threading.Thread(target=open_browser, daemon=True).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n服务已停止，感谢使用 BioNova。")
        httpd.server_close()


if __name__ == "__main__":
    main()
