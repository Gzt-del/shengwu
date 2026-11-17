/* ═══════════════════════════════════════════════════════════
   BioNova · 前端交互与可视化（零第三方库 · 原生 SVG 图表）
   ═══════════════════════════════════════════════════════════ */
"use strict";

const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];

const PALETTE = ["#38bdf8", "#a78bfa", "#34d399", "#fbbf24", "#fb7185",
                 "#22d3ee", "#c084fc", "#f472b6"];
const BASE_COLOR = { A: "#4ade80", T: "#fb7185", U: "#fb7185", G: "#fbbf24",
                     C: "#38bdf8", N: "#94a3b8" };

const PAGE_META = {
  stats:   ["序列概览", "序列组成、GC 含量与滑窗变异分析"],
  align:   ["序列比对", "Needleman-Wunsch 全局比对 / Smith-Waterman 局部比对"],
  orf:     ["ORF 与翻译", "双链六框翻译、开放阅读框与蛋白质理化性质"],
  enzyme:  ["酶切图谱", "限制性内切酶位点搜索与片段分析"],
  phylo:   ["系统发育", "遗传距离矩阵与邻接法（NJ）进化树"],
  cluster: ["序列聚类", "k-mer 距离 → PCoA 降维 → k-means 聚类"],
};

/* ───────────────────────── 基础设施 ───────────────────────── */

function toast(kind, title, msg = "", ms = 3400) {
  const host = $("#toastHost");
  const el = document.createElement("div");
  el.className = `toast ${kind}`;
  el.innerHTML = `<div><b></b><span></span></div>`;
  el.querySelector("b").textContent = title;
  el.querySelector("span").textContent = msg;
  host.appendChild(el);
  setTimeout(() => { el.classList.add("out"); setTimeout(() => el.remove(), 350); }, ms);
}

function setLoading(on) { $("#loading").hidden = !on; }

async function api(path, payload) {
  setLoading(true);
  try {
    const res = await fetch(path, {
      method: payload === undefined ? "GET" : "POST",
      headers: payload === undefined ? {} : { "Content-Type": "application/json" },
      body: payload === undefined ? undefined : JSON.stringify(payload),
    });
    const data = await res.json();
    if (data.error) throw new Error(data.error);
    return data;
  } catch (e) {
    toast("err", "分析失败", e.message || String(e), 4600);
    throw e;
  } finally {
    setLoading(false);
  }
}

function esc(s) {
  return String(s).replace(/[&<>"']/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}
function fmt(n) { return Number(n).toLocaleString("en-US"); }
function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }

/* ───────────────────────── SVG 工具 ───────────────────────── */

const NS = "http://www.w3.org/2000/svg";
function svg(tag, attrs = {}, parent = null) {
  const el = document.createElementNS(NS, tag);
  for (const k in attrs) el.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(el);
  return el;
}
function clearSvg(el) { while (el.firstChild) el.removeChild(el.firstChild); }

function niceTicks(min, max, n = 5) {
  if (!isFinite(min) || !isFinite(max)) return [0, 1];
  if (min === max) { min -= 1; max += 1; }
  const span = (max - min) / n;
  const mag = Math.pow(10, Math.floor(Math.log10(span)));
  const norm = span / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step)
    out.push(+v.toFixed(10));
  return out.length ? out : [min, max];
}

/* 通用折线/面积图 */
function lineChart(host, data, opts = {}) {
  host.innerHTML = "";
  if (!data || !data.length) { host.innerHTML = "<p class='input-hint'>暂无数据</p>"; return; }
  const W = 900, H = 300, ml = 52, mr = 18, mt = 16, mb = 40;
  const iw = W - ml - mr, ih = H - mt - mb;
  const xk = opts.x || "pos", yk = opts.y || "gc";
  const xs = data.map(d => d[xk]), ys = data.map(d => d[yk]);
  let x0 = Math.min(...xs), x1 = Math.max(...xs);
  let y0 = Math.min(...ys, 0), y1 = Math.max(...ys);
  if (opts.yMax !== undefined) y1 = opts.yMax;
  const pad = (y1 - y0) * 0.08 || 1;
  y1 += pad; if (opts.yMin === undefined) y0 -= pad; else y0 = opts.yMin;
  const sx = v => ml + (x1 === x0 ? iw / 2 : (v - x0) / (x1 - x0) * iw);
  const sy = v => mt + ih - (v - y0) / (y1 - y0) * ih;

  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` }, host);
  const defs = svg("defs", {}, s);
  const gid = "lg" + Math.random().toString(36).slice(2, 7);
  const grad = svg("linearGradient", { id: gid, x1: 0, y1: 0, x2: 0, y2: 1 }, defs);
  svg("stop", { offset: "0%", "stop-color": opts.color || "#38bdf8", "stop-opacity": .38 }, grad);
  svg("stop", { offset: "100%", "stop-color": opts.color || "#38bdf8", "stop-opacity": .02 }, grad);

  // 网格 + y 轴
  const yt = niceTicks(y0, y1, 5);
  yt.forEach(v => {
    svg("line", { x1: ml, x2: W - mr, y1: sy(v), y2: sy(v),
      stroke: "rgba(148,163,214,.12)", "stroke-width": 1 }, s);
    const t = svg("text", { x: ml - 8, y: sy(v) + 4, "text-anchor": "end",
      fill: "#6b7aa3", "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = (Math.abs(v) >= 1000 ? v.toFixed(0) : +v.toFixed(2));
  });
  // x 轴
  const xt = niceTicks(x0, x1, 6);
  xt.forEach(v => {
    const t = svg("text", { x: sx(v), y: H - 14, "text-anchor": "middle",
      fill: "#6b7aa3", "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = fmt(v);
    svg("line", { x1: sx(v), x2: sx(v), y1: mt, y2: mt + ih,
      stroke: "rgba(148,163,214,.07)" }, s);
  });
  svg("line", { x1: ml, x2: W - mr, y1: mt + ih, y2: mt + ih,
    stroke: "rgba(148,163,214,.35)", "stroke-width": 1.2 }, s);

  // 面积 + 折线
  let dArea = `M ${sx(xs[0])} ${sy(y0)}`;
  let dLine = "";
  data.forEach((d, i) => {
    dArea += ` L ${sx(d[xk])} ${sy(d[yk])}`;
    dLine += (i ? " L " : "M ") + sx(d[xk]) + " " + sy(d[yk]);
  });
  dArea += ` L ${sx(xs[xs.length - 1])} ${sy(y0)} Z`;
  svg("path", { d: dArea, fill: `url(#${gid})` }, s);
  svg("path", { d: dLine, fill: "none", stroke: opts.color || "#38bdf8",
    "stroke-width": 2.2, "stroke-linejoin": "round" }, s);

  // 零参考线
  if (y0 < 0 && y1 > 0)
    svg("line", { x1: ml, x2: W - mr, y1: sy(0), y2: sy(0),
      stroke: "rgba(251,191,36,.5)", "stroke-dasharray": "5 4" }, s);

  // 采样点 + 原生 tooltip
  const step = Math.max(1, Math.ceil(data.length / 160));
  for (let i = 0; i < data.length; i += step) {
    const d = data[i];
    const c = svg("circle", { cx: sx(d[xk]), cy: sy(d[yk]), r: 2.6,
      fill: opts.color || "#38bdf8", opacity: .85 }, s);
    const ti = svg("title", {}, c);
    ti.textContent = `${opts.xLabel || "位置"} ${d[xk]} → ${opts.yLabel || yk} ${d[yk]}`;
  }

  // 轴标题
  const ax = svg("text", { x: W - mr, y: H - 14, "text-anchor": "end",
    fill: "#46547c", "font-size": 11 }, s);
  ax.textContent = opts.xLabel || "";
  const ay = svg("text", { x: ml - 8, y: mt + 4, "text-anchor": "end",
    fill: "#46547c", "font-size": 11 }, s);
  ay.textContent = opts.yLabel || "";
}

/* 通用柱状图 */
function barChart(host, data, opts = {}) {
  host.innerHTML = "";
  if (!data || !data.length) { host.innerHTML = "<p class='input-hint'>暂无数据</p>"; return; }
  const W = 900, H = opts.height || 280, ml = 52, mr = 16, mt = 16, mb = opts.rotateX ? 74 : 42;
  const iw = W - ml - mr, ih = H - mt - mb;
  const xk = opts.x || "x", yk = opts.y || "y";
  const ys = data.map(d => +d[yk]);
  let y1 = Math.max(...ys, 0.0001), y0 = Math.min(...ys, 0);
  y1 += (y1 - y0) * 0.1 || 1;
  const sy = v => mt + ih - (v - y0) / (y1 - y0) * ih;
  const bw = iw / data.length;

  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` }, host);
  niceTicks(y0, y1, 5).forEach(v => {
    svg("line", { x1: ml, x2: W - mr, y1: sy(v), y2: sy(v),
      stroke: "rgba(148,163,214,.12)" }, s);
    const t = svg("text", { x: ml - 8, y: sy(v) + 4, "text-anchor": "end",
      fill: "#6b7aa3", "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = +v.toFixed(2);
  });
  svg("line", { x1: ml, x2: W - mr, y1: sy(Math.max(y0, 0)), y2: sy(Math.max(y0, 0)),
    stroke: "rgba(148,163,214,.35)" }, s);

  data.forEach((d, i) => {
    const v = +d[yk];
    const x = ml + i * bw + bw * 0.16;
    const w = bw * 0.68;
    const y = sy(Math.max(v, 0));
    const h = Math.abs(sy(v) - sy(Math.max(y0, 0)));
    const color = d.color || opts.color || PALETTE[i % PALETTE.length];
    const rect = svg("rect", { x, y, width: Math.max(w, 1.5), height: Math.max(h, 0.5),
      rx: 3.5, fill: color, opacity: .88 }, s);
    rect.style.transition = "opacity .15s";
    rect.addEventListener("mouseenter", () => rect.setAttribute("opacity", 1));
    rect.addEventListener("mouseleave", () => rect.setAttribute("opacity", .88));
    const ti = svg("title", {}, rect);
    ti.textContent = `${d[xk]}: ${v}`;

    const showLabel = data.length <= 30 || opts.alwaysLabel;
    if (showLabel) {
      const t = svg("text", {
        x: x + w / 2, y: opts.rotateX ? mt + ih + 14 : H - 22,
        "text-anchor": opts.rotateX ? "end" : "middle",
        transform: opts.rotateX ? `rotate(-45 ${x + w / 2} ${mt + ih + 14})` : "",
        fill: "#8ea0c8", "font-size": opts.rotateX ? 10.5 : 11,
        "font-family": "Consolas,monospace",
      }, s);
      t.textContent = d[xk];
    }
  });
}

/* 甜甜圈图 */
function donutChart(svgEl, items, opts = {}) {
  clearSvg(svgEl);
  const cx = 120, cy = 120, r = 86, sw = 30;
  const total = items.reduce((a, b) => a + b.count, 0) || 1;
  let angle = -90;
  svg("circle", { cx, cy, r, fill: "none", stroke: "rgba(148,163,214,.1)", "stroke-width": sw }, svgEl);
  items.forEach((it, i) => {
    const frac = it.count / total;
    const a = frac * 360;
    const c = it.color || PALETTE[i % PALETTE.length];
    const circ = 2 * Math.PI * r;
    const arc = svg("circle", {
      cx, cy, r, fill: "none", stroke: c, "stroke-width": sw,
      "stroke-dasharray": `${(a / 360) * circ} ${circ}`,
      "stroke-dashoffset": -(angle / 360) * circ,
      transform: `rotate(90 ${cx} ${cy})`,
      "stroke-linecap": "butt",
    }, svgEl);
    const ti = svg("title", {}, arc);
    ti.textContent = `${it.symbol}: ${it.count} (${it.pct}%)`;
    angle += a;
  });
  // 内圈
  svg("circle", { cx, cy, r: r - sw / 2 - 4, fill: "#0c1330", opacity: .0 }, svgEl);
}

/* 散点图（PCoA） */
function scatterChart(host, points, opts = {}) {
  host.innerHTML = "";
  const W = 900, H = 420, ml = 62, mr = 30, mt = 24, mb = 52;
  const iw = W - ml - mr, ih = H - mt - mb;
  if (!points.length) { host.innerHTML = "<p class='input-hint'>暂无数据</p>"; return; }
  const xs = points.map(p => p.x), ys = points.map(p => p.y);
  const allX = xs.concat((opts.centers || []).map(c => c[0]));
  const allY = ys.concat((opts.centers || []).map(c => c[1]));
  let x0 = Math.min(...allX), x1 = Math.max(...allX);
  let y0 = Math.min(...allY), y1 = Math.max(...allY);
  const px = (x1 - x0) * .14 || 1, py = (y1 - y0) * .14 || 1;
  x0 -= px; x1 += px; y0 -= py; y1 += py;
  const sx = v => ml + (v - x0) / (x1 - x0) * iw;
  const sy = v => mt + ih - (v - y0) / (y1 - y0) * ih;

  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` }, host);
  niceTicks(x0, x1, 6).forEach(v => {
    svg("line", { x1: sx(v), x2: sx(v), y1: mt, y2: mt + ih,
      stroke: "rgba(148,163,214,.09)" }, s);
    const t = svg("text", { x: sx(v), y: H - 26, "text-anchor": "middle",
      fill: "#6b7aa3", "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = +v.toFixed(2);
  });
  niceTicks(y0, y1, 5).forEach(v => {
    svg("line", { x1: ml, x2: W - mr, y1: sy(v), y2: sy(v),
      stroke: "rgba(148,163,214,.09)" }, s);
    const t = svg("text", { x: ml - 8, y: sy(v) + 4, "text-anchor": "end",
      fill: "#6b7aa3", "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = +v.toFixed(2);
  });
  // 坐标轴标题
  const lx = svg("text", { x: W - mr, y: H - 6, "text-anchor": "end",
    fill: "#8ea0c8", "font-size": 12 }, s);
  lx.textContent = opts.xLabel || "PC1";
  const ly = svg("text", { x: 14, y: mt + 6, fill: "#8ea0c8", "font-size": 12 }, s);
  ly.textContent = opts.yLabel || "PC2";

  points.forEach(p => {
    const c = PALETTE[(p.cluster || 0) % PALETTE.length];
    const g = svg("g", {}, s);
    const halo = svg("circle", { cx: sx(p.x), cy: sy(p.y), r: 13,
      fill: c, opacity: .14 }, g);
    const dot = svg("circle", { cx: sx(p.x), cy: sy(p.y), r: 6.5,
      fill: c, stroke: "#0b1224", "stroke-width": 1.6 }, g);
    const ti = svg("title", {}, dot);
    ti.textContent = `${p.label}  (PC1 ${p.x}, PC2 ${p.y})  簇 ${p.cluster + 1}`;
    const t = svg("text", { x: sx(p.x) + 11, y: sy(p.y) + 4,
      fill: "#cbd7f5", "font-size": 11 }, g);
    t.textContent = p.label.length > 22 ? p.label.slice(0, 21) + "…" : p.label;
    g.style.cursor = "default";
    dot.addEventListener("mouseenter", () => halo.setAttribute("r", 18));
    dot.addEventListener("mouseleave", () => halo.setAttribute("r", 13));
  });
  // 聚类中心
  (opts.centers || []).forEach((c, i) => {
    const g = svg("g", {}, s);
    const col = PALETTE[i % PALETTE.length];
    svg("path", {
      d: `M ${sx(c[0]) - 8} ${sy(c[1])} H ${sx(c[0]) + 8} M ${sx(c[0])} ${sy(c[1]) - 8} V ${sy(c[1]) + 8}`,
      stroke: col, "stroke-width": 2.4, opacity: .95,
    }, g);
    const ti = svg("title", {}, g);
    ti.textContent = `簇 ${i + 1} 中心 (${+c[0].toFixed(3)}, ${+c[1].toFixed(3)})`;
  });
}

/* ───────────────────────── 页面导航 ───────────────────────── */

function gotoPage(name) {
  $$(".nav-item").forEach(b => b.classList.toggle("active", b.dataset.page === name));
  $$(".page").forEach(p => p.classList.toggle("active", p.id === "page-" + name));
  const [t, d] = PAGE_META[name];
  $("#pageTitle").textContent = t;
  $("#pageDesc").textContent = d;
  window.scrollTo({ top: 0, behavior: "smooth" });
}
$$(".nav-item").forEach(b => b.addEventListener("click", () => gotoPage(b.dataset.page)));

function currentPage() {
  const el = $(".page.active");
  return el ? el.id.replace("page-", "") : "stats";
}

/* ───────────────────────── 示例数据 ───────────────────────── */

const SAMPLE_MAP = {
  stats: "statsInput", orf: "orfInput", enzyme: "enzInput",
  phylo: "phyloInput", cluster: "cluInput",
};

async function loadSamples() {
  try {
    const data = await api("/api/samples");
    const sel = $("#globalSample");
    data.samples.forEach(s => {
      const o = document.createElement("option");
      o.value = s.id; o.textContent = s.name;
      o.title = s.desc;
      sel.appendChild(o);
    });
    sel.addEventListener("change", () => {
      const s = data.samples.find(x => x.id === sel.value);
      if (!s) return;
      applySample(currentPage(), s);
      sel.value = "";
    });
  } catch (e) { /* 静默 */ }
}

function applySample(page, s) {
  const recs = parseFastaLite(s.content);
  if (page === "align") {
    if (recs.length < 2) { toast("info", "该示例不足两条序列", "比对模块需要 2 条序列"); return; }
    $("#alignA").value = ">" + recs[0].name + "\n" + wrap60(recs[0].seq);
    $("#alignB").value = ">" + recs[1].name + "\n" + wrap60(recs[1].seq);
  } else if (page === "orf" || page === "enzyme") {
    $("#" + SAMPLE_MAP[page]).value = s.content;
  } else if (SAMPLE_MAP[page]) {
    $("#" + SAMPLE_MAP[page]).value = s.content;
  }
  toast("ok", "已加载示例", s.name + " → " + PAGE_META[page][0]);
}

function parseFastaLite(text) {
  const out = [];
  let name = "Seq", buf = [];
  text.split(/\r?\n/).forEach(line => {
    if (line.startsWith(">")) {
      if (buf.length) out.push({ name, seq: buf.join("") });
      name = line.slice(1).trim().split(/\s+/)[0] || "Seq";
      buf = [];
    } else if (line.trim()) buf.push(line.trim());
  });
  if (buf.length) out.push({ name, seq: buf.join("") });
  return out;
}
function wrap60(s) { return s.match(/.{1,60}/g)?.join("\n") || s; }

/* ───────────────────────── 1. 序列概览 ───────────────────────── */

function tilesHtml(items) {
  return items.map(t => `
    <div class="tile">
      <div class="k">${esc(t.k)}</div>
      <div class="v">${t.v}<span class="u"> ${esc(t.u || "")}</span></div>
    </div>`).join("");
}

async function runStats() {
  const text = $("#statsInput").value.trim();
  if (!text) { toast("info", "请先输入序列"); return; }
  const data = await api("/api/stats", {
    text, window: +$("#gcWindow").value || 100, step: +$("#gcStep").value || 50,
  });
  const s = data.summary;
  $("#statsTiles").hidden = false;
  $("#statsTiles").innerHTML = tilesHtml([
    { k: "序列条数", v: s.count, u: "条" },
    { k: "总长度", v: fmt(s.total_length), u: "bp/aa" },
    { k: "平均长度", v: fmt(s.avg_length), u: "" },
    { k: "平均 GC", v: s.avg_gc ?? "—", u: s.avg_gc != null ? "%" : "" },
    { k: "最长", v: fmt(s.max_length), u: "" },
    { k: "最短", v: fmt(s.min_length), u: "" },
  ]);

  // 组成甜甜圈（第一条序列）
  const comp = data.composition[0];
  const colors = comp.kind === "nucleotide"
    ? { A: "#4ade80", T: "#fb7185", U: "#fb7185", G: "#fbbf24", C: "#38bdf8", N: "#94a3b8" }
    : {};
  const items = comp.items.map((it, i) => ({
    ...it, color: colors[it.symbol] || PALETTE[i % PALETTE.length],
  }));
  $("#compCard").hidden = false;
  donutChart($("#compDonut"), items);
  const top = items[0];
  $("#donutCenter").innerHTML =
    `<div class="big">${comp.kind === "nucleotide" ? (data.summary.avg_gc ?? 0) + "%" : comp.items.length}</div>
     <div class="small">${comp.kind === "nucleotide" ? "GC 含量" : "氨基酸种类"}<br>${esc(comp.name)}</div>`;
  $("#compLegend").innerHTML = items.map(it =>
    `<span><i style="background:${it.color}"></i>${it.symbol} <b>${it.count}</b> (${it.pct}%)</span>`
  ).join("");

  // 明细表
  $("#seqTableCard").hidden = false;
  $("#seqTable").innerHTML =
    `<thead><tr><th>#</th><th>名称</th><th>类型</th><th class="num">长度</th>
      <th class="num">GC%</th><th class="num">N</th><th class="num">缺口</th></tr></thead><tbody>` +
    s.detail.map((d, i) => `
      <tr><td>${i + 1}</td><td class="mono">${esc(d.name)}</td>
      <td><span class="tag ${d.type === "nucleotide" ? "dna" : "aa"}">${d.type === "nucleotide" ? "DNA/RNA" : "蛋白"}</span></td>
      <td class="num">${fmt(d.length)}</td>
      <td class="num">${d.gc ?? "—"}</td>
      <td class="num">${d.n_count ?? "—"}</td>
      <td class="num">${d.gap_count}</td></tr>`).join("") + "</tbody>";

  // GC 曲线 / Skew
  if (data.gc_curve && data.gc_curve.length) {
    $("#gcCard").hidden = false;
    $("#gcFocus").textContent =
      `序列 ${data.focus.name} · 窗口 ${data.focus.window} / 步长 ${data.focus.step}`;
    lineChart($("#gcChart"), data.gc_curve, {
      x: "pos", y: "gc", color: "#38bdf8", xLabel: "基因组位置 (bp)", yLabel: "GC%",
    });
  } else $("#gcCard").hidden = true;

  if (data.gc_skew && data.gc_skew.length > 1) {
    $("#skewCard").hidden = false;
    lineChart($("#skewChart"), data.gc_skew, {
      x: "pos", y: "skew", color: "#a78bfa", yMin: -1, yMax: 1,
      xLabel: "位置 (bp)", yLabel: "skew",
    });
  } else $("#skewCard").hidden = true;

  toast("ok", "分析完成", `${s.count} 条序列 · 用时 ${data._elapsed_ms} ms`);
}

/* ───────────────────────── 2. 序列比对 ───────────────────────── */

let ALIGN_MODE = "global";
$$("#alignMode button").forEach(b => b.addEventListener("click", () => {
  $$("#alignMode button").forEach(x => x.classList.remove("active"));
  b.classList.add("active");
  ALIGN_MODE = b.dataset.v;
}));

function renderAlignment(r) {
  const A = r.aligned_a, B = r.aligned_b;
  const MAX = 6000;                       // 展示上限，防止超长比对卡顿
  const shown = Math.min(A.length, MAX);
  const CHUNK = 60;
  let html = "";
  for (let s = 0; s < shown; s += CHUNK) {
    const e = Math.min(s + CHUNK, shown);
    let ra = "", rm = "", rb = "";
    for (let i = s; i < e; i++) {
      const a = A[i], b = B[i];
      const isGap = a === "-" || b === "-";
      const same = a === b && !isGap;
      const clsA = isGap ? "gap" : same ? "match" : "mismatch-a";
      const clsB = isGap ? "gap" : same ? "match" : "mismatch-b";
      ra += `<span class="base ${clsA}">${a === "-" ? "·" : esc(a)}</span>`;
      rb += `<span class="base ${clsB}">${b === "-" ? "·" : esc(b)}</span>`;
      rm += `<span class="base ${same ? "match" : "mid"}">${same ? "|" : isGap ? " " : "·"}</span>`;
    }
    html += `<div class="chunk-head">列 ${s + 1} – ${e}</div>`;
    html += `<div><span class="row-label">seq A</span>${ra}</div>`;
    html += `<div><span class="row-label"> </span>${rm}</div>`;
    html += `<div><span class="row-label">seq B</span>${rb}</div>`;
    if (e < shown) html += `<div style="color:#46547c;font-size:11px;margin:4px 0 8px 74px">— 第 ${e + 1} 列起省略（比对全长 ${fmt(A.length)} 列）—</div>`;
  }
  $("#alnView").innerHTML = html;
  $("#alnMeta").innerHTML = [
    `A 覆盖 ${r.start_a}–${r.end_a} / ${fmt(r.len_a)} (${r.coverage_a}%)`,
    `B 覆盖 ${r.start_b}–${r.end_b} / ${fmt(r.len_b)} (${r.coverage_b}%)`,
    `比对长度 ${fmt(r.aligned_len)} 列`,
    `模式 ${r.mode === "global" ? "全局 NW" : "局部 SW"}`,
  ].map(x => `<span>${x}</span>`).join("");
}

function renderDp(p) {
  const s = $("#dpHeat");
  clearSvg(s);
  const rows = p.grid.length, cols = p.grid[0].length;
  const W = 640, H = 420, ml = 44, mt = 40;
  const cw = Math.min(18, (W - ml - 10) / cols, (H - mt - 20) / rows);
  const lo = p.min, hi = p.max, span = (hi - lo) || 1;
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      const t = (p.grid[i][j] - lo) / span;
      const h = 214 + t * 66, sat = 45 + t * 45, light = 8 + t * 46;
      const rect = svg("rect", {
        x: ml + j * cw, y: mt + i * cw, width: cw - 1.2, height: cw - 1.2, rx: 2,
        fill: `hsl(${h} ${sat}% ${light}%)`,
      }, s);
      const ti = svg("title", {}, rect);
      ti.textContent = `H[${i}][${j}] = ${p.grid[i][j]}`;
    }
  }
  // 序列标签
  for (let j = 0; j < cols - 1; j++) {
    const t = svg("text", { x: ml + (j + .5) * cw, y: mt - 8, "text-anchor": "middle",
      fill: BASE_COLOR[p.b[j]] || "#8ea0c8", "font-size": 11,
      "font-family": "Consolas,monospace" }, s);
    t.textContent = p.b[j];
  }
  for (let i = 0; i < rows - 1; i++) {
    const t = svg("text", { x: ml - 7, y: mt + (i + .7) * cw, "text-anchor": "end",
      fill: BASE_COLOR[p.a[i]] || "#8ea0c8", "font-size": 11,
      "font-family": "Consolas,monospace" }, s);
    t.textContent = p.a[i];
  }
  const cap = svg("text", { x: ml, y: 18, fill: "#6b7aa3", "font-size": 11.5 }, s);
  cap.textContent = `DP 矩阵 ${rows}×${cols} 预览（分数范围 ${lo} ~ ${hi}）`;
  s.setAttribute("viewBox", `0 0 ${Math.max(ml + cols * cw + 12, 300)} ${Math.max(mt + rows * cw + 16, 160)}`);
}

async function runAlign() {
  const a = $("#alignA").value.trim(), b = $("#alignB").value.trim();
  if (!a || !b) { toast("info", "请输入两条序列"); return; }
  const r = await api("/api/align", {
    a, b, mode: ALIGN_MODE,
    match: +$("#alMatch").value || 0, mismatch: +$("#alMismatch").value || 0,
    gap: +$("#alGap").value || 0,
  });
  $("#alignTiles").hidden = false;
  $("#alignTiles").innerHTML = tilesHtml([
    { k: "比对得分", v: r.score, u: "" },
    { k: "一致性", v: r.identity, u: "%" },
    { k: "匹配", v: r.matches, u: "" },
    { k: "错配", v: r.mismatches, u: "" },
    { k: "空位", v: r.gaps, u: "bp" },
    { k: "比对长度", v: fmt(r.aligned_len), u: "列" },
  ]);
  $("#alnViewCard").hidden = false;
  renderAlignment(r);
  $("#dpCard").hidden = false;
  renderDp(r.matrix_preview);
  toast("ok", "比对完成", `得分 ${r.score} · 一致性 ${r.identity}%`);
  $("#copyAln").onclick = () => {
    copyText(`${$("#alnView").innerText}`);
  };
}

function copyText(t) {
  navigator.clipboard?.writeText(t).then(
    () => toast("ok", "已复制到剪贴板"),
    () => toast("err", "复制失败", "浏览器拒绝了剪贴板访问"),
  );
}

/* ───────────────────────── 3. ORF 与翻译 ───────────────────────── */

let ORF_DATA = null;

async function runOrf() {
  const text = $("#orfInput").value.trim();
  if (!text) { toast("info", "请先输入核酸序列"); return; }
  const starts = $$(".start-codon:checked").map(x => x.value);
  if (!starts.length) { toast("info", "请至少选择一个起始密码子"); return; }
  const data = await api("/api/orf", {
    text, min_aa: +$("#orfMin").value || 30, starts,
    both_strands: $("#orfBoth").checked,
  });
  ORF_DATA = data;
  $("#orfTiles").hidden = false;
  const tp = data.top_properties;
  $("#orfTiles").innerHTML = tilesHtml([
    { k: "序列长度", v: fmt(data.length), u: "bp" },
    { k: "ORF 数量", v: data.orf_count, u: "个" },
    { k: "最长 ORF", v: data.orfs[0]?.length_aa ?? "—", u: data.orfs[0] ? "aa" : "" },
    { k: "代表蛋白 MW", v: tp ? fmt(tp.molecular_weight) : "—", u: tp ? "Da" : "" },
    { k: "代表蛋白 pI", v: tp ? tp.isoelectric_point : "—", u: "" },
  ]);

  $("#orfTableCard").hidden = false;
  const rows = data.orfs;
  $("#orfTable").innerHTML =
    `<thead><tr><th>#</th><th>读框</th><th class="num">起点</th><th class="num">终点</th>
      <th class="num">长度(nt)</th><th class="num">肽长(aa)</th><th>起始</th>
      <th>蛋白质（前 36 aa）</th></tr></thead><tbody>` +
    (rows.length ? rows.map((o, i) => `
      <tr class="clickable ${i === 0 ? "selected" : ""}" data-i="${i}">
        <td>${i + 1}</td>
        <td><span class="tag ${o.strand === "+" ? "c1" : "c5"}">${o.frame}</span></td>
        <td class="num">${fmt(o.start)}</td><td class="num">${fmt(o.end)}</td>
        <td class="num">${fmt(o.length_nt)}</td><td class="num">${o.length_aa}</td>
        <td class="mono">${o.start_codon}</td>
        <td class="mono">${esc(o.protein.slice(0, 36))}${o.protein.length > 36 ? "…" : ""}</td>
      </tr>`).join("")
      : `<tr><td colspan="8" style="text-align:center;color:#6b7aa3">未找到满足条件的 ORF，可尝试降低最小肽长</td></tr>`) +
    "</tbody>";

  $$("#orfTable tr.clickable").forEach(tr =>
    tr.addEventListener("click", () => {
      $$("#orfTable tr").forEach(x => x.classList.remove("selected"));
      tr.classList.add("selected");
      showOrfDetail(+tr.dataset.i);
    }));
  if (rows.length) showOrfDetail(0);

  // 六框
  $("#orfDetailRow").hidden = false;
  $("#framesView").innerHTML = data.frames.map(f => `
    <div class="frame-row" title="${esc(f.protein)}">
      <div class="fl">${f.frame}</div>
      <div class="fp">${esc(f.protein.slice(0, 400)) || "（空）"}</div>
    </div>`).join("");

  // 密码子使用
  if (data.codon_usage.length) {
    $("#codonCard").hidden = false;
    barChart($("#codonChart"), data.codon_usage.slice(0, 20).map(c => ({
      x: c.codon, y: c.pct, color: "#818cf8",
    })), { x: "x", y: "y", rotateX: true, yLabel: "频率 %",
           height: 300, alwaysLabel: true });
  }
  toast("ok", "ORF 分析完成", `找到 ${data.orf_count} 个开放阅读框`);
}

function showOrfDetail(i) {
  const o = ORF_DATA.orfs[i];
  if (!o) return;
  $("#protTitle").innerHTML =
    `蛋白质 · ORF ${i + 1} <small>${o.frame} 链 · ${o.start}–${o.end} · ${o.length_aa} aa</small>`;
  const p = o.protein;
  $("#protView").innerHTML = esc(p).replace(/\*/g, '<span class="stop">*</span>')
    .replace(/^M/, '<span class="start-m">M</span>');
  // 服务端仅对最长 ORF 计算理化性质
  const t = ORF_DATA.top_properties;
  if (i === 0 && t) {
    $("#protTiles").innerHTML = tilesHtml([
      { k: "分子量", v: fmt(t.molecular_weight), u: "Da" },
      { k: "等电点", v: t.isoelectric_point, u: "pI" },
      { k: "GRAVY 疏水", v: t.gravy, u: "" },
      { k: "芳香族占比", v: t.aromaticity, u: "%" },
      { k: "pH7 净电荷", v: t.charge_ph7, u: "" },
    ]);
  } else {
    $("#protTiles").innerHTML =
      `<div class="tile"><div class="k">提示</div><div class="v" style="font-size:14px">理化性质仅对最长 ORF 计算</div></div>`;
  }
}

/* ───────────────────────── 4. 酶切图谱 ───────────────────────── */

let ENZ_SEL = new Set();
let ENZ_COLORS = {};

async function initEnzymes() {
  const data = await api("/api/enzyme_catalog");
  const grid = $("#enzGrid");
  data.catalog.forEach((e, i) => {
    ENZ_COLORS[e.name] = PALETTE[i % PALETTE.length];
    const chip = document.createElement("div");
    chip.className = "enz-chip";
    chip.dataset.name = e.name;
    chip.dataset.site = e.site.toLowerCase();
    chip.innerHTML = `<span class="nm">${esc(e.name)}</span>
      <span class="st">${e.site.replace("^", "<b>^</b>")}</span>`;
    chip.addEventListener("click", () => {
      if (ENZ_SEL.has(e.name)) ENZ_SEL.delete(e.name); else ENZ_SEL.add(e.name);
      chip.classList.toggle("on", ENZ_SEL.has(e.name));
      $("#enzCount").textContent = `已选 ${ENZ_SEL.size} 种`;
    });
    grid.appendChild(chip);
  });
  $("#enzSearch").addEventListener("input", ev => {
    const q = ev.target.value.trim().toLowerCase();
    $$(".enz-chip", grid).forEach(c => {
      c.style.display = (!q || c.dataset.name.toLowerCase().includes(q) ||
        c.dataset.site.includes(q)) ? "" : "none";
    });
  });
  $("#enzAll").addEventListener("click", () => {
    $$(".enz-chip", grid).forEach(c => {
      if (c.style.display !== "none") { ENZ_SEL.add(c.dataset.name); c.classList.add("on"); }
    });
    $("#enzCount").textContent = `已选 ${ENZ_SEL.size} 种`;
  });
  $("#enzNone").addEventListener("click", () => {
    ENZ_SEL.clear();
    $$(".enz-chip", grid).forEach(c => c.classList.remove("on"));
    $("#enzCount").textContent = "已选 0 种";
  });
  // 默认预选常用 4 种
  ["EcoRI", "BamHI", "HindIII", "PstI"].forEach(n => {
    ENZ_SEL.add(n);
    $(`.enz-chip[data-name="${n}"]`)?.classList.add("on");
  });
  $("#enzCount").textContent = `已选 ${ENZ_SEL.size} 种`;
}

async function runEnzyme() {
  const text = $("#enzInput").value.trim();
  if (!text) { toast("info", "请先输入 DNA 序列"); return; }
  if (!ENZ_SEL.size) { toast("info", "请至少选择一种酶"); return; }
  const data = await api("/api/enzyme", { text, enzymes: [...ENZ_SEL] });
  const dg = data.digest;

  renderEnzMap(dg);

  // 明细表
  $("#enzTableCard").hidden = false;
  $("#enzTable").innerHTML =
    `<thead><tr><th>酶</th><th>识别位点</th><th class="num">位点数</th><th>切点位置（左侧碱基数）</th>
      <th class="num">片段数</th></tr></thead><tbody>` +
    dg.enzymes.map(e => `
      <tr>
        <td><span style="color:${ENZ_COLORS[e.name]};font-weight:700">■</span> ${esc(e.name)}</td>
        <td class="mono">${esc(e.site)}</td>
        <td class="num">${e.count}</td>
        <td class="mono" style="max-width:260px;white-space:normal">${e.positions.join(", ") || "—"}</td>
        <td class="num">${e.fragments.length}</td>
      </tr>`).join("") + "</tbody>";

  // 片段分布
  $("#fragCard").hidden = false;
  barChart($("#fragChart"), dg.fragment_hist.map(f => ({
    x: f.len + "bp", y: f.count, color: "#34d399",
  })), { x: "x", y: "y", color: "#34d399", yLabel: "片段数", height: 280, alwaysLabel: true });

  toast("ok", "酶切完成",
    `共 ${dg.total_cuts} 个切点 · ${dg.combined_fragments.length} 个片段`);
}

function renderEnzMap(dg) {
  $("#mapCard").hidden = false;
  const s = $("#enzMap");
  clearSvg(s);
  const W = 960, ml = 46, mr = 46, y = 132;
  const iw = W - ml - mr;
  const L = dg.length || 1;
  const sx = p => ml + (p / L) * iw;

  const defs = svg("defs", {}, s);
  const g = svg("linearGradient", { id: "dnaGrad", x1: 0, x2: 1 }, defs);
  svg("stop", { offset: "0%", "stop-color": "#38bdf8" }, g);
  svg("stop", { offset: "50%", "stop-color": "#a78bfa" }, g);
  svg("stop", { offset: "100%", "stop-color": "#34d399" }, g);

  // 主链
  svg("rect", { x: ml, y: y - 5, width: iw, height: 10, rx: 5,
    fill: "url(#dnaGrad)", opacity: .9 }, s);
  // 端点标注
  [["5'", ml - 8, "end"], ["3'", ml + iw + 8, "start"]].forEach(([t, x, anch]) => {
    const el = svg("text", { x, y: y + 4, "text-anchor": anch, fill: "#6b7aa3",
      "font-size": 12, "font-family": "Consolas,monospace" }, s);
    el.textContent = t;
  });
  const lenT = svg("text", { x: ml + iw / 2, y: 22, "text-anchor": "middle",
    fill: "#8ea0c8", "font-size": 12.5 }, s);
  lenT.textContent = `序列长度 ${fmt(L)} bp · ${dg.enzymes.length} 种酶 · ${dg.total_cuts} 个切点`;

  // 切点与标签
  let row = 0;
  dg.enzymes.forEach(e => {
    const col = ENZ_COLORS[e.name] || "#38bdf8";
    e.positions.forEach((p, idx) => {
      const x = sx(p);
      svg("line", { x1: x, x2: x, y1: y - 34, y2: y + 14,
        stroke: col, "stroke-width": 2, "stroke-linecap": "round" }, s);
      svg("circle", { cx: x, cy: y, r: 3.4, fill: "#0b1224",
        stroke: col, "stroke-width": 2 }, s);
      const ly = 46 + ((row + idx) % 3) * 17;
      const tx = svg("text", { x: clamp(x, ml + 6, ml + iw - 6), y: ly,
        "text-anchor": "middle", fill: col, "font-size": 10.5,
        "font-family": "Consolas,monospace", opacity: .95 }, s);
      tx.textContent = `${e.name}(${p})`;
      svg("line", { x1: x, x2: x, y1: ly + 4, y2: y - 34,
        stroke: col, "stroke-width": .7, "stroke-dasharray": "2 3", opacity: .4 }, s);
    });
    row += Math.max(1, Math.ceil(e.positions.length / 3));
  });

  // 片段长度（全部酶共同作用）
  const frags = dg.combined_fragments;
  let acc = 0;
  const maxF = Math.max(...frags, 1);
  frags.forEach((f, i) => {
    const x1 = sx(acc), x2 = sx(acc + f);
    acc += f;
    const cy = y + 46;
    if (i > 0) svg("line", { x1, x2, y1: cy, y2: cy, stroke: "rgba(148,163,214,.4)",
      "stroke-dasharray": "3 3" }, s);
    if (f / maxF > 0.035 || frags.length <= 14) {
      const t = svg("text", { x: (x1 + x2) / 2, y: cy + 16, "text-anchor": "middle",
        fill: "#cbd7f5", "font-size": 11, "font-family": "Consolas,monospace" }, s);
      t.textContent = f;
    }
    // 刻度
    if (i > 0) svg("line", { x1, x2: x1, y1: y + 6, y2: cy + 4,
      stroke: "rgba(148,163,214,.35)" }, s);
  });
  const ft = svg("text", { x: ml, y: y + 86, fill: "#6b7aa3", "font-size": 11 }, s);
  ft.textContent = "↑ 数字为共同酶切后各片段长度 (bp)";

  // 图例
  let lx = ml;
  dg.enzymes.forEach(e => {
    const col = ENZ_COLORS[e.name] || "#38bdf8";
    const gEl = svg("g", {}, s);
    svg("rect", { x: lx, y: 208, width: 11, height: 11, rx: 3, fill: col }, gEl);
    const t = svg("text", { x: lx + 16, y: 217.5, fill: "#9aa8cc", "font-size": 11.5 }, gEl);
    t.textContent = `${e.name} (${e.count})`;
    lx += 26 + e.name.length * 7.4 + 24;
  });

  s.setAttribute("viewBox", `0 0 ${W} 232`);
}

/* ───────────────────────── 5. 系统发育 ───────────────────────── */

let TREE_DATA = null;

function heatColor(t) {
  const h = 214 + t * 66, sat = 45 + t * 45, light = 9 + t * 44;
  return `hsl(${h} ${sat}% ${light}%)`;
}

async function runPhylo() {
  const text = $("#phyloInput").value.trim();
  if (!text) { toast("info", "请先输入多条序列"); return; }
  const data = await api("/api/phylo", {
    text, model: $("#phyloModel").value,
  });
  if (data.error) throw new Error(data.error);
  TREE_DATA = data;

  // 距离矩阵
  $("#dmCard").hidden = false;
  $("#dmModel").textContent = data.model === "jc69" ? "Jukes-Cantor 校正" : "p-distance";
  const n = data.names.length;
  let flat = [];
  data.dist_matrix.forEach(r => flat.push(...r));
  const dmax = Math.max(...flat) || 1;
  let html = `<thead><tr><th></th>${data.names.map(x => `<th>${esc(x)}</th>`).join("")}</tr></thead><tbody>`;
  for (let i = 0; i < n; i++) {
    html += `<tr><th>${esc(data.names[i])}</th>`;
    for (let j = 0; j < n; j++) {
      const v = data.dist_matrix[i][j];
      const t = v / dmax;
      html += `<td class="num" style="background:${heatColor(t)};color:${t > .55 ? "#e8eefc" : "#9aa8cc"}">${v.toFixed(3)}</td>`;
    }
    html += "</tr>";
  }
  $("#dmTable").innerHTML = html + "</tbody>";

  // Newick
  $("#nwCard").hidden = false;
  $("#newickOut").textContent = data.newick;
  $("#copyNewick").onclick = () => copyText(data.newick);

  // 树
  $("#treeCard").hidden = false;
  drawTree();
  toast("ok", "建树完成", `${n} 条序列 · 邻接法 NJ`);
}

function drawTree() {
  if (!TREE_DATA || !TREE_DATA.tree) return;
  const useLen = $("#treeScale").checked;
  const s = $("#treeSvg");
  clearSvg(s);

  // 布局：节点 x = 根到该节点的累计枝长（cladogram 模式每层 +1）
  let leafY = 0;
  (function walk(nd, x0, isRoot) {
    const dx = isRoot ? 0 : (useLen ? (nd.length || 0) : 1);
    nd._x = x0 + dx;
    if (nd.leaf) {
      nd._y = leafY++;
    } else {
      nd.children.forEach(c => walk(c, nd._x, false));
      nd._y = (nd.children[0]._y + nd.children[nd.children.length - 1]._y) / 2;
    }
  })(TREE_DATA.tree, 0, true);

  const leaves = [];
  (function collect(nd) { if (nd.leaf) leaves.push(nd); else nd.children.forEach(collect); })(TREE_DATA.tree);

  const W = 960, H = Math.max(300, 90 + leaves.length * 40);
  const ml = 30, mr = 210, mt = 34, mb = 34;
  const iw = W - ml - mr, ih = H - mt - mb;
  const maxY = Math.max(leafY - 1, 1);
  let maxX = 0;
  (function mx(nd) { maxX = Math.max(maxX, nd._x); if (!nd.leaf) nd.children.forEach(mx); })(TREE_DATA.tree);
  const xSpan = maxX || 1;
  const sx = v => ml + (v / xSpan) * iw;
  const sy = v => mt + (maxY ? v / maxY : 0.5) * ih;

  s.setAttribute("viewBox", `0 0 ${W} ${H}`);

  // 背景网格（横）
  for (let i = 0; i <= leaves.length; i++) {
    svg("line", { x1: ml, x2: W - mr + 150, y1: sy(i - .5), y2: sy(i - .5),
      stroke: "rgba(148,163,214,.05)" }, s);
  }

  const stroke = "url(#treeGrad)";
  const defs = svg("defs", {}, s);
  const tg = svg("linearGradient", { id: "treeGrad", x1: 0, x2: 1 }, defs);
  svg("stop", { offset: "0%", "stop-color": "#38bdf8" }, tg);
  svg("stop", { offset: "100%", "stop-color": "#a78bfa" }, tg);

  // 枝条
  function drawBranches(nd) {
    const px = sx(nd._x), py = sy(nd._y);
    if (!nd.leaf) {
      nd.children.forEach(c => {
        const cx = sx(c._x), cy = sy(c._y);
        svg("path", {
          d: `M ${px} ${py} L ${cx} ${py} L ${cx} ${cy}`,
          fill: "none", stroke, "stroke-width": 2.1,
          "stroke-linecap": "round", "stroke-linejoin": "round",
        }, s);
        const ti = svg("title", {}, s.lastChild);
        ti.textContent = `${c.name || "内部节点"} · 枝长 ${c.length}`;
        drawBranches(c);
      });
      svg("circle", { cx: px, cy: py, r: 3, fill: "#0b1224",
        stroke: "#818cf8", "stroke-width": 1.8 }, s);
    } else {
      // 叶标签
      const t = svg("text", { x: cx2(nd), y: py + 4.5, fill: "#e8eefc",
        "font-size": 13, "font-weight": 600 }, s);
      t.textContent = nd.name || "Leaf";
      const d = svg("text", { x: cx2(nd) + 8 + (nd.name || "").length * 7.6, y: py + 4.5,
        fill: "#6b7aa3", "font-size": 10.5, "font-family": "Consolas,monospace" }, s);
      d.textContent = useLen ? (nd.length || 0).toFixed(4) : "";
      svg("circle", { cx: sx(nd._x), cy: py, r: 4.5,
        fill: "#38bdf8", stroke: "#0b1224", "stroke-width": 1.6 }, s);
    }
  }
  function cx2(nd) { return sx(nd._x) + 10; }

  // 根
  svg("circle", { cx: sx(0), cy: sy(TREE_DATA.tree._y), r: 5,
    fill: "#0b1224", stroke: "#a78bfa", "stroke-width": 2.2 }, s);
  const rt = svg("text", { x: sx(0), y: sy(TREE_DATA.tree._y) - 14,
    "text-anchor": "middle", fill: "#a78bfa", "font-size": 11 }, s);
  rt.textContent = "root";

  drawBranches(TREE_DATA.tree);

  // 比例尺
  if (useLen && maxX > 0) {
    const barVal = niceScale(xSpan);
    const barPx = (barVal / xSpan) * iw;
    const y0 = H - 16;
    svg("line", { x1: ml, x2: ml + barPx, y1: y0, y2: y0,
      stroke: "#8ea0c8", "stroke-width": 2.4 }, s);
    svg("line", { x1: ml, x2: ml, y1: y0 - 5, y2: y0 + 5, stroke: "#8ea0c8", "stroke-width": 2 }, s);
    svg("line", { x1: ml + barPx, x2: ml + barPx, y1: y0 - 5, y2: y0 + 5,
      stroke: "#8ea0c8", "stroke-width": 2 }, s);
    const t = svg("text", { x: ml + barPx + 8, y: y0 + 4, fill: "#8ea0c8",
      "font-size": 11, "font-family": "Consolas,monospace" }, s);
    t.textContent = `${barVal} 替换/位点`;
  }
}
function niceScale(span) {
  const mag = Math.pow(10, Math.floor(Math.log10(span)));
  const n = span / mag;
  return +(n < 1.5 ? 0.1 : n < 3 ? 0.2 : n < 7 ? 0.5 : 1) * mag;
}

$("#treeScale")?.addEventListener("change", drawTree);

/* ───────────────────────── 6. 序列聚类 ───────────────────────── */

async function runCluster() {
  const text = $("#cluInput").value.trim();
  if (!text) { toast("info", "请先输入多条序列"); return; }
  const k = clamp(+$("#cluK").value || 3, 2, 8);
  const data = await api("/api/cluster", { text, k });
  if (data.error) throw new Error(data.error);

  // PCoA 散点
  $("#pcoaCard").hidden = false;
  $("#pcoaExplained").textContent =
    `PC1 解释 ${data.explained[0]}% · PC2 解释 ${data.explained[1]}% · ${data.method}`;
  const points = data.coords.map((c, i) => ({
    x: c[0], y: c[1], label: data.names[i], cluster: data.labels[i],
  }));
  scatterChart($("#pcoaChart"), points, {
    centers: data.centers, xLabel: `PC1 (${data.explained[0]}%)`,
    yLabel: `PC2 (${data.explained[1]}%)`,
  });
  $("#cluLegend").innerHTML = Array.from({ length: data.k }, (_, i) =>
    `<span><i style="display:inline-block;width:11px;height:11px;border-radius:50%;background:${PALETTE[i % PALETTE.length]};margin-right:5px"></i>簇 ${i + 1}</span>`
  ).join("");

  // 结果表
  $("#cluTableCard").hidden = false;
  $("#cluTable").innerHTML =
    `<thead><tr><th>#</th><th>序列名称</th><th>簇</th>
      <th class="num">PC1</th><th class="num">PC2</th></tr></thead><tbody>` +
    data.names.map((nm, i) => `
      <tr><td>${i + 1}</td><td class="mono">${esc(nm)}</td>
      <td><span class="tag c${(data.labels[i] % 5) + 1}">簇 ${data.labels[i] + 1}</span></td>
      <td class="num">${data.coords[i][0].toFixed(4)}</td>
      <td class="num">${data.coords[i][1].toFixed(4)}</td></tr>`).join("") +
    "</tbody>";

  // 距离矩阵
  $("#cluDmCard").hidden = false;
  const n = data.names.length;
  let flat = [];
  data.dist_matrix.forEach(r => flat.push(...r));
  const dmax = Math.max(...flat) || 1;
  let html = `<thead><tr><th></th>${data.names.map(x => `<th>${esc(x.slice(0, 10))}</th>`).join("")}</tr></thead><tbody>`;
  for (let i = 0; i < n; i++) {
    html += `<tr><th>${esc(data.names[i].slice(0, 12))}</th>`;
    for (let j = 0; j < n; j++) {
      const v = data.dist_matrix[i][j];
      html += `<td class="num" style="background:${heatColor(v / dmax)}">${v.toFixed(3)}</td>`;
    }
    html += "</tr>";
  }
  $("#cluDm").innerHTML = html + "</tbody>";

  toast("ok", "聚类完成", `${n} 条序列 → ${data.k} 个簇`);
}

/* ───────────────────────── 事件绑定 ───────────────────────── */

$("#runStats").addEventListener("click", runStats);
$("#runAlign").addEventListener("click", runAlign);
$("#runOrf").addEventListener("click", runOrf);
$("#runEnzyme").addEventListener("click", runEnzyme);
$("#runPhylo").addEventListener("click", runPhylo);
$("#runCluster").addEventListener("click", runCluster);

$$("[data-clear]").forEach(b => b.addEventListener("click", () => {
  const el = $("#" + b.dataset.clear);
  if (el) { el.value = ""; el.focus(); toast("info", "已清空输入"); }
}));

// 使用说明弹窗
$("#btnHelp").addEventListener("click", () => $("#helpModal").hidden = false);
$("#helpClose").addEventListener("click", () => $("#helpModal").hidden = true);
$("#helpModal").addEventListener("click", e => {
  if (e.target === $("#helpModal")) $("#helpModal").hidden = true;
});
document.addEventListener("keydown", e => {
  if (e.key === "Escape") $("#helpModal").hidden = true;
  if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) {
    const map = { stats: runStats, align: runAlign, orf: runOrf,
                  enzyme: runEnzyme, phylo: runPhylo, cluster: runCluster };
    map[currentPage()]?.();
  }
});

/* ───────────────────────── 启动 ───────────────────────── */

(async function boot() {
  loadSamples();
  initEnzymes();
  // 健康检查
  try {
    const h = await fetch("/api/health").then(r => r.json());
    if (h.ok) console.log(`[BioNova] 服务就绪 v${h.version}`);
  } catch (e) { console.warn("[BioNova] 健康检查失败", e); }
})();
