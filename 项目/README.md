# BioNova · 生物信息学智能分析平台

> **方向：序列分析与分子进化（Sequence Analysis & Molecular Evolution）**
> 福州大学 · 人工智能专业 · 生物信息学课程项目

BioNova 是一个**开箱即用的本地 Web 分析平台**，围绕"序列分析与分子进化"方向集成了六大学习与研究中最常用的分析模块。整个项目**只依赖 Python 3.8+ 标准库**（零第三方包、零联网、数据不出本机），前端图表全部为自研原生 SVG 实现，双击一个 `.bat` 文件即可启动。

---

## 一、快速开始

### 运行环境
- Windows / macOS / Linux
- Python **3.8 及以上**（本项目在 Python 3.11 下开发测试）
- 任意现代浏览器（Chrome / Edge / Firefox）

### 启动方式

| 方式 | 操作 |
|------|------|
| **方式一（推荐）** | 双击 `启动项目.bat`，浏览器会自动打开 `http://127.0.0.1:8765` |
| **方式二** | 命令行执行 `python -X utf8 server.py` |

停止服务：在命令行窗口按 `Ctrl+C`（或直接关闭窗口）。

> 端口被占用时会自动尝试 8765 ~ 8790 区间内的空闲端口。

### 三步上手
1. 点击顶部 **「加载示例数据…」**，选择一个内置示例；
2. 点击 **「开始分析 / 运行比对 / 构建进化树」** 等按钮；
3. 查看统计瓦片、图表、表格与可视化结果。

---

## 二、功能模块

| 模块 | 核心功能 | 关键算法 |
|------|----------|----------|
| **1 序列概览** | 序列组成统计、GC 含量、滑窗 GC 曲线、GC-Skew、多序列明细 | 前缀和滑窗 O(n) |
| **2 序列比对** | 全局/局部双序列比对、打分参数可调、比对视图、DP 矩阵热图 | **Needleman-Wunsch / Smith-Waterman**（线性 gap） |
| **3 ORF 与翻译** | 双链六框翻译、ORF 查找（可选起始密码子）、蛋白质分子量/等电点/GRAVY、密码子使用频率 | 标准遗传密码表、二分法求 pI（Bjellqvist pKa） |
| **4 酶切图谱** | 26 种常用限制酶、位点搜索（支持 IUPAC 简并碱基）、彩色图谱、片段长度分布 | 位点表达式 `G^AATTC` 解析与匹配 |
| **5 系统发育** | 遗传距离矩阵、进化树可视化、Newick 导出 | p-distance / **Jukes-Cantor** 校正 + **邻接法 Neighbor-Joining** |
| **6 序列聚类** | k-mer 距离矩阵、PCoA 二维坐标图、k-means 聚类结果 | **Jaccard 距离 → PCoA（Jacobi 特征分解）→ k-means++** |

机器学习视角：模块 6 把每条序列编码为 3-mer 特征，用无监督方法回答"哪些序列相似"，是生物学与人工智能的交叉点。

---

## 三、项目结构

```
生物信息学/项目/
├── server.py                 # HTTP 服务器 + 全部 API 路由（标准库 http.server）
├── 启动项目.bat               # 一键启动（自动选端口、自动开浏览器）
├── requirements.txt          # 依赖说明（无第三方依赖）
├── README.md                 # 本文件
├── bio/                      # 核心算法包（纯 Python）
│   ├── __init__.py
│   ├── fasta.py              # FASTA / 裸序列解析、序列类型判定、反向互补
│   ├── stats.py              # 组成统计、滑窗 GC、GC-skew、密码子频率
│   ├── align.py              # NW / SW 动态规划比对 + 回溯
│   ├── genetic.py             # 密码子表、六框翻译、ORF 查找、蛋白质理化性质
│   ├── enzyme.py              # 限制酶位点搜索与酶切片段计算
│   └── phylo.py               # 距离矩阵、NJ 建树、PCoA、k-means、Newick
├── static/                   # 前端（无框架、无第三方 JS 库）
│   ├── index.html            # 单页应用结构
│   ├── css/style.css         # 深空蓝玻璃拟态主题
│   └── js/app.js             # 交互逻辑 + 全部自研 SVG 图表
├── samples/                  # 内置示例数据
│   ├── 示例A_人β珠蛋白HBB片段.fas      # 真实序列（前 93 bp）
│   ├── 示例B_模拟编码基因ORF.fas       # 模拟（ORF / 酶切演示）
│   ├── 示例C_模拟病毒演化树.fas         # 模拟（NJ 建树演示）
│   ├── 示例D_模拟同源基因对.fas         # 模拟（比对演示）
│   └── 示例E_模拟微生物群落聚类.fas     # 模拟（聚类演示）
└── tools/
    ├── generate_samples.py   # 示例数据生成脚本（固定随机种子，可复现）
    └── test_smoke.py         # 后端算法冒烟测试
```

---

## 四、API 接口

服务器同时提供 REST 风格 JSON 接口，方便二次开发或写报告时调用：

| 方法 | 路径 | 说明 |
|------|------|------|
| GET | `/api/health` | 健康检查 |
| GET | `/api/samples` | 获取全部内置示例 |
| GET | `/api/enzyme_catalog` | 限制酶目录 |
| POST | `/api/stats` | `{text, window, step}` 序列统计 |
| POST | `/api/align` | `{a, b, mode, match, mismatch, gap}` 序列比对 |
| POST | `/api/orf` | `{text, min_aa, starts, both_strands}` ORF 查找 |
| POST | `/api/translate` | `{text}` 六框翻译 + 蛋白质性质 |
| POST | `/api/enzyme` | `{text, enzymes[]}` 酶切分析 |
| POST | `/api/phylo` | `{text, model}` 距离矩阵 + NJ 建树 |
| POST | `/api/cluster` | `{text, k}` PCoA + k-means 聚类 |

示例：

```bash
curl -X POST http://127.0.0.1:8765/api/align \
  -H "Content-Type: application/json" \
  -d "{\"a\":\"GATTACA\",\"b\":\"GCATCG\",\"mode\":\"global\"}"
```

---

## 五、算法原理摘要（可用于实验报告）

### 1. Needleman-Wunsch 全局比对
递推式（线性 gap 罚分 `g`）：

```
H(i,j) = max { H(i-1,j-1) + s(aᵢ,bⱼ)    匹配/错配
              { H(i-1,j)   - g            上方引入 gap
              { H(i,j-1)   - g            左侧引入 gap
```

边界 `H(0,j) = -j·g`，`H(i,0) = -i·g`；自右下角沿最大分方向回溯得到最优比对。
本实现时间/空间复杂度 **O(mn)**，并通过暴力枚举验证了小规模用例的最优性。

### 2. Smith-Waterman 局部比对
与 NW 的两点差异：① 边界全为 0；② 递推结果若 ≤0 则置 0；③ 自矩阵最大值处回溯，遇到 0 终止。

### 3. Jukes-Cantor 距离校正
观测差异比例 `p`，模型假设四种碱基替换速率相等：

```
d = -¾ · ln(1 - 4p/3)
```

当 `p ≥ 0.75` 时视为饱和，赋一个大值避免 `ln` 无定义。

### 4. 邻接法（Neighbor-Joining, Saitou & Nei 1987）
1. 计算净分歧率 `rᵢ = Σⱼ dᵢⱼ / (m-2)`；
2. 找 `Q(i,j) = (m-2)·dᵢⱼ - rᵢ - rⱼ` 最小的节点对；
3. 枝长 `Lᵢ = ½dᵢⱼ + ½(rᵢ-rⱼ)`，`Lⱼ = dᵢⱼ - Lᵢ`（负值截断为 0）；
4. 合并新节点，距离 `d(u,k) = ½(dᵢₖ + dⱼₖ - dᵢⱼ)`，回到第 1 步；
5. 剩余两节点等分相连，输出 Newick 格式。

### 5. PCoA（经典多维标度）
1. 由 3-mer 计数向量计算加权 Jaccard（Ruzicka）距离 `d = 1 - Σmin(fₐ,f_b) / Σmax(fₐ,f_b)`（短序列上集合版 Jaccard 会因 k-mer 空间饱和而失效，故采用计数版）；
2. 双中心化 `B = -½ · J D² J`；
3. 对实对称矩阵 `B` 做 **Jacobi 旋转法**特征分解，取前两个特征向量作为二维坐标；
4. 解释方差 = λᵢ / Σλ。

### 6. k-means++
随机种子固定（seed=42），k-means++ 概率初始化 + Lloyd 迭代，保证结果可复现。

---

## 六、测试

```bash
python -X utf8 tools/test_smoke.py
```

覆盖：FASTA 解析、统计、NW/SW（含暴力枚举最优性校验）、翻译、ORF 坐标、酶切片段长度、NJ 建树、PCoA/k-means、全部示例数据的端到端链路。

---

## 七、数据说明

- **示例 A**：人 β-珠蛋白基因（HBB）编码区前 93 bp —— 真实经典教学片段，翻译产物与教材一致（`MVHLTPEEKSAVTALWGKVNV…`）；
- **示例 B–E**：由 `tools/generate_samples.py` 以固定随机种子生成的**教学用模拟序列**，不对应任何真实生物，仅用于演示算法流程。

---

## 八、可扩展方向

- [ ] 多序列比对（Clustal 风格渐进比对 + 引导树）
- [ ] BLAST 风格的本地快速相似检索（k-mer 索引 + 打分）
- [ ] 系统发育 bootstrap 支持率评估
- [ ] RNA 二级结构预测（最小自由能 / 动态规划）
- [ ] 基于机器学习的序列分类（如 GC 含量 + k-mer 特征 → 随机森林）
- [ ] 导出 PNG/SVG 报告与分析历史记录

---

## 九、许可与声明

本项目为课程学习用途开源示例，算法均按教科书标准实现并附测试。
界面字体与配色为自研设计；示例数据来源见第七节。
