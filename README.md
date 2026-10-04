# eDNA 固定参考树短序列放置与群落距离后端

基于 FastAPI + Biopython 的环境 DNA 短序列放置（phylogenetic placement）服务：
在**固定参考拓扑和原始枝长**的前提下，把每条待测序列独立地枚举到参考树的
每一条枝上，用 JC69 模型通过 Felsenstein 递推最大化整树对齐的对数似然，
并以 jplace version 3 格式交付结果。在此之上提供 `POST /compare` 群落比较：
把不同采样地点的读数计数映射为树上的质量分布，计算两两之间
**p=1 的 Kantorovich–Rubinstein（树 Wasserstein-1 / 推土机）距离**。

## 范围

- 参考序列 3–20 条，待测序列 1–5 条；所有序列等长且 ≤ 500 列。
- 允许 ACGT 与 IUPAC 歧义码（R/Y/S/W/K/M/B/D/H/V/N）；DNA 不接受 RNA 碱基
  `U`（会定位到列并整单拒绝），缺口（`-` `.` `?`）视为未知、不提供碱基证据。
- 空序列、非法字符、重名、非正/非有限枝长、树叶与参考 ID 不一一对应等
  问题会被**定位**（记录名、列号或树节点）并**整单拒绝**（HTTP 422）。
- 对每条待测序列：枚举所有枝，插入节点把原枝拆成两段（均非负、总长不变），
  新叶枝长在 [0, 2] 内优化（粗网格 + 模式搜索细化），返回全部枝的最优位置、
  新枝长、对数似然与归一化的相对似然权重（like_weight_ratio）。
- 权重是跨枝的**相对支持度**，不是物种归属正确的概率。
- 完全无有效碱基证据（全缺口/全 N）的待测序列不强制归属，返回无法放置原因。
- 计算采用逐列缩放避免长序列下溢；不使用汉明距离或最近叶启发式。

## 模块协作

- `app/validation.py` — FASTA/Newick 解析与定位校验（整单拒绝）
- `app/tree.py` — 无根二叉树图、稳定枝号、带 `{edge_num}` 的 jplace 树序列化
- `app/likelihood.py` — JC69 转移矩阵、Felsenstein 双向消息（rerooting）与逐列缩放
- `app/placement.py` — 逐枝优化（拆分点 + 新枝长）、权重归一化与排序
- `app/jplacefmt.py` — jplace v3 文档组装（五个标准放置字段）
- `app/jplacetree.py` — 从 jplace 带枝号 Newick 还原有向树几何（叶名含逗号时
  按引号语法解析）
- `app/community.py` — 读数→插入点质量分配、按可放置读数归一化、p=1 KR 距离
  与逐枝贡献
- `app/main.py` — FastAPI 入口：`POST /place`、`GET /place/{job_id}/jplace`、
  `POST /compare`

叶 ID 若含 `, ( ) : ; [ ] =`、空白等特殊字符，导出 Newick 时自动按
Newick 引号规则包裹（内部单引号翻倍转义），例如 `sp,A` 导出为 `'sp,A'`，
不会再在逗号处被切成两个叶。

## 群落距离 `POST /compare`

请求体为 2–10 个唯一样本（保持输入顺序），每样本：

- `sample_id`：样本（地点）标识，同请求内不可重复；
- `job_id`：已有放置任务 ID，**不会重新放置或修改任务**，只读取已存结果；
- `counts`：待测 ID → 非负整数读数映射；遗漏的待测 ID 视为 0。

整单拒绝（HTTP 422，定位到样本）的情形：样本数越界、`sample_id` 重复、
未知 `job_id`、映射中出现该任务未知的待测 ID、非法（负/非整数）计数、
引用任务的带枝号 Newick 文本不一致（只允许在同一棵参考树上比较）、
某样本归一化后的有效总数为 0（全部读数都无法放置）。

### 质量构造

对每条读数，按其全部候选放置的 `like_weight_ratio` 把计数质量分配到
对应参考枝的插入位置——**不是只取最佳枝，也不平铺/移动到枝端**。
插入位置按 jplace 的 `distal_length` 解释（相对远端子节点），
`pendant_length` 不参与几何。无法放置的 ID（如全缺口序列）连同
`id`、`count` 与原因列入 `excluded` 并排除；随后每个样本按
**可放置读数总和**（`effective_totals`，恒为正）归一化为概率分布。

### p=1 KR 距离的含义

对两棵归一化质量分布 μ、ν，p=1 Kantorovich–Rubinstein 距离等于把 μ 的
质量沿枝搬运成 ν 的最小总工作量（单位质量 × 枝长）。按割边流解释：

`distance = Σ_枝 ∫_枝 |该点远端子树内的累计净质量差| dx`

每条枝在其内部插入位置处被分段，分段内累计质量差为常数，逐段积分后
跨枝求和。响应里每对样本给出**全部枝号**的非负 `contribution`
（无流量的枝为 0），且每对的枝贡献之和等于该对距离。

### 响应格式

- `samples`：按输入顺序排列的样本 ID；
- `distance_matrix`：同序对称矩阵，对角为 0；
- `pairs`：按输入顺序的无序样本对（i<j），含 `distance` 与
  `edge_contributions`（每个元素为 `{edge_num, contribution}`，覆盖所有枝）；
- `excluded`：每样本被排除的 `{id, count, reason}`；
- `effective_totals`：每样本用于归一化的可放置读数总和。

请求样例见 `examples/compare_template.json`（把 `<PASTE_JOB_ID>` 换成
`/place` 返回的任务 ID；同一任务可被多样本复用）。

## 运行

```bash
PYTHONPATH=. .venv/bin/python -m uvicorn app.main:app --port 8169
```

## 自测

```bash
.venv/bin/python -m compileall -q app scripts
PYTHONPATH=. .venv/bin/python scripts/selftest.py
```

自测覆盖：示例放置、权重归一化、似然降序排序、jplace 下载结构、
非法字符与 `U` 定位整单拒绝、含逗号叶 ID 的引号导出往返、群落距离矩阵的
对称性/对角零/枝贡献非负且加总等于距离，以及未知任务、未知待测 ID、
非法计数、重复样本、零有效总数与参考树不一致等比较拒绝路径。

## curl 示例

```bash
# 组装请求（examples/ 内含 5 参考序列、3 待测序列、示例树）
.venv/bin/python - <<'EOF' > /tmp/payload.json
import json, pathlib
ex = pathlib.Path('examples')
print(json.dumps({
  "reference_fasta": ex.joinpath('reference.fasta').read_text(),
  "query_fasta": ex.joinpath('queries.fasta').read_text(),
  "newick": ex.joinpath('tree.nwk').read_text()}))
EOF

# 提交放置
curl -s -X POST http://127.0.0.1:8169/place \
  -H 'Content-Type: application/json' --data-binary @/tmp/payload.json

# 下载 jplace（job_id 取自上一步响应）
curl -OJ http://127.0.0.1:8169/place/<job_id>/jplace

# 群落比较（同一 job 表示这些地点共用同一次放置任务）
sed "s/<PASTE_JOB_ID>/<job_id>/g" examples/compare_template.json \
  | curl -s -X POST http://127.0.0.1:8169/compare \
      -H 'Content-Type: application/json' --data-binary @-
```

## jplace 输出

`tree` 字段为带 `{edge_num}` 注释的 Newick 树，枝号与结果中的
`edge_num` 一致；`fields` 为五个标准放置字段：
`edge_num, likelihood, like_weight_ratio, distal_length, pendant_length`。
同一查询的放置按似然降序排列，同值按枝号升序。
