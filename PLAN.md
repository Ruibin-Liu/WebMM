# Plan: RGD(R-基团分解)+ 骨架频率分析(纯 app 端,反应 SMARTS 路线)

## 探测结论(已实测)

- vendored `get_rgd` 绑定恒 null(疑似未实现)——弃用
- **反应路线全通**:`get_rxn('[*:1]c1ccccc1[*:2]>>[*:1].[*:2]')` +
  `run_reactants(MolList)` → 产品集(多匹配=多集)→ 每集 MolList 的片段
  `get_smiles()`(苯乙酮例:R1=CCO、R2=O=CO 正确;产品序=反应式标签序,
  确定性);MCS 可用(`get_mcs_as_json(MolList)` → SMARTS,备后续)
- RGD 核心 v1 需显式 `[*:n]` 标记(库内无 mol 编辑 API,自动标记不可行,
  诚实文档化;DataWarrior 同惯例)

## 设计(Search 标签页下两面板,复用库基建)

### T1 R-基团分解面板

- 输入:核心 SMILES(占位示例 `[*:1]c1ccccc1[*:2]`,标签数 k 由正则
  解析)+「Example」按钮一键填 aspirin-苯环双取代示例
- 分解:每库条目 → MolList → run_reactants → 取**第一产品集**(多匹配
  只取首集,文档化;RDKit 完整 RGD 会跨匹配打分)
- 结果表:Name | SMILES | R1..Rk(片段 canonical SMILES)+ 匹配数状态行
  (matched/total,不含核心者不计)+ CSV 导出
- Python 奇偶:`scripts/gen_rgd_refs.py` 用 rdChemReactions 同反应
  SMARTS,产物片段剥除 `[n*]` 虚原子后 canonical 对拍(冻结参考 +
  Playwright 全量断言)

### T2 骨架频率面板

- 「Analyze scaffolds」:对库逐条 murckoScaffold(现有 JS 函数)→
  频率表:骨架 SMILES | 条数 | 占比(降序);CSV 导出
- 行点击 → 回载该骨架(单分子视图,复用 loadSearchRow 式流程)

### T3 验收

1. RGD 奇偶:参考集(核心=对位双取代苯;库=demo 55)片段逐一相等
2. 交互:Example 一键、分解表 R1..Rk 正确、无核心条目被排除且计数
   正确、CSV 列齐;骨架表计数总和=库数、行点击回载
3. 390px、零 page error、三模式检索回归;门禁 sanity(293 测试,纯 JS)
