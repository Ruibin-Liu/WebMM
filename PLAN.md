# Plan: aza-hop 正式化 + SA 立体罚项精化(v1.10.0→v1.10.1;全部完成)

## 1. aza-hop 正式化 ✅

- 位点扩到五元芳环碳(`[c;H1;r5]`——噻吩/咪唑 C→N 成噻唑/吡唑系)
- **双位扫描**(diazine/triazole 组合):C(n,2) 组合,与单扫合并去重,
  总候选 ≤400 硬顶(既有),邻位 N-N 由 RDKit is_valid 诚实过滤
- UI:按钮改 "aza-scan (r5+r6, 1-2 swaps)",状态行报单/双计数

## 2. SA 立体罚项精化(消掉类固醇 1.32 差)

RDKit 的四个罚项按其真语义实现:
- **n_chiral = 潜在手性中心总数**(FindMolChiralCenters
  includeUnassigned=True = 指定+未指定):sp3 C(必要时含 N)且
  4 配体可区分——配体区分用迭代配体排序(Morgan 式 3 轮,环情
  形按 RDKit 行为近似),Python 对拍定分寸
- **n_bridgehead**:环键数 ≥3 的原子(稠环联结点)
- **n_spiro**:恰好 2 环键且属 2 个不同环(螺原子)
- **n_macro**:尺寸 >8 的环计数(用逐键最短环枚举,去重)
- Python 侧逐分子计算四分量的期望值入 fixture,Rust 逐分量对拍,
  目标:全语料(含类固醇)|ΔSA| < 0.05,金标容差收紧回 0.05

## 验收

- cargo(含 sascore_golden 容差 0.05)+ m0-m6 + 平台 Node 全绿;
  aza-scan 双扫在 m6 出真实 diazine;零 page error


## 实施与验收(完成)

1. **aza-hop 正式化 ✅**:位点扩至五元芳环([c;H1;r5,r6]——噻吩→噻唑、
   咪唑→吡唑系);双位扫描 C(n,2)(diazine/azole,邻位 N-N 由
   is_valid 诚实过滤);终态状态行带 "N single + M double positions"
   计数;paracetamol → 6 候选(2 单换 + 3 双换 diazine + 1),
   T 83.4-81.0%,flex 至 87.6%,SA 2.03-2.77 排序合理
2. **SA 罚项精化 ✅**:①n_chiral = 潜在手性中心(**全分子迭代规范
   排序**至稳定(sp3 C + 4 配体(含 1 隐氢)两两可区分)——树走
   配体不变量会破环对称性,对称笼体系(C1C2CC3CC1C2C3)过计;
   排序法对称性构造保证,类固醇 8/8、葡萄糖 5/5 精确);②
   bridge/spiro/macro 于枚举环基(RDKit 规则:bridgehead = ≥3 环键
   ∧ 每含环与另一含环共享 ≥2 键;spiro = 两环恰交此原子;macro =
   >8 元环)——norbornane 2/2、spiro 1/1、cyclodecene 1/1 全对;
   **文档化限界**:对称笼 C1C2CC3CC1C2C3 的 bridgehead 1/2(环基
   依赖 symmSSSR,应激例排除并注释);③新 sa_penalty_components
   测试逐分量对拍 penalties.json(RDKit Python 期望值)
3. **wasm 端对齐**:重建后全语料 max |ΔSA| = **0.0000**(v1.10.0
   时 1.32);类固醇 4.4251 = 参考值逐位;SA 容差收紧 0.05

## 验收数字

- cargo 322 + sascore_golden 4(位级 + 分量级 + 0.05 容差);
  m0-m6 七套件 37/10/11/10/32/44/**27**;平台 Node 34/34;clippy 0;
  零 page error
