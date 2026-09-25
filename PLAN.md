# Plan: ETKDG 4D 阶段换真 L-BFGS — 替换固定步长梯度下降 → v1.2.3

## 背景

minimize_4d_first(400 迭代)与 minimize_4d_collapse(200 迭代)的
默认路径使用**固定步长梯度下降**(step = 0.1/max_g)——收敛效率低
(审计:4D 需要 ~60-90 次迭代才停滞)。而 rdkit_all() 分支已用
`lbfgs_minimize`(共享 L-BFGS);3D 的 minimize_etkdg 也是 L-BFGS
(1-46/300 收敛)。**把默认路径的 4D 从固定步长换成同一个
lbfgs_minimize**,预计迭代数从 ~60-90 降到 ~10-30(与 3D 的收敛
率一致),每次迭代的 O(n²) 梯度/能量求值次数同比下降。

## 任务

1. **minimize_4d_first 默认路径重写**:删除手动固定步长循环 +
   best_coords 追踪,改为与 rdkit_all 分支相同的
   `lbfgs_minimize(&mut x, n, 4, &energy_at, &gradient_at, max_iter,
   force_tol)` 调用(用 energy_4d/gradient_4d + 阶段权重)。
   两个分支合并为一个(不再需要 rdkit_all 分支区分)。
2. **minimize_4d_collapse 同理**:合并为单一 lbfgs_minimize 调用
   (权重 FOURTH_MIN_WEIGHT_CHIRAL / FOURTH_MIN_WEIGHT_FOURTH)。
3. **保留 stagnation 中断**:lbfgs_minimize 已有 force_tol 收敛;
   外层 max_iter 不变(400/200),但 L-BFGS 会提前退出。
4. **嵌入质量裁决**:ensemble_stats_vs_rdkit 6/6 门禁(min/median
   不回归)+ ETKDG 既有种子回归测试;能量序列**会变**(收敛到
   不同精度的 4D 起点)——不是 bug,是更收敛的起点,下游 3D
   精修兜底。
5. **发布**:1.2.2→1.2.3;CODE_STATUS/PLAN;wasm;commit+tag。

## 验收(实施后实测记录)

- `cargo test` 275/275 全绿;ensemble 6/6;benchmark 230/230 逐位一致;
  clippy 0;fmt;wasm(node 冒烟 1.2.3);API 零变化
- **诚实速度结论:4D 换 L-BFGS 在药物分子尺度无可测墙钟收益**
  (交错 A/B:55-57 vs 53-61 ms/embed,噪声级)。原因:L-BFGS 每迭代
  收敛更快但线搜索需要更多函数评估,总 O(n²) 求值次数相当。
- **结构收益**(非速度):统一优化器(所有阶段同一个 lbfgs_minimize,
  消除固定步长分支与 rdkit_all 双路径);4D 最小更收敛(更高质量的
  3D 起点);代码减少 ~120 行。
- 嵌入速度的诚实总结(五轮迭代后):**55-60 ms/embed @33 原子是
  当前架构的本质成本**,由分散的 O(n²) 梯度/能量/检查项构成,
  无单一 >1.5× 杠杆;进一步提速需 SIMD(受逐位对拍约束)或算法
  层重设计。
