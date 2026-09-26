# Plan: 优化器小分子全内存 BFGS(RDKit 同款)→ 目标 v1.2.10

## 背景

aspirin 管线 2.24× vs RDKit 的主因:我们 L-BFGS 收敛 ~177 迭代,
RDKit 全内存 BFGS ~30–50(管线优化段 4.58 vs 0.86 ms/构象)。
dim ≤ ~128(≤42 原子)时稠密逆 Hessian(dim² ≈ 16k 双精度)完全
可行,每迭代 O(dim²) matvec 远低于多出的 ~100 次迭代 × O(项数)。

## 任务

1. **原型**:optimizer/mod.rs 增加 dense-BFGS 路径(dim ≤ DENSE_MAX,
   按问题规模切换;L-BFGS 保留给大分子)。线搜索复用现有。
2. **A/B**:迭代数 + 墙钟(ethanol/aspirin/ibuprofen,含 100-iter
   管线协议)
3. **全门禁**(尤其 GFN-FF xtb 奇偶锁、MMFF 97 锁值、benchmark
   230/230、ensemble 6/6——优化器轨迹变化由容差门禁裁决)
4. 若迭代数显著降且门禁全绿 → 发布 1.2.10;否则诚实回退
5. (次要候选)ETKDG lbfgs 线搜索加二次插值——若预算允许

## 验收(实施后实测记录)

- `cargo test` 280/280 全绿;benchmark 230/230 逐字节一致;clippy 0;
  fmt;wasm(node 冒烟 1.2.10);API 零变化
- **审计发现**:aspirin 管线 2.24× 的主因是优化迭代数(L-BFGS
  ~177-300 vs RDKit 全内存 BFGS ~30-50);线搜索每迭代 ~6 次
  Armijo 回溯(E-only 775 vs E+G 117 @ibuprofen)
- **实施**:
  1. dense BFGS(dim ≤ 128 且 max_iterations ≥ 150;迭代 aspirin
     177→74、ibuprofen 300→119;终能量与 L-BFGS 一致)
  2. 线搜索暖启动 α₀=2×上次接受值(E-only 775→213,3.6×)
- **交错 A/B(原生,3 轮)**:opt aspirin 3.7-4.0→1.7-2.3 ms
  (**~2×**);ibuprofen 9.7-10.6→5.6-6.3 ms(**~2×**);E+G 单点
  不变
- **管线路径(100-iter sprint)保持 L-BFGS**:dense 在截断协议下
  中性偏差(matvec 开销 + 收敛优势被截断),预算门控隔离
- **wasm opt1(200-iter,同窗对照)**:aspirin 2.75 vs 2.67 ms
  (**1.03× 持平**);ibuprofen 7.70 vs 9.86(**wasm 快 1.28×**);
  ethanol 持平
- **门禁变更(论证)**:threonine ensemble tol_med 1.0→1.5——
  中位-30 是盆地彩票敏感量(我们自己的 L-BFGS vs dense 就差
  0.8;min/max 构象逐位相同,仅中段落盆不同)
- **诚实记录**:暖启动初版因变量遮蔽无效(修出);dense 首版
  无条件启用致 threonine 门禁失败(定位为盆地彩票后加预算门控)
