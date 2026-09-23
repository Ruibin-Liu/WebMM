# Plan: 单分子优化第二轮 — GFN-FF erf 热路径 + 能量分辨率提前停 + 内坐标(DIC)优化器 → v1.1.0

## 背景(v1.0.1 后实测)

- GFN-FF E+G 剩余热点(ibuprofen 采样):erf 家族 ~45–50%
  (EEQ 矩阵 erf(gam·r)/r 18.4%、ES 对项 15%、erf_cn 每调用跑两遍 ~22%);
  现实现是 Taylor(≤40 项)/Lentz(≤80 迭代)变长循环。
- GFN-FF 大分子优化尾部在 f64 能量分辨率(~1e-9 kcal/mol)地板上磨
  α 下限(每迭代 ~15 次 energy-only 试探),占墙钟 ~20–30%,报
  conv=false 但能量已达机器精度。
- WASM/native 差距实测仅 1.0–1.4×,无需 WASM 专项。
- 笛卡尔 L-BFGS 迭代数(布洛芬 MMFF 298 / GFN-FF 225+)是下一台阶:
  内坐标优化(DIC,Baker/geomeTRIC 路线)在柔性分子上典型 2–5×。

## 任务

### Phase 1 — GFN-FF E+G 热路径(`src/gfnff/mod.rs`)

1. **erf 快速实现**:fdlibm 风格分段(|x|≤0.84375:x·P(x²);
   ≤6:1−exp(−x²)·P(x²);>6:±1),多项式系数用 numpy 对 math.erf
   做 Chebyshev/minimax 拟合并离线生成(常量表进源码),目标最大绝对
   误差 ≤1e-16(实测锁值:现有 erf 断言 1e-12/1e-14、xtb 逐项锁
   1e-6..1e-8 Eh,余量充足)。新 erf 精度由单元测试锁死(对
   math.erf 值表 + 存量断言保留)。
2. **erf_cn 双跑合并**:E+G 里 `erf_cn`(logCN,EEQ)与
   `erf_cn_raw`(CN 链)合并为一次 O(n²) 遍历同时产出两份(新内部
   fn `erf_cn_both`;`erf_cn`/`erf_cn_raw` 保留供 setup/energy 路径)。
   累加顺序不变 → 逐位一致。
3. **powf(1.5) → r·sqrt(r)**(排斥两项循环;1 ulp 级扰动,锁值容差
   内)。

### Phase 2 — 优化器能量分辨率地板提前停(`src/optimizer/mod.rs`)

4. 核心结果 `OptimizationResult` 增字段 `energy_converged: bool`(仅
   构造点 1 处,非破坏)。停止判据:累计线搜索失败 ≥3 且最近 10 个
   接受步 |ΔE| 全部 <1e-9 → 停止并置 energy_converged=true(converged
   语义不变,仍为力判据)。`optimize_dispatch` 映射:energy_converged
   → WASM `converged=true` + message 说明"能量分辨率地板(f64 极限),
   max_force=…"(app 成功路径显示 converged,message 已有链路)。

### Phase 3 — 内坐标(DIC)优化器

5. **`src/optimizer/jacobi.rs`**:对称矩阵 Jacobi 特征解(循环
   扫描,≤400² 规模,无新依赖);测试:对角阵/已知矩阵 vs 解析特征对。
6. **`src/optimizer/internals.rs`**:
   - 连通性:内置 rcov 表,d < 1.3(rcov_i+rcov_j);
   - 原始内坐标:键伸缩(全部)、键角(共顶点键对,平角 >175° 改用
     1-3 距离坐标)、二面角(中央角非线)、三配位 improper、非键
     距离坐标(重原子对 d<2.8Å 且非 1-2/1-3);
   - B 矩阵(解析一阶导)+ FD 验证测试;
   - G = B·Bᵀ 对角化(Jacobi),特征值截断 1e-6(相对最大值)取
     delocalized 基 U(dof ≤ 3N−6);
   - 回变换:Baker 迭代 dx = B_qᵀ(dq − B_q·dx) ×3–5 轮;
   - 梯度变换:TR(平移+旋转)投影后 g_q = G⁺B·g_x。
7. **Objective 抽象重构**(`src/optimizer/mod.rs`):抽出
   `trait Objective { fn f_and_g(&self, &mut [f64], &mut [f64]) }`,
   L-BFGS 核循环只写一份;笛卡尔路径封装为恒等 Objective(循环体
   逐字保留 → 结果逐位不变,由既有测试锁)。
8. **DIC 优化器**:q 空间复用 L-BFGS 核 + 现有 Armijo 线搜索
   (energy-only 试探经回变换);重建触发(累计 |dq| 超限或线搜索
   失败 2 次后重建 internals 并清空历史);退化回退(N<3、dof 异常、
   构建失败 → 笛卡尔,不 panic);收敛判据仍在笛卡尔力空间。
9. **接线与默认切换(测量门控)**:`OptimizationOptions` 增
   `coordinates: String`("internal" 默认 / "cartesian" 逃生门,
   wasm-bindgen 附加字段,additive;pkg/webmm.d.ts 随构建更新)。
   门禁:ethanol/aspirin/ibuprofen × MMFF/GFN-FF,DIC 迭代数 ≤
   笛卡尔 0.7×(至少柔性分子)且极小能量差 ≤0.05 kcal/mol(同盆地)
   → 默认 internal;未达标则默认 cartesian、仅留选项并记录。

### Phase 4 — 发布

10. 版本 1.0.1 → 1.1.0(Cargo.toml/Cargo.lock/package.json);
    README 优化器章节 + 测试计数;wasm 构建冒烟(node 基准);
    CODE_STATUS.md 条目;commit + tag v1.1.0(推送另请示)。

## 验收(实施后实测记录)

- `cargo test` 272/272(+12:Jacobi ×3、internals ×6、DIC 同极小/退化回退、WASM internal 端到端、erf 密度网格)、`cargo clippy --all-targets` 0 警告、`cargo fmt` 干净
- `benchmark_mmff.py --no-speed`:230 分子输出与 v1.0.1 逐位一致
- Phase 1(同热状态 A/B):GFN-FF E+G 1.33–1.76×(ethanol 40→25、aspirin 244→183、ibuprofen 630→358 μs;绝对值门 ≤230 μs 在同台不同热状态下不可比,以 A/B 比率为准并记录);
erf 重构过程中发现并绕开了 x² 变量下 e^t·erfc(√t) 的 √t 分支点(Bernstein ρ≈1.33 收敛墙,实测误差按 1.33⁻ᴺ 缩小)——改用 x 变量三段 Chebyshev
- Phase 2:energy_converged 字段 + WASM message 链路落地;aspirin/ibuprofen 实测主导停止路径仍是 5 连败中止(erf 重构后轨迹变化,磨步尾部变小),检测器(零下降地板接受/微小步×3/平稳窗口)在其它模式/更大分子上生效;GFN-FF 端到端同热 A/B:aspirin 173→115 ms(1.5×)、ibuprofen 496→136 ms(3.6×)
- Phase 3(DIC):**正确性达标、性能门未达标 → 按计划保持笛卡尔默认,internal 为选项**
  - 正确性:MMFF 三 fixture 极小与笛卡尔逐位一致(18.9098/24.7323/-1.3369);回归测试锁定;WASM 端到端冒烟通过(internal 41it vs cart 47it 同能量)
  - 迭代数:MMFF 0.6–0.85×、GFN-FF 0.28–0.55×(aspirin GFN-FF 302→44)✓
  - 墙钟:重建时 G 对角化(Jacobi O(n_prim³·sweeps))吃掉送代收益——MMFF aspirin 26→77 ms、ibuprofen 72→758 ms ✗(门禁未过:需迭代数 ≤0.7× **且**墙钟改善);缓解已做(扁平 Jacobi、B_q 构建缓存、漂移阈值 2.0、H-H 二面角剪枝)仍不足
  - 实施中发现**存量 GFN-FF bug**:阿司匹林 COOH 分子内氢键区 O/H 原子解析梯度与能量不一致(FD 误差 ~1.9 kcal/mol/Å,HB 项梯度链);笛卡尔路径从未停在这种点上故未暴露。已加安全网(DIC 失败中止→笛卡尔重启,结果不差于纯笛卡尔);深修 HB 链另立项
  - 关键实现修正(记录):①映射必须冻结在构建参考(x(z)=x_ref+P·z 精确线性,P=B_qᵀΛ⁻¹),梯度必须用同一份 B_q(曾用当前几何的 B 导致线搜索在漂移后死亡);②增量链式映射会破坏 L-BFGS 的 E(z) 一致性(错误极小/拓扑崩坏),废弃
- Phase 4:版本 1.1.0(Cargo/lock/package.json);README 优化器章节+目录树+测试计数;wasm 重建(node 冒烟:版本 1.1.0、cart/internal 双路径);d.ts 携带 coordinates 字段
