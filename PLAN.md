# Plan: 遗留项推进 — GFN-FF 氢键梯度链深修 + DIC 重建成本削减 → v1.1.1

## 背景

- v1.1.0 期间发现存量 GFN-FF bug:阿司匹林 COOH 分子内氢键区 O/H 原子
  解析梯度与能量不一致(FD 逐原子扫描,误差 ~1.9 kcal/mol/Å)。笛卡尔
  优化路径恰好不停在这种点上故未暴露;DIC 安全网(笛卡尔重启)只是绕过。
- DIC(opt-in)迭代数 0.28–0.85× 但墙钟被重建时的 G 对角化吃掉
  (ibuprofen MMFF 72→758 ms)。漂移阈值 2.0 下仍有 ~10–15 次重建/运行。

## 任务

### 1 — GFN-FF 氢键/O-H 梯度链深修(正确性,优先)

1. **触发几何捕获**:`internal_opt.rs` 加 env 门导出(DIC_DUMP=路径,
   仅失败中止路径写 XYZ),复现不一致几何入库
   `tests/fixtures/gfnff/aspirin_hb_stall.xyz`。
2. **逐项隔离**:`Gfnff::energy_and_gradient` 增 `pub(crate)` 可选
   per-term 梯度输出参数(调用方不传零开销);诊断/测试用它对
   bond/angle/torsion/rep/es/disp/hb/xb/batm 逐项 FD 对拍,定位错误项。
3. **对拍 xtb 定分歧**:触发几何跑 `~/.local/xtb-gxtb/bin/xtb --gfnff`
   (能量+梯度),与 xtb 参考源(`gfnff_eg.f90`,已拉取)逐链比对,
   修正移植缺口(候选:egbond_hb 的 hb_cn 链、eg2/eg3 氢键项、
   gfnffdampa、bATM)。
4. **验证与锁值**:触发几何逐原子 FD 一致性回归测试(与其它原子同量级,
   ≤1e-4 相对);全量测试、既有 xtb 锁值、benchmark_mmff 230/230 不回归;
   aspirin GFN-FF 笛卡尔端到端复测(能量/迭代)。

### 2 — DIC 重建成本削减(性能,测量门控)

5. **漂移阈值实验矩阵**:DRIFT_REBUILD ∈ {2.0, 5.0, 1e9(仅拓扑应变)}
   × 3 fixtures × 2 引擎,记录迭代数/墙钟/重建数;采纳最优。
6. **H 侧翼二面角剪枝实验**(若 5 不足):保留重原子侧翼;完整性由
   现有 FD/roundtrip 测试 + 同极小能量门把关。
7. **决策门**:GFN-FC 侧任一 fixture 墙钟较笛卡尔改善 ≥1.2× 且 MMFF
   不劣于 1.5× → 重新评估默认切换;否则维持 opt-in,数据入
   CODE_STATUS。

### 3 — 发布

8. 版本 1.1.0 → 1.1.1(bugfix + DIC 调优,无 API 变化);
   CODE_STATUS 条目;README 测试计数如有变化同步;commit + tag
   v1.1.1(推送另请示)。

## 验收(实施后实测记录)

- `cargo test` 273/273(+1:aspirin hb-stall 逐原子 FD 一致性锁值)、clippy 0、fmt 干净、benchmark_mmff 输出与 v1.1.0 逐位一致(phosphirane/cyclobutene 既有离群点)
- 任务 1(GFN-FF 梯度链深修)完成,根因与修复:
  - 逐项隔离(TermGradients 新诊断通道 + CN 链按 bond/es/disp 拆分)定位到 bond 项;
  - 根因:egbond_hb 的 hb_cn 梯度链——xtb `dncoord_erf` 的 dtmp 是 tmp=0.5(1+erf) 导数的**一半**(上游怪癖),我们逐字移植并在 B 侧还顺了符号,合计 **−½×真值**;hb 高斯(kn=27.5)仅在 r≈rc 处活着——阿司匹林 COOH 分子内氢键恰好命中;此前误判“高斯已死”系诊断脚本 rcov 索引差一
  - 修复:d(hb_cn)/dr 全导数(×2)+ B 侧符号修正 → 触发几何 FD==解析逐位一致(全部 9 项、全原子);既有 xtb 锁值全不回归(所有验证几何处高斯指数级死亡,water-dimer hb 链贡献 ~5.8e-11)
  - 端到端效果:GFN-FC aspirin/ibuprofen 优化首次全部 conv=true(此前 maxf 停在 0.3 的“f64 地板”实为本 bug 的幻影力),极小更深(aspirin −2537.763、ibuprofen −3570.803),aspirin 34 ms;v1.1.0 的 energy_converged 机制保留作真噪声的安全网
- 任务 2(DIC 重建成本):漂移阈值矩阵 {2.0, 5.0, 1e9} 实测单调改善 → 采纳 1e9(不再漂移重建;映射对任意漂移精确线性,重建只为条件数)
  - 修后 DIC 迭代数 0.44–0.66×;墙钟:aspirin GFN-FC DIC ≈笛卡尔(0.75–1.3×,波动),ibuprofen GFN-FC 1.5–1.8× 慢、MMFF 1.7–1.8× 慢(单次变换开销 ~0.5–1 ms/迭代)
  - 决策门未全过(GFN-FC 需 ≥1.2× 改善且 MMFF ≤1.5× 劣化)→ **默认保持笛卡尔,internal 维持 opt-in**;后续方向:变换增量化(缓存 b_matrix 增量)/更廉价的重入
- 发布:1.1.0→1.1.1(纯 bugfix+调优,无 API 变化);wasm 重建(node 冒烟:1.1.1、cart/internal 双路径);README 测试计数 273
