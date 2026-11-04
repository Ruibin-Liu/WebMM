# 引擎侧能力参考(Engine Capabilities)

> WebMM v1.13.0 的 Rust/WASM 引擎面能力总览。数字全部取自实测档
> (examples/perf_probe.rs、examples/climb_molecules.rs、README 性能表、
> 各轮 PLAN/CODE_STATUS 记录);协议面见 docs/task-contract.md 与
> docs/worker-protocol.md,接口冻结政策见 docs/STABILITY.md。

## 0. 总则

- **确定性**:全引擎位等确定(bit-exact)——嵌入/搜索/优化/MD/NEB 全部由
  显式 seed 驱动,复读逐位一致(NEB 攀岩有位等复读测试守卫);worker 边界
  f64 走 17 位指数字符串(Rust str::parse 正确舍入,JS 只搬运)。
- **诚实原则**:不伪造能力——GFN-FF 拒绝 NEB(见 §9)、粗带攀岩回退报告
  带峰(climb_started=0)、相位 1 硬边诚实 conv=false、失败任务给出真实
  取消原因、本地确定性引擎不做 QSAR/ADMET 预测。
- **接口冻结**:v1.12.0 起接口面(wasm 导出/类方法集/worker 消息/任务
  KINDS)由 m22 golden 钉死;变更纪律见 STABILITY.md。

## 1. 分子输入与解析(`src/molecule`)

- SDF/MOL(V2000)解析;加氢(拓扑 + 3D 附着);Kabsch 对齐导出。
- SMILES 不经引擎——前端 RDKit-js 转 molblock 后进引擎。
- 金属键合体系 MMFF typing 拒绝(返回 NULL,与 RDKit 一致,绝不乱配型);
  此类体系走 GFN-FF。

## 2. 3D 嵌入(`src/etkdg`)

- ETKDG v3 距离几何:立体化学感知、seed 化;etkdg 组合规则同 RDKit。
- 验证:30 seeds × 6 分子,与 RDKit 全局最低能位等一致(6/6)。

## 3. 力场(`src/mmff`,`src/gfnff`)

| | MMFF94 / MMFF94s | GFN-FF |
|---|---|---|
| 覆盖 | 有机常规(H/C/N/O/F/S/P/Cl/Br/I + 电荷) | 广覆盖元素表(含金属配位) |
| 精度 | 230/230 vs RDKit < 0.01 kcal/mol;97 单点锁 0.00000 | 32/32 有机 < 4×10⁻⁶ Eh vs xtb;5 金属络合物(最好 Ni(CO)₄ 逐位,最差 Zn(NH₃)₄²⁺ 5.5×10⁻² Eh) |
| 梯度 | 解析梯度(230 分子 + GFN-FF fixtures FD 审计 0 失败,最差 2×10⁻⁵ rel) | 同左 |
| 用途 | 全管线默认(MMFF94s) | MMFF 拒配型体系的搜索/优化/MD;**不参与 NEB**(键感知依赖几何,见 §9) |
| 性能 | 乙醇 9 原子 1.7–3.6 µs/评估(核内零分配) | 淬灭缓存:精确坐标命中,建一次 535 µs |

## 4. 几何优化(`src/optimizer`)

- L-BFGS(默认)/BFGS;收敛判据可调(迭代/力/能量)。
- RestraintSet:冻结原子 + 距离/角度/二面角约束(受限扫描、柔性精修用)。
- 乙醇位移构型 300 步上限淬灭 ~115 µs;ibuprofen 全优化 7.7 ms(浏览器)。

## 5. 分子动力学(`src/md`)

- NVE(Velocity-Verlet)与 NVT(BAOAB Langevin,friction 可调);可分段
  步进(浏览器轨迹动画)。
- 隐式溶剂 GBSA(`src/solvation`):OBC2 Born 半径(HCT ψ + OBC2 α/β/γ)
  + LCPO SASA;默认气相,可选附加。
- MMFF NVT 步 ~27 µs(33 原子,浏览器,≈37k steps/s)。

## 6. 元动力学(`src/metad`)

- Well-tempered(γ 可调);CV:二面角(弧度,−π..π)与距离(Å)。
- 山丘偏置 4σ 截断——每步代价不随沉积山丘数增长。
- FES 重构(fes_s/fes_f 网格)、hill centers 全量导出。
- 底层引擎(well-tempered + 底层 MD)自 v1.x 稳定;iMTD 调度在其外层。

## 7. iMTD-GC 构象搜索(`src/metad/imtd.rs`)

CREST iMTD-GC 谱系:**多轮短 metad + 遗传交叉**,非单条长 metad。

- **调度器在引擎内、driver 在前端**(`app/tasks/imtd.js` scatter-gather
  并行,worker farm 属前端);单轮 = runs_per_round 条短 metad(偏置强度
  与 CV 轮换)+ crossing_pairs 次遗传交叉。
- **CV 池**:可旋转二面角(自动识别)+ 距离补足;逐 run 轮换。
- **两层去重**:坐标预筛 + 淬灭后能量/结构判层;energy_window / rmsd_tol /
  max_conformers 三旋钮;ensemble 上限 250。
- **确定性**:seed 派生 seed+round*1000+run;浏览器并行结果与串行逐位一致
  (m15 守卫)。
- **verify 相**:搜索后 unbiased MD 梯子(可多轮)产出转迁计数/占据——
  景观的无偏统计源。
- 搜索成本:乙醇 2×3×800 ≈ 17 ms;ibuprofen 2×3×800 ≈ 220–300 ms(原生)。

## 8. 能量景观(`ImtdSearch::landscape` → 前端 SVG)

- **梯子**(metad + MD)折叠统计:帧数、访问、转迁对(pairs + count)。
- **Disconnectivity 森林**:Kruskal 最小生成树;阈值用精修鞍点(可得时)
  而非带峰——树随精修更新。
- **盆地网络图**:确定性力导向(220 迭代,布局按盆地数|边表|视口|梯源缓存,
  选中点击只重绘);边宽 ∝ 转迁计数,节点色 ∝ 能量,半径 ∝ 占据。
- **盆地表**:E、G(QH)、Boltzmann.E、Boltzmann.G、简并度 g、占据;
  **转迁表**:a/b/count。
- 刚性体系(单盆地)正确产出空景观(aspirin/naphthalene 验证)。

## 9. NEB 与鞍点精修(`src/neb.rs`)

- CI-NEB:端点 Kabsch 对齐、improved tangents(Henkelman 2000)、端点
  带内重淬灭;**默认攀岩开启**(v1.13.0)。
- **积分器**:默认 FIRE(全带全局自适应,Bitzek 2006),QuickMin 可选。
  丁烷边 108.8 ms/1937 迭代(QuickMin)→ 41.5 ms/803+24 迭代(FIRE+攀岩);
  FIRE 精度经受限扭转扫描独立仲裁(势垒 3.92 vs 参考 3.95,QuickMin
  下垂 0.6)。
- **收敛像冻结**:力 < tol/2 连续 40 迭代跳过评估与积分;**峰值像(鞍点
  载体)永不冻结**;诊断字段 evals_skipped/n_frozen。
- **攀岩(钉住 CI)**:相位 1 全收敛后激活;攀岩像一次选定永不重指(旧
  极限环根因是角色跳换,非积分器);严格内极大守卫(峰不与端点相邻且
  高于双端 0.1 kcal);**相位 1 快照回退**——攀岩结果低于带峰 0.1 以上时
  诚实回退带峰(climb_started=0);能量失速收敛(脊两侧 1e-4 kcal 定位)。
- 粗带(≤10 像)攀岩实测修复离散化下垂:8 像 -1.325→-1.155、10 像
  -1.262→-1.132(扫描参考 -1.09)。
- 大分子验证(ibuprofen/threonine 全链路):攀岩激活边全部收敛或诚实
  回退;预存相位 1 硬边(复合扭转 20000 迭代)与攀岩无关。
- **GFN-FF 拒绝 NEB**:键感知依赖几何,端点对齐后感知失真——报错
  "NEB refinement supports MMFF engines only"。
- 成本:丁烷边 41.5 ms;大分子边 ~1.2–2 s(原生串行;浏览器 worker 按边
  并行,墙钟取最慢边)。

## 10. 振动与热力学(`src/vib.rs` + imtd thermo)

- 有限差分 Hessian → 符号频率;**准谐(QH)地板 100 cm⁻¹**
  (QH_CUTOFF_CM);harmonic_correction 输出 (ZPE+热校正, S, 模式计数);
  单位常数由甲烷振动带锁定。
- **Boltzmann 居住**:`boltzmann_populations(_g)` 按简并度加权。
- **简并度 = 观察到的井数**:变体几何计数统一覆盖镜像孪生/置换旋转/
  C3 转子(g ≤ 对称数);乙醇 σ=3、丁烷 3–5 井有测试锚定。
- 盆地表中 Boltzmann.E / Boltzmann.G 与占据由上述链条导出。
- 性能:乙醇 frequencies ~290–390 µs。

## 11. 3D 形状(`src/shape`)

- Grant–Pickup 高斯重叠(shape-it 常数系);shape Tanimoto =
  V_AB/(V_A+V_B−V_AB);多起点 + L-BFGS 刚体对齐;色版(color)对齐可选。
- 验证:自对齐 1.000000;与 RDKit 网格法相关 ρ=0.96。
- ECFP4 指纹 + Tanimoto、Ertl SA score(RDKit sascorer 位等移植)在
  同一 wasm 面(shape_align_wasm / ecfp4_fingerprint_wasm / sa_score_wasm …)。

## 12. 性能速查(量尺:examples/perf_probe.rs,M 系 release)

| 路径 | 数字 |
|---|---|
| MMFF energy_and_gradient(乙醇 9 原子) | 1.7–3.6 µs(核内零分配) |
| L-BFGS 淬灭(300 步上限) | ~110–120 µs |
| iMTD 搜索(乙醇 2×3×800) | 17–22 ms |
| iMTD 搜索(ibuprofen 2×3×800) | 220–300 ms(原生) |
| NEB 单边(丁烷,FIRE+攀岩默认) | 41.5 ms / 827 迭代 |
| NEB 单边(大分子,ibuprofen 级) | ~1.2–2 s(原生串行) |
| vib::frequencies(乙醇) | ~290–390 µs |
| MD 步(MMFF NVT,33 原子) | ~27 µs(浏览器) |
| ETKDG 嵌入(ibuprofen) | ~11 ms(浏览器) |

## 13. 限制与诚实降级

- 不做:蛋白质、周期体系、GFN-FF 之外的 QM、显式溶剂、QSAR/ADMET。
- GFN-FF:不进 NEB(拒配);搜索/优化/MD 可用。
- NEB 攀岩:粗带(≤10 像)或复合扭转硬边可能诚实回退/不收敛——输出
  始终带 converged/climb_started 语义,消费方不猜。
- dimer 单端鞍点搜索:**关闭**(无单端场景;若立动力学/速率功能再评)。
- m1_3d_pipeline 依赖外部 xtb 二进制(本机环境项,非引擎缺陷)。

## 14. 验证状态汇总

| 面 | 守卫 |
|---|---|
| 引擎单测 | cargo test 372+7(含 6 NEB:势垒带/端点固定/扫描参考收敛/粗带回退/退化拒绝/确定性) |
| 属性测试 | 3 性质(proptest,seed 固定) |
| 浏览器 CDP | 36 套件 + platform 153(m22 golden 24 项钉接口面) |
| RDKit/xtb 对拍 | MMFF 230/230 门禁(scripts/benchmark_mmff.py)、GFN-FF 32/32 |

## 附:文档地图

- 接口冻结与版本纪律 → docs/STABILITY.md
- 任务信封协议 → docs/task-contract.md;执行器设计 → docs/task-executor-design.md
- worker 消息面 → docs/worker-protocol.md
- 原子类型覆盖 → docs/atom-type-coverage.md;GFN-FF 移植 → docs/gfnff-porting-notes.md
- 能量验证分析 → docs/validation-energy-analysis.md
- 开发轮次日志(排雷存档)→ CODE_STATUS.md + docs/CODE_STATUS-archive.md + 各轮 PLAN.md
