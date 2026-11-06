# 引擎侧能力参考(Engine Capabilities)

> WebMM v1.15.0 的 Rust/WASM 引擎面能力总览。所有性能与验证数字均取自实测基准(`examples/perf_probe.rs`、`examples/climb_molecules.rs`、`README` 性能表、各轮 `PLAN` / `CODE_STATUS` 记录)。协议规范见 `docs/task-contract.md` 与 `docs/worker-protocol.md`,接口冻结政策见 `docs/STABILITY.md`。

---

## 0. 总则

- **确定性(Determinism)**:全引擎实现位等确定(bit-exact)。嵌入、构象搜索、几何优化、分子动力学(MD)与微动弹性带(NEB)全部由显式随机数种子(`seed`)驱动,多次复读结果逐位一致(NEB 攀岩设有位等复读测试守卫)。Worker 边界的 `f64` 浮点数统一采用 17 位指数字符串传输(由 Rust `str::parse` 执行正确舍入,JavaScript 仅作透传搬运)。
- **诚实原则**:绝不伪造计算能力。GFN-FF 明确拒绝 NEB(见 §9);粗带攀岩触发回退时如实报告带峰(`climb_started=0`);相位 1 硬边未收敛时诚实返回 `conv=false`;失败任务给出真实取消原因;作为本地确定性物理引擎,不做 QSAR/ADMET 经验预测。
- **接口冻结**:自 v1.12.0 起,引擎接口面(WASM 导出函数、类方法集、Worker 消息协议、任务 `KINDS`)由 `m22` golden 测试钉死;变更纪律详见 `docs/STABILITY.md`。

---

## 1. 分子输入与解析(`src/molecule`)

- **格式与预处理**:支持 SDF/MOL(V2000)解析、加氢(拓扑构建 + 3D 坐标附着)以及 Kabsch 刚体对齐导出。
- **SMILES 边界**:SMILES 字符串不直接经过引擎,由前端 `RDKit-js` 转换为 `molblock` 后再输入引擎。
- **金属体系路由**:含金属键合体系在 MMFF 原子配型(typing)阶段会被明确拒绝(返回 `NULL`,行为与 RDKit 一致,绝不强行乱配型);此类体系路由至 GFN-FF 处理。

---

## 2. 3D 嵌入(`src/etkdg`)

- **算法实现**:采用 ETKDG v3(Experimental-Torsion Knowledge Distance Geometry)距离几何算法,具备立体化学感知与显式 `seed` 控制;ETKDG 组合规则与 RDKit 保持一致。
- **基准验证**:在 30 seeds × 6 分子的测试集中,与 RDKit 找到的全局最低能构象实现位等一致(6/6)。

---

## 3. 力场(`src/mmff`, `src/gfnff`)

| 维度 | MMFF94 / MMFF94s | GFN-FF |
| :--- | :--- | :--- |
| **元素覆盖** | 常规有机元素(H/C/N/O/F/S/P/Cl/Br/I + 形式电荷) | 广覆盖元素周期表(含过渡金属配位体系) |
| **能量精度** | 230/230 分子对比 RDKit 误差 < 0.01 kcal/mol;97 个单点锁定 0.00000 | 32/32 有机分子对比 `xtb` 误差 < 4×10⁻⁶ Eh;5 个金属络合物(最优 Ni(CO)₄ 逐位一致,最差 Zn(NH₃)₄²⁺ 误差 5.5×10⁻² Eh) |
| **梯度实现** | 解析梯度(230 个 MMFF 分子 + GFN-FF 测试集经有限差分 FD 审计 0 失败,最大相对误差 2×10⁻⁵ rel) | 同左(解析梯度,同享 FD 审计标准) |
| **管线定位** | 全管线默认力场(MMFF94s) | 用于 MMFF 拒配型体系的构象搜索、优化与 MD;**不参与 NEB**(因键感知依赖几何坐标,见 §9) |
| **评估性能** | 乙醇(9 原子)1.7–3.6 µs/次评估(核内零内存分配) | 淬灭缓存机制:精确坐标命中复用,单次拓扑构建耗时 535 µs |

---

## 4. 几何优化(`src/optimizer`)

- **优化器**:支持 L-BFGS(默认)与 BFGS,收敛判据(最大迭代步数、力阈值、能量阈值)均可调。
- **约束系统(`RestraintSet`)**:支持冻结特定原子,以及施加距离、键角、二面角的平底简谐(flat-bottom harmonic)约束(服务于受限扭转扫描与柔性精修)。
- **性能**:乙醇位移构型在 300 步上限内完成淬灭耗时 ~115 µs;布洛芬(ibuprofen)全结构优化在浏览器端耗时 7.7 ms。

---

## 5. 分子动力学(`src/md`)

- **系综与积分器**:支持 NVE(Velocity-Verlet 积分器)与 NVT(BAOAB Langevin 恒温器,摩擦系数 `friction` 可调);支持分段步进执行(适配浏览器端实时轨迹动画渲染)。
- **隐式溶剂模型(`src/solvation`)**:实现 GBSA 模型,包含 OBC2 Born 半径计算(HCT ψ 积分 + OBC2 α/β/γ 参数)与 LCPO 溶剂可及表面积(SASA)非极性项;默认采用气相,可按需附加。
- **性能**:MMFF NVT 单步耗时 ~27 µs(33 原子体系,浏览器环境,折合 ≈37k steps/s)。

---

## 6. 元动力学(`src/metad`)

- **算法与集体变量(CV)**:实现 well-tempered 元动力学(偏置因子 γ 可调);集体变量(Collective Variables, CV)支持二面角(弧度制,−π..π)与原子间距离(Å)。
- **计算复杂度控制**:高斯山丘(hills)偏置势采用 4σ 截断策略,确保单步评估代价不随已沉积山丘总数的增长而膨胀。
- **数据导出与架构**:支持自由能面(Free Energy Surface, FES)网格重构(`fes_s` / `fes_f`)及山丘中心(hill centers)全量导出。底层引擎(well-tempered + 基础 MD)自 v1.x 起保持稳定,iMTD 调度层构建于其外围(见 §7)。

---

## 7. iMTD-GC 构象搜索(`src/metad/imtd.rs`)

承袭 CREST iMTD-GC 算法谱系:**采用「多轮短时元动力学 + 遗传交叉(Genetic Crossing, GC)」策略**,而非依赖单条长时元动力学轨迹。

- **混合架构与单轮结构**:调度器(scheduler)内置于 Rust 引擎,驱动器(driver)位于前端(`app/tasks/imtd.js` 采用 scatter-gather 模式并行调度前端 Worker farm)。单轮搜索包含 `runs_per_round` 条短时元动力学轨迹(交替轮换偏置强度与 CV 组合)以及 `crossing_pairs` 次结构遗传交叉。
- **CV 池构建**:自动识别可旋转二面角,不足时以原子对距离补足,并在每个独立 run 间轮换分配。
- **两层去重机制**:由「坐标预筛」与「淬灭后能量/结构精细判别」组成;提供 `energy_window`、`rmsd_tol`、`max_conformers` 三个控制旋钮,系综(ensemble)容量上限为 250。
- **并行确定性**:子任务种子按 `seed + round*1000 + run` 公式派生;浏览器端多线程并行结果与单线程串行执行结果逐位一致(受 `m15` 测试守卫保护)。
- **验证阶段(verify 相)**:构象搜索结束后,运行无偏(unbiased)MD 梯子(ladders,可多轮),直接产出盆地间的转迁计数与热力学占据统计——作为能量景观的无偏统计数据源(见 §8)。
- **性能**:乙醇(2×3×800 规格)耗时 ≈ 17 ms;ibuprofen(2×3×800 规格,原生环境)耗时 ≈ 220–300 ms。
- **构象子集契约(WebMM↔WebQM)**:`conformer_subset_json` 提供三档确定性采样——ETKDG 快速档 `lowest_per_cluster`(默认)与 `energy_stratified`,以及 CREST iMTD-GC 谱系档 `imtd_lowest_per_cluster`(固定 2×3×800 小规模预算)。统一输出能量升序的 `energy_kcal_mol` / `coords_angstrom` / `cluster_id`,provenance 记录 `webmm_version` / `seed` / `energy_engine` / `strategy`;`n_target`、能量窗和 RMSD 容差均设硬上限,超限显式报错。iMTD 无结果时不静默回退 ETKDG;任何档位均不声称全局覆盖或全局最低。契约全文见 `docs/webqm-conformer-subset-schema.md`。

---

## 8. 能量景观(`ImtdSearch::landscape` → 前端 SVG)

- **轨迹梯子(Ladders)折叠统计**:汇总元动力学与验证期 MD 轨迹,统计总帧数、盆地访问频次及转迁对(`pairs` + `count`)。
- **不连通性树(Disconnectivity Forest)**:基于 Kruskal 最小生成树算法构建;能垒阈值优先采用 NEB 精修后的鞍点能量(若已就绪)而非原始路径带峰,树结构随鞍点精修进度动态更新。
- **盆地网络图(Basin Network)**:采用确定性力导向布局(固定 220 次迭代;布局坐标按「盆地数 | 边表 | 视口 | 梯源」联合键缓存,点击选中节点仅触发重绘而不重算布局)。视觉编码:边宽 ∝ 转迁计数,节点颜色 ∝ 能量,节点半径 ∝ 占据率。
- **数据表输出**:
  - **盆地表**:输出势能 `E`、准谐自由能 `G(QH)`、`Boltzmann.E`、`Boltzmann.G`、简并度 `g` 及轨迹占据率。
  - **转迁表**:输出稀疏转迁矩阵三元组 `a` / `b` / `count`。
- **边界退化处理**:对单盆地刚性分子体系(经 aspirin、naphthalene 验证)能够正确退化并产出空景观。

---

## 9. NEB 与鞍点精修(`src/neb.rs`)

- **CI-NEB 基础管线**:支持端点 Kabsch 对齐、改进切线估计(Improved Tangents, Henkelman 2000)以及端点带内重淬灭;自 v1.13.0 起**默认开启攀岩像(Climbing Image, CI)**。
- **积分器选择**:默认采用全局自适应 **FIRE** 积分器(Bitzek 2006),提供 QuickMin 作为备选。以丁烷扭转边为例,由 QuickMin 的 108.8 ms / 1937 迭代提速至 FIRE + 攀岩的 41.5 ms / 803+24 迭代;FIRE 的鞍点精度经受限扭转扫描独立仲裁(测得势垒 3.92 vs 扫描参考值 3.95,而 QuickMin 存在 0.6 的离散化下垂)。
- **收敛像冻结优化**:当单个像的受力连续 40 次迭代满足 < tol/2 时,跳过其力场评估与坐标积分;**峰值像(鞍点载体)永不冻结**。结果中导出 `evals_skipped` 与 `n_frozen` 诊断字段。
- **攀岩策略(钉住 CI)**:
  - **激活与角色锁定**:待相位 1 全带收敛后激活攀岩;攀岩像一经选定即永久钉住、绝不中途重指(彻底根除旧版因角色跳换引发的极限环振荡,已证实非积分器缺陷)。
  - **严格内极大守卫**:要求候选峰不得与两端点直接相邻,且能量必须同时高于双端至少 0.1 kcal。
  - **相位 1 快照回退**:若攀岩收敛结果反常低于相位 1 带峰达 0.1 以上,则触发安全机制,诚实回退至相位 1 带峰快照并标记 `climb_started=0`。
  - **能量失速判据**:在鞍点脊两侧实现 1e-4 kcal 级能量定位收敛。
- **粗带下垂修复与大分子验证**:
  - 在粗离散带(≤10 像)实测中显著修复离散化下垂:8 像体系由 -1.325 改善至 -1.155,10 像体系由 -1.262 改善至 -1.132(受限扫描参考值为 -1.09)。
  - 大分子全链路验证(ibuprofen / threonine):所有触发攀岩的边均正常收敛或诚实回退;个别预存的相位 1 硬边(涉及复合扭转、耗尽 20000 次迭代)属路径拓扑本身难度,与攀岩模块无关。
- **力场限制**:**GFN-FF 明确拒绝 NEB**。因 GFN-FF 的成键拓扑感知强依赖初始几何,端点对齐与插值会导致拓扑感知失真,调用时直接报错:`"NEB refinement supports MMFF engines only"`。
- **计算成本**:丁烷单边耗时 41.5 ms;大分子单边耗时 ~1.2–2 s(原生单线程串行测试;浏览器端通过 Worker 按边并行调度,总墙钟时间取决于最慢边)。

---

## 10. 振动分析与热力学(`src/vib.rs` + imtd thermo)

- **频率与准谐近似(QH)**:通过有限差分(FD)构造质量加权 Hessian 矩阵并求解带符号振动频率;设立**准谐(Quasi-Harmonic, QH)低频地板 100 cm⁻¹**(常量 `QH_CUTOFF_CM`)。`harmonic_correction` 函数统一输出 (ZPE+热校正, S, 模式计数),底层物理常数与单位换算由甲烷标准振动带严格锁定。
- **玻尔兹曼分布与简并度**:
  - `boltzmann_populations`(及自由能版本 `boltzmann_populations_g`)在计算热力学占据率时自动按构象简并度 g 加权。
  - **简并度 g = 实际观察到的对称等价势能井数**:通过变体几何计数算法,统一覆盖镜像孪生、原子置换旋转及 C₃ 甲基转子等情形(满足 g ≤ 点群对称数)。乙醇(σ=3)与丁烷(3–5 井)均设有专项测试锚定。
  - 能量景观盆地表(见 §8)中的 `Boltzmann.E`、`Boltzmann.G` 及最终占据率均由上述计算链条导出。
- **性能**:乙醇 `frequencies` 完整振动分析耗时 ~290–390 µs。

---

## 11. 3D 分子形状与指纹(`src/shape`)

- **3D 形状对齐**:基于 Grant–Pickup 高斯体积重叠模型(采用 `shape-it` 常数体系),定义形状相似度 Shape Tanimoto = V_AB / (V_A + V_B − V_AB)。采用多起点初始化结合 L-BFGS 刚体对齐优化,支持可选的药效团色版(color)对齐模式。
- **精度验证**:相同分子自对齐得分为 1.000000;与 RDKit 网格积分法计算结果的秩相关系数达 ρ=0.96。
- **2D 指纹与合成可及性**:同属该 WASM 导出面(`shape_align_wasm` / `ecfp4_fingerprint_wasm` / `sa_score_wasm` 等)还集成了 ECFP4 摩根指纹(含 Tanimoto 相似度)以及与 RDKit `sascorer` 位等移植的 Ertl 合成可及性评分(SA score)。

---

## 12. 性能速查表

> 评测基准:`examples/perf_probe.rs`,Apple Silicon M 系 release 编译档(注有「浏览器」者为 WASM 实测)。

| 评测路径 | 性能指标 |
| :--- | :--- |
| `MMFF energy_and_gradient`(乙醇,9 原子) | 1.7–3.6 µs(核内零分配) |
| L-BFGS 几何淬灭(300 步上限) | ~110–120 µs |
| iMTD 构象搜索(乙醇,2×3×800) | 17–22 ms |
| iMTD 构象搜索(ibuprofen,2×3×800) | 220–300 ms(原生) |
| NEB 单边精修(丁烷,FIRE + 攀岩默认) | 41.5 ms / 827 迭代 |
| NEB 单边精修(大分子,ibuprofen 级) | ~1.2–2 s(原生串行) |
| `vib::frequencies` 振动分析(乙醇) | ~290–390 µs |
| MD 单步积分(MMFF NVT,33 原子) | ~27 µs(浏览器) |
| ETKDG 3D 嵌入(ibuprofen) | ~11 ms(浏览器) |

---

## 13. 边界限制与诚实降级

- **明确不涉足的领域**:生物大分子/蛋白质体系、周期性边界条件(PBC)体系、GFN-FF 以外的量子力学(QM)计算、显式水分子溶剂化模型、QSAR/ADMET 数据驱动预测。
- **GFN-FF 能力边界**:禁止进入 NEB 管线(主动拒配报错);开放用于构象搜索、几何优化与分子动力学。
- **NEB 攀岩降级语义**:面对粗离散带(≤10 像)或复合扭转硬边时,算法可能触发快照回退或未能收敛——返回结构始终携带明确的 `converged` 与 `climb_started` 状态标志,供上层消费方直接读取而无需猜测。
- **Dimer 单端鞍点搜索**:**当前关闭**(现有管线均为双端已知的盆地间连边场景,无单端盲搜需求;未来若立项反应动力学/速率常数功能再行评估)。
- **外部环境依赖项**:`m1_3d_pipeline` 验证套件依赖宿主机的外部 `xtb` 二进制可执行文件(属本机测试环境配置项,非 WebMM 引擎自身缺陷)。

---

## 14. 验证状态汇总

| 验证维度 | 守卫规模与覆盖详情 |
| :--- | :--- |
| **引擎单元测试** | `cargo test` 共 372+7 项(含 6 项 NEB 核心守卫:势垒带/端点固定/扫描参考收敛/粗带回退/退化拒绝/确定性) |
| **属性测试(Property Testing)** | 3 项核心物理不变量性质(基于 `proptest`,固定 seed) |
| **浏览器端到端(CDP)** | 36 个测试套件 + platform 153 项断言(含 `m22` golden 24 项接口面锁定守卫) |
| **RDKit / `xtb` 基准对拍** | MMFF 230/230 门禁测试(`scripts/benchmark_mmff.py`)、GFN-FF 32/32 门禁测试 |

---

## 附:相关文档地图

- **接口冻结与版本纪律** → `docs/STABILITY.md`
- **任务信封协议与执行器设计** → `docs/task-contract.md`、`docs/task-executor-design.md`
- **Worker 消息协议面** → `docs/worker-protocol.md`
- **力场参数与移植笔记** → `docs/atom-type-coverage.md`、`docs/gfnff-porting-notes.md`
- **能量验证与误差分析** → `docs/validation-energy-analysis.md`
- **开发轮次日志(历史排雷存档)** → `CODE_STATUS.md`、`docs/CODE_STATUS-archive.md` 及各轮次 `PLAN.md`

<!-- optimized by gemini-4-argon, facts frozen v1.13.0 -->
