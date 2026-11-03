# Plan: 简并度恢复 + 谐振自由能校正(占据科学性一轮交付)(已完成)

用户裁决「按流程来」= 设计页建议项 ①+② 连做:恢复被两级去重折叠的镜像简并度(Boltzmann 乘 g_i),并把布居从能量升级为自由能(G = E + 谐振 ZPE + S_vib 项)。

## 交付
1. **简并度恢复(src/metad/imtd.rs)**
   - 能量判据分支细分:折叠发生时做**镜面叠合检验**(固定反射 z→−z 后再做真旋转 Kabsch=无叠合)——反射叠合 < tol ⇒ 镜像对 → `degeneracy[i] += 1`(折叠保留 + 记 g);否则 ⇒ 甲基置换 → 普通折叠(g 不变)。C3 类置换相关筒并(methanol 甲基转动)不计——v1 诚实边界,记录。
   - Ensemble 增 `degeneracy: Vec<u64>` 平行向量(默认 1,逐出同步收缩);landscape basins 增 "g";Boltzmann 升级 `boltzmann_populations(energies, degeneracies, t)`(权重 g·exp(−ΔE/kT));旧签名 wasm 导出保留(g=1 包装),新增 `imtd_boltzmann_populations_g_json`。
   - 语义:镜像帧仍是该折叠盆地的"访问"(占据归入折叠条目);g=该盆地观测到的**不同井数**(≤2)。
   - 测试:手工构造乙醇 OH 转动异构三元组(H8 绕 C–O 轴 0°/±120° 淬灭)经 gather_verify 折叠 → 2 成员、basins g=[1,2];丁烷镜像对(重原子镜像,不折叠)各自 g=1 且布居求和天然等价 g=2(断言)。
2. **谐振自由能校正(新模块 src/vib.rs)**
   - `frequencies(ff, coords, masses) -> Vec<f64>`(**带符号** cm⁻¹):解析梯度中心差分 Hessian(h=0.005 Å,3N×2 次梯度),质量加权,`sym_jacobi` 对角化;单位换算常数从基本常数推出,甲烷/水频带测试锁定(最高频 2800–3300 等)。
   - `harmonic_correction(freqs, t) -> (g_corr, zpe, n_imag)`:<25 cm⁻¹ 模式剔除(平动/转动残迹);虚频计 n_imag 并按 |ν| 计入(诚实旗标);G_i(T) = E_i + Σ[½k θ + kT ln(1−e^{−θ/T})],θ=1.43877·ν̃ K;同分子转动/平动熵自动相消故不计。
   - 测试:H 对称;淬灭甲烷 6 个近零模 < 25;频带;T→0 ⇒ G_corr→ZPE;确定性。
3. **协调器 + WASM(加法)**
   - `ImtdSearch.thermo_json(t)`(原生+wasm 方法):对最终系综逐盆地算频率与 G(引擎同系综,GFN-FF 亦可——FD 梯度可得);输出 {i, G, gCorr, zpe, nImag, nuMax, nuMin}。
   - d.ts:ImtdSearch +1 方法 + 新自由函数;既有导出零变化。
4. **UI(index.html)**
   - 盆地表:增「简并度 g」「G(T)」「Boltzmann(G)」列;Boltzmann(E) 列改乘 g;虚频盆地 ⚠ 标注;T 输入 onchange 触发 thermo 重算(协调器缓存按 (系综版本,T) 失效)。
   - 驱动 handle.thermo(t) 或直接协调器方法(结果已含)——走后者,thermo 属视图派生非批次。
5. **测试**:原生(vib 模块 + 简并度构造 + thermo 乙醇 n_imag=0、ZPE 量带);CDP 新套件 **m20_degeneracy_thermo**(确定性种子钉实测:g 列存在且至少一盆地 g=2;Boltzmann(E) 与手算 g·exp 一致;G 列有限、Boltzmann(G) 和=100%;换 T 更新;零 page error);run.sh 挂 m20。
6. **不做(记录)**:C3+ 置换筒并、准谐振(QT)低频处理、转动/平动熵显式项(同分子相消)、简并度对 metad/MD 采样占据的分解。

## 验收(实测)
cargo fmt/clippy 双目标 0(redundant_clone/len>1 三处);cargo test 360(lib,+3:镜像孪生折叠构造/简并玻尔兹曼/热力学层;vib 模块另 4)+ 7 全绿;wasm 重编 d.ts 增 imtd_boltzmann_populations_g_json(+此前轮的 neb/ boltzmann 基线)+ ImtdSearch.thermo_json 方法;run.sh 32/33 CDP 全绿(m1 xtb 环境项)+ m20 9/9;m18/m19 因表列变更同步断言(排雷);m6_platform 出现两次 import-pins flake,stash 二分 + 复跑证实非本轮回归(既有脆弱家族)。

## 排雷存档(本轮)
- **镜像判据被标签噪声污染**:固定标签下甲基近自由旋转使 proper/improper RMSD 同涨(实测孪生 vs 成员 0.411/1.178 与 1.178/0.411 互换)——绝对 improper 阈值失效;改**边际判据** proper−improper > 0.15(置换帧两者相等,天然排除),探针实证 g=2。
- **孪生构造两课**:±120° 旋转起点会被再淬灭把 C–C 扭转拖进另一转子(能量差恰为 anti/gauche 0.18);镜像孪生必须由**重原子平面反射**构造(反射极小=极小)。gather_verify 前忘 plan_verify 会静默丢弃(守卫返回空)。
- **表头/行单元格错位**:MD 占据列条件渲染表头但行恒发 11 格——m19 列索引读到 '–' 得 NaN;m20 硬编码索引恰逢错位巧合通过。统一恒发(占位列未验证时显示 –),五套件列索引重核。表头文案改名 Boltzmann→Boltz.E 也牵动 m19 正则。
- **ln(1−e^{−θ/T}) 符号**:exp_m1 少取负 → ln(负数) = NaN(T→0 暴露);直接写 1−exp。振动熵项恒负,G_corr(298) < ZPE(首版断言方向写反)。
- **平转残迹含小负本征值**:淬灭极小的 FD Hessian 六个近零模正负兼有——"真极小"断言须容 cutoff(谐波求和已 <25 cm⁻¹ 剔除,虚频只在 |ν|≥25 计)。
- 单位换算常数(108.6 cm⁻¹/√(kcal·mol⁻¹·Å⁻²·amu⁻¹))由甲烷 C–H 顶带 2800–3300 端到端锁定;水/乙醇频带与 ZPE≈48.8 量带交叉验证。
