# Plan: Metad 构象搜索 iMTD-GC 迭代-交叉调度(引擎 + WASM)(已完成)

用户指令「先1」:按已记录的设计基线实施——CREST iMTD-GC 迭代路线:多轮短程 metad(轮换偏置强度与 CV 组合)收集构象,低能构象遗传交叉产新种子,迭代至系综收敛;不做单次长程 metad;现有 well-tempered metad 引擎原样保留,迭代-交叉调度加在外层。

## 交付
1. **原生调度器 `src/metad/imtd.rs`(新文件;`src/metad/mod.rs` 仅加 `pub mod imtd;`——well-tempered 引擎零语义改动)**
   - **CV 池(确定性)**:可旋转二面体 CV——无环重原子键 j–k(两端各存在邻居),i/l 取端点邻居(优先重原子、最小下标);池不足(<2)时补距离 CV(图距 ≥4 的重原子对,确定性排序截断);池空(刚性分子)→ 仅淬灭输入为唯一构象,结果如实记录(rounds=0/1 构象)。
   - **轮次调度**:每轮 `runs_per_round` 条短程 metad 顺序执行;CV 按全局 run 计数轮换池,偏置强度沿确定性梯 [1.0, 0.5, 0.25]×base 轮换;MD 种子 = base_seed + r·1000 + m;起点结构:round 0 = 输入(优化后),后续轮 = 交叉新种子优先 + 既有低能极小轮换。
   - **收集**:每条 run 每 `quench_interval` 步淬灭快照(共享 MMFF + L-BFGS `optimize`,迭代封顶)得局部极小;能量窗(默认 6 kcal/mol,相对系综当前最小)+ **对齐 RMSD 去重**(重原子;Horn 四元数 Kabsch,复用 `sym_jacobi`,坐标级私有实现,只允许真旋转)入系综。
   - **遗传交叉**:每轮末对低能极小配对,沿可旋转键 (j,k) 切片(图连通分片,H 随母重原子):B 按 F1 重原子对齐到 A,child = A 的 F1 坐标 + B 的 F2 坐标,淬灭后按同规则去重入系综;新极小即下轮种子。
   - **收敛**:某轮新增唯一极小 ≤ `convergence_new_minima`(默认 0)→ 收敛停止;硬上限 `max_rounds`;系综容量 `max_conformers`(保最低能)。全程确定性:同 seed 同输入逐位可复现,无墙钟。
   - API:`ImtdConfig`(Default 电池参数)+ `imtd_conformer_search(mol, mmff_variant, config) -> ImtdResult`(coords/energies/round_found/每轮新增计数/rounds_executed/converged/pool 大小)。
2. **WASM 导出(加法,不破坏既有契约)**:`ImtdOptions`(setter 风格镜像 MetaDOptions)+ `run_imtd_conformer_search(sdf, options) -> ImtdResult`(扁平数组 conformer-major,镜像 OptimizedConformers 形状,success/error 通道同款)。构建后核对 pkg/webmm.d.ts 仅新增符号。
3. **测试**:imtd.rs 单元测试(CV 池:ethanol 两条可旋转二面体;对齐 RMSD 旋转恢复/镜像拒绝;交叉子代有效性;端到端 n-butane/ethanol 小预算找到 ≥2 极小;methane 刚性 1 极小;确定性逐位复现;max_rounds 封顶)+ lib.rs WASM 测试(数组一致性/确定性/坏输入 error 通道)。
4. **不做(记录,后续轮次)**:前端 workbench/任务层接线(把 IMTD 接进 conformers kind 或独立 kind)、GFN-FF 支持(v1 仅 MMFF,镜像现有 metad WASM 面)、并行多 metad(顺序执行;worker farm 归前端层)。

## 验收
cargo fmt clean;cargo clippy --all-targets 双目标(native + wasm32)0 警告;cargo test 全绿(含新测试);wasm 重编且 webmm.d.ts 仅增量;run.sh 既有套件零回归(无新前端套件)。

## 排雷存档(本轮)
- **Kabsch 方向约定**:lib.rs `kabsch_align` 的 Horn 实现是 R·p+t≈q(p=动、 q=参考);照抄数学时把残差方向对调(p↔q 互换)会静默给出错旋转——测试用已知旋转+平移抓出(rmsd 2.59 而非 ~0)。坐标级实现必须同约定。
- **fragment_split 点标记陷阱**:用 seen[k]=true 封堵制会把割边端点 k 排除在两侧之外(丁烷 13≠14)——必须按边封堵(只断 j→k 遍历,k 归 F2)。
- **可旋转键规则**:若要求两端各≥2 个重邻居会丢掉羟基因子(H-O-C-C);放宽为“任一邻居”后丁烷得 3 条(端基 C-C 的 H 扭转也是合法 CV),交叉测试须显式取中央键 (1,2)。
- **断言阈值**:120° 片段旋转下 4 重原子对齐 RMSD 仅 ~0.36(对齐会抵消单端位移),去重/差异断言阈值取 0.2 而非 0.5。
- **收敛阈值语义**:convergence_new_minima=usize::MAX 会让第 1 轮即判收敛(每轮新增≤MAX 恒真)——“预算耗尽”路径无法用阈值构造,测试改验“达阈即停”(rounds=1、converged)。
- **环境项(非本轮四引起)**:m1_3d_pipeline 的 xtb 交叉核对依赖 `~/.local/xtb-gxtb/bin/xtb`,本机缺失即 exec 127(与 wasm 内容无关);m5 amine_smarts 全量跑偶发 flake(单跑 49/49 过,与 m12 引擎加载偶发同族)。

## 验收(实测)
cargo fmt clean;clippy native(--all-targets)与 wasm32(--lib)均 0 警告;cargo test 338(lib,+9 新)+ 7(sascore_golden)全绿;wasm 重编后 webmm.d.ts 增量纯净(仅 ImtdOptions/ImtdSearchResult/run_imtd_conformer_search 三符号);run.sh 26/27 CDP 套件全绿(m1 为上述 xtb 环境缺失,exec 127 不达页面)+ Node platform 153/153。
