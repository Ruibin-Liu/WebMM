# Plan: WebMM↔WebQM 构象子集接口(方案 B 实现,additive minor → v1.14.0)(已完成)

用户批准启动(此前 t8↔t7 契约互锁,docs/webqm-conformer-subset-schema.md 定稿)。
性质:**additive 变更**(冻结面只增不改,m22 golden 显式更新 + 版本 bump 三处 + STABILITY 增条)。

## 交付
1. **Rust 端 `conformer_subset_json`(src/lib.rs 新增导出)**
   - 输入:molblock/SDF + n_target + energy_window_kcal + rmsd_tol + seed + strategy('lowest_per_cluster' | 'energy_stratified')。
   - 管线:ETKDG(seed 派生)→ 几何加氢 → MMFF94s 优化 → 能量窗过滤 → 全对重原子 Kabsch-RMSD 簇判定(贪心,按能量升序代表)→ strategy 子集选择;坐标扁平化 [x,y,z,…] 序列。
   - **硬上限 + 显式报错(Err,不静默钳制)**:n_target ≤ 200 且 >0;energy_window_kcal ∈ (0,50];rmsd_tol ∈ (0,2];strategy 仅两枚举值;嵌入失败/空窗 → Err 诚实报错。
   - 输出 JSON 契约:provenance{webmm_version, seed, energy_engine:'MMFF94s'};conformers[{rank, energy_kcal_mol, coords_angstrom, cluster_id}];坐标嵌套数组;能量排序升序 rank=1 起。
   - energy_stratified:唯一能量序列(1e-3 kcal 容差合并)等距取代表,按窗口幅度自适应。
2. **单测(src/lib.rs tests 段)**:乙醇 n_target=3 命中 3(最低能 gauche?实测确定);同 seed 位等;空窗报错;超限报错(n_target=0/201,窗口=0/100,rmsd=0/5);坏 strategy 报错;两 strategy 结果差异;naphthalene 刚性分子簇数=1;三策略下 cluster_id 自洽。
3. **m22 golden 更新**:导出函数清单 + conformer_subset_json(唯一新增);版本三处(Cargo/package/CODE_STATUS)+ STABILITY v1.14.0 条目;契约文档标注“已实现”。
4. **完成通知 t7**(hc):接口面、示例输出、测试状态。

## 验收(实测)
cargo test **377+7** 全绿(+5 subset 测试:硬上限拒绝/确定性+排序+坐标维数/萘单簇/双 strategy 差异/能量窗过滤);clippy 双 0(needless_question_mark 一处);m22 golden 更新后 24/24(导出清单 +conformer_subset_json,版本 1.14.0 三处一致);重建 wasm 后 golden 过;契约文档标"已实现 v1.14.0"。
实现定稿:conformer_subset_json(wasm 壳)→ conformer_subset_core(原生内核,JsValue 在 non-wasm32 panic 故测试走内核);oversample 6×(pool ≤400);贪心重原子 Kabsch-RMSD 簇代表;两 strategy(lowest_per_cluster 默认 / energy_stratified 能量梯等距);全部超限显式 Err。

## 排雷存档(本轮)
- **JsValue 在 non-wasm32 目标上 panic**(wasm-bindgen 0.2.106 拦截)——wasm 导出函数不能被原生单测直接调用;拆 壳(JsValue 边界)+ 内核(纯 String)双层,原生测试只调内核。所有未来新增 wasm 导出沿用此模式。
- aligned_rmsd_over 原为 imtd 私有 fn——提升为 pub(crate) 复用(簇判定与两层去重同一 Kabsch 语义)。
- 乙醇 MMFF94s 全局最低能实测 **-1.517**(非印象中的 -3.5)——测试常量必须跑完实测再写断言,凭记忆写必返工。
- golden 更新走纪律:导出清单 +1 行 + 注释标 additive 版本 + 三处版本号 + STABILITY 条目,一次到位。
