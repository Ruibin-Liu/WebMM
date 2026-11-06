# Plan: conformer_subset 新增 iMTD-GC 采样策略(additive → v1.15.0)(已完成)

用户批准在既有 WebMM↔WebQM 契约中增加 `imtd_lowest_per_cluster`。保留 `lowest_per_cluster`(ETKDG 快速档)与 `energy_stratified` 的行为和默认不变;函数签名、导出清单与 JSON 字段均不变,仅扩展 strategy 合法值与采样来源。

## 交付
1. **采样策略接入(src/lib.rs)**:`strategy='imtd_lowest_per_cluster'` 走既有 `imtd_conformer_search`(SearchEngine::MMFF94s),固定小规模预算 `max_rounds=2, runs_per_round=3, steps_per_run=800`;seed、energy_window、rmsd_tol、n_target 沿用契约参数,系综容量按 `min(6*n_target,250)` 封顶。输出仍按能量升序、RMSD 簇代表、rank/cluster_id/coords_angstrom 契约生成;不声称全局覆盖或全局最低。
2. **错误语义**:新 strategy 纳入合法值校验;iMTD 无可用构象时显式 Err,不静默回退 ETKDG(调用方选择 CREST 谱系就必须得到该谱系结果或明确失败)。
3. **测试**:乙醇新策略确定性、坐标维数、能量排序与最低能锚;萘刚性单簇;新策略与 ETKDG 默认分离(provenance 标记 sampler 便于消费者区分);既有 5 项 subset 测试不回归。
4. **文档/版本**:契约文档补采样档位与固定预算;ENGINE-CAPABILITIES 同步;STABILITY 记录 additive 语义扩展,版本升 v1.15.0(Cargo/package 对齐;m22 导出清单零变化)。
5. 完成后 hc 通知 t7;commit+push。

## 验收(实测)
cargo test 379+7 全绿(+2:iMTD 乙醇确定性/能量序/坐标维数/provenance 与萘刚性单簇);clippy 双 0;m22 24/24(导出清单零变化,版本 1.15.0 对齐);run.sh 36/37(唯 m1 外部 xtb);新策略逐位确定,无 ETKDG 静默回退。

## 排雷存档(本轮)
- iMTD 路径先对输入重原子分子统一加氢再调用 imtd_conformer_search,保证输出坐标维数与 ETKDG 档一致(乙醇 9 原子);重原子索引在 H 追加后保持不变,簇判定可共用。
- 采样档与选择策略合在一个 strategy 枚举中,避免新参数和导出面变化;provenance.strategy 是 additive 字段,消费者可区分采样来源。
- 固定 2×3×800 是明确的小规模预算而非“完整 CREST 搜索”;文档同步强调不保证全局覆盖/全局最低,空结果必须报错而非换采样器。
