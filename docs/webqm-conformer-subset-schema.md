# WebMM↔WebQM 构象子集接口契约(提案已双方确认,待 WebMM 用户批准排期)

> 通道:herdr-app hc(t8 WebMM-dev ↔ t7 WebQM-dev),2026-11-05。
> 性质:**已实现 v1.14.0**(additive minor:m22 golden 导出清单 +1;导出 conformer_subset_json,
> 按 docs/STABILITY.md 纪律;内核 conformer_subset_core 为原生测试入口)。WebQM 不是本仓构建依赖。

## 分工边界(已确认)
- WebMM 提供小分子构象子集(生成/去重/力场能序);不做 QM 溶剂/梯度
  (WebMM 的 GFN-FF 是力场移植,非 QM 引擎)。
- WebQM 做片段级 QM(局部性质);不做整分子构象搜索。
- 片段切割/封端/局部电荷/质子化由 WebQM 侧另行定义——**切割不是坐标子集**,
  WebMM 只提供整分子构象,不预设切割语义。
- CREST/构象集合不是 WebQM 当前工作的前置阻塞项。

## 端点(方案 B,双方确认)
`conformer_subset(molblock, { n_target, energy_window_kcal, rmsd_tol, seed, strategy })`
→ JSON:

```json
{
  "provenance": {
    "webmm_version": "…",      // 引擎版本锚点
    "seed": 42,                 // 同 seed 同输入逐位一致
    "energy_engine": "MMFF94s" // 能序的力场标注
  },
  "conformers": [
    { "rank": 1, "energy_kcal_mol": -12.34, "coords_angstrom": [[x,y,z], …], "cluster_id": 0 }
  ]
}
```

## 字段与语义(按 t7 回执锁定)
- 字段名固定:`coords_angstrom` / `energy_kcal_mol` / `energy_engine`
  (显式单位命名,不用泛 `coords`)。
- `strategy` 默认 `'lowest_per_cluster'`(ETKDG 嵌入池,每 RMSD 簇取最低能);备选
  `'energy_stratified'`(ETKDG 池能量分层抽样);高采样档
  `'imtd_lowest_per_cluster'`(v1.15.0,CREST iMTD-GC 谱系,固定 2×3×800
  小规模预算后按簇取最低能)。iMTD 无结果时显式报错,不静默回退 ETKDG;
  固定预算不声称全局覆盖或全局最低。provenance 增 `strategy` 记录采样档。
- **硬上限(实现必须提供)**:`n_target`、`energy_window_kcal`、`rmsd_tol`
  三参数各自的合法上限与拒绝语义(超限报错,不静默钳制)。
- consumer(WebQM)侧可自行记录 `webqm_version` 于其 provenance——双向版本锚点。
- energy 为 MMFF94s 力场序,供 WebQM 片段级 QM 重排;非 QM 能量。

## 实现入口(排期后)
- `src/lib.rs` 新增 wasm 导出 `conformer_subset`(冻面前为 0 新增,故 additive);
- 复用既有管线:ETKDG 嵌入 → MMFF94s 优化 → 两层去重(imtd.rs 同款语义);
- 测试:确定性(同 seed 位等)、硬上限拒绝语义、簇代表正确性、
  strategy 两模式差异;浏览器侧 m25 级 CDP 套件。
