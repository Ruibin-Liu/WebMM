# Plan: 两个提速杠杆 —— 起点预筛重打分(rescore_top)+ 初筛构象降档(v1.6.1)

## 测量结论(先行)

- 预制成本分解:ibuprofen iter=250 → 105.7ms,iter=1 → 64ms(嵌入+建场
  地板),iter=50 → 57.5ms——**ETKDG 嵌入占 55-60%,iter=50 的松弛近乎
  免费**(aspirin 同构:36.5 / 19.7 / 31.7)
- 对齐成本(上轮已知):每起点 ~24ms 全量 IE 重打分占 ~95%

## 杠杆一:起点预筛后重打分(引擎)

- `AlignOptions.rescore_top: usize`(默认 **3**,opts_json
  "rescore_top";usize::MAX = 旧行为全重打)
- align_colored 两段:先跑全部起点收集代理分(便宜),按代理排序后仅
  对 **top-K 姿态**做全量 IE + 颜色重打分;polish 从最优全量分姿态出发
- 等价性测试:fixtures 8 分子两两配对(K=3 vs K=MAX),断言 combo 差
  < 1e-3 并报告实际最大差;既有测试(自对齐 1.000000 等)不回归

## 杠杆二:初筛构象降档(页面)

- 两段式的 phase 1 预制改 `max_iter: 50`(SCREEN_PREP_ITER;嵌入价买
  半松弛几何);phase 2 的 top-50 仍用完整 250 预制(e.sdf3d 与
  e.sdf3dScreen 双缓存,ensureEntry3D 增 tier 参数)
- 小库单段路径不变(全 250)

## 验收

1. Rust:等价性测试 + 293→294+;clippy 1.98(-D warnings)、fmt;版本
   1.6.0→1.6.1(wasm 重建,node 冒烟)
2. 性能断言(node + Playwright):ibuprofen 单对全质量对齐 472ms →
   预期 ≤150ms;55 库两段式 10-11s → 预期 ≤6s
3. 质量断言(Playwright):两段式(降档预制 + top-3 重打)top-20 vs
   全质量单段参考,召回 ≥ 18/20(旧口径 19/20,允许小幅让步并如实
   记录)
4. 390px、零 page error、三模式回归
