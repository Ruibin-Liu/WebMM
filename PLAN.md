# Plan: Batch 标签 E2E 补全——m4_batch 扩为全功能套件(描述符奇偶/排序/过滤/错误行/取消/SDF 输入/行回载/模式隔离)

## 背景

m4_batch 现仅 7 项(跑通+CSV 头+3D SDF provenance),不断言数值。
Batch 标签实际功能面远大于此:描述符/QED/PAINS 列、三态排序
(aria-sort)、五档过滤、错误行、Cancel、SDF 输入、行点击回载、
以及新落地的模式隔离。冻结参考 tests/fixtures/lbdd/refs.json
(51 分子描述符/QED/PAINS/Murcko,与页面端口 51/51 逐位)一直没被
batch 套件使用——参照 m5 模式把数值断言打上去。

## 实施(扩 tests/cdp/m4_batch.test.js,存量 7 项不动,追加新节;零产品代码改动)

1. **描述符/警示奇偶**:9 分子(caffeine/aspirin/ibuprofen/naproxen/
   paracetamol/rhodanine_nmethyl[PAINS=1]/metformin[QED<0.5]/
   cholesterol[Lipinski 失败 logP 7.39, QED<0.5]/hexane[TPSA=0,
   无骨架])+ 1 行垃圾输入 → 逐行对 refs 断言 MW/cLogP/TPSA(2dp)、
   HBD/HBA/RotB、QED(2dp)、PAINS 计数+badge class
2. **错误行**:错误文本渲染;filter=all 可见、lipinski 下隐藏
3. **过滤语义由 refs 现算**(非硬编码):lipinski/veber/qed≥0.5/nopains
   的期望可见名集合 vs 页面实际(覆盖胆固醇出局、metformin 出局、
   rhodanine 出局三个判别点)
4. **排序三态**:MW asc(己烷首)→ desc(胆固醇首),aria-sort/
   ▲▼ 同步;第三击 key 清空回 none(行序不复原是既有语义,不断言行序)
5. **CSV 内容**:逐行 QED mean 6dp、PAINS 计数、Murcko 字符串
   (含 metformin/hexane 空)与 refs 精确一致
6. **行点击回载 single**
7. **多记录 SDF 输入**(页内现造两条 molblock)按名解析成两行
8. **Cancel**:3D 批量中途取消 → status 'Cancelled'、batchRun 空、
   Run 按钮复位
9. **模式隔离**:batch 跑出结果后切 single(panel 隐、输入不串)、
   切回(批量文本+panel 恢复)
10. README 的 m4 描述行同步

## 验收

- `node m4_batch.test.js` 全绿;m5/m3 抽查回归

## 验收结果(实施后)

- `node m4_batch.test.js` **32/32,三连跑稳定**:存量 6 + 新 26——
  9 分子逐字段描述符/警示奇偶(MW/cLogP/TPSA 2dp、HBD/HBA/RotB、
  QED 2dp、PAINS 计数+badge class)、错误行渲染与过滤行为、五档过滤
  (期望集由 refs 现算,胆固醇/metformin/rhodanine 三判别点全中)、
  MW 三态排序+aria-sort、CSV 逐行 QED 6dp/PAINS/Murcko 精确、行点击
  回载、SDF 双记录按题名解析、Cancel 立即取消确定性断言、模式隔离
  往返
- **过程中的诚实发现**:refs.json 的 props.HBA 是 QED 口径(咖啡因 3)
  与 batch 表 NumHBA 语义不同;且 **RDKit 2026.03 把 NumHBA 从 N+O
  计数改为严格受体口径(咖啡因 6→3)**——vendored wasm 2026.03.6 与
  Python 2026.03.6 的 CalcNumHBA 对 9 分子逐值一致,本地 homebrew
  Python 3.14 的 2025.09.3 只在咖啡因上发散。故新设
  scripts/gen_batch_refs.py(强制 RDKit≥2026.03 解释器门禁)→
  tests/fixtures/lbdd/batch_refs.json 作为 batch 专用金标
- 回归:m5 26/26;m3 9/10(axe nested-interactive 存量)
