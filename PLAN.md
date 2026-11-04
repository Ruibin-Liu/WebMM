# Plan: 设计循环闭环——骨架探索结果一键推送到分子工作台(已完成)

用户点题:骨架探索(neighborhood explorer / aza-scan)设计出的分子,在结果列表里可一键推送到分子工作台(single 视图管线),形成 工作台→探索→新分子→回工作台 的循环。
**性质:用户明示的功能需求(冻结豁免);纯前端(app/index.html),零接口面改动(m22 golden 不涉 app DOM)。**

## 现状依据
- 探索结果:`renderAnalogRows`(表 #analogRows,行含 r.idx/r.smiles/r.mw/tanimoto/rotb/sa/flexT + ★pin 列);`analogLastResults()`/`__analogResults` 保底;aza-scan 复用同一表同一渲染 → 同列自动受益。
- 工作台入口:`#input`(SMILES/MolBlock)+ `process(true)`(RDKit 解析→2D 渲染→历史/URL);`switchMode('single')` 切单分子视图。
- 测试:m6 系只查 rows>0/fact 字段/thumbnail,不锁列数;fact 层 regex 在项目态不在 DOM。

## 交付
1. **UI(app/index.html)**:analog 表头 + 列「⇗」(title:推送到工作台);每行按钮 `analogToBench(idx)`——从 analogLastResults() 取 smiles,写 `#input`,`switchMode('single')`,`process(true)`,`#analogStatus` 反馈「analog #N → 工作台」;窗口函数挂 window(analog 同款模式)。
2. **测试(tests/cdp/m24_bench_loop.test.js)**:mock `__analogResults`(两条假记录)+ 直调 `renderAnalogRows` → 点 ⇗ → 断言:input.value === 该行 SMILES、模式切到 single、analogStatus 含反馈、无 page error;真跑一条小探索再推一条作 E2E(沿用 m6 模式,小分子快速档)。
3. **文档**:README 能力表 Ligand workbench 行补"design loop"短语;CODE_STATUS 记用户豁免的功能例外;ENGINE-CAPABILITIES 不涉(引擎面无改动)。

## 验收(实测)
m24_bench_loop 11/11(双页接线 8 项 + E2E 3 项:真探索→⇗→single 渲染);run.sh 36/37(唯 m1 xtb 环境项,m6 系全绿不受新列影响);m22 golden 零差异(零 Rust/wasm diff);run.sh 套件注册表加 m24。

## 排雷存档(本轮)
- **真实行 idx 是数字**(worker 索引,results.push({...cands[msg.idx], idx: msg.idx})),而 onclick 模板插值永远是字符串——严格相等静默失配(mock 测试用字符串 idx 完全测不出,E2E 才暴露);跨 DOM 边界的 id 比较一律 String() 两侧。
- **switchMode 会按模式暂存/恢复 #input**:先填输入再切模式会被恢复覆盖——填框必须在 switchMode 之后;mock 接线测试(没经过模式切换)对此不敏感,顺序语义靠 E2E 兜底。
- 工作台(index.html)的 analog 面板在 single 模式**内部**(switchMode('analog') 实际落在 single),平台页(platform.html)的 analog 是独立搜索模式段——两页模式机制不同但同名函数,dual-page 补丁必须各自核对语义。
- 两页内联脚本是手工同步的副本(无共享模块)——改一处必须同步另一处,m24 双页各测一遍。
