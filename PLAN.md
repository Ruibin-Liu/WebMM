# Plan: Workbench file:// 打开时的裸 TypeError 修复(RDKit 未加载防御 + 可操作提示)

## 背景

用户报告:`TypeError: Cannot read properties of null (reading 'get_mol')`。

诊断(Playwright 实测):
- gh-pages 部署实例健康(RDKit 2026.03.6 正常加载);
- 本地 http 服务健康;
- **file:// 直开 app/index.html**:vendored `RDKit_minimal.wasm` 的
  fetch 被 file:// 阻断 → `initRDKitModule()` reject(RuntimeError:
  Aborted — both async and sync fetching of the wasm failed)→
  `rdkitModule` 保持 null → 任何交互(process() 第 1144 行
  `rdkitModule.get_mol(input)`)抛出所报 TypeError。error 框其实已有
  "Failed to load RDKit: …" 但交互路径给出的是裸 TypeError,且无
  file:// 归因与修复指引。

site 两页在 file:// 下同样不可用(webmm.js 动态导入失败),属同类
环境限制;本计划只修 app 的报错质量,不改架构、不加 file:// 支持
(wasm 需 http)。

## 任务(app/index.html,纯 JS 防御层)

1. `initRDKit()` catch 增强:`location.protocol === 'file:'` 时给出
   归因 + 指引("serve over HTTP, e.g. `python3 -m http.server` from
   the repository root, then open http://localhost:8000/app/");
   非 file:// 保留原始错误文本。
2. 新增 `rdkitAvailable()`:rdkitModule 就绪返回 true;否则把上述
   可操作提示写入 #error 并返回 false。
3. 用户入口防御:`process()`、`runBatch()` 顶部调用
   `rdkitAvailable()`,未就绪直接 return(取代裸 TypeError)。
   其余按钮(Embed/Opt/Conformers 等)本来就 disabled直到分子
   加载,恢复路径经 process(),无需逐个加。

不做:file:// 下强行可用(wasm 体积/构建不可行);site 两页同类
提示(另行立项);任何引擎/vendor 改动。

## 验收(实施后实测记录)

- Playwright 8/8:file:// 打开——输入 SMILES + Enter 与 Batch 入口均不再
  抛 TypeError,#error 显示 file:// 归因 + `python3 -m http.server`
  指引;http 打开——行为零变化(版本 2026.03.6、CCO 正常解析、
  error 框空、零 page error);
- file:// 下残留的唯一 pageerror 是设计内的 webmm.js 双布局探测
  (既有行为,非本次范围);
- `cargo test` 281/281、clippy 0、fmt 干净(零 Rust 改动)。
