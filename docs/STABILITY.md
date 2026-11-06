# 接口稳定性政策(自 v1.12.0 起)

WebMM 自 v1.12.0 起处于 feature freeze。本文件定义哪些面是稳定接口、如何
安全变更。守卫:`tests/cdp/m22_api_surface.test.js` 对下述各面做**精确相等
快照比对**,增删改任何条目都会红灯。

## 稳定域

1. **WASM 导出面**(`pkg/webmm.d.ts`):全部 `export function` 与
   `export class` 的方法集(wasm-bindgen 内部符号 `init`/`initSync`/
   `free`/`toJSON` 等除外)。
2. **Worker 消息协议**(见 docs/worker-protocol.md):conf/shape/analog
   三 worker 的 `msg.type` 集合与载荷形状。
3. **任务契约**(docs/task-contract.md):envelope 版本与 kind 注册表
   (`analogScan`/`conformers`/`imtd`/`optimize`/`shapeSearch`),各 kind
   input v1 冻结(imtd 为 v1.1 含可选 phase/sourceHash)。
4. **JS 全局面**:`window.webmm`(wasm 模块)、`window.IMTD`(驱动:
   start/runParallel/optsToWasm/DEFAULTS;handle: promise/verify/refine/
   thermoJson)、`window.Tasks.*`(contract/scheduler/pools/funnel/store)。
5. **文件布局**:`app/`(index.html、platform.html、三 worker、tasks/、
   platform/、vendor/)、`site/`(演示页)、`pkg/` 构建产物布局
   (scripts/stage-pkg.sh 镜像 pages.yml)。
6. **引擎原生 pub 面**:`src/lib.rs` 顶层 pub 项、`metad::imtd` 的
   pub 类型与函数、`neb`、`vib` 模块公开签名。

## 变更纪律

- **patch(1.12.x)**:仅缺陷修复与性能优化,任何稳定域零差异
  (m22 必须保持全绿且 golden 不动)。
- **minor(1.13+)**:允许**纯增量**(新导出/新消息类型/新 kind),流程:
  1) 更新本文件或 worker-protocol.md;2) 更新 m22 golden(显式);
  3) 实现与测试;4) 版本号三处对齐(Cargo.toml/package.json/
  webmm_version)。
- **major(2.0)**:任何移除/改名/语义破坏。当前无计划。

## 已知非稳定(冻结线以下,不承诺)

- 内部函数签名(非 pub)、私有结构体字段。
- `pkg/` 内部生成物(随构建变化,仅导出面受保护)。
- 性能特征(迭代次数、耗时)——不属于接口。

## 冻结时间线

- 2026-11:feature freeze 宣布(a5a0b42 收官轮)。
- v1.12.0:接口 golden 固化(本文件 + m22)。
- v1.14.0:additive 新增——wasm `conformer_subset_json`(WebMM↔WebQM 构象子集契约,docs/webqm-conformer-subset-schema.md;硬上限显式报错;确定性 seed;Rust 内核 conformer_subset_core 供原生测试)。golden 更新:导出清单 +1。
- v1.13.0:数值语义升级(接口面零差异,m22 golden 不变):①NEB 积分器默认 QuickMin→FIRE(2.7×,精度经受限扭转扫描独立仲裁更优);②攀岩相默认开启——refine 输出从"带峰上界"升级为"CI 精修鞍点"(钉住攀岩像+相位 1 快照回退;粗带爬坏时诚实回退带峰,climb_started=0)。NebConfig 增 Rust 侧内部字段 integrator(不入 wasm 面)。
