# Plan: Phase 1 — WebMM Workbench（对齐 Molecule Clipboard 的轻量 ligand-based CADD 平台）

## 产品定位

在 `app/` 子目录新增静态单页应用（工作名 **WebMM Workbench**，UI 级名字、零 repo 影响）。
定位：**轻量 ligand-based CADD 工作台** —— 输入结构 → 2D 查看 → 理化性质 → 3D 构象 →
力场优化/应变评估 → 构象系综 → 导出。基本风格、输入输出、属性计算完整对齐
Molecule Clipboard（下称 MC），以同样设计语言长出 3D 面板。

**交互式 MD 与 MetaD 不进 app**：可交互/可视化的动力学过程与轻量 CADD 流程目的不同，
属于探索型工具——归 Playground（交互 MD、原子拖拽/扭键，需持续打磨升级）和
Demo（MD/MetaD + FES）所有。app 顶部导航互链两者，"深究细节去 Playground"。

职责划分复刻 MC 的 "RDKit.js + JSME" 模式：

| 职责 | 工具 | 说明 |
|---|---|---|
| SMILES 解析、Canonical SMILES、InChI/InChIKey、2D SVG、描述符、Lipinski/Veber、构象间 RMSD | **RDKit.js**（MC 同款官方 minimal 构建） | 与 MC 逐项一致 |
| 2D 绘制/编辑 | **JSME**（MC 同款） | molblock 进出 |
| ETKDG 3D 构象（批量）、MMFF94/MMFF94s + GFN-FF 优化（L-BFGS）、能量逐项分解、构象系综能级排序、导出 | **WebMM WASM**（本仓库引擎） | MC 明确列出的 "No 3D conformer generation" 限制，正是我们补上的 |

隐私主张与 MC 完全一致：全本地计算、无上传、可离线（第三方库全部 vendor，不用 CDN）。

## 设计系统（直接采用 MC 的 token）

```css
--bg:#f8fafc; --panel:#fff; --border:#e2e8f0; --text:#0f172a;
--muted:#64748b; --accent:#2563eb; --good:#166534; --bad:#b91c1c;
```
- 复用 MC 组件骨架：`.container`、`.panel(.grid)`、`.action-buttons`、`#error` 条、
  `.modal`（header/body/close）、`.drop` 拖放高亮、`.result-badge`（pass=good/fail=bad）、
  `#historyModal`、`sourceModal`、footer `version + copyright`。
- 按钮：`.primary`（主操作）、`.action-btn`（次级）+ `.copy-btn`，同 hover/disabled 语义。
- 徽章行：`<div class="badges">` —— "230/230 match RDKit"、"GFN-FF exact vs xtb"、
  "runs locally"。

## 页面结构（保持 MC 的 input→output 阅读顺序）

```
┌ header: 标题 + badges + [Playground] [Demo] 互链（深究动力学去那边）──┐
│ Row1  #input 面板(SMILES/MOL/SDF textarea + drop + Draw/Edit(JSME弹窗) │
│       + Clear/History)          │  "2D Structure" 面板(RDKit.js SVG;    │
│                                 │  即 MC 原版 "Structure" 仅改名;       │
│                                 │  Download SVG/PNG·Copy SVG·Share link │
│                                 │ 按钮组原样保留)                        │
│ Row2  "3D Structure" 面板(3Dmol,同款 <strong> 标题命名) ——             │
│       Embed 3D · Optimize[MMFF94s|MMFF94|GFN-FF] ·                     │
│       Conformers[N ▸ 批量构象 + 能量排序 + RMSD 剪枝] ·                 │
│       Export(SDF 多构象/XYZ/PNG)                                        │
│       └ 子条: 构象能量表(rank, E(MMFF)/E(GFN-FF), ΔE, 选中即显示)       │
│ Row3  #props "Properties" 表: RDKit 段(MW/ExactMW/cLogP/TPSA/HBD/HBA/  │
│       RotB + Lipinski/Veber 徽章 + InChI/InChIKey, 逐行=MC, 含 "More   │
│       Properties" <details> 折叠段)                                     │
│       + WebMM 段(E(MMFF)/E(GFN-FF) + 逐项分解 bond/angle/torsion/rep/  │
│       es/disp/hb/xb/batm, 收敛信息, 耗时)                               │
└ footer: 隐私声明 · version · copyright(同 MC 版式) ───────────────────┘
```

## 里程碑（每个带验收与门禁）

### M0 — 骨架 + RDKit 对齐基线（"MC 克隆"）
`app/index.html`（+vendored rdkit.js/3Dmol/JSME），design tokens、三面板静态布局；
input→2D→props 链路跑通（SMILES/MOL/SDF、拖放、JSME 弹窗、描述符表、Lipinski/Veber、
InChI/InChIKey、错误条）。
验收：对 MC 逐项核对描述符数值一致（caffeine/ibuprofen/ascorbic acid 三例）；
面板命名 = MC 原版仅 "Structure"→"2D Structure"，Download SVG/PNG·Copy SVG·
Share link 按钮组原样可用；布局截图 1440/820px 两档；零 console error。
门禁：cargo test 231/231、clippy 0（无 Rust 改动）。

### M1 — 3D 面板：Embed + Optimize + 能量分解 + 导出
Embed 3D（ETKDG）→ 3Dmol 渲染；Optimize 按所选引擎跑 L-BFGS，几何更新；
`#props` WebMM 段展示总能量 + 逐项分解表（新 Rust 导出，见下）；
Export SDF/XYZ/PNG 下载。
验收：caffeine MMFF 能量与 demo 页一致（-123.49 kcal/mol 量级）；GFN-FF 能量与
cargo test 基准一致；导出文件可被 RDKit.js 回读。
门禁：CDP headless 脚本（load→input→embed→optimize→export 全链路，零 console error）。

### M2 — 构象系综（替代原"交互 MD"：CADD 核心需求）
批量 ETKDG（N 默认 50，上限 200，进度条）→ 每构象双力场优化 → 能量排序表
（E(MMFF)/E(GFN-FF)/ΔE 相对最低）→ RMSD 剪枝（RDKit.js 或 JS 内 Kabsch+重原子，
阈值可调）→ 点击行切换 3D 视图 → 导出多构象 SDF（能量的 SDF 属性字段写入）。
实现注意：`generate_initial_coordinates_wasm` 目前硬编码 `random_seed: 42`
（每次同构象）——需新增带 seed 的批量导出（见 Rust 侧）。
验收：ibuprofen 50 构象：能量跨度合理、RMSD 剪枝后簇数 < N、多构象 SDF 可被
RDKit.js 读回且构象数一致；批量进度可取消。
门禁：CDP 全链路 + 零 console error；cargo test 231/231、clippy 0。

### M3 — 状态、历史与打磨
URL 分享（query 编码 SMILES+engine+N，MC parity）；history 弹窗（localStorage，
同 MC 语义）；键盘可达性、响应式、a11y、隐私/source 弹窗、README 更新。
验收：URL 载入恢复完整状态；history 增删可用；lighthouse a11y ≥ 95；
CDP 全链路 + 截图归档。

## Rust 侧小幅补充（M1/M2，不破坏任何现有门禁）

1. `OptimizationResult` 增加 `energy_components: Vec<f64>`（EnergyComponents 九项），
   供 props 分解表；GFN-FF 路径已有 EnergyComponents，透出即可。（M1）
2. 批量构象导出：`generate_conformers_wasm(sdf, n, seed_base)` —— 内部循环 ETKDG
   （seed = seed_base+i，配置已支持 `random_seed`，现 wasm 入口硬编码 42 需解开），
   避免从 JS 调 N 次、每次重解析 SDF。（M2）
3. （可选，M3）qa 电荷导出 + 3D 电荷着色，放 v1.1。

## 部署

同 repo GitHub Pages：`app/` 静态目录，wasm 产物按 site 现有 staging 流程进 `pkg/`；
RDKit.js/JSME/3Dmol 全部 vendor 进 `app/vendor/`（MC 单文件哲学 → 我们"全 vendor 可
离线"，不用 CDN，与隐私声明一致）。

## Out of scope（v1）

- **交互式 MD 与 MetaD/FES**（归 Playground/Demo；app 只保留互链）；
- 大规模构象农场（Workers 并行，>500 构象级，v1.5）；
- WebGPU（Phase 2，触发条件见前期讨论）；
- 多分子/批量模式、账号/后端、PBC/金属体系、docking。

## 与 MC 的功能对照表

| MC 功能 | app/ 处理 |
|---|---|
| SMILES/MOL/SDF 输入、拖放 | 原样复刻 |
| Canonical SMILES / InChI / InChIKey | 原样复刻（RDKit.js） |
| 2D SVG 渲染（"Structure"→"2D Structure"） | 原样复刻，仅改名 |
| MW/cLogP/TPSA/HBD/HBA/RotB + Lipinski/Veber | 原样复刻（同阈值同配色） |
| JSME 编辑弹窗 | 原样复刻 |
| URL 分享 / 历史 | 复刻并扩展（engine/N 入 state） |
| 单文件、本地、无上传 | 同主张；vendor 依赖，可离线 |
| ——（MC 限制项）—— | **3D 构象与批量系综、双力场优化、能量分解、导出**（WebMM 新增） |
| ——（我们主动不做）—— | 交互 MD、MetaD/FES → Playground/Demo |
