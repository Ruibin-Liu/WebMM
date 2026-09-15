# 金属配位复合物 + 带电物种 typing — v1.0 实施计划

> 承接 v0.7.0（两引擎有机体系逐位对拍完成）。前置依赖：构象系综端到端对拍
> （排在本计划之前）。方法论沿用已验证的审计管线：xtb 6.7.1 二进制为裁判，
> `--gfnff --verbose` 逐项能量 + `gfnff_topo` 二进制解析（chieeq/gameeq/alpeeq/qa
> 逐位对比）定位 typing 分歧。

## 范围界定（明确过的研究结论）

| 对象 | 处置 |
|---|---|
| **金属配位复合物**（TM/主族金属中心 + 配体，单分子，含带电配合物） | ✅ 本计划 |
| 金属晶体 / PBC 体系 | ❌ 永久 out of scope（另一层物理；力场工具不碰） |
| 裸金属团簇（Mₙ） | ❌ out of scope（GFN-FF 无 M–M 成键参数；出口是 GFN1/2-xTB） |

- **GFN-FF 侧**：xtb 官方元素覆盖 H–Zn + Br/I，过渡金属有专门扩展参数；
  gxtb 源码全部逻辑在（`gfnff_ini2.F90` hyb/etacoord、`gfnff_ini.f90` 670-716
  EEQ 三阶、`gfnff_param.f90` metal(103)/mchishift=−0.09），是移植工作。
- **MMFF 侧**：MMFF94 元素集 H,B,C,N,O,F,Na,Mg,Al,Si,P,S,Cl,K,Ca,Zn,Br,I，
  金属–配体键无参数 → RDKit `MMFFGetMoleculeProperties` 返回 NULL →
  我们要做的是**失败模式对拍**（同样的"不支持"路径，而非编造参数算出垃圾能量）。

## 现状盘点（2026-09-14）

- ✅ `data/gfnff_params.json` 全部 103 元素表（rad/repan/repz/chi/gam/alp/cnf/
  en/normcn/group/**metal**），metal: Sc–Zn=2(TM), Ga/Sn..=1(主族)
- ❌ `metal_is()` 桩返回 false（"organic subset: no metals"，mod.rs ~2994）
- ❌ `mchishift` 不在 json/Params 中
- ❌ imetal 分支（dgam ff=−0.08/−0.9、ff_alp +0.3/−0.1、TM chieeq 阻尼）未移植
- ❌ 金属键半径 ×2（get_nb f1/f2=fq·2）、nbf/nb/nbm 三邻居表制度、etacoord
  （η²/η⁵）、金属 hyb 规则、M-NC 腈、M-O-X 共轭、Sn/Pb/Bi 低 CN 去金属化
- ✅ rep 侧 H–M 0.85 因子已实现（唯一在位的金属分支）
- ⚠️ 带电物种：MMFF 侧 90 分子对拍集基本中性；阴/阳离子/两性离子未系统覆盖

---

## Task Group 1 — 参数管线与分类（小，先行）

**Files:** `data/gfnff_params.json`, `src/gfnff/mod.rs` (Params)

1. `mchishift: -0.09` 加入 json + `Params`（gfnff_param.f90:768）。
2. `metal_class(z) -> u8`（读 metal 数组 0/1/2）替换 `metal_is` 桩。
3. **Sn/Pb/Bi 去金属化**：邻居数 ≤4 且 group>3 → 按 0 处理（ini 297-300，
   "The number of neighbors can only decrease from first to second qloop"）。
   注意它依赖 qloop 两遍的邻居数单调性 —— 我们的 2 次 detect_bonds 循环同构。

**验收:** 单测：metal_class(30)=Zn→? (查表), Sc–Zn=2, Sn(cn=4)=0 vs Sn(cn=6)=1。

## Task Group 2 — 键检测与邻居表制度（结构性前置）

**Files:** `src/gfnff/mod.rs` (`detect_bonds` ~3094, `assign_hyb` ~4278)

**Problem:** xtb 的 get_nb 分三档（nbf 全量 / nb 常规 / nbm 去金属），金属对
半径阈值 ×2（ini2 100-115: `if param%metal>0: f1=fq*2`），且金属"键"判定
完全不同。我们的 detect_bonds 单表、无金属放大。

**Steps:**
1. detect_bonds 支持金属对阈值 ×2（用 metal_class）。
2. 引入三表语义：先精确弄清 get_nb 三档阈值差异（读 mctc get_nb 源码；
   预计是 rtmp 缩放/截断不同）→ 决定用"全量 nbf + 过滤子集"模拟还是真三表。
3. nbmdiff = |nbf|−|nbm|、nbdiff = |nbf|−|nb| 暴露给 hyb（Group 3 依赖）。

**验收:** Zn(acac)₂ / [Cu(NH₃)₄]²⁺ 的键列表与 xtb `--verbose` bond 表逐条一致
（btyp/pibo/fqq 列）。

## Task Group 3 — hyb 与特殊态（最大块）

**Files:** `src/gfnff/mod.rs` (`assign_hyb`), 新 `etacoord.rs`（如独立更清晰）

按 gfnff_ini2.F90 逐分支移植（行号以审计快照为准）：
1. **金属 hyb**：H 桥连（nb20i=2→1, >2→3, >4→0）；M⁺ 四配位；>4 且 ati>10
   且 nbdiff=0 → hyb=5（hypervalent）；B/Be 同构分支。
2. **etacoord**（ini2 155-205）：C nb20i≥4 且 nbm=3（Cp-η⁵）；C nb20i=3 且
   nbm=2（炔烃 η²）；邻接金属判定（nm=0 则 etacoord 无效；nm=1 时 ncm 区分
   σ-烷基 vs η-配位）；itag=−1 标记；nbdum 切换到 **nbm**（去金属邻居表）供
   hyb 判定 —— 这是 η-配体 sp³/sp² 正确判定的关键。
3. **N 族金属相关规则**：M-NC 腈（idxdum 为金属 → hyb=1）；R-N=C / R-N=N /
   N=N=N；NO₂/B-N/N-SO₂（已实现，复验不回归）；**pyridine-N 配位金属**
   （nbmdiff>0 且 nn=0 → hyb=2, ini2 304）。
4. **O 族**：M-O-X 共轭（nn_nearest_noM CN=3→hyb2 / 4→hyb3, ini2 346-352）。
5. **qloop 第二遍语义**：金属体系下邻居数只减不增的假设依赖两遍循环
   （我们的 `for _iter in 0..2` 已同构，补齐中间状态传递）。

**验收:** 测试集上 xtb 原子表（`atom neighbors erfCN metchar sp-hybrid imet pi
qest`）逐原子一致 —— sp-hybrid 负值语义（amide=−hyb、carbene=−hyb）一并核对。

## Task Group 4 — EEQ 参数金属分支

**Files:** `src/gfnff/mod.rs`（EEQ-param 段 ~526-560, solve_eeq ~3281）

1. dgam ff 覆盖序修正：元素分支之后 `imetal=1 → ff=−0.08`、`imetal=2 →
   ff=−0.9`（ini 694-695，注意它们**覆盖**元素值，包括 >10 的 −0.02）。
2. ff_alp：imetal=1 → +0.3、imetal=2 → −0.1（ini 713-714）。
3. **mchishift**：TM 的估计 chieeq −= 0.09（ini 459，估计阶段；
   最终 chieeq 段不受影响）。
4. dxi TM：group 7 多价 Cl/Br.. 且邻居有金属 → +nn·0.05 替代 −nn·0.021
   （ini 440-446, nm 计数）。
5. rabd（floyd 拓扑距离）在金属键下的 rad 来源核对（param%rad 对金属的值）。

**验收:** gfnff_topo 解析对比：金属测试集 chieeq/gameeq/alpeeq/qa 逐位一致
（方法学已验证 —— 乙酰胺审计即此法）。注意二进制 vs 源码可能再有分歧
（amide −0.16 教训）：以二进制为准，分歧记录进 docs/gfnff-porting-notes.md。

## Task Group 5 — 能量项金属规则

**Files:** `src/gfnff/mod.rs` (energy/energy_and_gradient)

1. 键拉伸：M–X 的 bstren/参考半径路径（egbond M 分支）；M–M 键的行为
   （预期：不生成键项或极弱 —— 用 xtb 逐项能量验证）。
2. rep：金属 CN 进入 fn 因子的路径核对（nrepscal 已实现，验证金属 CN 值）；
   H–M 0.85 已在；桥连 H 的 rep 特例。
3. HB/XB：XH···M 配位键计入 HB 体系的规则（hb setup ini 896+，
   "make list of HB H atoms but only if they have a positive charge" 的金属侧）；
   M 作为 BATM 中心的行为。
4. 扭转：金属参与键的 nrot/φ0/fc 路径（btyp=3 sp-X 无扭转等）。

**验收:** 金属测试集逐项能量（9 项）≤1e-4 Eh 对拍 xtb。

## Task Group 6 — MMFF 失败模式对拍 + 带电物种集

**Files:** `src/mmff/mod.rs`（typing 失败路径）, 对拍 harness

1. **金属复合物**：RDKit `MMFFGetMoleculeProperties` → NULL 的体系集合
   （Fe/Co/Ni/Cu/Zn 配合物数例），我们的 MMFF 入口必须走同样的
   "typing not supported" 错误路径（信息可不同，行为必须一致：无能量、
   无编造参数）。核对现有行为，必要时加守卫。
2. **带电物种扩展集**（MMFF + GFN-FF 双侧）：乙酸根、甘氨酸两性离子、
   铵根、硝酸根、硫酸根二负、[Fe(CN)₆]³⁻；并入 90 分子套件扩成常驻回归。

**验收:** 扩展集 MMFF 全对拍（预期多数直接过 —— 电荷增量路径已有）；
金属体系双侧行为 = RDKit/xtb 的"不支持"或精确能量。

## Task Group 7 — 测试集、夹具与 CDP

1. 金属测试集构造（≥10 分子）：Zn/Ni/Cu/Co/Fe 的乙酰丙酮盐、六氨/四氨水
   合配合物、[Fe(CN)₆]³⁻、二茂铁（η⁵ 压力测试）、金属卟啉片段（大体系压轴）。
   几何来源：xtb `--opt` 或 RDKit ETKDG+MMFF（金属侧 ETKDG 不支持 → 用 xtb）。
2. fixtures 入 `tests/fixtures/gfnff/metals/`；cargo 回归按家族锁值。
3. CDP：Workbench 输入金属 SDF → GFN-FF 单点/优化链路冒烟（若 MMFF 侧
   不支持则 UI 错误提示正确显示）。
4. docs/gfnff-porting-notes.md 增补金属章节（分歧记录处）。

---

## 里程碑与门禁

| 阶段 | 内容 | 门禁 |
|---|---|---|
| MG-1 | Group 1-2（分类+键检测） | xtb bond 表逐条一致；既有 32/32 不回归 |
| MG-2 | Group 3-4（hyb+EEQ） | 原子表逐原子一致 + gfnff_topo 逐位一致 |
| MG-3 | Group 5（能量项） | 金属集逐项 ≤1e-4 Eh |
| MG-4 | Group 6-7（MMFF 失败模式+集成为回归） | 全套件（90+扩展+金属）常驻绿；CDP 全绿 |

## 风险

- **nbf/nb/nbm 三表语义**是最大不确定点（Group 2）：若阈值差异复杂，
  可能需要真移植三表；预算先行 spike。
- **二进制 vs 源码分歧**（amide −0.16 前科）：每条金属规则落地时优先用
  gfnff_topo 验证二进制实际行为，不轻信源码注释。
- **金属测试集几何**：ETKDG 不支持金属 → 依赖 xtb 优化几何，夹具生成
  流程要写清（复现性）。
- 带电 GFN-FF 侧 qfrag 分配（>2 fragment 电荷放置）此前是近似实现，
  金属配合物若多 fragment 可能触发 —— Group 4 一并核对（ini 570-588
  的 try-both 逻辑我们已实现，验证金属集）。
