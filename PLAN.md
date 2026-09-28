# Plan: 3D 形状相似性(ROCS 谱系高斯形状对齐)——引擎新模块 + Search 集成

## 背景

LBDD 第二层主打差异化:浏览器内 3D 形状检索。纯静态无服务器场景下,竞品
(LigandScout/ROCS/docking 前过滤)均不可得;自家引擎有 ETKDG 构象、MMFF、
L-BFGS——形状对齐是天然拼图。用户确认:不做内置大库,3D 形状主打。

## 参考语义与常数(已探明)

- **参考家族**:Grant & Pickup 1996 原子高斯形状(shape-it 参数化谱系,与
  ROCS 同族);score = 形状 Tanimoto = V_AB/(V_A+V_B−V_AB)
- **原子高斯**:g(r)=GCI·exp(−α|r−r_i|²),**GCI=2√2**(shape-it config.h),
  α 为逐元素表(=κ/r_Bondi²,κ≈2.418;H 1.679158/C 0.836674/N 1.006447/
  O 1.046567/…,取自 shape-it GAlpha;冷门元素回退 κ/r_Bondi²)
- **解析锚点(免费单测)**:单原子体积 = GCI·(π/α)^{3/2} ≡ 4/3·π·r³
  (硬球体积,常数自洽,如 H: 7.238 Å³)
- **乘积高斯**:α_c=Σα,center 加权,C=ΠC·exp(−Σ_{i<j}αiαj/α_c·d_ij²),
  V=C·(π/α_c)^{3/2};分子体积 = 包含-排斥展开(按 C 衰减剪枝)
- **RDKit ShapeTanimotoDist 是硬球栅格模型**(EncodeShape 球占用 +
  grid Tanimoto,默认 0.5 Å 间距)——与我们的高斯是**同族不同模型**:
  交叉验证为相关性级(Spearman + 偏差界),非逐位奇偶;如实文档化
  (先例:GFN-FF vs xtb 残差惯例)

## 设计决策

1. **优化目标 = 成对高斯重叠代理**(Σ_ij V(g_i·h_j),解析梯度,6-DOF:
   旋转矢量+平移,复用引擎 L-BFGS m=20);**报告分数 = 全包含-排斥体积
   Tanimoto**(逐起点全量重打分选优)。代理与全量的 argmax 轻微失配由
   多起点 + 全量选优兜底;shape-it 自身用四元数梯度爬升 + 全量目标,
   全量梯度化为后续优化项(记录,不做)
2. **多起点**:主轴对齐 4 起点(惯性主轴正负组合)+ N 随机旋转
   (默认 8+4,种子可复现);每起点 L-BFGS 收敛即停
3. **氢原子计入体积**(shape-it 谱系;RDKit 栅格默认忽略 H——相关性
   脚本两口径都测)

## 任务

### T1 引擎模块 `src/shape/mod.rs`

- ShapeAtom{center, alpha};α 表(H..Br 常用 + Bondi 回退);GCI
- self_volume(全包含-排斥,递归乘积枚举 + C 剪枝)、overlap_full(A,B
  交叉项)、overlap_pairwise + 6-DOF 解析梯度(ShapeObjective impl
  optimizer::Objective)
- align(query, target, opts{starts, seed, max_iter}) → {tanimoto, transform,
  iterations}
- 单测:单原子/双原子重合/远离解析恒等式;FD 梯度 proptest(种子锁定);
  自对齐随机姿态 ≥0.999;Tanimoto 对称性/值域

### T2 WASM 导出(additive)

- `shape_align_wasm(query_sdf, target_sdf, opts_json) -> JSON`(d.ts 同步)
- 版本 1.3.2 → **1.4.0**(引擎功能新增);wasm 重建 + node 冒烟

### T3 交叉验证 `scripts/shape_correlation.py`

- ~20 对分子固定姿态:Rust 侧高斯 Tanimoto(经 native 测试例程导出)
  vs RDKit ShapeTanimotoDist(栅格,含/不含 H 两口径)——Spearman、
  偏差界,写入 CODE_STATUS(相关性级,非奇偶)

### T4 Search 集成「Shape (3D)」模式

- 模式下拉增 Shape (3D);查询:Use current(有 3D 用当前 3D SDF,否则
  引擎单构象 ETKDG+MMFF94s 现场生成);库:has_coords 的 SDF 直用,SMILES
  条目逐条生成 1 构象(分块异步 + 进度状态)
- 结果列:Shape Tanimoto(降序);行点击回载**对齐后**目标坐标
  (transform 施回 molblock,复用固定列手术工具)
- 性能预算:对齐 ≤10ms/对(wasm)→ 2000 库全扫 <30s + 3D 预制 ~1-2min
  分块进度

### T5 验收

1. Rust 单测全绿(解析恒等式 + FD + 自对齐);clippy 0/fmt
2. node 冒烟:版本 1.4.0;benzene↔benzene 对齐 Tanimoto≈1
3. Playwright(仓库根布局服务,引擎路径可达):Shape 模式全流程
   (aspirin 查询 → 自身 1.0 居首/水杨酸次席)、SMILES 库 3D 预制进度、
   行回载对齐坐标、390px、零 page error
4. 既有回归:LBDD parity 51/51、Search 奇偶、281→N 测试

## 边界与不做

- 不做:颜色力场(color Tanimoto/combo)、全量目标的解析梯度(shape-it
  式,后续)、构象系综逐构象对齐(每分子 1 构象)、对齐叠合 3D 可视化
- 不引入依赖;不复制 shape-it 代码(方程与参数实现 clean-room,引用
  Grant & Pickup 1996;常数表为物理参数)
