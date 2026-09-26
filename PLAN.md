# Plan: 准确性与正确性加强 → 目标 v1.3.0

## 背景

性能目标已达(RDKit 持平/反超)。本会话引入大量新解析梯度(MMFF
vdW/torsion/oop/融合键合、ETKDG hb/dihedral/improper/linear)+
优化器变更(dense BFGS/暖启动)——每个改动有针对性 FD 测试,但
缺**全语料系统性审计**。另有已知精度缺口未根因。

## 任务

1. **全语料梯度审计**(最高优先):examples/grad_audit.rs——
   MMFF 230 分子语料 + GFN-FF 32 分子集,每分子 3 个随机扰动
   几何,逐原子中心差分 vs 解析梯度;容差:相对 1e-5(FD 截断
   ~1e-6 起步,按数据定);报告最差分子/原子/项。
   附加不变性检查:平移/旋转能量不变(~1e-10 相对)。
2. **MMFF 既有离群根因**:phosphirane、cyclobutene(benchmark
   2 个 >0.01 偏差)——逐项能量分解 vs RDKit 定位差异项,
   修 typing/参数或定性为上游分歧。
3. **GFN-FF ferricyanide kb/pibo**(+0.4% 键项,Hückel pibo
   0.971 vs 0.959)——定位 pibo 来源差异。
4. **门禁**:发现即修 + 回归锁值;cargo test、benchmark、
   ensemble、clippy、fmt;无性能回归(E+G 微基准抽查)
5. **发布**:1.2.10→1.3.0;CODE_STATUS/PLAN;commit+tag;冒烟
   必须命中引擎输出行

## 验收(实施后实测记录)

- `cargo test` 281/281 全绿(+1 金属梯度 FD 回归锁);**benchmark
  首次 230/230 全过且 0 能量差 > 0.01**(原 228/230 带 2 离群);
  clippy 0;fmt;wasm(node 冒烟 1.3.0,引擎输出行确认);API 零变化
- **新常驻审计工具**(examples/grad_audit.rs):MMFF 230 分子语料 ×
  3 扰动几何 + GFN-FF 夹具 × 2,逐原子 FD vs 解析梯度
  - MMFF:0 失败,最差 2.1e-5(FD 截断噪声内)——v1.2.5–v1.2.9
    的全部梯度工作经全语料验证
  - GFN-FF:初测 6 失败(最差 7.7e-2)→ 修复后 0 失败(1.0e-5)
- **修复 1(GFN-FF)**:rep 梯度路径漏 H-金属 ff=0.85 规则(能量
  路径有)——失败集恰好全是含氢金属配合物(Co 氨、二茂铁、锌氨;
  Ni(CO)₄/铁氰化物无氢故过)。含氢金属配合物的 GFN-FF 优化此前
  受到微妙错误的力
- **修复 2(MMFF 扭转估计遮蔽 bug)**:`let v1 = 0.0`(扭转参数
  V1)遮蔽 Table VI 元素 V 值 → 规则 d/e/f/h 的 (v0·v1).sqrt()
  全为 0 → phosphirane 的 C-P 扭转族(V3=0.3759)被整体丢弃
- **修复 3(MMFF 扭转表查找过度通配)**:阶段内额外的 i/l→0 回退
  替换命中 l 特异行 (0,i=0,j=20,k=30,l=30,v3=−0.5),而 RDKit
  只做阶段精确查键 → cyclobutene 的 H-C(sp2)-C(sp3)-H 扭转多出
  −0.5×2。移除回退后逐阶段精确查找
- **逐项根因法**:SetMMFF*Term(bool) 差分得 RDKit 逐项能量 →
  两离群 100% 是扭转项;GetMMFFTorsionParams 逐四元组对照定位
  参数来源;数值复现(sqrt(2.12·2.40)/6 = 0.3759 精确命中)
- **phosphirane 精确修复(追加)**:根因两层——①手工 C-P 键行被
  截断到 4 位小数(2.7618/1.8331 vs RDKit 全精度 2.7618357/
  1.8330689),r0 顺带污染角 ka 链;②四个 P 角表行(环×2 + 非
  环×2)ka 为 3 位小数表值,而 RDKit 表缺这些行→经验规则全精度
  计算(ka=0 标记走既有"规则 ka + 表 θ0"路径)。修后逐参数
  逐位一致(ka 0.144579/0.142591/0.660651/0.426428 = RDKit 精确),
  **phosphirane 28.39248 = RDKit 28.3925,cyclobutene 9.31478 =
  9.3148,benchmark 最差残差降至 +0.0011(thietane,既有)**
- **ferricyanide kb/pibo(链路完全映射,修复延期)**:xtb 的
  gfnff_charges 逐配体响应几何不对称(Fe-C 1.8617/1.8631);
  我们 setup EEQ(键数 CN + Floyd 拓扑距离)完全对称(Fe +0.32
  vs xtb +0.19)→ Hückel 对角线偏差 → pibo 0.971 vs 0.959 →
  kb → bond −0.019;es +0.001 同源。数值验证:xtb setup qa ≈
  我们 runtime 式求解(二程参数 + erf-CN + 实际距离,0.1908 ≈
  0.1899)。但 setup 重解触发反馈环(topo_q → 键检测 → 二程
  参数 → qa)破坏 21 个测试(nicarbonyl 0.000000→0.0044)——
  正确修复需移植 xtb setup-qloop 的精确语义,独立立项
