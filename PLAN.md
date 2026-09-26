# Plan: xtb setup-qloop 精确移植(ferricyanide pibo 修复)→ 目标 v1.3.1

## 背景

ferricyanide 链路已映射:我们 setup EEQ(键数-CN + Floyd,对称)
vs xtb gfnff_charges(逐配体几何响应)→ Hückel 对角 → pibo
(0.971 vs 0.959)→ kb → bond −0.019。naive 重解触发反馈环破坏
21 测试。需精确移植 xtb gfnff_ini 的 setup-qloop 语义。

## 任务

1. **源码定位**:找到 gxtb/xtb 的 gfnff_ini.f90(或其重写源)——
   setup 的电荷迭代循环(qloop):迭代次数、每轮的求解参数
   (一程/二程)、CN 来源(键数 vs erf)、距离矩阵(Floyd vs 几何)、
   最终 qa 落点。
2. **逐行对照我们的 setup 流**(Gfnff::new 605-810 区):标出每处
   语义差异;解释 Fe +0.32 vs +0.19 与剩余 runtime ~1e-3 差。
3. **实施移植**(保持 32/32 有机锁 4e-6 不回归 + 金属锁全绿为
   硬门禁;nicarbonyl 必须保持 0.000000)
4. **验收**:ferricyanide 逐项对拍(bond/es/total 收敛至参考)、
   全门禁、grad_audit
5. **发布**:1.3.0→1.3.1(若移植成功);失败则文档化精确差异清单

## 验收(实施后实测记录)

**结论:qloop 语义对照完成,setup 侧无分歧;残差定位到 Hückel
重试矩阵的二进制-vs-源码漂移,不可再分(无 5 月构建树),文档化
关闭。零行为变更。**

- **qloop 对照**(xtb 源 gfnff_ini.f90 400-705 vs 我们 Gfnff::new):
  - xtb 最终 topo%qa 也是 Floyd 拓扑求解(qloop 2 轮,末轮
    goedeckera(rtmp))——**与我们一致**;gfnff_topo 重启文件解析
    证实:xtb topo qa Fe +0.3177/C −0.0137/N −0.5393(对称),
    我们 +0.3208/−0.0141/−0.5394(Δ≤3e-3)——此前"setup 电荷
    分歧"假设被推翻(gfnff_charges 是 runtime q,非 setup qa)
  - gfnff_charges = chk%nlist%q(calculator.f90 315)= 逐几何求
    解的 runtime 电荷——其不对称性与 setup 无关
- **残差真身**(verbose 键表):0.959 是 **C≡N** pibo(Fe-C
  piBO=0,Fe 不入 π 系;6 个 C≡N 是独立 2 轨道 π 系统)
- **逐数值追踪**:每系统 ipis=−1 → nel=3 → HOMO>0.4 → nel−1=2
  重试 → 重试矩阵 B:我们 0.6631(Pold=0.4855,收敛定点)vs
  二进制 0.595(隐含 Pold≈0.39)→ pibo 0.971 vs 0.959
- **我们与 9 月源码逐行一致**(occu 奇电子 na=⌈nel/2⌉ 同、
  fermismear 4000K kT=0.3447 eV 同、迭代/断点同);二进制为
  5 月构建(30c6303),参考树 9 月浅克隆(fa035bc)无该提交
  ——Pold 轨迹差异疑似版本漂移(候选:htriple 1.45 vs 注释
  "1.4"、Pold 初始化),无源树不可裁决
- **旁证**:nicarbonyl(电中性,C≡N nel=2 偶)精确 0.000000
  ——偶 nel 路径两实现逐位同;有机 32/32 同理;唯奇 nel
  (荷电 CN⁻ 络合物)触雷——解释了为何只有 ferricyanide 残差
- 残差量级:total −0.018 Eh(0.28%),测试容差 5e-2 内;按
  "不可裁决上游分歧"惯例文档化,拒绝单夹具魔法常数
