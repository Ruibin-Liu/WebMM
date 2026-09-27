# Plan: Demo 页 Metadynamics 经典实验第一批 — 丁烷/甲基环己烷/水杨酸(接入现有 recipe 框架)

## 背景

环己烷翻环示例上线后的候选案例筛选(原生验证,临时例程验收后删除)。
选定第一批三个,加上已有环己烷共四个"一键经典实验";trans-1,2-DMCH
与联苯列为第二档(暂不做);丙氨酸二肽需 2D CV,另立项。

已验证的事实基础(引擎零改动):

- **丁烷**(14 原子,CV=dihedral 0,1,2,3,60k 步,hillW 0.25):
  三阱 FES——anti ±171°(简并基态)、gauche ±66°(+0.9–2.0,真值 0.78)、
  eclipsed 势垒 0°(+5.7)/±126°(+4.1)。优化参考(引擎↔RDKit 逐位):
  anti E=−5.0760、gauche−anti ΔE=**+0.78**(实验 ≈0.9)。gauche 帧优化
  落在 −4.29 = anti+0.78(Use frame 闭环可验证)。
- **甲基环己烷**(21 原子,CV=环二面角 5,4,3,2,100k 步):非对称双阱
  (eq +54.4° 起始侧 / ax −57°)。优化参考:ax−eq ΔE=**+1.37**
  (MMFF94s;实验 A 值 1.74)。**FES 盆差不报数字**:demo 预算内驻留偏置
  噪声主导(100k/200k/deposit25 实测 +3.6/+4.2/+5.0,物理真值 1.37),
  evidence 明示"lobe balance is run-dependent,以优化参考为准"。
- **水杨酸**(16 原子,CV=**distance** 9,2(酚 O···羰基 O),60k 步,
  hillW 0.3):两态 FES——H-bonded 2.75 Å(全局)vs open 4.15 Å
  (+0.85–1.14)。起始几何为氢键构象(O···O=2.66 Å,X-ray ≈2.65);
  优化参考 E=8.8294。展示 distance CV 类型(此前 demo 只演示 dihedral)。
- 三个新 preset 的冻结奇偶参考(RDKit 单点,引擎对拍 Δ=0.0000):
  butane −5.0760 / methylcyclohexane 0.6983(eq 椅)/ salicylic 8.8294。

诚实边界(沿用环己兰先例):FES 盆差半定量(1D 投影 + 有限采样);
定量锚点是优化参考值,全部可经 "scrub → Use frame → Optimize" 闭环验证。

零引擎/WASM/API/Rust 改动;纯 site/index.html。

## 任务

1. **Preset 扩充**:MOLS 增 butane(anti 构象)/ methylcyclohexane
   (eq 椅)/ salicylic(H-bonded 构象)三个 SDF;RDKIT_REF 增对应冻结值;
   preset 顺序:caffeine, ethanol, butane, benzene, cyclohexane,
   methylcyclohexane, salicylic, dmso(.presets 已 flex-wrap,无需改 CSS)。
2. **实验表驱动重构**:新增 `METADEXP` 表(4 条:ringflip/butane/mch/
   salicylic),每条含 label/title/mol/cvType/cvAtoms/steps/hillW/marks/
   gap(或 fesLine)/ref/verify/msg;`runExperiment(key)` 载入分子并钉死
   全部 MD+Metad 表单字段(dt 1.0/T 300/摩擦 1.0/seed 42/步数/hill
   0.3+hillW/每 50 步/γ10)保证复现性;替代原 metadRecipe/isRingFlip
   单例代码。命中判定保持数据驱动:(preset key, cvType, cvAtoms 规范化)。
3. **UI**:步骤 5 的单按钮换为 4 个 step-btn link 实验按钮
   (cyclohexane: chair flip / butane: anti/gauche / methylcyclohexane:
   eq/ax / salicylic acid: H-bond),title 由 METADEXP 注入;init 解禁、
   disableActions 纳管(.btn-exp 类)。claim 文案改为四实验通用引导。
4. **evidence/status**:命中实验时——有 gap 的报盆差行(半定量标注);
   MCH 报 fesLine(run-dependent 声明)+ 优化参考;全部追加 reference 行
   与 verify 行;状态栏消息按实验定制。drawFES 的 marks 机制不变
   (距离 CV 的 mark 单位为 Å,复用同一绘制路径)。

不做:引擎/WASM/Rust 改动;trans-1,2-DMCH、联苯(第二档);丙氨酸
二肽(需 2D CV,另立项);通用坐标轴刻度(既有延期项)。

## 验收(实施后实测记录)

Playwright 23/23 全绿,要点:

- 四实验各一键跑通,字段钉死正确;evidence 各含 reference 行
  (+5.93 / +0.78 / +1.37 / 2.66 Å)与 verify 行;状态栏定制消息;
- MCH 与水杨酸无 'FES basins' 数字盆差行(实施中实测:水杨酸浏览器
  运行盆差符号翻转为 open −3.1——两盆近简并 + 驻留偏置,与 MCH 同治,
  改报 run-dependent 声明 + 几何锚点);
- 丁烷 C 闭环逐位命中:gauche 帧(CV=−40°)→ Use frame → Optimize →
  E_final = −4.2938 = anti+0.78(与原生一致);
- 新 preset:8 按钮出现;butane 奇偶 MATCH(vs 冻结 −5.0760);
- 环己烷实验回归一致;caffeine 默认 CV 无实验行;
- 水杨酸 FES 截图目检:H-bonded/open 标签 + 虚线标注正常;
- 布局:390px 无溢出、桌面 nav 54;
- `cargo test` 281/281、clippy 0、fmt 干净;临时例程
  (metad_cand/metad_cand2/opt_e/sp_e)已删除。
