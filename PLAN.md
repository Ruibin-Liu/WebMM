# Plan: Demo 页 Metadynamics 经典实验第二批 — trans-1,2-DMCH + 联苯

## 背景

第一批(环己烷/丁烷/甲基环己烷/水杨酸)上线后的第二档两个案例,
接入同一 METADEXP 框架;引擎零改动。丙氨酸二肽(2D CV)仍另立项。

已验证的事实基础(本轮原生复测,临时例程验收后删除):

- **trans-1,2-二甲基环己烷**(24 原子,ee 椅起始,RDKit 单点/引擎
  对拍 6.2046,CV=环二面角 5,4,3,2,100k 步):翻环成功(CV ±89/95°),
  非对称双阱;起始侧 +55.0° = ee。优化参考(引擎↔RDKit 对拍):
  diaxial − diequatorial = **+1.80**(b1 aa 优化 8.0068 − b0 ee 6.2046;
  实验 ≈1.86——2×A 值减 Me–Me gauche 抵消)。FES 盆差不报数字:
  驻留偏置(100k 本种子 camping 在 ee 侧 → ax +3.92;200k 上轮
  camping 在对侧——符号都会翻),与 MCH 同治。
- **联苯**(22 原子,RDKit 单点/引擎对拍 39.3292,CV=环间二面角
  4,5,6,7,100k 步):退化 ±54° 双阱(MMFF94s 优化扭转角 −53.9°;
  气相实验 ≈45°)+ 平面位阻壁(0° +3.3 / ±180° +4.2);90° 肩部
  近乎平坦——MMFF94s 的共轭项弱于教科书 ~1.6 kcal/mol(FF 局限,
  诚实标注)。两阱简并(对映),gap 行无意义。

零引擎/WASM/API/Rust 改动;纯 site/index.html。

## 任务

1. **Preset 扩充**:MOLS 增 dmch(trans-1,2-Dimethylcyclohexane,ee 椅
   几何)/ biphenyl(twist 几何);RDKIT_REF 增 6.2046 / 39.3292;
   preset 顺序:…, cyclohexane, methylcyclohexane, dmch, salicylic,
   biphenyl, dmso(共 10 个,.presets 已 flex-wrap)。
2. **METADEXP 增两条**:
   - dmch:CV 5,4,3,2、100k、hillW 0.2;marks 'ee chair' +0.96 /
     'aa chair' −1.00;gap=null + run-dependent fesLine;ref 行
     (+1.80,实验 1.86,gauche 抵消注记);verify 行(−55° 帧 →
     Optimize → E = ee + 1.80);定制 msg。
   - biphenyl:CV 4,5,6,7、100k、hillW 0.2;marks 'twist' ±0.94 +
     'planar wall' 0;gap=null(简并对映阱)+ fesLine(平面壁来源 +
     MMFF 90° 肩部近零的诚实注记);ref 行(54° MMFF94s vs 实验 45°);
     verify 行(twist 架任意帧 → Optimize → E = 39.33);定制 msg。
3. **claim 文案**:补两实验枚举(ee/aa 平衡、联苯扭转)。
4. 实验按钮增至 6 个(.btn-exp,机制不变;命中判定仍是数据驱动的
   (mol, cvType, cvAtoms))。

不做:引擎/WASM/Rust 改动;丙氨酸二肽 2D CV(另立项);FES 坐标轴
刻度(既有延期项);第二档之外的扩容。

## 验收(实施后实测记录)

Playwright 23/23 全绿,要点:

- 两实验一键跑通,字段钉死(5,4,3,2 / 4,5,6,7、100k、hillW 0.2);
  evidence 含 ref 行(+1.80 / twist minimum)与 verify 行;无 'FES basins'
  数字行;状态栏定制消息;preset 10 按钮;两分子奇偶 MATCH(6.2046 /
  39.3292);butane 回归无污染;390px 无溢出、nav 54;
- biphenyl C 闭环:任意帧 → Optimize → E = 39.3292(±54° 极小,逐位);
- **实施中两个诚实发现并修正文案**:
  1. **引擎 DihedralCV 符号与 RDKit GetDihedralDeg 相反**(同坐标同原子
     序:引擎 −55.0 vs RDKit +55.0)——MCH/DMCH 的 marks 与 CV 提示符号
     全部翻正(batch-1 上线的 MCH 标注同步修正);
  2. **取代环己烷的 ax/aa 椅在轨迹中从未被沉积**:原生逐帧核查
     (100k/200k、种子 42/7/123)远叶帧优化 100% 落扭船架(DMCH
     12.5136 = ee+6.3;MCH 6.9/7.386 = eq+6.1~6.7),ax/aa 椅是 1D 投影
     内的窄子盆地——MCH/DMCH 的远叶标注改为 'twist-boat shelf',
     verify/ref 改为承诺可验证的扭船值 + 注明取代椅"藏在其内"
     (batch-1 MCH 的误导性 verify 提示一并修正);
  浏览器实测一致:DMCH 远叶 9 帧全部优化到 12.514。
- `cargo test` 281/281、clippy 0、fmt 干净;临时例程
  (metad_cand/sp_e/dih_sign/dmch_check/mch_check)已删除。
