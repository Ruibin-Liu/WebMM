# Plan: Top-2 —— color 深化 + 排名可解释(收敛路线图第二项)

## 范围
1. **引擎小扩展**:`AlignOptions.color_scale`(默认 1.0)—— 联合姿态优
   化中 color 项整体缩放(双侧站点 w×scale);事后 Tanimoto 保持不变
   (自重叠与交叉重叠同缩,比例不变);`shape_align_color_wasm` opts 解析
   (clamp 0–20);Rust 测试 color_scale=0 逐位还原纯 shape 叠合。
2. **Shape↔Color 滑杆**(0–1,默认 0.5):
   - 拖动即时**固定姿势重排**(状态行明示"固定姿势,未重新叠合"——
     Astra 红线:改权重≠重新叠合);排序用保序形式 score+k·colorT
     (k=(1-w)/w;w=0.5 ⇒ k=1 恰为历史 score+colorT,行序与 m6 lead-hop
     契约不变——首版用归一化混合打破契约被全量回归抓住);列显示用
     归一化混合(0–1);null color 行按 shape 排(无证据≠零分)。
   - w≠0.5 显式"按此权重重新叠合检索"按钮(color_scale 经 worker 的
     colorOpts 通道进联合优化)。
3. **特征贡献 tooltip**:pharmMatch 返回逐查询位点 types;Pharm 列
   tooltip 列 per-feature ✓/✗ + 匹配率(≥0.5 判匹配)。
4. **系综公平性(Astra 盲点)**:Conf 列 tooltip 载明获胜构象号/有效构象
   数/采样优势警示;状态行在有效构象数不均时追加 ⚠ 提示。

## 验收
- cargo test 329(+1 color_scale)、clippy 0、fmt clean;WASM 重编重暂存。
- 新 CDP 套件 m8_combo 9 断言(滑杆/重排标签/纯 shape 混合/25-75 混合
  数值/特征 tooltip/conf tooltip);run.sh 全量 21 套件。
