# Plan: Demo 页 Metadynamics 经典示例 — cyclohexane 椅式↔船式翻环(A+B+C)

## 背景

外部评审方向:Demo 页第 5 步(metadynamics)缺一个"有意义"的示例。
原生可行性实验(examples/metad_cyclohexane.rs,临时验证例程,验收后删除)
已确认当前引擎零改动即可演示教科书级物理:

- cyclohexane preset(已有)优化后 = 椅式,E = −3.5609 kcal/mol,
  环二面角 C1-C2-C3-C4 = −54.3°(atoms 0,1,2,3,0 基)。
- CV = 该二面角,页面默认参数(100k 步、dt 1 fs、300 K、摩擦 1/ps、
  hill 0.3 kcal/mol、σ 0.2 rad、每 50 步沉积、γ=10):原生 1.2–2.4 s,
  种子 42/7/123 全部多次翻环,CV 扫过 ±86–92°,~25–35% 帧进入船区
  (|s|<0.5 rad);从**未优化的 preset 几何直接跑**(一键路径)同样成立。
- FES 形状三种子一致:±55° 两个椅式深盆 + ±30° 扭船浅盆 + 0° 船式鞍。
- 船区轨迹帧去偏置优化 → 扭船局部极小 E = +2.3688 = chair + 5.93
  kcal/mol(实验值 ≈ 5.5)。
- 诚实边界:FES 盆差数值(0.25–3.3)系统低于 5.93 — 1D 投影熵 +
  100k 步欠收敛,属半定量;页面文案按"半定量 FES + 定量极小值差"措辞,
  不宣称 ΔF 定量。

零引擎/WASM/API/Rust 改动;纯 site/index.html(JS/HTML/CSS)。

## 任务

1. **A — 一键示例按钮(步骤 5)**:
   - step-actions 增 `cyclohexane ring flip` 按钮(step-btn link 样式,
     title 说明:载入 cyclohexane、CV=C1-C2-C3-C4 环二面角、100k 步);
     init() 解禁、disableActions 纳管。
   - 行为:loadMol('cyclohexane') → CV type=dihedral、CV atoms=0,1,2,3、
     md-steps=100000 → runMetad()。
   - 步骤 5 claim 追加一句经典示例引导。
2. **B — FES 标注 + 盆差证据(仅环翻示例命中时)**:
   - 判定(数据驱动,非标志位):currentPresetKey==='cyclohexane' &&
     cvType==='dihedral' && CV atoms 规范化为 '0,1,2,3'。
   - drawFES 增可选 marks 参数:±54.4°/±29.2° 处画虚线竖标 + 'chair'/
     'twist-boat' 文字(示例命中时传入;非示例运行绘图逐字节不变)。
   - evidence 增三行:FES 盆差(chair |s|>0.9 rad 最小值 vs 船区
     |s|<0.6 rad 最小值,半定量标注)、参考值(扭船极小 +5.93 vs 实验
     ≈5.5)、验证指引(scrub→Use frame→Optimize)。成功状态栏消息同步
     定制。不加通用坐标轴刻度(能量图坐标轴已立项为后续设计候选)。
3. **C — "Use frame" 验证闭环(播放条)**:
   - 播放条增 `Use frame` 文本按钮(width auto,默认禁用,
     enable/disablePlayback 纳管,title 说明用途)。
   - 行为:当前帧坐标 buildSdfFromCoords → 写回 sdf-input、currentSdf、
     currentPresetKey=null、renderSdf、markStep('parity', false)(与手输
     SDF 行为一致,奇偶参考失效)→ 状态栏提示去跑 Optimize 看落盆。
   - metad 轨迹坐标是 COM 中心化/去旋转的(平移/旋转不变性保证优化
     能量不受影响)。

不做:引擎/WASM/Rust 改动;通用 FES/能量图坐标轴;Playground/app 页;
图标化按钮(沿用既有文本按钮惯例)。

## 验收(实施后实测记录)

- 一键(Playwright 实测):点击后自动载入 cyclohexane、表单字段
  正确(dihedral / 0,1,2,3 / 100000)、步骤 5 卡标记 done、metad 完整
  跑完(2000 山丘,headless 端到端 ~17s,含 120 帧动画);
- FES 标注:截图目检确认 chair/twist-boat ×2 虚线竖标;
  evidence 含盆差行(该次运行 boat region +2.5 kcal/mol)、
  +5.93 参考行、verify 指引行;状态栏定制消息;
- C 闭环:scrub 到 frame 25(CV=35°)→ Use frame → textarea 更新、
  parity 卡取消 → Optimize → **E_final = 2.3688 kcal/mol,与原生
  扭船参考逐位一致**(65 迭代,0.01s);
- 非示例回归:caffeine 默认 CV(8,2,1,0)跑 metad,evidence 无环翻行、
  纯净证据保持;
- 布局:初版 .playback 被 Use frame 按钮撑到 472px → 补 flex-wrap 后
  390px scrollW==clientWidth==390,桌面 nav 高 54 不变;
- Playwright 断言 21/21;`cargo test` 281/281、clippy 0、fmt 干净
  (零 Rust 改动);examples/metad_cyclohexane.rs 已删除。
