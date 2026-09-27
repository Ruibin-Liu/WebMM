# Plan: Demo 步骤 5 按钮分行 + 能量图/FES 坐标与释义 + 三页 3Dmol 视图重置按钮

## 背景

用户三项反馈:
1. Demo 步骤 5 的第一个实验按钮与 'Run MetaD' 没有分行(7 个按钮
   挤在一行 flex-wrap 里,首行混排不清晰);
2. MD/MetaD 下方的能量曲线无坐标(此前"能量图坐标轴"曾立项延期,
   现用户点名要);MetaD 的 FES 图"意义不明确"(无轴、无解释);
3. 三个 3Dmol 框(Demo/Playground/Workbench)都缺视图重置按钮,
   且要求**在框内**。

零引擎/WASM/Rust 改动;site/index.html、site/playground.html、
app/index.html。

## 任务

1. **按钮分行(site/index.html)**:实验按钮移入独立容器
   `.step-actions.experiments`(复用 flex-wrap,margin-top 6px),
   'Run MetaD' 独占一行;按钮 id/类/接线不变。
2. **能量图坐标(drawChart,site/index.html)**:
   - 左边距 40px 画 PE 纵轴刻度(min/mid/max 三条水平网格线 +
     kcal/mol 数值);有温度数据时右边距 36px 画 T 刻度(红色系);
   - 底部 14px 画时间轴三刻度(0/中点/末点,ps);
   - 折线/游标映射改用新绘图区;图例行保持。
3. **FES 坐标与释义(drawFES,site/index.html)**:
   - signature 增 cvType 参数;绘图区 pad L34/R10/T16/B16;
   - 纵轴:ΔF(相对全局最小,0 起)三刻度;横轴:5 刻度,二面角
     显示度数(±180/±90/0)、距离显示 Å;
   - marks 虚线/标签适配新绘图区;
   - #legend-fes 释义文案:`ΔF along CV (lower = favored) ·
     F = −γ/(γ−1)·Σhills · N hills ●`(讲清"低=更有利"与构造式)。
4. **视图重置按钮(三页)**:
   - 各 viewer 容器内右下角 `.view-reset`(30px,⟲,title/aria-label
     "Reset view",z-index 高于 canvas;三页 CSS 同构);
   - 点击 = `viewer.zoomTo(); viewer.render();`(Demo/Playground
     直接调;Workbench 以 sdf3d 存在为前置,无模型时 no-op);
     轨迹播放不受影响(updateViewerFrame 本就不动相机,重置只影响
     当前视图)。

不做:Playground 能量图坐标(用户所指为 Demo 的 MD/MetaD 图;
Playground 图另行跟进);引擎/布局级改动。

## 验收(实施后实测记录)

- Playwright 12/12 + 截图目检:
- 分行:Run MetaD 底边 1438 < 首实验按钮顶边 1444;
- 坐标:能量图左轴 kcal/mol 刻度 + 网格线 + 时间轴(ps,metad 无
  温度时右轴正确缺席);FES 图 ΔF 纵轴(0/中/顶)+ 度数横轴
  (−180..180),实验 marks 正常;#legend-fes 释义
  "ΔF along CV — lower = favored (F = −γ/(γ−1)·Σ hills) · N hills ●";
- 重置按钮:三页均在框内(Demo zoomTo spy 命中;Playground 旧框外
  按钮已删;Workbench 全链路 CCO→Embed 3D→in-frame→点击零错误);
- 实施中排掉的两个关联坑:①Workbench initViewer/clear3DViewer 会
  重写 viewer3d 容器(createViewer 与占位文本恢复都清子节点)——
  按钮改为清空后重新 appendChild(clear3DViewer 须在 innerHTML
  重写**前**取引用);②**既有 bug 顺修**:3Dmol canvas 的 intrinsic
  宽度把 grid/flex 列撑到 753px("390 直接打开+跑实验"流程从未测过,
  以前 1440→390 缩窗触发 resize 自愈)——三页容器加 min-width: 0,
  fresh-390 流程实测 390/390;
- 回归:实验 evidence(+0.78/verify)不变、零 page error;
- `cargo test` 281/281、clippy 0、fmt 干净(零 Rust 改动)。
