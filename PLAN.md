# Plan: Workbench 单分子处理态 390px 溢出修复(grid 内在轨道 + 按钮行无折行)

## 背景

LBDD 收口轮验收发现(存量,HEAD 复现):单分子处理态在 390px 视口下
`scrollWidth = 552px`。本轮 Playwright 精确定位根因(注入式 CSS A/B 实测):

- **驱动元素 = 2D 导出按钮行**(884 行 `.panel.grid` 第一个 grid 子项内的
  `<div style="display:flex; gap:0.75rem; justify-content:center;">`,行内
  SVG/PNG/MOL/SDF + Copy 系列 + Link 按钮无折行):其 min-content ≈ 511px
- grid 子项默认 `min-width: auto` → 511px 作为轨道基底尺寸传播 →
  `.grid` 单列轨道(620/900px 断点后 1fr)被撑到 511.328px → 文档 552px
- 隐藏该行实测轨道塌回 308px、文档 390px,确证唯一驱动源
- 与 v1.3.2 已修的 3Dmol intrinsic 宽度撑破 grid 同类(grid/flex 子项
  auto 最小尺寸);`.action-buttons` 类本就带 `flex-wrap: wrap`(3D 导出行
  用的正是它,故无恙)——2D 行是漏网的内联样式

## 修复(最小外科,纯 CSS/标记)

1. `.grid` 定义后追加 `.grid > div { min-width: 0; }`(仓库既有模式,
   两个 grid 子项:2D 列 + 属性列)
2. 888 行内联 flex 行改用 `class="action-buttons"`(该类 = display:flex +
   gap:0.75rem + flex-wrap:wrap + justify-content:center,语义完全等价 + 补上
   缺失的 wrap)——不引入新 CSS 规则

注入式实测(C 方案):轨道 511.328px → 308px,scrollW 552 → 390,按钮行
高 31px(nowrap 溢出)→ 75px(两行整齐折行)。

## 验收(Playwright)

1. 390px:aspirin/rhodanine/cholesterol 处理态 scrollW == 390 == clientW;
   Embed 3D 后仍 390(3D 行复用 .action-buttons 本就带 wrap,顺带回归);
   batch 流(上轮已干净)不回归
2. 桌面 1440:双列布局、面板高度、按钮单行(9 按钮总宽 < 列宽)不变
3. 截图目检 390 处理态(按钮折行整齐、无裁剪/重叠)
4. 零 page error;CSS-only 改动,parity/LBDD 不受影响(不重跑全量,门禁
   sanity:281 cargo test、clippy 0、fmt)

## 边界与不做

- 不动 site/ 两页(v1.3.2 已各自修复且无此结构)
- 不重构其它内联样式;不加新组件
