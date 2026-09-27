# Plan: Playground 能量图坐标补齐(上一轮 scope 外的跟进项)

## 背景

上一轮图表坐标只做了 Demo 页(MD/MetaD 所指),Playground 的实时能量图
(900×100,300 样本滚动窗,PE/T 双线各自归一)记为"另行跟进";用户
现在点名要。

## 任务(site/playground.html)

- timeBuf 缓冲(随 peBuf/tBuf push/splice/clear 三处同步)提供滚动窗
  的 MD 时间跨度;
- drawChart 重写为带轴版本(与 Demo 同构):左 PE 轴(kcal/mol,蓝,
  min/mid/max 三网格线)、右 T 轴(红)、底部时间轴三刻度(ps,
  按窗口实际 time_fs 跨度标注);双线保持各自归一。

## 验收(实施后实测记录)

- Playwright 实跑 12s:图有绘制(11.5k 像素)、零 page error;截图
  目检——左蓝轴 −98/−104/−111、右红轴、网格线、时间轴 0.0/4.5/8.0 ps
  (与 overlay 的 8.7 ps 一致),无裁剪/重叠;
- `cargo test` 281/281、clippy 0、fmt(零 Rust 改动)。
