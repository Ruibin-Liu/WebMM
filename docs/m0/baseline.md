# M0 基线会话 — 当前 app(平台化前)的可证伪基线

> 复跑:`python3 -m http.server 8901 &` + `node scripts/baseline/m0_baseline.js`
> 机器相关性:墙钟绝对值随机器浮动;**动作数是机器无关主指标**。
> 平台化后(M1c/M2)同脚本复测,对照 §成功判据。

## 数字(2026-10-03,Chromium headless,Apple Silicon)

| 会话 | 动作数 | 墙钟 (s) | 到短清单 (s) | 短清单 |
|---|---|---|---|---|
| A 命中发现 | **11** | 3.7 | 3.7 | 10 |
| B 批量甄别+SAR | **8** | 0.9 | 0.3(批量段) | 10 |
| C Lead hopping 一轮 | **11** | 9.0 | 9.0 | top5 重叠 3 |

## 会话与观察

**A 命中发现**(库装载→相似性 0.5 阈→shape 3D→pharm≥60% 过滤):
11 个动作里有 4 个纯模式切换/阈值重跑(onSearchModeChange+runSearch
成对出现 ×2)——正是 B4(模式割裂)的直接计量。短清单 10 条
(paracetamol/nicotine/fragments 级)。

**B 批量甄别+SAR**(批量 10 分子→QED 排序→Lipinski 过滤→行检视
→切 Search 装库):批量段 0.3s 到短清单(2D 描述符很快,真实痛点
不在速度在**跨标签搬运**——行检视后要再切 Search 重装库才能做
SAR,B1 的计量)。注:三态排序首击为升序,真实用户到降序需两击,
基线按一击记录(确定性优先)。

**C Lead hopping**(shape 检索→取命中#2→作新查询→再检索):
**一轮 hop = 11 个动作**——含行点击回载、切模式、清空查询框、
手输 SMILES(库内无 hit→SMILES 快捷映射,真实用户要另开工具查)。
B5 断点的量化:平台化目标 = 单轮 hop ≤ 3 动作(hit-as-query)。
top5 重叠 3/5(ibuprofen/naproxen/anthracene 复现,naproxen 视角
引入 nicotine/naphthalene)。

## 成功判据(M1c/M2 出口复测)

- A:动作数 11 → ≤ 6(递进下钻免模式切换);墙钟不劣化
- B:跨标签动作消除后 ≤ 5;批量段墙钟持平
- C:单轮 hop 11 → ≤ 3;墙钟 ≤ 6s
- 新增能力不得使任何会话动作数上升(零学习成本下限的计量形式)
