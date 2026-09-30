# Plan: Search 标签 E2E 固化入库——tests/cdp/m5_search_e2e.test.js

## 背景

上轮全功能审查(Search 重点)用的是 /tmp 一次性脚本,审完即弃。用户问
"能 e2e 测试吗"——能,且应固化为可重复跑的入库套件。tests/cdp 已有
m0-m4 惯例(独立 node 脚本、:8901 repo-root 伺服、硬编码 Chromium 路径、
RESULT 行汇总),Search 标签(相似性/子结构/shape/RGD/骨架/药效团 +
库持久化)是唯一无 E2E 覆盖的大块——正好补上。

## 实施

1. 新增 `tests/cdp/m5_search_e2e.test.js`,沿用 m0-m4 惯例;冻结参考直接
   读 `tests/fixtures/lbdd/search_refs.json`(2 sim 查询 × 5 指纹 × 55
   全精度)与 `rgd_refs.json`(6 命中名+片段)
2. 覆盖(53 项审查的 Search 子集,断言修正两处脚本伪差的教训):
   - demo 库 55 装载 → reload 自动恢复 → 末尾 Clear(含 localStorage)
   - 相似性奇偶 550 值 === 精确;UI 路径自命中居首/降序/阈值
   - 子结构 4 查询命中集精确
   - shape:自命中 100% 居首、Color/Pharm/Combo 列、Pharm ≥60% 恰余
     {aspirin, salicylic acid}、行点击回载自动 embed+Features
   - RGD:Example 6/55 命中名+片段全精确;Auto core 值断言
   - 骨架:计数和=34、行点击回载
   - 药效团:默认 4 勾选+距离矩阵、±1.5/N=1 命中恰 {aspirin,
     salicylic acid}、Conf=1、行点击回载
   - 全程零 page error;exit code 反映成败
3. tests/cdp/README.md 补 m5 行(含运行命令)
4. 纯新增测试文件 + README,不碰 app/src(无产品代码改动)

## 验收

- `node m5_search_e2e.test.js` 全绿(本地 :8901 伺服)
- m0-m4 不回归(抽查 m4 + m3)
- README 与套件一致

## 验收结果(实施后)

- `node tests/cdp/m5_search_e2e.test.js` **26/26 一次通过**:库装载/
  reload 恢复/Clear;相似性奇偶 550 值 === 精确;UI 路径自命中居首/
  降序/阈值;子结构 4 查询命中集精确;shape 自命中 100%+三列+
  Pharm ≥60% 恰余 {aspirin, salicylic acid}+行回载自动 embed/
  Features;RGD 6/55 命中名+片段精确、Auto core 值断言;骨架计数和
  =34+行回载;药效团默认 4 勾选+距离矩阵、±1.5/N=1 命中恰
  {aspirin, salicylic acid}、Conf=1、行回载;零 page error
- 回归:m4 7/7;m3 9/10(axe nested-interactive 存量,文献在案)
- README 运行命令与套件清单已含 m5;零产品代码改动
