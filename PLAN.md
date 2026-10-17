# Plan: aza-scan 双向化(aza C→N + 可选 deaza N→CH)

## 设计
- aza 卡加复选框「含反向 deaza(N→CH)」(id=azaDeaza,默认关)—— 按
  "默认单向、勾选反向"决策;组合数翻倍故不默认全开。
- deaza 枚举:SMARTS [n;r5,r6;!H1](吡啶型芳环氮;[nH] 型由 is_valid
  兜底过滤)在 m0 上匹配(与 aza 同款索引口径),元素列手术 N→C
  ('C   '),重原子 surgeryBase 上隐式氢由价态自动补;单+双取代、
  规范 SMILES 去重、is_valid 过滤全部复用 addCand。
- 计数与契约:__azaCounts 保持 {n, single, double}(m6 契约:Done 行
  /single \+ \d+ double/),新增 deazaSingle/deazaDouble 仅在勾选时
  上报并追加到状态行 ' · deaza X+Y'。
- 不做:跨方向组合双换(一 C→N + 一 N→C)—— 组合空间再翻倍,无明确
  需求;[nH] 型 deaza(吡咯→环戊二烯非芳香,由过滤自然排除)。
- 卡副题 tooltip 与按钮 title 同步提及双向。

## 改动(app/platform.html 单文件)
1. aza 卡 markup:复选框 + 说明 tooltip。
2. runAzaScan:ringN 匹配(勾选时)、deaza 单/双枚举入同一 cands、
   计数与状态行、早退条件改为两集合均空。

## 验收(m6_analog_ui 追加段)
- 4-甲基吡啶父体:默认跑 → 所有行含 n(纯 aza);勾选 deaza →
  出现无 n 行(甲苯,Cc1ccccc1)、__azaCounts.deazaSingle===1、
  状态行含 'deaza';零 page error。
- 回归:m6_platform(aza 段契约不变)+ run.sh 全量。
