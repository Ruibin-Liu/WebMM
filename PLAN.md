# Plan: platform 页无障碍收口(键盘行 + 模态焦点圈 + m6 挂 axe)

## 现状盘点(实测)
- 已达标:search 行/hit 行(tabIndex+role+Enter/Space)、侧栏 role=tab、
  drop role=button、confChart 方向键、模态 Escape/背景关闭。
- 缺口:①RGD 结果行仅 click;②scaffold 行有 tabIndex/role 但无 keydown;
  ③queryHistory 条目为 div[onclick],纯鼠标;④三个模态(jsme/history/source)
  无 Tab 焦点圈、关闭后焦点不归还("focus returns" 注释是未兑现的承诺)。

## 改动(app/platform.html)
1. RGD 行补 tabIndex=0 + role=button + Enter/Space(镜像 hit 行模式)。
2. scaffold 行点击体提取为局部 act(),click 与 keydown 共用。
3. queryHistory 条目 div 加 tabindex=0 + role=button;容器一次性委托
   keydown(Enter/Space → 目标条目 .click())。
4. 模态:共享 focusModal(id)(打开时聚焦首个可聚焦元素、记录外部焦点)
   与 closeModal(id)(隐藏 + 焦点归还);document 级 Tab 圈定在最上层
   打开模态内;Escape/背景路径改走 closeModal;display 逻辑与既有测试
   断言零变化。

## 测试(tests/cdp/m6_platform.test.js 追加)
- axe-core 扫描 0 serious/critical(m3 同款注入)。
- 键盘检查:search 行 tabindex/role;queryHistory 条目 Enter 触发重跑
  (searchQuery 被写入);模态 Tab 圈定(焦点不出模态)+ Escape 归还焦点。

## 验收
- m6 全绿(新检查含);loop/csv/sa_warmup/platform Node 复跑绿;
  零 page error。index.html 不动(m0-m5 不受影响)。
