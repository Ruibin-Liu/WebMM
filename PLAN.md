# Plan: SA 表空闲预热失效修复(预热早于 WASM 就绪)

## 根因
- ensureSaTable 的 requestIdleCallback 预热在 window.webmm 尚为 undefined 时
  触发:fetch fpscores.bin 成功后 window.webmm.sa_load_table_wasm(...) 抛
  TypeError → catch 告警并重置缓存。缓存重置使正式路径(analog 4104/4309)
  重试成功,功能无损 —— 但预热白做(5.6MB 白拉),首次 analog 运行仍吃
  ~1.4s 同步解析,且每次冷加载都出一条误导性告警。

## 改动(app/platform.html,单点)
- 预热改为轮询守卫:webmm 未就绪时每 500ms 重试,就绪后调用 ensureSaTable
  (ensureSaTable 自身的失败-重置语义保留不动)。

## 验收
- 新增确定性 CDP:按住引擎脚本让 webmm 保持 undefined 越过 idle 触发窗
  —— 旧代码此时必出告警(复现证据),修复后零告警;释放引擎、webmm
  就绪后预热自动补上(sa_score_wasm 直调返回有限值 = 表已由预热装载,
  正式 SA 路径从未运行)。
- m6 + loop 复跑全绿;常规加载零新增告警。

## 不做
- 不动 ensureSaTable 本体与正式调用点;不暴露新的页面 API。
