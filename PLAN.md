# Plan: 收尾三件(gitignore / og:url / worker 生命周期审计修复)

## 审计结论
- applyDashLayout:仅移动节点(appendChild/insertBefore),从不 clone——
  监听器随节点走,无泄漏;恢复顺序 reverse 正确(后兄弟先归位)。
  结论:干净,不改代码,记录于 CODE_STATUS。
- analog.worker 生命周期:两个真 bug——
  ①init 失败挂死:ensure() 30ms 轮询 __ready 无超时,worker 内 import(wasm)
  失败只回 done 不回 ready → Promise.all 永久 pending,≤8 个 interval 永续
  泄漏 + worker 永不 terminate + 取消按钮无效;
  ②陈旧 token 只 return 不 terminate:B 流式/C Promise.all/D flex 三个分支,
  取消或重跑后 NW 个 worker(各持整份 WASM)存活到页面关闭。
- shape 农场同款②(无 ready 轮询故无①):onmessage 陈旧分支同样不 terminate。
  同族 bug 一并修(审计中发现),两行同款。

## 改动
1. .gitignore 追加 nohup.out。
2. app/platform.html og:url → app/platform.html(当前指向 index.html,
   平台页分享预览落错页;index.html 自身正确不动)。
3. runAnalogFromCands:ensure() 加 30s 超时 + w.onerror → 失败即全 terminate
   + analogStatus 报错;提取 killAll(),B/C/D 陈旧 token 分支统一 terminate。
4. shape 农场 onmessage 陈旧分支补 terminate。

## 验收
- 新增 m6_analog_cancel.test.js:①wasm 动态 fetch 拦截 → 真实 UI 触发探索 →
  状态行报 init 失败(旧代码此处永挂);②运行中取消 → 状态行 settle、
  后续可再跑(terminate 本身无法从页面观测,以不挂死+可重跑为行为断言)。
- m6/loop/csv/sa_warmup/platform Node 全绿;m5(shape 路径)复跑绿。
