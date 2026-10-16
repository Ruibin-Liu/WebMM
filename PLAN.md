# Plan: M3 worker 池规范完形(持久池 + 中毒回收 + 堆监控 + 重派)

## 规范条款 → 实施决策
- ✅ 无状态小作业幂等/队列重派:现有 score 作业本就无状态;新增
  pool.submit 队列 + 中毒重派(每作业 1 次重试)。
- ✅ 取消 ≠ 丢堆:池进程级持久(首次运行惰性 spawn,跨运行复用),
  取消/完成只停消费不 terminate —— 规范核心痛点("丢 300MB 堆")根治。
- ✅ 异常即中毒回收:worker 侧区分 Result 错误(化学失败,done+failed
  照旧)与 JS 异常/wasm trap(poison 上报);主线程 terminate+重生+重派。
- ⚠️ HEAP 监控:协议与阈值回收路径已就位(worker 应答 'heap'、主线程
  25 作业一探、>256MB 回收)——但 pkg/webmm.js 未导出 wasm.memory
  (闭包内部),worker 读得 0,监控暂空转。不为一个 getter 重编 WASM:
  引擎将来导出 memory 即自动激活,测试已注明。
- ✅ 池规模:规范 min(核数−1,4)(现为 8,收敛到规范)+ deviceMemory
  ≤2GB 时减半。
- ❌ transferable:结果是小对象/字符串(分数+SDF 串),结构化克隆成本
  可忽略;规范该条针对大坐标数组场景,不适用,记录不做的理由。
- ❌ 构建哈希戳:仓库无 vendored 哈希基建,记录不做。

## 改动
1. app/analog.worker.js:ready 附带 wasmHeap;'heap' 查询应答;done 回
   __job 路由戳;queue catch 改报 poison(化学失败早退不变)。
2. app/platform.html:新增 analogPool(持久池:惰性 spawn、就绪屏障含
  失败上报与重生上限、submit 队列、空闲分发、poison 回收+重派、heap
  轮询回收、stats 调试口);runAnalogFromCands 改用池 —— 删除
  spawn/ensure/waiters/killAll,B/C/D 陈旧分支只停消费不 terminate;
  init 失败状态行语义保留(现测试契约)。
3. 范围:仅 analog 农场(规范条文明列 M3 生成式);shape/conf 农场
  维持现状并注明。

## 验收(新增 tests/cdp/m6_analog_pool.test.js)
- 连续两次真实探索均 Done,且第二次池零重生(复用实证,stats 断言)。
- 人造毒作业(molblock 非字符串 → worker JS 异常)触发 poison:重派一次
  后按 failed 结算、worker 重生计数 +1,随后真实运行仍 Done。
- heap 上报为有限数;取消后再跑不重建池。
- 回归:m6_analog_cancel 两场景不变绿(init 失败秒报 + 重跑吃掉旧运行)、
  m6_platform 39/39、run.sh 全量 13/13。
