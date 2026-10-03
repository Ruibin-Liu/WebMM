# Plan: M3b — 漏斗 C 级(扩构象精化)+ D 级(flex 精修 top-10)(explore/lbdd-platform,M3 漏斗收口)

> M3a 已完成;本 PLAN 为 M3 第二轮。SA score Rust 移植独立立项
> (需 Rust 侧展开型 Morgan + fpscores 数据表 + Python 对拍,~2-3 天,
> 不塞本轮;合理性启发式维持 RotB+ΔMW 列)。

## 实施

1. worker 增 returnBest:回传最优构象 SDF + 变换(向后兼容)
2. C 级:B 完成后 top-30 以 30 构象重评(seed 偏移 +10000),两轮
   系综取 max,更新 bestSdf
3. D 级:flexAnalogCandidate——确定性重跑胜者对齐取姿态变换 →
   MCS 配对(≥6)→ Kabsch 核心入位 → 位置约束(fc25/tol0.5)
   MMFF 松弛 → 重打分;Flex T 独立列(与刚体并陈,升降都诚实
   显示);自叠合 ≥0.9995 跳过守卫沿用
4. 状态行叙述漏斗阶段(B done→C refining→D→Done with N flexed)

## 验收结果

- paracetamol 全类别:**35 计分 + top-10 flex 全成,10s 总墙钟**
  (规范预算 45-60s 内);flex 化学有信息量:羧酸 83.8→86.9、
  母体 83.5→88.6、硫代氨基甲酸 82.3→88.0 向上(扭转适应),
  Br 90.5→87.1 向下(MMFF 松弛偏离最优刚体姿态)——两列并陈
- m5 **65/65**(+1 漏斗全通断言);六套件 37/10/11/10/32/65;
  平台 34/34;零 page error
- 教训:测试文件里 python 转义把 \d 写成字面反斜杠-d(正则永不
  匹配→超时);修批量的服务器 kill 误杀后续批次

> M0/M1/M2 已完成;本 PLAN 为 M3 第一轮(M3a)。M3b(flex 精修 top-K、
> C 精化扩构象、SA score Rust 移植)另轮。

## 范围

单取代点基团替换探索:位点检测(6 类末端基 SMARTS)→ 核心 molblock
手术(R# 哑原子封口)→ 内置片段库(44 个七类)经 canonical SMILES
拼装 + 化学合理性过滤(价键/重原子窗口/转子≤10/≤400 产物)→
worker 农场(analog.worker.js,逐候选构象系综 + 逐构象刚性叠合,
流式渲染)→ 结果表(shape T/MW/RotB/★pin=合成库身份)。硬顶全落:
单位点/单类别≤44(全选)/产品≤400/运行令牌取消。

## 实施记录(五个真 bug 的排雷链)

1. 位点原子 0 基(get_substruct_matches 先例)——我误减 1
2. 附着原子被组内 H 抢走——附着必须重原子
3. 键级被手术硬编码为 1——芳香环全饱和、羰基消失、假手性出现;
   保留原键级(含 type-4 芳香)
4. get_molblock() 重写丢芳香性——手术直吃 sdf3d 原文
5. 移除组不带走其 H→断连 [HH] 组分→引擎构象入口 trap;组 H 随组走
   + [H] 字面量/空括号清理(引擎构象入口从未被喂过含氢 molblock,
   全管线皆无氢输入——入库知识)
6. Pin/Exclude 校验放宽:合成分子与孤儿 id 是合法甄别对象(规范:
   甄别挂 molId,库籍贯非前提);空 molId 仍拒

## 验收结果

- paracetamol 全类别:35/44 计分 4s(8-worker);烷基类 7/8 计分 1s;
  化学全对(羟乙酰胺 85.0% 居首、乙基/环丙基/丙基同系物序合理、
  母体自换 83.3%)
- 类似物 pin ★(合成库身份)→工作集计数/撤销/持久化全通
- m5 **64/64**(+4);六套件 37/10/11/10/32/64;平台 34/34;
  390px;零 page error

> M0/M1/M2a 已完成;本 PLAN 为 M2 第二轮(M2b,收口)。

## 范围

项目 JSON 导出(命令日志本体)/导入(校验 schemaVersion→锁内清库
重写→显式 reload);Web Locks 序列化追加(cmdSeq 取自店非本地计数,
双标签发散内存态不再碰撞命令 id)+ BroadcastChannel 通知只读标签
(stale 标志);血缘条目带命中计数。

## 实施记录

1. db.js:withLock(navigator.locks 回退直执行)+ createStoreChannel
2. project.js:applyInternal 锁内取店 cmdSeq 定 id(双标签安全)+
   广播 applied;replaceProject(校验/清 commands+snapshots+facts/
   顺序重写/meta 项目戳)
3. compat.js:activeStorage/activeLock/projectModule 记录,import/
   exportProjectData 门面,channel.onmessage→needsReinit+stale()
4. app:⬇project/⬆import 按钮(confirm 后 reload);血缘条目 (N) 计数
5. 测试教训:Playwright 新 context=隔离存储(广播/IDB 不跨)——
   多标签测试须同 context,且默认 context 不许 newPage→m5 启动改
   显式 newContext;段落 pin 污染后续段→节末清理

## 验收结果

- 浏览器实测:导出 3 命令(ImportLibrary/CreateRound/Pin)→清库→
  导入→reload 后 **55 库+1 轮次+pin+备注全部回来**;双标签:A 写
  →B(同 context)stale=true 广播达;A/B 各写 log 无 id 碰撞
- m5 **60/60**(+3);六套件 37/10/11/10/32/60;平台 33/33;
  零 page error

## 范围

四件套:轮次分数入 L3 事实(facts 库,DB v2)→ 共识=全轮平均秩
显式融合列(公式进 tooltip)→ 溯源 CSV(逐轮 rank+consensus+甄别
态)→ RGD/骨架以命中集为输入(scope 下拉)。

## 实施记录

1. db.js v2 增 facts 库;compat 事实门面(内存缓存+持久化带构建
   溯源);project 暴露 storagePut/Get
2. recordSearchRound 记分(rows→hits 过滤 molId/score);consensusFor
   = 全部已记轮平均秩(<2 轮显 — );Cons. 列插 Combo 后(行索引
   ≤8 不扰动)
3. exportResultsCSV:name/smiles/MW + 逐轮 rank:sim:query 列 +
   consensus + pinned/note/excluded
4. sarScopeEntries + rgdScope/scaffoldScope 下拉;runRGD/
   analyzeScaffolds 主循环与分母文案切换
5. 实施三教训:多步补丁整批失败时步骤间必须逐段落盘(两次全批
   回退);ternary 中间不可插语句(行级手术修复);m5 断言勿写死
   数值(consensus 聚合全历史轮、RGD 命中数取决于取代模式化学)

## 验收结果

- 浏览器实测:两轮重叠查询(aspirin/salicylic)→ cons 1.5/1.5
  公式 tooltip 正确;CSV 表头 rank:sim:×2+consensus+甄别四列;
  scaffold(hits) 19 scaffolds/39 环状;RGD(hits) 6/50 分母正确
- m5 **57/57**(+4);六套件 37/10/11/10/32/57;平台 33/33;
  零 page error

## 范围

每次检索 = CreateRound 命令入 DAG(spec 从 DOM 读,含 mode/query/
threshold/confs/flex);行级 ⇄ hit-as-query 一次动作建子轮(基线 C
的 11 动作 → 1);血缘面板(缩进树,点击重执行,reload 持久);
批量表/pharm 面板不入轮(文档化)。分数列为 L3 事实持久化留给 M2
溯源列。

## 实施记录

1. recordSearchRound(四个检索完成点接线:sim/sub/shape×2)
2. hitAsQuery(molId):设父边+查询+runSearch;⇄ 图标挂甄别列
   (stopPropagation 防行点击)
3. renderQueryHistory:深度缩进血缘,↳ 子轮,当前轮 ◂ 标记,
   末 12 条;rerunRound 保 DAG 位置重执行
4. 库装载/恢复后渲染血缘(与 renderWorkSet 同点)

## 验收结果

- 浏览器实测:round1 入史+⇄ 在位;⇄ 一次点击 → ↳ 子轮 + 查询
  自动设为命中;**reload 后 DAG 持久**(命令日志);rerun 重执行
- m5 **53/53**(+4);六套件 37/10/11/10/32/53 全绿;平台 33/33;
  零 page error
- **基线 C' 复测**:同会话 6 动作(旧 11),hop 本身 = **1 动作**,
  墙钟 8.6→8.3s,top5 重叠不变(3/5)——成功判据 C 达成 ✓
- 实施教训:探针方法错误一次(chromium.launch 每次全新临时 profile,
  IndexedDB 不跨浏览器进程——持久化只能同浏览器内 reload 验证)
> 在表格固化语义前先验证甄别语义(v0.5 §5)。

## 范围

检索结果表三态操作列(末列,不扰动既有列索引)+ 工作集面板
(计数/固定列表+备注编辑/unhide-excluded/undo)+ compat 甄别门面;
批量表不在此轮(文档化)。轮次/DAG 传播留 M1c。

## 实施记录

1. compat.js 甄别门面:pin/unpin/exclude/include/undo/getPins/
   getExcludes/projectApi;saveLibraryInputs 返回富化条目(molId
   附着进内存库)
2. app:★/✕ 末列(☆→★ 翻转、排除行 0.45 透明)、工作集面板
   (pinned elsewhere 孤儿计数)、renderSearchResultsCurrent 重渲染
3. **undo/redo 双栈修复**:初版 undo 的逆命令又压回撤销栈 →
   undo/redo 永久乒乓栈走不下去;正解 = undo 应用逆命令不压栈入
   redo 栈、新用户命令清空 redo(命令日志仍追加逆命令,replay
   不变);拒绝性 undo 不消费栈顶
4. **耐久性竞态修复(真回归)**:旧 localStorage 写同步、IDB 异步
   ——demo 装载后立刻 reload 丢库;修 = 装载/清除状态行等持久化
   落地再报(loadDemoLibrary/loadSearchLibrary/clearSearchLibrary
   转 async + await persisted)

## 验收结果

- 平台 Node 测试 **33/33**(+9:门面 6 含 undo 栈下行不乒乓/redo
  LIFO 重钉带最后备注/拒绝不消费/先值捕获往返)
- 浏览器实测:pin→★+面板 1 pinned、exclude→0.45 透明+计数、
  undo→"Undid: Exclude" 行恢复、**reload 后 pin 存活**
- m5 **49/49**(+5);六套件 37/10/11/10/32/49 全绿;零 page error

## 范围

app/platform/ 五模块 + 旧 UI 库路径改造 + Node 测试 + 金标 #1。
零用户可见变化(持久化 bug 在表格出现前暴露)。

## 实施

1. identity.js(hash64/molId/合成库钩子/结构键注入/墓碑分配器)
2. state.js(命令集含 RemoveLibrary 增补 + 纯 reducer + 逆命令 +
   先值捕获 + 依赖拒绝)
3. db.js(IDB/mem 双适配器 + 持久性话术 + DB 三态消失检测 +
   迁移钩子注册位)
4. project.js(init=快照+尾重放/apply/undo=逆命令/混合快照节拍/
   export=facts)
5. compat.js(旧 UI 桥:saveLibraryInputs/restoreInputs(legacy 只读
   迁移)/clearLibrary)
6. app/index.html:五 script 标签 + saveSearchLibrary/
   restoreSearchLibrary/clearSearchLibrary 三处经 store 写
7. tests/platform/run.js:单元+性质(200 种子命令流不变式)+金标 #1

## 验收

- Node 24+ 平台测试全绿;浏览器内 IDB 全生命周期(装载 55→持久化
  →**删 legacy key 后 reload 仍恢复**→clear 后不复活);六套件
  m0-m5 全绿;零 page error;node --check 全部

## 验收结果(实施后)

- Node 平台测试 **24/24**(身份 5/reducer 6/性质 3[200 种子重放
  等价、逆命令往返、非法命令不抛]/project 4/compat 3/金标 3)
- 浏览器实测:persisted 55 → 删 wb-searchLibrary 后 reload 恢复
  55(证 IDB 为持久层)→ clear+reload 0;persist() 未授予正确告警;
  零 page error
- 六套件全绿 37/10/11/10/32/44;金标 m1a-project.json 已铸
  (3 命令投影断言);spec 增补三条(RemoveLibrary/history=KV/
  先值捕获)

## 范围

纯设计 + 两个探针 + 基线录制,零产品代码。产出 = M1a 可直接实施的
规范与基线数字。

## 任务

1. **探针 A:身份/评分可用性**(页面内 evaluate)
   - vendored RDKit-minimal 是否有 InChI/InChIKey(get_inchi /
     get_inchikey_for_inchi 或等价导出)→ structureKey 主选或回退定案
   - 展开型 Morgan 环境 ID 是否可获取(SA score 移植的 M3 前置)
2. **基线会话录制**(当前 main UI,Playwright 脚本入库可复跑):
   三个脚本化真实会话,记录动作数/墙钟/到短清单耗时——v0.5 §8
   的可证伪基线,UI 变更前最后窗口
   - A 命中发现:库装载→相似性→下钻 shape→药效团过滤→短清单
   - B 甄别与 SAR:批量处理→排序过滤→行检视→RGD/骨架→短清单
   - C 迭代(lead hopping):shape 检索→取命中作新查询→再检索→短清单
3. **M0 规范文档** docs/m0/spec.md(探针结果并入):
   - 部署支持矩阵(§0 事实确认)
   - 身份规范:molId/structureKey(层级定死)/行分配器墓碑/对账命令/
     合成库钩子
   - 存储与单元格:五态+computedAgainst 派生 Stale/整集存分/四层
     持久化/有序自逐出/列式布局
   - 命令日志:命令类型清单/快照混合节拍/撤销=逆命令/带依赖拒绝+
     归档子树/GC 语义
   - schema 版本化:首写即有/迁移钩子注册位/金标夹具政策
   - 视图域规则:行=节点输出/列=祖先路径/全局工作集/排除传播/DAG 修剪
   - M1a 五项入口条件核验清单
4. 基线报告 docs/m0/baseline.md(数字+复跑说明)

## 验收

- 两探针有结论;三基线会话有数字与可复跑脚本;spec.md 覆盖上述全部
  小节且无 TODO 悬空;M1a 入口条件逐条可核
