# Plan: M1a — 隐形店(身份 + IndexedDB + 命令日志,旧 UI 经新 store 写)(explore/lbdd-platform)

> M0 已完成(见 git 历史);本 PLAN 为 M1a 实施轮。

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
