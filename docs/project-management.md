# 项目管理与数据持久化设计(M4)— 三轮评审收敛稿

> 评审:Opus 5.5 三轮(v1 三阻塞 → v2 两阻塞五应修 → v3 Converged)。
> 本文档为冻结契约;标注「冻结」的条款变更需重新评审。

## 0. 背景与痛点

现状实测:命令日志存储(webmm-platform 单例 IDB)已持久化检索库/轮次
DAG/pin·exclude·备注/轮次分数,但:

- 项目是隐式单例,概览页无创建/管理能力;
- 骨架探索结果不保存(仅 pin 入库);
- 分子工作台历史在 localStorage(wb-compoundHistory),与 IDB 割裂;
- **既有缺陷(评审中实测核实)**:`exportProject()` 只导出 commands,
  facts(轮次分数)导出即丢。

原设计(lbdd-platform-redesign §4.3)的持久化四层:L1 事实始终存 /
L2 派生确定重算 / L3 随机结果带谱系存储、版本不一致显式警告 /
L3' 构象缓存可驱逐。多项目当时未 specced,本设计补齐。

## 1. 项目模型(冻结)

**每项目一个 IDB + 注册库。**

- 项目库:`webmm-project-<uuid>`;**注册库** `webmm-projects` 存
  `{id, dbName, name, createdAt, lastOpened, state, approxBytes,
  schemaVersion, meta{createdWithVersions}}`。
- 选型理由(修正后仅三条):删除干净(drop 即释放)、锁/广播天然
  分域、**commands+facts 同库单事务提交**。
- `dbName` 为不透明 uuid,**name 仅显示字段**(IDB 不能改名库,
  重命名只改注册库)。
- **commands+facts 同库同 upgrade 创建是锁定决策**:未来若 facts
  需独立驱逐域,接受退化两步提交 + pending 补偿,不静默。

### 1.1 注册库状态机(P1 首件交付物)

```
creating ──写入完成──▶ ready
ready ──恢复流程──▶ repointing ──指针切换+旧库 trashed──▶ ready(新 dbName)
ready ──删除──▶ trashed ──保留期/配额清理──▶ deleting ──drop──▶ (无)
(对账派生态)orphan:库存在但无注册条目;missing:有条目无库
```

状态全集:`creating / ready / repointing / trashed / deleting / orphan / missing`。
所有迁移只经上表路径,实现不得自造状态。

### 1.2 `indexedDB.databases()` 对账规则(启动时)

| 观察 | 判定 | 处置 |
|---|---|---|
| 库名匹配 `webmm-project-<uuid>` 且不在注册库 | orphan | 不自动删;UI 让用户选"恢复为项目/丢弃" |
| 库名不匹配模式 | 未知 | 忽略,不碰 |
| 注册库有条目但无同名库 | missing | 禁用打开,显示恢复引导 |
| Safari `databases()` stale 返回 | — | 对账仅作提示;以注册库为准 |

### 1.3 迁移(冻结:零拷贝)

现有 `webmm-platform` 单例 → 注册条目
`{id:'default', dbName:'webmm-platform', name:'默认项目', state:'ready'}`。
**不移动数据、幂等**。老测试兼容:无活跃项目参数时打开 default。

### 1.4 锁与广播(冻结)

- **锁名与 DB 无关、origin 全局唯一** —— 显著文档化,勿改为
  `project:default` 以为能分域(会与老代码失去互斥)。
- `lock(projectId) = projectId === 'default' ? 'project' : 'project:'+id`
  —— **默认项目永久沿用旧名**,保证新旧代码标签互斥。
- 注册库操作走 `'registry'` 锁 + 独立广播频道。

### 1.5 活跃项目(每标签)

- URL `?project=<id>`,两页互跳携带;每次打开从注册库解析 dbName
  (不缓存绑定)。
- `lastOpened` 仅作无参默认,节流写。

## 2. 导出/导入(冻结:format.version 1)

### 2.1 导出格式

```json
{ "format": "webmm-project", "version": 1,
  "projectMeta": {…}, "sourceProjectId": "…", "logHash": "…",
  "commands": […], "facts": [仅 layer ∈ {fact, stochastic}],
  "exportedAt": …, "engineVersions": {webmm, rdkit, fragmentSet} }
```

- **facts 分层(冻结)**:每条 facts 记录带
  `layer: 'fact'|'derived'|'stochastic'`;upcaster 为历史记录补标
  (轮次分数 → stochastic);L2 与缓存不导出。
- 已逐出截断的 run 导出 `truncated: true` + 谱系,不导空对象。

### 2.2 logHash 规范化算法(V3-1 修正版,冻结)

- **命令 payload 内浮点一律以十进制字符串或定点整数在写入时编码**;
  logHash 只规范化结构、不规范化数值。
- 结构规范化:递归排序 key;排除字段表(各命令中的时间戳等 volatile
  字段)显式列出;数字按 ECMAScript `Number.prototype.toString()` 的
  结果串,**禁止 toFixed 舍入**。
- hash 覆盖规范化后 commands 数组;导出同时带 rawCommands 与 logHash,
  导入以规范化重算校验。

### 2.3 导入谱系与去重

- 同 `sourceProjectId` 且 logHash 相同 → 提示"已存在,打开它 / 仍作
  为副本";新项目必生成新 id 并记 `importedFrom{sourceProjectId,
  logHash, at}`;"复制"走同路径。
- "恢复到当前项目" = 非破坏指针切换(§3)。

## 3. 指针重指协议(N2,v2 定案)

1. 全流程在 `'registry'` 锁内。
2. 置旧库条目 `state:'repointing'` → 持有该库的**写路径全部拒绝**
   (写闸)→ 批量写入新库(`creating`,分批事务)→ 注册库原子切
   `dbName` 指针 → 广播 `projectRepointed{projectId, oldDb, newDb}` →
   旧库入回收站。
3. 各标签收到广播:close 连接、清投影、模态提示"项目已恢复到新版本,
   点击重载"。
4. 写闸三缝隙(实现清单):
   - 注册库自身轻写(lastOpened)不闸,无害,文档承认;
   - 快照写入复用同一闸函数;
   - **闸读注册库必须直连 IDB,不得读缓存副本**。

## 4. 骨架探索持久化(P2)

- 有随机成分(ETKDG 种子,B/C 阶段),按 L3 实现。
- 命令 `CreateAnalogRun {runId, parentMolId, site, fragCat|aza, deaza,
  fragmentSetVersion, engineVersion, rdkitVersion}`(runId 创建时生成
  入 payload,重放不派生;版本字段显式拆分)。
- 结果 `facts['analog:<runId>']`:仅排行最小字段(SMILES/T/flexT/SA/
  MW/rotb),单 run top-N 截断 + 完整表标记可驱逐。
- 提交时序:计算全程不持锁;完成后在项目锁内**单事务跨 commands+facts**
  提交。**取消语义:commit 前取消 = 零写入;commit 已发起 = 视同成功**。
- undo = 追加 `HideAnalogRun{runId}`;facts 永不删;GC 口径 live/hidden,
  hidden 进 P3 逐出序列,逐出后 UI"已回收,可按谱系重跑(新 run)"。
- Pin 自包含:Pin 命令携带 productCanonical + 生成谱系,不依赖 facts;
  跨 run 归并键 = InChIKey(不可得回退 canonical SMILES)。
- 跨图谱系边:轮次命令 seed 支持引用 `{kind:'analogProduct', runId,
  molId}`;列表分开展示,谱系一张图。
- 过期判断按 run 自带谱系扫描(项目 meta 版本只记"创建于")。

## 5. 草稿板(替代 localStorage)

- 独立全局库 `webmm-scratch`,不属于任何项目(解决容量/同步/竞态/
  导出);一次性迁移 `wb-compoundHistory`。
- "发送到项目"动作生成命令,带 `source:'scratch'` 谱系 —— 痛点③的
  真修法。
- 导入文件中的草稿默认合并进全局 scratch(canonical SMILES 去重,
  带 importedFrom);导出可选附带。

## 6. 概览页(P1 UI)

新建/重命名/复制/切换/删除(软)+ 导入为新项目 + 导出 + 用量展示
(写入自累加字节,UI 标"约";删库后 estimate() 滞后不承诺即时回落)
+ persist() 状态与"上次导出 N 天前"引导。

## 7. GC 与配额(P3)

- 跨项目顺序:非活跃项目确定性缓存 → 各项目归档随机结果(hidden
  run / 归档轮)→ 旧快照;**命令日志与前沿轮永不**。
- 软删除保留期;配额紧时先清回收站。
- `QuotaExceededError`:整事务回滚 + UI 报错,禁止半提交。

## 8. 其它冻结契约

- 重放走主键游标顺序不排序;facts 用 key 前缀区间。
- schema-per-DB 代价:升级 on-open,upcaster 长期保留;注册库缓存
  schemaVersion 供概览统计。
- 注册库至少保留一个项目;全删后进入空态(创建/导入 CTA),老测试
  显式创建项目。
- 多标签命令广播(含 HideAnalogRun)沿用现有 applied → stale 机制,
  不新增协议。

## 9. 分期

| 期 | 内容 |
|---|---|
| P1 | 注册库 + 每项目一库 + 零拷贝迁移 + URL 活跃项目 + 概览管理 UI + 导入=新建 + 导出格式 v1(含 facts/layer/logHash)+ persist 引导 |
| P2 | analog run 持久化(§4 全部) |
| P3 | 跨项目逐出、回收站清理、批量输入入 L1、(可选)草稿板迁移 |

实现 checklist(评审遗留):状态迁移表为首件交付物;写闸三缝隙;
logHash 穷举测试(含 -0/NaN/Infinity 拒绝路径)。
