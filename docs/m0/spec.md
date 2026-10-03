# M0 规范 — LBDD 平台 v0.5 的实施契约(explore/lbdd-platform)

> 本文档是 M1a/b/c 的直接输入。所有决策已定稿;探针结果已并入。
> 上游:`docs/lbdd-platform-redesign.md`(v0.5,四轮评审)。

## 1. 部署支持矩阵(按仓库事实定稿)

| 目标 | Workers | IndexedDB | 支持级别 |
|---|---|---|---|
| GitHub Pages(现役) | ✓ | ✓ | 一级 |
| localhost 静态服务器 | ✓ | ✓ | 一级(开发/企业内网) |
| file:// | ✗(vendored RDKit fetch 失败,既有已知问题) | 不适用 | **不支持**(不变更现状) |
| PWA 安装 | ✓ | ✓(+ITP 豁免) | M1b 选项(须带 §8 SW 守卫) |

共享内存:默认**无**(多 worker 独立实例化,`WebAssembly.Module`
主线程编译一次 postMessage 分发)。CI 矩阵:Chromium + WebKit。

## 2. 身份规范(探针 A 已验证)

### structureKey = InChIKey(链路实测可用)

`key(mol) = rdkitModule.get_inchikey_for_inchi(mol.get_inchi())`

实测层语义(caffeine/aspirin/顺反/盐/互变探针):

| 维度 | 实测 | 决策 |
|---|---|---|
| 立体 | E≠Z 键不同(NSCUHMNNSA vs UHFFFAOYSA) | 立体未定≠已定 ✓ 直接用 |
| 盐 | 游离酸 -N ≠ 钠盐 -M | **默认不去重**(视图选项,后开);无归一化(盐剥离不在 minimal 构建) |
| 互变异构 | 2-吡啶酮两式同键(InChI 移动 H 归一) | 自动归并 ✓ 免费获得 |
| 规范化 | 无(原始输入直算) | **按项目不可变**;改层级=一次身份迁移命令 |

### molId 与行分配器

- molId = `hash(libraryId, recordOrdinal, rawRecord)`(同库重复行不碰撞)
- **行 id:追加式分配器 + 墓碑**。Merge→旧行墓碑化(永不物理删),
  Split→两侧新 id;物理压缩仅发生在快照时且仅限无列引用的行
- 对账命令:`IdentityChanged | IdentityMerged | IdentitySplit`(进
  命令日志);Split 两侧都不继承旧格,溯源指前键;Merge 两份溯源
  留史显新
- 合成库钩子(M3):生成分子 molId =
  `hash(parentMolId, swapSiteIdx, groupId, groupSetVersion, canonical)`

## 3. 存储与单元格

### 单元格(存储五态 + 派生 Stale)

存储:`NotInInput | Scored(v) | BelowThreshold(v) | Failed(reason) |
Pending` + 每格 `computedAgainst: {structureKey, buildHash}`。
**Stale 派生**:`computedAgainst ≠ (行当前键, 当前构建哈希)`。
性质不变式(入 CI):身份命令后,行内所有旧键格必读作 Stale。
Stale 默认排除于排序/融合;无"接受现状"出口。

### 列式布局

- 每 (轮, 分数型) 一 Float32Array + 状态位掩码(Uint8);
  排序/过滤索引惰性物化(未物化列排序 → 列头 "materializing…")
- 2D 缩略:worker 内 RDKit 画 → `ImageBitmap`(可转移)→ 主线程
  canvas LRU;SVG 字符串不进行模型;行高固定(UI 取舍已认可:
  无内联 3D 缩略,长名截断)
- **构建戳混合**:列内禁止,列间允许 + 融合列头警告

### 持久化四层 + 有序自逐出

| 层 | 内容 | 逐出序 |
|---|---|---|
| L1 事实 | 输入/查询/参数/pin/备注/molId/DAG/命令日志 | **永不** |
| L3 前沿 | 前沿轮的 3D 结果(姿态+分数+seed+戳) | **永不** |
| L2 确定缓存 | 指纹/描述符/QED/警示(带版本戳,重算=戳失配回退) | 1st |
| L3′ 可逐出 | 归档轮 3D 结果、构象系综(内容寻址)、旧快照(留 3) | 2nd/3rd |

耐久性话术:`persist()` 调用+结果呈现(被拒且无 PWA=常驻横幅);
`storage.estimate()` 阈值告警;命令日志导出 **opt-in**(用户选目录,
命名 `project-schemaVersion-date`,at-rest 不加密——立场显式声明);
导出在 Web Lock 内;启动 "DB 消失" 三态区分(从未有/已迁移/被逐出)。

## 4. 命令日志

### 命令类型(初版清单,M1a 冻结)

`ImportLibrary{libraryId, version} · ReconcileIdentities{…→
IdentityChanged/Merged/Split} · CreateRound{roundId, parentId, inputRef,
querySpec, params} · SetThreshold{roundId, value} · Pin{molId, note?} ·
Unpin · Exclude{molId, provenanceReason} · Include(override) ·
ArchiveSubtree{roundId} · RenameRound · Snapshot{…}(系统命令)`

**M1a 实施增补**(实施中发现的规范缺口,如实记录):
- `RemoveLibrary{libraryId}` 加入冻结清单(原 11 型漏了库移除)
- 单视图 history 等非项目语义状态 = IndexedDB **键值事实**,不入命令
  日志(命令日志只承载需要撤销/血缘/回放的语义变更);M2 再评估迁入
- `Unpin`/`Include` 的逆命令需**先值捕获**(apply 前存下旧 pin/exclude)

- **甄别规则**:分子级覆盖层 + 溯源 + DAG 前向传播 + 轮内覆写;
  融合轮多亲冲突 = **最新命令胜**(在命令内,不在视图);传播值
  永不存储(每轮派生列)
- **撤销 = 追加逆命令,永不截断日志**。甄别恒可撤销;带依赖结构
  命令 = 拒绝并列依赖 + 一键"归档子树";无级联
- **快照节拍(混合触发)**:≥500 命令 ∨ 尾部重放成本>1s(200k 档
  实测定数)∨ 检查点(导入完成/轮完成/对账完成)。快照 = 列式数组
  本体 + 行分配器 + DAG + 覆盖层;worker 内、Web Lock 守卫、永不
  在交互路径;留最近 3 份
- GC:墓碑行在快照时压缩(无列引用者);L2/L3′ 按上表序自逐出

## 5. schema 版本化

- `schemaVersion` 字段**从 M1a 第一个写入起存在**
- 迁移钩子注册位:`MIGRATIONS: Map<fromVersion, (db) => toVersion>`
  (M1a 只注册 identity hook,不写实际迁移)
- **金标夹具政策**:每个发布版本导出一个真实小项目入库;
  CI 加载投影状态 diff;RDKit 金标输出(canonical SMILES/指纹位/
  种子构象)按 vendored 构建哈希钉死。**金标 #1 于 M1a 出口铸造**
- 一次性迁移:M1a 读现有 localStorage 库(只读迁移,不回写)

## 6. 视图域规则

- 行 = 当前选中 DAG 节点输出集;列 = 祖先路径分数列(列选择器可
  加兄弟分支);**前沿视图**只显叶轮;死分支归档(可逆命令)+
  归档列惰性卸载
- 工作集全局:他分支 pin 折叠 "pinned elsewhere (N)";排除全局
  传播,常驻 "N excluded" 徽标
- 融合列(共识/组合)= 显式命名列,公式进 tooltip,混模式警告
- 规模验收:5k/50k/200k 三档各定加载/排序/滚动预算;**200k 档
  单列 WebKit 内存预算**

## 7. Worker 池(M1c/M3 定稿,M1a 只留接口)

无共享内存;池 = min(硬件并发−1, 4) 受 deviceMemory 调节;
作业无状态幂等可重派;作业尺寸定数:2D 每 200–500 分子、3D 每
10–20;RDKit 异常 = 中毒即回收;HEAP 滞回绝对阈值回收;结果走
transferable;每 worker 实例化后**校验模块哈希**与主线程预期一致
方可计算。

## 8. PWA/SW(M1b 选项,守卫为前置)

SW 缓存按构建哈希键控;启动校验实例化模块哈希;Web Lock 守卫的
schema 版本检查强制滞后标签重载。做不到 → 不上 PWA,只发 ITP
常驻横幅。

## 9. M1a 入口条件核验(Fable 终裁五项)

1. ✅ structureKey 归一化层级定死(§2:InChIKey 原始输入,按项目
   不可变)——本规范即为交付
2. ✅ Stale 派生不存储(§3)
3. ✅ 行分配器墓碑制(§2)
4. ⬜ schemaVersion 首写即有 + 金标 #1 于 M1a 出口(实施项)
5. ⬜ PWA 守卫落地或明示推迟 M1b(决策项,默认推迟)

## 10. 探针记录

- **A1 InChIKey**:`get_inchi`(mol 方法)+ `get_inchikey_for_inchi`
  (模块级,吃字符串)链路可用;三层语义实测如 §2 表 ✓
- **A2 Morgan 展开型**:minimal 构建仅有折叠指纹
  (get_morgan_fp / _as_uint8array);SA score 移植需 Rust 侧计算
  展开型 Morgan 环境 ID(M3 项,已记)
