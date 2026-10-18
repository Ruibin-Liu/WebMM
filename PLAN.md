# Plan: P1 — 项目注册库 + 每项目一库 + 概览管理 UI(docs/project-management.md 收敛稿)

## 范围
- app/platform/db.js:createIdbStorage(dbName) 参数化(默认 'webmm-platform'
  零行为变化);新增注册库 createRegistryStorage('webmm-projects')。
- 新增 app/platform/registry.js:状态机(creating/ready/repointing/trashed/
  deleting + 对账 orphan/missing)、list/create/rename/trash/setBytes、
  'registry' 锁与广播、零拷贝迁移(registerLegacy)、databases() 对账。
- app/platform/compat.js:init({projectId}) —— lock 名映射
  (default→'project',其余 'project:<id>')、频道按项目分域。
- app/platform.html:URL ?project=<id> 每标签活跃项目;概览页管理 UI
  (列表/新建/重命名/复制/切换/软删除/导入=新建/导出 v1 格式/用量"约"/
  persist 引导);导出升级 {format:'webmm-project', version:1, facts(layer
  补标), sourceProjectId, logHash(canonical 序列化,V3-1 修正:结构排序 +
  volatile 排除 + toString 数字、禁 toFixed)};导入接受 v1 老格式
  (commands-only)+ 新格式,导入=新建项目(importedFrom 谱系,同 id 同
  hash 提示已存在)。
- facts 层标记:saveRoundScores 写 layer:'stochastic';导出过滤
  layer∈{fact,stochastic}。

## 不做(P2/P3)
- analog run 持久化、指针重指(恢复到当前)完整流程(P1 先退位破坏性
  替换为"导入=新建",重指协议 P2 随 analog 一起)、草稿板迁移、逐出。

## 安全
- 迁移幂等零拷贝;老路径(无 ?project)打开 default,15 套件零改动为
  验收线;锁名 default 沿用 'project'(新老代码互斥)。

## 验收
- Node 新增:registry 状态机/迁移幂等/导入去重谱系/logHash 规范化稳定
  性(含 volatile 排除)/导出含 facts。
- m6:概览页出现项目列表与操作;新建项目切换后存储隔离(库互不可见);
  导入文件→新项目;默认项目迁移后原 55 库无损。
- run.sh 全量 15/15 复跑绿。
