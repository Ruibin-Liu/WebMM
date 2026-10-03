# Plan: M0 — LBDD 平台重设计的设计定稿轮(explore/lbdd-platform)

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
