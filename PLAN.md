# Plan: Top-5 阶段 4 —— 剩余三项收口(已完成)

用户指令「先完成剩余的」,承接上一轮列出的三项。

## 交付
1. **conformers kind(工作台薄壳,设计 D1 / 契约 §9.4 冻结)**
   - input:molblock(重原子,内容即身份)/n≤500/rmsdCut/engine/maxIter/workers 快照/seedBases 数组(42+w·500,快照即真相——W 硬件派生,跨机不可复现故入档供审计);inputRevision=structure 哈希。
   - 工作台侧表独立 IDB `webmm-workbench-tasks`(Platform.DB.createIdbStorage 参数化;两页 DB 所有权契约保持);boot 无条件对账(无 facts⇒abandoned)。
   - runConformers 任务化:提交即入旁表、farm 令牌 min(W,调度器总量)、完成记 nConfs/kept/duration;取消诚实化:process() 换分子→cancelled+原因 'molecule changed'(epoch 守卫留作纵深防御)。
2. **Ketcher 画布高亮(反向桥)**:模态打开→setMolecule→highlightPicked;3D 拾取原子经 editor.highlightAtomById 防御式点亮;状态行注明 molfile 序。
3. **CODE_STATUS 归档**:315 条巨段切前 10,其余 305 条+三处历史重复段+Fixed Issues 迁 docs/CODE_STATUS-archive.md;主文件 1212→46 行,Current Focus/In Progress/Upcoming 刷为真实现状。

## 验收
Node 153(+2);m14_wb_tasks 10/10(含 Worker 子类扣消息法确定性触发换分子取消);m11 扩高亮段 8/8;m2 零回归;run.sh 27 套件(26 CDP + Node)全绿。

## 排雷存档(本轮)
- aspirin 100 构象 <600ms 完成——睡眠竞态必输,测试须扣 worker 消息;
- 新 context 空库无需 deleteDatabase(boot 对账已开库,删除即竞态);
- `picked` 词法全局,window.picked 第三次踩雷。
