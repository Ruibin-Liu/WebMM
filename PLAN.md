# Plan: 收官轮——攀岩像两阶段 / C3 置换筒并 / 准谐振 / rerun UI / 任务层记账 / GFNFF 缓存度量(已完成)

用户指令「剩余的都做,然后冻结 feature」。dimer 按已批准建议保持关闭。本轮后仓库进入 feature freeze(仅修复与回归)。

## 交付
1. **筒并统一:变体几何计数(src/metad/imtd.rs,替代镜像特判)**
   - 能量分支帧统一走 `note_variant(i, coords)`:与该盆地已存**变体代表几何**(含成员自身)真旋转全原子 RMSD < tol ⇒ 重访;否则新变体入档、degeneracy+1。镜像对/甲基置换/C3 转子统一为"观测到的不同井数"(≤对称数,诚实);边际镜像判据与 improper_rmsd 降级为测试工具。
   - 测试:乙醇 OH 镜像孪生 g=2(回归);甲醇甲基转子三井 g=3(新);重访不涨。
2. **攀岩像两阶段(src/neb.rs)**
   - Phase1 标准带松弛至力阈(或 60% 预算);Phase2 激活攀岩像继续同判据——此前失败根因是攀岩在未收敛带上启动,预算被带松弛耗尽。默认 climbing=true 恢复;`NebResult` 增 phase 拐点信息。
   - 测试:丁烷 anti↔gauche:CI 收敛、鞍点落在 [带峰−0.1, 带峰+0.8] 且 barrier∈[2,6];景观 refine 路径自动受益(默认配置即 CI)。
3. **准谐振低频处理(src/vib.rs + thermo)**
   - `harmonic_correction` 增 QH 变体:ν̃' = max(ν̃, 100 cm⁻¹)(Grimme 式,缓解软扭转谐振熵高估);thermo 同算 plain 与 QH 双 G,JSON 增 "Gqh"/"qhCut";UI G 列用 QH,tooltip 注明并给 plain 值。
   - 测试:全频 ≥ 截止时与 plain 逐位相等;软模系统 QH G > plain G;T 极限不变式。
4. **任务层记账(契约 §9.5 v1.1,可选 phase 字段)+ rerun UI**
   - contract.js:imtd input 增可选 `phase:'verify'|'refine'`(默认缺省=搜索;verify 记 runs/steps/workers,refine 记 jobs)+ effectiveInput 钳制;docs 契约注记。
   - index.html:runMdVerify/runNebRefine 提交 typed task(令牌已持有,补 spec/旁表/托盘/终态 facts:mdFrames/refinedCount);boot 对账沿用薄壳(无 facts⇒abandoned 自动覆盖)。
   - **rerun**:托盘终态记录(conformers/imtd)增 ↻;重跑用**快照参数**(绝不取当前默认),守卫=快照 molblock 哈希 vs 当前结构,不符诚实拒绝并提示(镜像 optimize kind 语义);conformers 回填 confN/confRmsd/engineSel 后调 runConformers,imtd(含 phase 记录)回填轮/步/T 后重放对应动作。
   - 测试:CDP m21(verify/refine 各留 complete 记录含 facts;托盘 ↻ 重跑产生新 complete 记录;结构变更后 ↻ 拒绝且状态行明示)。
5. **GFNFF quencher 缓存:度量裁决**
   - 基准:乙醇/阿司匹林量级 GfnffForceField::new 计时;若 <100μs/次(淬灭占比可忽略)→ 记录度量与"不做"结论(诚实);否则实现**零风险精确坐标缓存**(同坐标重复淬灭场景:NEB 端点、refine 重算)。结果落入 CODE_STATUS。
6. **冻结准备**:CODE_STATUS 记录 freeze 姿态(后续仅修复/回归/性能,无新面)。

## 验收(实测)
cargo fmt/clippy 双目标 0(死码测试内迁/type alias/冗余克隆六处);cargo test 362(lib:筒并统一 2 新+攀岩断言改带峰语义+QH 1 新)+ 7 全绿;wasm 重编 d.ts 导出集合与上轮一致(本轮无新 wasm 面);run.sh 33/34 CDP 全绿(m1 xtb 环境项)+ m21 11/11;m20 的 g 断言随语义升级同步(≤2→观测井数合理域);GFNFF 度量入档(乙醇 9 原子 535µs/次构建)。**feature freeze 自本提交生效。**

## 排雷存档(本轮)
- **攀岩像两阶段:机制就绪、默认仍关**——QuickMin 在攀岩相存在极限环(实测 maxF 1.66/2.13 交替、鞍点冻在 -1.999 低于真鞍点 -1.68),需阻尼投影 LBFGS 类积分器,越过冻结线;带峰=上界语义保持,两阶段代码与 climb_started 字段保留待后续。
- **筒并统一揭示更深化学**:变体计数自动捕获转子对称数——乙醇实测单盆地观测 3~5 个等能井(甲基 σ3 × OH 双位组合),镜像对只是其特例;边际镜像判据随之退役(降级为测试工具)。
- **甲醇孪生构造过滤器 bug**:`src[i]==base[o1]` 位置比较 intending 跳过 OH 氢,实际只匹配 O 自身——OH 氢被误旋转,三井退化为两井;改邻接 BFS 圈定甲基侧。
- **imtd 记录身份双轨**:搜索记录 molblock=嵌入全 H 3D,而 rerun 守卫比较的是 sourceMolblock——苹果比橙子;§9.5 v1.1 增可选 sourceHash(源结构身份),守卫按 kind 分流(conformers 直比 molblock,imtd 比 sourceHash)。
- **脚本断言中止=未写盘再犯两次**:runMdVerify 头部改动与 lastImtdMolblock 声明分落两次脚本,一次断言中止丢头部、一次丢声明——运行期 ReferenceError 才暴露;多段拼接后必须 node --check + 全符号 grep。
- GFNFF 度量:535µs/构建(9 原子)——精确坐标缓存只惠及同坐标重算(thermo 换 T、NEB 端点);拓扑哈希缓存会重跑它想省的感知,不做。
