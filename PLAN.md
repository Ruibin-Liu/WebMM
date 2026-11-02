# Plan: iMTD 收口四项——GFN-FF 引擎 / 镜像构象语义 / 取消语义 / workbench 任务层接线 + UI(已完成)

用户指令「今天这些都做完」= 上轮未做清单四项。远期能量景观视图仍按用户既定指示「仅记录,不排期」。

## 交付
1. **GFN-FF 引擎支持(src/metad/imtd.rs)** ✅
   - `SearchEngine{MMFF94, MMFF94s, GFNFF}`(`from_label` 解析)取代 MMFFVariant 参数;`Quencher` 公开枚举:MMFF 共享实例(Box 尺寸纪律)/ GFNFF{at,charge} **逐次重建**(键感知读几何,与 analog 管线同纪律,数值地板 max_force 0.05/rms 0.005 同款);execute_run 的 metad 包裹力场按引擎构建;ImtdOptions.mmff_variant 兼作引擎选择器(取值 MMFF94s/MMFF94/GFNFF,既有取值语义不变)。
   - 测试:ethanol GFNFF 与 MMFF 面不同、确定性、能量有限。
2. **镜像构象语义(明确并锁定)** ✅ 模块文档「Conformer-equivalence semantics (v1)」:重原子镜像对(真旋转不可叠合)保留;同重原子构架+等能(甲基 H 置换/氢排布镜像)折叠。测试:丁烷密采样得 3 极小(anti+gauche± 等能差 1e-9)。
3. **取消语义(app/tasks/imtd.js)** ✅ `IMTD.start()->{promise, cancel(reason)}`;cancel latch(Promise.race)解围轮内未决的 Promise.all——**原设计在轮内散射未决时 cancel 死等(排雷)**;取消丢弃未完成轮、系综=最后完整轮快照,promise resolve 带 cancelled/cancelReason;runParallel 变糖衣。
4. **workbench 任务层接线 + UI(契约 §9.5 冻结)** ✅
   - docs/task-contract.md §9.5 + 注册表行;contract.js KINDS/冻结表/effectiveInput(钳制)/inputRevisionFor(structure 哈希)。
   - index.html:Spatial 组 btnImtd + imtdRounds/imtdSteps 字段;runImtd() 薄壳(spec 快照/旁表/令牌 min(W,总量)/托盘/epoch);**process() 换分子→cancelImtd('molecule changed')**;结果喂 confEnsemble(全 H 模板逐构象 SDF,复用图/能量面板;**无 ΔE 再优化——极小已全淬灭**);按钮使能挂 showSdfIn3D(sdf3d 守卫,与 btnOpt 同组);taskTrayCancel 挂钩。
5. **测试** ✅ 原生 +2(GFNFF 路径;对映体对保留);CDP 新套件 **m16_imtd_ui** 15 断言(按钮使能/系综 ≥2 排序/面板+逐构象 SDF/任务 complete 含 nConfs+rounds+duration/spec 快照钳制/令牌全释/GFNFF 冒烟/换分子中途取消(m14 扣消息法)旁表 cancelled+原因/迟到消息不复活/零 page error);m0 按钮集合期望 3→4 同步;run.sh 挂 m16。
6. **不做(记录)**:rerun UI(conformers kind 同无,一致)、GFNFF quencher 缓存优化、能量景观视图(远期)。

## 验收(实测)
cargo fmt/clippy 双目标(native --all-targets + wasm32 --lib)0 警告;cargo test 345(lib)+ 7(sascore_golden)全绿;wasm 重编 d.ts **零增量**(GFN-FF 为引擎内部重构,WASM 导出面不变);run.sh 28/29 CDP 全绿(m1 为 xtb 外部二进制缺失环境项)+ m16 15/15 + m15 13/13 + Node 153/153。

## 排雷存档(本轮)
- **python str.replace 全局命中**:contract.js 中 `if (kind === 'conformers')` 在 effectiveInput 与 inputRevisionFor 两处出现——参数校验整块误插进 revision 函数;Node 冒烟(makeSpec+inputRevisionFor)当场暴露修正。脚本锚点断言失败即退出(写盘在末尾)避免了半写状态。
- **取消死锁**:轮内 `await Promise.all` 未决时 cancel() 只置旗无法解围——须 cancel latch 竞速(Promise.race);取消丢弃未完成轮;Worker 扣消息法(m14 同款)实证取消路径可达且迟到消息不复活。
- **按钮使能双路径**:reset3DState 的 for-循环是**禁用**列表,真正解锁在 enable3dActions(解析后)与 showSdfIn3D(3D 存在后);iMTD 属后者(sdf3d 守卫),挂错处按钮永禁(m16 抓出)。
- **m0_core 硬编码 .action-btn.spatial 数量**(3→4):UI 增量牵动核心套件,期望集合同步。
- **clippy 三连**:from_str 混淆 std trait(改 from_label)、Quencher 变体尺寸差(Box MMFF + as_ref())、conv.clone 冗余(闭包参数化)。
- **镜像语义探针**:丁烷默认预算只到 2 极小(anti+单 gauche)是**覆盖不足**非去重折叠——密采样(steps 2500/quench 125)得 3;语义测试锁定等能对差 1e-9。
- **手动起测试服务器 cwd 必须仓库根**:tests/cdp 下起会 404→webmm 永不 ready(120s 超时假阴性)。
