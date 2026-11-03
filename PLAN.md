# Plan: 鞍点精搜(CI-NEB)+ 占据玻尔兹曼重构(已完成)

用户指令:实施剩余两项「鞍点精搜、占据玻尔兹曼重构」。观测势垒上界升级为 NEB 精修鞍点能;占据列旁加 Boltzmann 重构布居(能量→平衡布居,与采样占据三角对照)。

## 交付
1. **NEB 引擎(新模块 src/neb.rs,顶层 pub mod,与 md/metad 平级)**
   - `neb(ff, a, b, config) -> NebResult`:Kabsch 对齐 b→a(复用 imtd 坐标级实现,提为 pub(crate));线性插值初始路径;改进切向(Henkelman 能量加权);垂直核力+切向弹簧;**攀岩像**(最高能内像,F 反转切向分量);FIRE 积分器(确定性,无 RNG);收敛=内像最大力 < tol;端点固定。
   - 输出:完整路径能量、saddle_energy(最高内像)、barrier(相对低端点)、converged/iterations。语义:TransitionStats.barrier 一直存**绝对穿越能**——观测上界=窗口 PE 峰,NEB 精修=鞍点能,可比。
   - 测试:丁烷 anti↔gauche(鞍点在合理带 2–6 kcal/mol 且 >双端点、≤线性插值峰);乙醇 OH 扭转(小势垒带);端点不动;同构象平凡路径;逐位确定性。
2. **协调器精修相位(src/metad/imtd.rs)**
   - 盆地稳定 ID(u64 单调计数器入 Conformer 旁路向量)——修指标漂移隐患:pair 统计键改盆地 ID(位置可因后续验证折叠变动,ID 不变),Landscape 映射按 ID→最终序。
   - `plan_refine() -> Vec<RefineJob>`:对每条有计数的转移边出一 job{aId,bId,xyzA,xyzB}(坐标取当前 finalize 映射);`gather_refine(aId,bId,saddleJson)`:TransitionStats 增 `refined: Option<f64>`,断连树在有精修值时用之;仅 finished 且无在途批次;批次可重复。
   - JSON:job/结果编解码(f64 字符串纪律);transitions 增 "refined"。
3. **WASM(加法)**:`run_imtd_neb_json(sdf, job_json, ImtdOptions)`(worker 端纯函数;**GFNFF 拒绝并诚实报错**——拓扑逐几何,逐像逐迭代重建不可承受,v1 仅 MMFF);`imtd_boltzmann_populations_json(energies, temperature_k)`;`ImtdSearch.plan_refine_json()/gather_refine(...)`。
4. **玻尔兹曼重构(pub fn,imtd.rs)**
   - `boltzmann_populations(energies, temperature_k) -> Vec<f64>`:Emin 平移 + log-sum-exp 稳定归一;T 钳制合理域。**诚实标注**:found-minima 理想布居,忽略简并度(镜像对已折叠)与熵差/谐振 ZPE。
5. **驱动 + UI**
   - handle.refine():NEB job 散射(worker 池同 verify 模式,共享取消闩);MMFF 引擎前置守卫。
   - 景观头增「势垒精修(NEB)」按钮与温度 mini-input(默认 298.15);盆地表增「Boltzmann%」列(随 T 重算);精修后盆地表下附**转移边表**(a–b | 计数 | 观测上界 | NEB 精修),断连图阈值自动切精修值;元行更新;GFNFF 时按钮诚实禁用提示。
6. **测试**:原生(NEB 物理带/端点/确定性;refine 集成:refined 入档、树用精修、ID 稳定不受后续 verify 折叠影响;boltzmann:和为 1、T→极低归地面态、大能隙无溢出);CDP 新套件 **m19_saddle_boltzmann**(精修按钮:转移表出现且 refined 有限并 ≥双端点能;Boltzmann 列和=100%、升温 flattened(高能盆布居升);换分子清理;零 page error);run.sh 挂 m19。
7. **不做(记录)**:NEB 的 GFN-FF 支持、dimer 法、简并度恢复、谐振自由能校正、精修的任务层记账。

## 验收(实测)
cargo fmt/clippy 双目标 0(loop-index/into_iter/type-complexity 四处清理);cargo test 353(lib,+2:refine 集成与 ID 稳定/boltzmann 稳定性;另 neb 模块 3)+ 7 全绿;wasm 重编 d.ts 导出恰增 imtd_boltzmann_populations_json 与 run_imtd_neb_json 两函数;run.sh 31/32 CDP 全绿(m1 xtb 环境项)+ m19 10/10 + 既有零回归 + Node 153/153;丁烷 NEB 实测:anti↔gauche 带峰 -1.75 kcal(观测上界 15.6/19.3,物理带内),7173 步收敛。

## 排雷存档(本轮)
- **弹簧单位陷阱**:文献 NEB 弹簧 0.1–0.5 **eV/Å²** = 2.3–11.5 kcal/mol/Å²——首版 0.35 kcal/mol/Å² 弱 ~20 倍,像滑下山脊伪造低鞍点(-4.29≈端点能);spring=5.0 后带子压住山脊。
- **dist() 语义 bug**:首版弹簧间距用逐原子距离**求和**(≈n_atoms×集合范数),等效弹簧又强 ~14 倍——先弱后强两头错的根因都是它;改集合位移范数。
- **端点欠淬灭钉死带子**:搜索系综坐标带 250 步淬灭截断的残梯度(实测 ~9 kcal/mol/Å)——固定端点旁的像镜像该应变,力范数永远 ~8.99;neb() 内部先 3000 步紧淬灭两端点(实测 Epath0 不变但残差消除)。
- **QuickMin 过冲振荡**:刚性键模式平衡偏移 ~0.002 Å,dt 0.05/vmax 0.12 的步长 0.006 Å 过冲一个量级→梯度方向翻转→v 清零→永振;dt 0.015/vmax 0.06 后 7173 步收敛、鞍点 -1.75 落物理带。FIRE 在非保守 NEB 力上黏滞弃用。
- **攀岩像 v1 关闭**:CI 爬升在本预算(2 万步)内不收敛(爬升速率 ~0.005 Å/步);v1 报**带峰**=鞍点能上界(随像数下降),climbing 字段保留待后续。
- **JSON 双处发射**:transitions 在顶层与 md 块两处生成,fmt 后文本不同——replace 只中一处,"refined" 字段静默缺失(原生探针抓出);两处都已补。
- 500 kcal 能隙在 f64 下 exp 下溢为精确 0(诚实的零布居),测试断言 ≥0 而非 >0。
- 收敛判据用每像力**范数**(分量最大值不现实地严);tol 0.5 kcal/mol/Å(文献 0.05 eV/Å≈1.2 更松)。
