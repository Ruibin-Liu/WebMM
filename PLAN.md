# Plan: 景观 v2——纯 MD 盆地验证(无偏占据)+ 力导向网络图(已完成)

用户指令「继续」= 按 CODE_STATUS v2 队列:设计页 v2 包(纯 MD 盆地模式 + 网络图)。队列其余(鞍点精搜、占据玻尔兹曼重构、rerun UI、GFNFF quencher 缓存)仍记录不排期。

## 交付
1. **引擎:纯 MD 验证阶段(src/metad/imtd.rs)**
   - `plan_verify(n)->Vec<RunSpec>`:hill_height=0(零偏置=纯 MD;CV 仍评估但零贡献),起点=低能盆轮换,seed=config.seed+1e6+i 与 metad 种子空间不相交;仅 finished 且无在途轮/验证批次可规划;验证批次可多次规划(「再来一批」)。
   - `gather_verify(idx, minima)`:缓冲收齐后**只做分类折叠**(访问/转移/驻留/pruned 断链)进 **MD 统计梯**,新极小照常入系综(诚实);无交叉、无收敛判定、rounds 不动。步长/淬灭间隔由驱动侧 ImtdOptions 覆盖(RunSpec 不载步数,worker 用 opts 执行)。
   - Ensemble 统计双梯化:`Ladder{Metad,Md}`,note_visit/residence/transition/pruned 带梯参数;try_add_quenched 创建改推零、由调用方梯计一次访问(排雷语义变化:交叉产盆不再虚计 1 访问——占据 0% 直到被轨迹访问,更诚实);逐出/重排同步双梯。
   - `Landscape` 增 md 块(md_basins/md_transitions/md_frames/md_pruned)+ dgraph_md(MD 边的断连森林,同 dgraph 纯函数复用);landscape_json 增 "md" 子对象(无验证数据时省略)。
2. **WASM(加法)**:`ImtdSearch.plan_verify_json(n)`、`gather_verify(idx, minima_json)->bool`(镜像 gather 错误通道);d.ts 导出集合不变(类 +2 方法)。
3. **驱动(app/tasks/imtd.js)**:start() 返回的 handle 增 `verify(nRuns, stepsOverride)->Promise<landscapeJson>`——主循环后协调器不弃,验证批次自建临时 worker 池(同 ready 门/轮内散射),与主循环共享 cancelFlag(取消即断);opts 覆盖 steps_per_run/quench_interval。
4. **UI(index.html 景观区)**
   - 视图切换「断连图 | 网络图」;**网络图**:确定性力导向(圆周初值+固定迭代斥力/弹簧/向心,无 RNG),节点=盆地(半径∝√占据%,能量热色,标注 #i),边宽∝转移计数,边源=MD 梯优先否则 metad;点节点→showEnsembleConf 联动。
   - 「MD 占据验证」按钮:调存活 handle.verify(W, max(2000, steps×2));期间持 min(W,调度器总量) interactive 令牌;完成刷新 lastLandscape 重绘;盆地表增「MD 占据」列(有数据时);元行追加 MD 帧数;reset3DState/新运行清 handle。
   - 断连图在 MD 数据在手时可切 dgraph_md(切换跟随边源选择)。
5. **测试**:原生(plan_verify 规格:零偏置/种子不相交/守卫;MD 梯折叠:md_frames 记账、Σin==Σout、md 占据和≤100%、确定性;交叉产盆 0 访问语义;JSON md 块);CDP 新套件 **m18_landscape_v2**(断连图↔网络图切换:节点==盆地数、点节点联动;MD 验证按钮:md 帧>0、MD 占列出现、元行更新、令牌全释;零 page error);run.sh 挂 m18。
6. **不做(记录)**:鞍点精搜、占据玻尔兹曼重构、rerun UI、GFNFF quencher 缓存、MD 验证的任务层记账(视图精化动作,interactive 令牌即可;后续可补)。

## 验收(实测)
cargo fmt/clippy 双目标 0(unused_parens 九处清理);cargo test 348(lib,+1 验证阶段)+ 7 全绿;wasm 重编 d.ts 导出集合零变化(ImtdSearch 类 +2 方法 plan_verify_json/gather_verify);run.sh 30/31 CDP 全绿(m1 xtb 环境项)+ m18 14/14 + m15/m16/m17 零回归 + Node 153/153。

## 排雷存档(本轮)
- **JS 闭包域**:verify 定义在 start 的 promise IIFE 之外,引用不到 IIFE 内的 search/wasmUrl——"search is not defined" 运行时才爆;修法:search 提升到 start 作用域由 IIFE 赋值,wasmUrl 在 verify 内自取。排雷手法:page.evaluate 直调 handle.verify 返回错误串,比等状态行超时快得多。
- **wasm 包滞后**:lib.rs 加了方法但没重编 pkg —— 浏览器里 "plan_verify_json is not a function";凡动 WASM 面,先 npm run build 再跑 CDP。
- **创建访问语义变化**:try_add_quenched 改推零、调用方梯计一次访问——交叉产盆从虚计 1 访问变为 0% 直到被轨迹访问(更诚实);frames 记账同步移到 note_visit/note_pruned。
- **Ladder 双梯重构**:六条平行向量收进 LadderStats 结构(on_new_basin/on_evict/note_*),Ensemble 持 metad+md 两梯——逐出/重排一处实现两梯共用,免六字段×2 的散弹式修改。
- **python 拼接残渣**:drawDgraph 抽取后旧 drawLandscape 体内留下引用未定义 W 的悬空变量——整块重写(6679..6922)比补丁可靠;node --check 内联脚本兜底。
- classify_frames 的 DRY:metad 折叠与验证折叠共用逐帧分类(同判据/同转移链语义),差异仅梯与 Origin;新极小 round_found 取 rounds_executed。
