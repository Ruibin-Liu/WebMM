# Plan: 并行多 metad——iMTD-GC 散射-收集并行化(引擎原语 + worker 通道 + 并行驱动)(已完成)

用户指令「并行多metad」。WASM 单线程、不引新依赖(无 rayon/threads),按既定架构「worker farm 归前端层」:把 iMTD-GC 调度器分解为**散射-收集原语**,轮内 k 条短程 metad 散射到浏览器 worker 池真并行,主线程协调器折叠去重/交叉/收敛判定。

## 交付
1. **引擎散射-收集原语(src/metad/imtd.rs;well-tempered 引擎仍零改动)**
   - `CvSpec` 公开;`RunSpec{round, run_index, start, cv, hill_height, seed}`(每 run 自包含、纯函数可执行);`Minimum{coords, energy}`(已淬灭极小)。
   - `execute_run(mol, variant, quench_ff, config, spec) -> Vec<Minimum>`:一条短程 metad + 周期淬灭(淬灭移入执行端)。
   - 协调器 `ImtdSearch`:`new`(淬灭输入/建池/fragments)→ `plan()`(轮首一次性定全部起点:交叉种子优先+低能池轮换;缓存幂等;**轮内起点只依赖上一轮末态 ⇒ 所有 run 相互独立**)→ `gather(idx, minima)` 缓冲,收齐按 idx 序折叠 + 遗传交叉 + 达阈判定(**折叠序=顺序路径序 ⇒ 与到达序无关,逐位确定**)→ `done()/result()`(非破坏快照)。
   - `imtd_conformer_search` 重写为同一原语的顺序组合(单代码路径,并行与顺序天然逐位一致)。
   - JSON 编解码(worker 边界):f64 全部走 `"{:.17e}"` 字符串(见排雷 1)。
   - **两级去重**(本轮语义修正):同重原子构架 AND(全原子对齐 RMSD < tol OR |ΔE| < 1e-3)——保留氢转动异构体(乙醇 anti/gauche),防甲基 H 置换假新构象;镜像同能对折叠为 v1 语义(记录)。
2. **WASM 加法导出(lib.rs)**:`run_imtd_metad_run_json(sdf, spec_json, ImtdOptions) -> Result<String,String>`(worker 端纯函数)+ `ImtdSearch` 类(new/success/error/plan_json/gather→bool/done/round/result);`imtd_config_from`/`imtd_result_to_wasm` 提取复用。d.ts 增量恰两符号。
3. **conf.worker 通道**:`imtdRun` 消息(sdf+specJson+[setter,value] optSets)→ `run_imtd_metad_run_json` → `imtdMinima`/`imtdError` 回帖(worker schema 无关)。
4. **并行驱动 app/tasks/imtd.js**(index.html 挂载;无 UI/不接任务层):`window.IMTD.runParallel(sdf, opts, onRound)`——W≤4 conf.worker 池(ready 门+轮内 Promise 散射),plan→散射→gather(任意到达序)→轮回直至 done;返回系综+并行统计(usedWorkers 等);暴露 `optsToWasm`/`DEFAULTS`。
5. **测试**:原生 +4(execute_run 纯性+JSON 逐位往返;plan 幂等;散射-收集(顺序+乱序到达)==顺序路径逐位;氢转动异构体保留且置换不膨胀);WASM +1(协调器类+worker 函数 == 一次性调用逐位,含坏 JSON gather 拒收);CDP 新套件 **m15_imtd_parallel** 13 断言(并行==顺序 energies 逐位、3 worker 全落位、onRound 轮报、坏 SDF 拒绝、零 page error);run.sh 挂 m15。
6. **不做(记录)**:任务层 kind 化、UI 按钮、GFN-FF、取消语义(驱动 terminate 池)、镜像构象区分。

## 验收(实测)
cargo fmt/clippy 双目标(native --all-targets + wasm32 --lib)0 警告;cargo test 343(lib)+ 7(sascore_golden)全绿;wasm 重编 d.ts 仅 ImtdSearch/run_imtd_metad_run_json 两符号增量;run.sh 27/28 CDP 全绿(m1 为 xtb 外部二进制缺失,exec 127 环境项)+ m15 13/13 + Node 153/153;乙醇小负载实测并行加速 1.37×(启动开销占主导,重负载趋近 min(runs, W)×)。

## 排雷存档(本轮)
- **serde_json 浮点解析 1 ULP 丢位**:发射端 ryu 最短表示正确,但解析器非正确舍入——实证 `0.10966905921216365` 往返后 bits 末位-1(原生测试逐位断言当场抓住)。worker 边界所有 f64 走 `"{:.17e}"` 字符串:Rust `str::parse` 正确舍入、17 位有效数字必往返,JS 侧只透传不解释。
- **重原子去重吞掉氢转动异构体**:乙醇 anti/gauche(ΔE=0.18)的重原子 C-C-O 三角形可被 Kabsch 完美叠合 → 全并成 1 个构象(m15 预检 n=1 引出)。两级判据:重原子构架相同 AND(全原子 RMSD 近 OR 能量差<1e-3);后者防甲基 H 置换假新构象(置换标签淬灭到同一能量)。副作用:镜像同能对(gauche±)也被能量判据折叠——v1 语义,区分需二面角符号,留后续。
- **round_found[0]==0 断言过强**:全局极小由 metad 第 1 轮发现是常态(输入恰在高能井),改 `round ≤ rounds_executed` 合法性断言。
- **wasm 对象不可跨 worker**(结构化克隆不含 wasm-bindgen 包装):spec/minima 以 JSON 字符串为载体;ImtdOptions 在 worker 内重建,驱动发 [setter, value] 对让 worker 保持 schema 无关。
- **起点选择语义细化**:原顺序实现轮内 run m 的起点在 run m-1 折叠后选取(与系综增长交织);并行要求轮首一次性定全部起点——顺序路径同此定义,故与上一轮顺序实现可有数值差异(无黄金数值契约,PLAN 已声明)。
- **JS DEFAULTS 表镜像 Rust 默认**是唯一漂移风险点:驱动总是发完整 optSets(两测侧同源),漂移只影响默认行为不改一致性;排雷记录,未锁测试。
- 环境项同前:m1 xtb 二进制缺失(exec 127);临时探针 example 诊断后即删(仓库 0 警告纪律)。
