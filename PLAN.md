# Plan: 单分子优化提速 — GFN-FF C6 参考表预计算 + 线搜索 warm-start

## 背景(实测诊断)

- GFN-FF E+G 比 MMFF 慢 ~17×(n=33:3939 μs vs 508 μs);`sample` 剖析
  4180 个采样点中 ~77% 落在 `alpha_ref`:D3/D4 色散的参考 C6 系数
  (23 点梯形 × refnᵢ×refnⱼ)对每对原子、每次力场调用重算,而它只依赖
  元素对 (Zᵢ,Zⱼ),是纯静态量(xtb/dftd4 均 setup 预计算)。
- Armijo 线搜索每次迭代从固定初始步长(1.5 Å 位移上限)开始、纯减半
  回溯,实测每迭代烧 10–26 次全梯度评估(MMFF 11–13×,GFN-FF 25–26×),
  且每次试探算完整梯度而 Armijo 判据只需能量。
- 放大器:aspirin/ibuprofen GFN-FF 1000 步不收敛(50 s/102 s)。

## 任务

1. **GFN-FF C6 参考表预计算**(`src/gfnff/mod.rs`)
   - `Gfnff` 新增字段 `c6ref: HashMap<(usize, usize), Vec<Vec<f64>>>`,
     `new()` 里对本分子出现的元素对(key 取 min/max 归一)按**与现行
     完全相同的累加顺序**构建:`s += tw[k]·alpha_ref(zi,a,k)·alpha_ref(zj,b,k)`
     (存梯形和 s 本身,不乘 thopi——保证 energy/E+G 两条路径各自保持
     逐位不变:energy 路径 `gw·gw·thopi·s`,梯度路径 `refc6=thopi·s`)。
   - `d4_dispersion` 与 `d4_dispersion_grad` 的逐对三重循环改查表。
2. **ForceField trait 增加 energy 默认方法**(`src/forces.rs`)
   - `fn energy(&self, coords) -> f64`,默认实现走 energy_and_gradient
     (丢梯度);非破坏性新增,MD/metad/solvation 实现方零改动。
   - MMFF 覆写:新增真正的 energy-only 路径(复用现有 *_energy 项函数,
     按与 compute_energy_and_gradient_into 相同的累加顺序求和——结果
     逐位不变,跳过全部梯度计算)。
   - GFN-FF 覆写:走既有 `Gfnff::energy()`(单点 Eh→kcal)。
3. **线搜索改造**(`src/optimizer/mod.rs`)
   - Armijo 试探一律用 `ff.energy(...)`(不再分配 dummy 梯度);slope
     (g·d)在循环外算一次传入。
   - 回溯由纯减半改为二次插值(用 f0、slope、f_trial),下限夹在
     [0.1·α, 0.5·α];退化时回退减半;min_alpha 地板与“有限且不增能量
     才接受”语义保持不变。
   - 初始试探步:L-BFGS 步改为标准单位步 α₀=1.0(两循环递归产生的
     方向自带逆曲率度量,这是 Nocedal 惯例;旧 1.5 Å 位移上限每迭代
     过冲 ~1000×烧 10–26 次试探)。最陡下降前几步保持 0.5 Å 位移缩放。
     (实现中试过 warm-start(上次接受位移),实测绘出 Armijo 弱判据
     下步长逐迭代坍缩的爬行,废弃;见 CODE_STATUS。)
4. **GFN-FF 收敛复测(测量门控,无预授权改动)**
   - 任务 1–3 落地后复测 aspirin/ibuprofen GFN-FF 优化:若 ≤1000 步
     收敛则记录即可;若仍不收敛,诊断失速点的 max_f 轨迹并在
     CODE_STATUS 记录结论(物理/阈值类修复另立项,不在本计划)。
5. **回归测试**
   - 新增:C6 表与逐对现算逐位相等(随机元素对抽样);MMFF
     energy-only 与 calculate_energy 逐位相等;优化器调用计数回归
     (ethanol MMFF 每迭代平均力场调用 ≤5)。
   - 既有 256 测试全绿(opt_compare 能量窗口、ensemble 统计容差、
     GFN-FF xtb 锁值、MMFF 230/230 单点奇偶)。

## 验收

- `cargo test` 全绿(实测 260/260,含 4 项新回归)、`cargo clippy --all-targets`
  0 警告、`cargo fmt` 干净
- `python3 scripts/benchmark_mmff.py --no-speed` 230/230 不回归(与改动前
  输出逐位一致;phosphirane/cyclobutene 两离群点为主分支既有)
- 性能门(native release,同机实测,3 次均值):
  - GFN-FF E+G:aspirin ~122 μs(原 2008,16×)、ibuprofen ~287 μs(原
    3939,14×)✓
  - 线搜索:MMFF 每迭代平均力场调用 2.1–2.3(原 11–13)✓;GFN-FF 实测
    ~7(偏差说明:尾部迭代在 f64 能量分辨率地板上磨到 α 下限拉高均值,
    与梯度 bug 修复后新发现的精度地板问题司源,见下)✗(宽松超出)
  - 端到端:aspirin MMFF 259→32 ms、GFN-FF 50,302→354 ms(142×);
    ibuprofen MMFF 2065→175 ms、GFN-FF 102,355→953 ms(107×)✓
- 实施中发现并修复存量 bug(超出原计划范围,属任务 4 诊断的直接产物):
  GFN-FF 适配层梯度双重换算 1186×(内侧 gxtb 移植已在尾部转 kcal/mol/Å,
  M1 时代包装层又乘 627.51/BOHR)——方向导数 FD 对拍证实;修复后
  GFN-FF 力阈值首次可达、极小更深(aspirin -2537.647→-2537.677)
- 任务 4 结论:梯度修复后 GFN-FF 大分子仍 conv=false(maxf 失速在
  0.2–0.3,分布(f64 能量分辨率 ~1e-9 之下的软扭转模式;FD 证实残余
  力真实、但其可兑现能量下降低于分辨率;xtb 自身阈值更宽)。能量已
  收敛到机器精度;物理/阈值类修复另立项
- API 契约:WASM 导出签名零变化;ForceField trait 仅增默认方法
- CODE_STATUS.md Recently Completed 顶部追加条目
