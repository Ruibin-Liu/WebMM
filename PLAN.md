# Plan: 性能主线第四轮——攀岩大分子验证(多 CV 体系 refine 全链路观察)(已完成)

CODE_STATUS 记录的下一候选。攀岩(钉住 CI-NEB)至今只在丁烷(单脊短路/双脊长路)验证;阈值(严格内极大守卫 +0.1 kcal、回退判据 −0.1 kcal、patience/失速参数)未在大分子上考察。

## 交付
1. **验证探针(examples/climb_molecules.rs)**:fixtures 四体系(ibuprofen 33 原/threonine/aspirin/naphthalene)跑原生 iMTD→landscape→plan_refine→逐边 NEB(默认配置,攀岩开)→gather_refine,统计:
   - 收敛率(conv)、攀岩激活率(climb_started>0)、回退率
   - 相位 2 迭代占比(丁烷 ~3%)
   - 鞍点物理性(高于两端、与 observed PE 峰的差距方向)
   - 每边耗时
2. **阈值评估**:若回退率过高(粗边/复合扭转边天然多)考虑放宽回退阈或按 n_images 缩放;若攀岩激活率异常低则查守卫。只在有实测证据时改,否则只记录。
3. **修复(若验证暴露缺陷)**:冻结期内 fix 范围;若涉及行为语义则并入 B 轮版本说明。
4. 文档:观察结论入 PLAN/CODE_STATUS;随后 B 轮版本收尾。

## 验收(实测,examples/climb_molecules.rs,2×3×800 搜索 + 2 轮 MD verify)
- **ibuprofen(33 原子,7 盆地,3 边)**:e0/e1 攀岩收敛(climb_from 954/542,Δ 势垒 +2.6/+0.2);e2 相位 1 硬边 20000 迭代不收敛、攀岩未激活——**climbing=false 对照逐位相同(27.265/20000),预存特性非回归**,诚实 conv=false。
- **threonine(17 原子,5 盆地,3 边)**:e0/e1 攀岩收敛(相位 2 占比 52%/2%);e2 攀岩跑 4486 迭代后回归→**诚实回退带峰 44.256**(climbing=false 时 3816 迭代直接收敛,差值即攀岩成本)。
- **aspirin/naphthalene**:刚性→1 盆地 0 边,管线正确无事可精修。
- 全部鞍点高于两端(Δ +0.2~+4.1)✓;无 panic;确定性;阈值(+0.1 守卫/-0.1 回退)无证据误伤,**不改**。
- 攀岩统计:6 边中 4 攀成/1 诚实回退/1 预存硬边;refine 成本大分子 ~1.2-2s/边(原生串行;浏览器 worker 并行取最慢边)。
cargo test 372 绿;clippy 双 0;m22 golden 零差异(无 src 改动,纯验证轮+诊断探针)。

## 排雷存档(本轮)
- refine 阈值(+0.1 守卫/-0.1 回退)在 33 原子多可旋转键体系上无证据误伤——验证轮的价值在于"不改"也是结论。
- climb_started=0 有两种语义(守卫拒绝激活 vs 攀岩后回归回退),对外不可区分;浏览器协议(refine_result_to_json 仅 saddle+barrier)不暴露此区别——记为冻结期后候选(worker 协议 v1.2 增 refined 字段),不在本轮做。
- 相位 1 硬边(复合扭转,20000 迭代不收敛)在攀岩开/关下逐位相同——攀岩轮的相位 1 改动确认零行为漂移;此类边的预算分配(如按 n_atoms 缩放 max_iterations)记为冻结期后候选。
- Quencher 在 webmm::metad::imtd 命名空间(非 metad 根);examples 诊断探针保留 WEBMM_NO_CLIMB 环境开关作攀岩开/关对照。
