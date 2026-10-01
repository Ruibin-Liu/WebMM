# Plan: GFN-FF 构象坍缩修复——引擎三层(ATM/HB/XB 平滑短程阻尼 + EEQ 自适应岭 + 线搜索能量地板) + app 两层(worker 重算守卫 + 系综 sanity 过滤)

## 背景

用户分子 C30H32ClN9O3 的 GFN-FF 100 构象任务实测:2/100 构象发散
(E=-1.9e4/-1.4e22,最近原子对 0.56/0.0009 Å),劫持系综下游。深挖
历经四轮假设排除(原生复现失败→平台分叉、EEQ 条件数"良态"、双路径
"一致"、硬门控),最终用"wasm 导出几何 + 原生以 wasm 嵌入建场复算"
闭环定位全部三个真根因:

1. **ATM 三体项**(batm)c9·(ang+1)/(r_ij·r_jk·r_ik)³ 的 1/r⁹ 当
   c9<0 时无下界(实测 0.02 Å 处 -1e4 kcal/mol);
2. **HB/XB 项** damp/r³ 的自带 damps 在 r→0 处正则化不足
   (实测 N-H 0.02 Å 处 -2.1e4);
3. **EEQ 运行时 KKT 求解**:两原子近乎重合时 erf 核行趋同 → 矩阵
   近奇异(实测 min pivot 4.3e-2→3.5e-4)→ 电荷爆到 |q|=1.1e4 →
   es 漏斗 -35329 kcal/mol(这是压死骆驼的最后一根,也是最深的一层)。

过程中证伪的假设与中间发现(全部存档):硬 0.8 Å 门控失败(项在边界
-150..-180,reopening 断崖造出假吸引子,线搜索收割数千 kcal);
GFN-FF 拓扑从建场几何读出(键检测),在别的几何上重建 FF 会得到不同
力场——多个"验证"因此失效;wasm/原生 ETKDG 嵌入有 ulp 级分叉被混沌
放大,决定哪个种子踩中奇异流形;线搜索 f_new=-inf 时 `<=` 为真会被
接受(存量隐患,顺修)。

## 引擎侧(v1.6.3→v1.6.4,src/gfnff + src/optimizer)

1. **clash_damp/clash_damp3**:f(r)=1-exp(-(r/0.66Å)¹²),三距离乘
   积 D 乘在 ATM/HB/XB 项上(能量与梯度严格链式 E=D·E₀、∇=D∇E₀+E₀∇D)。
   健康距离 f=1-10⁻¹⁸=**精确 1.0**(f64 逐位零漂移);0.8 Å 处
   f=0.99993(无断崖);<0.4 Å 项灭活。r→0 时 f,f'→0 完全良态。
   变量命名 cdamp/cgd 避开 HB 内部自带 damp 的遮蔽(曾致双乘 bug)。
2. **solve_eeq_kkt**(energy() 与 eg_core() 两处共用):|q|max≤2e 用
   精确解(逐位零漂移);爆炸时 ×3 粗定位 + 18 次二分调 Tikhonov 岭
   使 |q|max≈2——电荷与 es 有界且随几何连续,确定性保证两路径一致。
   岭区内梯度近似(文档化折衷)。
3. **armijo_line_search 能量地板** -1e8 kcal/mol:拒绝非有限与超地
   板试步(顺修 -inf 被接受存量隐患)。任何未来未知的奇异项都被兜住。
4. 新测试 5 个:ATM/HB 坍缩有界(修复前 +1.9e6 失败实证)、阻尼区
   连续性(0.005 Å 步长无 >25 kcal 跳变)、EEQ 岭坍缩有界 + 双路径
   一致、线搜索拒收 -1e22/-inf、该分子 3 构象全管线健康带断言
   (夹具 tests/fixtures/gfnff/c30h32cln9o3.mol 入库)。295→301 测试,
   GFN-FF 全部 xtb 奇偶夹具零漂移通过。
5. 版本 1.6.4;wasm 重建;node 冒烟:100 构象零 |E|>10000 离群
   (min -9258/中位 -9093,修复前 -1.4e22)。

## app 侧(app/index.html + app/conf.worker.js)

6. **worker E 重算守卫**:量化 SDF 上重建 FF 的键感知可能换拓扑致
   NaN→JSON null(实测 rep 项 sqrt 负数)——`Number.isFinite` 守卫,
   非有限回退批量 E(真正优化该构象的 FF 之值)。
7. **finishConformers sanity 过滤**:非有限 E 或 min-dist<0.3 Å 的
   构象剔除 + 状态行注明(任何引擎的最终兜底)。
8. m2_conformers 增 GFN-FF 节:该分子 100 构象 → 能量全有限且
   -10000<min、max<0、轴文本良构。

## 验收结果

- 原生:301 测试全绿、clippy -D warnings 0、fmt 干净;GFN-FF xtb
  夹具逐位零漂移;坍缩阶梯扫描 14→0 负爆;EEQ 岭 z_30 es
  -35329→-160.5
- wasm 冒烟:100 构象零离群(min -9258.4/中位 -9093.2/max -8214.4)、
  min-dist 0.551、全部分项 sane;版本 1.6.4
- CDP:m2 11/11(GFN-FF 节过)、m1 10/10、m4 32/32×3、m5 26/26;
  m0 33/37 与 m3 9/10 为文献在案存量(MC 参照缺 QED 行/axe 项);
  m1 的 /tmp/caff24.sdf 配方错误(README 写了茶碱)顺修
- 临时诊断 example 与探针 hook 全部清除;README 配方修正
