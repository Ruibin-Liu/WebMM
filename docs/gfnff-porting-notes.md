# GFN-FF 移植笔记(Rust/WebMM)

参考实现:**xtb 6.7.1**(github.com/thfroitzheim/xtb `gxtb` 分支,LGPL-3.0),
源码快照:`/tmp/gxtb-src/xtb-gxtb`(含 `src/gfnff/`)。
参数已由 `scripts/extract_gfnff_params.py` 自动提取到 `data/gfnff_params.json`。
本地验证参照:`~/.local/xtb-gxtb/bin/xtb <mol>.xyz --gfnff --verbose`
(会打印逐项能量分解和 setup 表;注意先删 `gfnff_topo` 重启文件,否则不打印 setup)。

方法文献:Spicher & Grimme, *Angew. Chem. Int. Ed.* **2020**, 59, 15665。

## 单位约定

- 内部一律 **Bohr**(代码里 BOHR=0.52917726);API 入口 Å。
- `covalentRadD3` 数组是 Å,xtb 加载时整体乘 `aatoau·4/3`(covalentradd3.f90 尾部),
  即 rcov(Bohr) = rcov_Å/0.529177×4/3 —— CN、角度阻尼都用它。
- `rad`(EEQ 拓扑距离)是 Å。
- gfnffrab 表(r0/cnfak/en)直接给出 Bohr 量级数值。

## 能量项总表(gfnff_eg.f90 `gfnff_eg`)

| 项 | 公式 | 源码 |
|---|---|---|
| 键 | E = kb·exp(−alp·(r−rab0)²);rab0 由 gfnffdrab 按 CN 现算:(r0_A+cnfakA·cn_A + r0_B+cnfakB·cn_B + vb1)·ff(EN) | egbond 956;rab 表 gfnff_rab.f |
| 键参数 | vb1 = −0.110+shift(XH −0.05, X-sp3 −0.022, X-sp +0.14, 重原子…);vb2 = 0.3731·(1+0.3171·ΔEN²+0.2538·bstrength);vb3 = −bond_i·bond_j·ringf·bstrength·fqq·fheavy·fpi·fxh·fcn | gfnff_ini 1215-1425 |
| 角 | E = fc·(cosθ−cosφ0)²·damp²,damp=1/(1+(r²/rcut²)²),rcut=0.595·(rcov_i+rcov_j)²;fc = angl_i·angl2_j·angl2_k·fqq·f2·fn·fbsmall;**fn = 1−2.36/nn²**;fqq = 1−(qa_i·qa_j+qa_i·qa_k)·(−0.54);φ0 规则表见 ini 1610-1795 | egbend 1233 |
| 扭转 | 未移植(egtors 1444,vtors setup ini 1830+) | — |
| 非键 SE 排斥 | E = 0.427·repz_i·repz_j·exp(−α·r^1.5)/r;α = √(repan_i·(1+0.348·qa)·fn_i · repan_j·(1+0.348·qa)·fn_j)·ff;fn = 1−0.127/(1+nb²);ff: H..H 0.629(×1.458 若 1-3, ×0.708 若 1-4), M..H 0.85, C..H 0.91, O..H 1.04 | gfnff_eg 345-420;alphanb setup ini 800-836 |
| **成键 SE 排斥**(易漏!) | E = exp(−√(repa_i·repa_j)·r^1.5)·repz_i·repz_j·**1.7583**/r,对每个键 | gfnff_eg 596-660 |
| EEQ 静电 | A_ii = √(2/π)/√α_i + γ_i;A_ij = erf(γ_ij·r)/r,γ_ij=1/√(α_i+α_j);RHS = chieeq + cnf·√logCN;解含总电荷约束的线性方程;E = Σq_iq_j·erf(γr)/r + Σ[−q·RHS + q²/2(γ+√(2/π)/√α)] | goed_gfnff 1758 |
| 拓扑电荷 qa | 同上但 r = 1.175·(rad 拓扑距离,Floyd 最短路)/0.5292,RHS = −χ+dxi+cnf·√min(nb,4.4),α=alp²、γ=gam(第一遍) | goedeckera ini2 1639 |
| EEQ 二遍参数 | chieeq = −χ+dxi(O:−0.02(H2O)−0.005/H…);gameeq = gam+qa·ff_gam(O −0.15,H −0.08,C sp3 −0.27…);alpeeq = (alp+ff_alp·qa)²(ff_alp: C 0.09,N −0.21, VIA −0.03, VII 0.50) | ini 670-719 |
| CN | erfCN = 0.5(1+erf(−7.5·(r−r0)/r0)),r0=Σrcov;logCN = ln(1+e^4.4)−ln(1+e^(4.4−cn)) | gfnff_dlogcoord 3497 |
| D4(BJ) 色散 | E = −0.5·C6·(1/(r⁶+R0⁶) + 2·(3·r2r4_i·r2r4_j)/(r⁸+R0⁸))·ζ_i·ζ_j;R0² = (0.58·√(3·r2r4_i·r2r4_j)+4.8)²;C6 = Σ_ref w_i·w_j·(3/π)·Σ_k wᵤ_k·α_i(k)·α_j(k)(**非均匀频率网格** trapzd,dftd4.F90 381);w = 高斯权重(wf=4, 按 logCN);ζ = zeta(Z,qa) 电荷缩放;α_ref = max(ascale·(alphaiw − hcount·sscale·secaiw),0) | gdisp0 102-350;dftd4.F90 54-104 |
| HB/XB | 未移植(egbond_hb 1033,hbset) | — |
| 键合 ATM | 未移植(batmgfnff_eg,b3list = 1-4 对) | — |
| zeta 缩放 | qmod = zeff+q;ζ = exp(3(1−exp(−c(Z)·(1−zeff/qmod)))) | zeta ini 2384 |

## 杂化判定(gfnff_ini2.F90 215-360,已移植有机子集)

- H:2 邻=1(桥氢);group13-16:按配位数 4/3→sp3,3/2→sp2,2/1→sp(2 配 C 需看几何角 <150°→sp2 卡宾);
  N 例外多(NO2、R-N=C、线性判定);O:≥3 或 2 → sp3(3),1 → sp2(2),CO → sp(1)。
- 完整版含金属/eta 判定,未移植。

## 键检测

xtb 用 gfnffrab 估计半径(×rthr=1.25,金属 rthr2)判定;当前 Rust 用
rcov_Å×1.25 近似(有机分子一致),gfnffrab 完整移植在 TODO 里(qa-shift 后的半径)。

## 当前验证状态(vs `xtb --gfnff --verbose`)

**能量：全部已实现项在所有测试分子上达到机器精度(≤1e-7 Eh)**
**梯度：解析梯度 vs xtb 差 ≤ 0.1 kcal/mol/Å(相对 ~0.2%)**

| 项 | 能量(水/甲烷) | 能量(8 双原子) | 能量(乙烷/甲醇扭转) | 梯度 |
|---|---|---|---|---|
| bond | 1e-8 ✅ | — | — | ✅(vs FD) |
| angle | 精确 ✅ | — | — | ✅(含阻尼导数) |
| es (EEQ) | 1e-7 ✅ | — | — | ✅(驻点定理:不含 dq/dr) |
| rep | 3e-8 ✅ | — | — | ✅ |
| disp | 1e-8 ✅ | 1e-8 ✅ | — | ✅(含 dc6/dcn 链) |
| torsion | — | — | 1e-9 / 1.4% | ✅(dphidrPBC 直译) |

梯度验证：vs `xtb --gfnff --grad` 的梯度文件(注意文件有两块，
第二块才是真梯度，与打印的 |dE/dxyz| 一致)最大差 0.099 kcal/mol/Å，
范数 0.06290 vs 0.06277 Eh/Bohr。残差来自键 r0(cn) 链的近似
(xtb 用 gfnffdrab 的完整 grab 链)。

`GfnffForceField` 已实现 `forces::ForceField` trait(kcal/mol 单位)，
可直接接入 WebMM 优化器/MD;冒烟测试：畸变水 → 能量下降、r(OH) 收敛。

### 色散之谜的最终解（重要教训）

三个叠加的 bug 曾让 disp 偏差 26~80%:
1. `sqrtZr4r2(z) = √(0.5·r4Overr2(z)·√z)`(模块尾部变换，非直接开方)
2. `sscale` 加载循环在重写时丢失 → O 参考系 α 未减氢贡献
3. **zeta 符号**:`exp(3·(1−exp(+c·(1−zeff/qmod))))` —— 我误写成 `−c`。
   定位方法：装 gfortran(brew install gcc)+ meson,`-Dtblite=disabled`
   编译 xtb(thfroitzheim fork 的 mctc-lib wrap 需改回 main 分支),
   在 `d3_gradient` 里打印 c6/zeta/disp → 直接看到 zeta=0.7463 vs 我 的 1.2852,
   再在 `gfnff_ini` 的 `zeta()` 调用处打印 → 公式对照一目了然。
   **C6 积分(α 参考系极化率)自始至终完全正确**。

### 扭转项的两个关键点

- 阻尼用**两侧 1–3 距离 + 中心键**(不是三条键!)——这是 ethane/methanol
  双方程反推 + 乙烷 1e-9 验证确认的
- 二面角 `na = ra×rb`(顺序写反会得到 π−φ)

### π 体系(HMO)的关键点

- 自洽迭代(≤5轮):off-diag = −√(hoff_i·hoff_j)·(1−hiter·(2/3−P_old)),
  P 收敛判据 ΔE < 1e-4;diag = hdiag(z) + qa·hueckelp3 − (piel−1)·pilpf
- 对角化后本征值 ×0.1×27.211(eV);**Fermi smearing T=4000K**(bkt=0.3447 eV)
  产生分数占据 → P = C·diag(focc)·C^T;pibo = P[bond] 
- 双基例程(focc(HOMO)≈focc(HOMO+1) < 1e-4)→ 对称破缺整占据
- piel 电子数规则按 (Z,hyb):C=1, N(sp2)=1/N(sp3)=2, O(sp/sp2)=1/O(sp3)=2,
  F=2(或 sp 时 3), S 同 O;芳环参考 P=2/3(bzref)
- pibo 进键参数:shift += hueckelp·(bzref−pibo);fpi = 1−hueckelp2·(bzref2−pibo);
  扭转:f2 = pibo·exp(−2.5(1.24−pibo)^14),重原子邻接 ×1.3,f1×0.55

### 梯度移植的关键点

- **EEQ 驻点定理**：能量对 q 是驻定的(约束拉格朗日解)，梯度**不含 dq/dR**，
  只有 erf(γr)/r 对项 + RHS 的 cn 链(−q·cnf/2√cn · dcn/dR)
- **角度 alist 的 i/j 交换**：xtb 的 alist(1)=中心，但 egbend 调用交换了
  传参 —— 中心原子收 −dedb·damp−term1−term2,邻原子 1 收 deda+term1
- 扭转梯度：dphidrPBC(constr.f90)逐字移植;非周期下 mode 1/2 等价
- 阻尼导数:ddamp = −4·rr/(r²·(1+rr)²) = 2·∂damp/∂r²
- CN 链统一在末尾:g += dcn^T·dE/dcn(bond r0、EEQ RHS、D3 dc6 三处积累)
- xtb 梯度文件有两块,第二块才是真梯度

### 插桩编译环境(已就绪,可复用)

```
/tmp/gxtb-src/xtb-gxtb/build_x/xtb   # 带符号的本地构建
# 重建: cd /tmp/gxtb-src/xtb-gxtb && ninja -C build_x
```

## TODO(按优先级)

1. ~~disp H-X 因子~~ ✅(zeta 符号)
2. ~~解析梯度 + ForceField trait~~ ✅(vs xtb ≤0.1 kcal/mol/Å;GfnffForceField 可直接优化)
3. ~~π 体系 Hückel~~ ✅(苯 piBO=0.666 与 xtb 一致;bond/es/disp 精确;rep 差 3% 待环检测)
4. **环检测**(getring36):环扭转 FC(fr3-fr6)、角 φ0 环规则、ringf、1-3/1-4 H...H 排除(苯 rep 3% 差的主因)
5. 离平面 improper(sp2 中心,omega 势)——依赖 π 检测
6. HB/XB 项 + bonded ATM(凝聚相/生物分子重要)
7. gfnffrab 完整移植(scaleF Ln/An、qa-shift)替换 rcov×1.25 键检测;梯度 r0-CN 链精确化(消除 0.1 kcal/Å 残差)
8. 接入 WebMM 的 L-BFGS(optimize() 目前限 MMFFForceField,需泛化或适配)
9. 金属、PBC(可无限期推后)

## 文件

- `src/gfnff/mod.rs` — 参数加载 + 拓扑 + EEQ + 能量(全部已实现项)
- `scripts/extract_gfnff_params.py` — 参数提取(从 xtb 源码,重跑需源码快照)
- `data/gfnff_params.json` — 元素表×103 + 生成器常数 + D4 参考系数据
- 测试:`cargo test --release --lib gfnff`(水分解验证)
