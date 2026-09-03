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
| 扭转+improper | E=(1+cos(n(φ−φ0)))·fc·∏dampt(两侧1-3距+中心键);**环键 lring 分支**(ini 1897-1920):顺式极小 φ0=0,中心键在 3-环时 rings4=3(f1=fr3=0.3),4/5/6 环 fr4/fr5/fr6(需最大公共环 + notpicon),无公共环且两外侧端基时 nrot=6 φ0=30° f1=0.3,5-环 C·N 酰胺 CB7 f1=5;sp3 特例(N-N/N-O/O-O)在两分支之后统一套用;额外 sp3-sp3 gauche 扭转仅限非环键(.not.lring);improper:sp2 中心 fc·(1−cosφ)·dampt(1st邻三条距),sat-N fc·(cosφ−cos80°)²;邻原子按到中心距离排序;fc=torsf3·(1−Σpibo·0.5)·(1+5qa) 羧基×38/卤代×10/硝基N ×10/f2 | egtors 1444/1520;ini 2040-2150 ✅ |
| 非键 SE 排斥 | E = 0.427·repz_i·repz_j·exp(−α·r^1.5)/r;α = √(repan_i·(1+0.348·qa)·fn_i · repan_j·(1+0.348·qa)·fn_j)·ff;fn = 1−0.127/(1+nb²);ff: H..H 0.629(×1.458 若 1-3, ×0.708 若 1-4), M..H 0.85, C..H 0.91, O..H 1.04 | gfnff_eg 345-420;alphanb setup ini 800-836 |
| **成键 SE 排斥**(易漏!) | E = exp(−√(repa_i·repa_j)·r^1.5)·repz_i·repz_j·**1.7583**/r,对每个键 | gfnff_eg 596-660 |
| EEQ 静电 | A_ii = √(2/π)/√α_i + γ_i;A_ij = erf(γ_ij·r)/r,γ_ij=1/√(α_i+α_j);RHS = chieeq + cnf·√logCN;解含总电荷约束的线性方程;E = Σq_iq_j·erf(γr)/r + Σ[−q·RHS + q²/2(γ+√(2/π)/√α)] | goed_gfnff 1758 |
| 拓扑电荷 qa | 同上但 r = 1.175·(rad 拓扑距离,Floyd 最短路)/0.5292,RHS = −χ+dxi+cnf·√min(nb,4.4),α=alp²、γ=gam(第一遍) | goedeckera ini2 1639 |
| EEQ 二遍参数 | chieeq = −χ+dxi(O:−0.02(H2O)−0.005/H…);gameeq = gam+qa·ff_gam(O −0.15,H −0.08,C sp3 −0.27…);alpeeq = (alp+ff_alp·qa)²(ff_alp: C 0.09,N −0.21, VIA −0.03, VII 0.50) | ini 670-719 |
| CN | erfCN = 0.5(1+erf(−7.5·(r−r0)/r0)),r0=Σrcov;logCN = ln(1+e^4.4)−ln(1+e^(4.4−cn)) | gfnff_dlogcoord 3497 |
| D4(BJ) 色散 | E = −0.5·C6·(1/(r⁶+R0⁶) + 2·(3·r2r4_i·r2r4_j)/(r⁸+R0⁸))·ζ_i·ζ_j;R0² = (0.58·√(3·r2r4_i·r2r4_j)+4.8)²;C6 = Σ_ref w_i·w_j·(3/π)·Σ_k wᵤ_k·α_i(k)·α_j(k)(**非均匀频率网格** trapzd,dftd4.F90 381);w = 高斯权重(wf=4, 按 logCN);ζ = zeta(Z,qa) 电荷缩放;α_ref = max(ascale·(alphaiw − hcount·sscale·secaiw),0) | gdisp0 102-350;dftd4.F90 54-104 |
| HB eg1(未键合H) | A···H···B:bas=(qa·ca·rAH⁴+qb·cb·rBH⁴)/(rAH⁴+rBH⁴),aci 同理;rdamp=damp/rAB³;out-of-line damp(hbacut/radab);完整梯度 | abhgfnff_eg1 2088 ✅ |
| 键合 ATM | b3list=无序1-4对×两端邻居(butane 90项);E=c9·(0.375·ijmk·imjk·mijk/r6+1)/r⁹,c9=∏clamp(1−3qa,±4)·zb3atm³(**zb3atm[i]=−(i+1)·∛batmscal**,Z=1 特例 −0.25);梯度除以长度非平方 | batmgfnff_eg 3348 ✅ |
| XB | A–X···B(X∈Cl/Br/I/S/Se/Te/P/As/Sb):expo=xbacut((rAX+rBX)/rAB−1) 无 radab,rdamp=damp/rBX³,**cb=1.0 硬编码**;列表:A-X键+bpair(X,B)>3+B的xhbas>0 | rbxgfnff_eg 3216 ✅ |
| zeta 缩放 | qmod = zeff+q;ζ = exp(3(1−exp(−c(Z)·(1−zeff/qmod)))) | zeta ini 2384 |

## 杂化判定(gfnff_ini2.F90 215-360,已移植有机子集)

- H:2 邻=1(桥氢);group13-16:按配位数 4/3→sp3,3/2→sp2,2/1→sp(2 配 C 需看几何角 <150°→sp2 卡宾);
  N 例外多(NO2、R-N=C、线性判定);O:≥3 或 2 → sp3(3),1 → sp2(2),CO → sp(1)。
- 完整版含金属/eta 判定,未移植。

## 键检测

xtb 用 gfnffrab 估计半径(×rthr=1.25,金属 rthr2)判定;当前 Rust 用
rcov_Å×1.25 近似(有机分子一致),gfnffrab 完整移植在 TODO 里(qa-shift 后的半径)。

## 2025-09 review 修复 · 第四轮(3-环角度修正 / 完整环扭转分支)

对拍上游 xtb 源码(github.com/grimme-lab/xtb gfnff_ini/ini2/eg)发现并修复两处真实保真度 bug:

1. **3-环角度 +4° 分支复活**:此前直译了 xtb 的字面量 `ringsj+ringsk.eq.102`,但
   xtb ringsatom 对"不在任何环里的原子"返回哨兵 **99**(循环不执行,rings 保持 99),
   而我们的 ring_size 用 0 表示无环 → 分支永假。102 = 99+3:中心在 3-环、
   两个角度邻居中一个在 3-环、另一个无环(例如甲基环丙烷的 H-C(ring)-C(ring) 角)。
   现在显式 0→99 映射;甲基环丙烷 angle +0.045553381 对 xtb +0.045553364714。
2. **环扭转分支补全**(ini 1897-1963,此前只有 acyclic 分支):
   - lring = 中心键在环里(ringsbond);环键扭转极小在顺式 φ0=0(btyp=2 时 nrot=2);
   - 中心键在 3-环 → rings4 直接取 3("the 3-ring is special"),f1=fr3=0.3;
   - rings4=4/5/6 需 ringstorl(最大公共环)==rings4 且 notpicon:
     fr4=1.0(nrot=6,φ0=30°)、fr5=1.5(nrot=6,30°)、fr6=5.7(nrot=3,60°);
   - rings4=0 且 kk/ll 均为端基:btyp=1 时 nrot=6 φ0=30° f1=0.30;
   - btyp=2 ∧ 环键在 5-环 ∧ C(6)·N(7):酰胺邻接时 f1=5(CB7 修正);
   - **sp3 特例(N-N/N-O/O-O)移到两分支之后统一套用**(xtb 顺序)——同时修正了
     pi-sp3 覆盖与特例的先后(如羟胺 N-O);
   - 额外 sp3-sp3-sp3-sp3 gauche 扭转加 .not.lring 门(环键不再重复加)。
   甲基环丙烷 torsion +0.007939443 对 xtb +0.007939436401(此前 +0.0335,偏大 4 倍);
   七项能量 + FD(3.7e-6 kcal/mol/Å)全部精确。新测试 tests_ring3::methylcyclopropane_vs_xtb。
3. **评审误报备案**:detect_bonds 的超配位过滤(逐对跳过)与 xtb icase-2 getnb 逐字一致
   ("if(nnfi.gt.hc_crit) cycle"按对跳过,无按距离裁剪);xtb 最终拓扑用的是完整 nbf,
   两表在正常有机价态(不过上限)下重合。注释已更新说明。
4. 清理:eg2_rnr 过期注释更新;pyr_h2o.xyz 从 /tmp 移入仓库(src/gfnff/pyr_h2o.xyz,
   include_str! 加载);删除死声明 rabd0、重复 nh/no、多余的 hb_ab.clone();
   ipis 的 Floyd-Warshall 提到 π 环循环外(每分子一次)。

## 2025-09 review 修复 · 第三轮(带电体系 / 全角度规则 / N-O 键扭转特例)

1. **ipis(π 体系电荷)**:对含 π 体系的片段做"中性化 EEQ → qheavy 凝聚 → 差分×1.1 取整",
   nelpi −= ipis(gfnff_ini 610-660)。Cp⁻ 的 Hückel Nel=5→6(芳香),七项能量全部精确。
   同时修复单片段时 qfrag 忘设总电荷的 bug(此前 Cp⁻ 电荷被约束成 0!)。
2. **双片段带电体系**:xtb 对两种电荷放置各解一次 EEQ 取 es 低者(ini 570-588),
   替换"最电正片段"启发式。
3. **完整角度 φ0/f2 规则**(ini 1610-1790 全子集):环 3/4/5(82/96/109°)、3-环中心 + 邻居一环一非环 +4(xtb ringsj+ringsk=102 即哨兵 99+3,并非"双邻居都在环")、
   N-sp2 环内 105°、O 的 Ph-O-Ph 109°、酰胺 N(115/1.2)与 π-N(113/1−0.7Σpibo)、
   重原子 aheavy3/4−5nh、卡宾 145/90、卤素 sp 90/linear+0.6/Z^0.15、CO2 f2=2、
   环成员判定改为 ringsbend(含三原子的最小环,收集 ≤6 元全部环)。
4. **完整 N/O 杂化**(ini2 280-357):硝基 N(3 配 + 端氧)→ sp2+itag(无自有 π 电子,
   nelpi N=2)、CO/OH 端氧 → sp、R-N=C / R-N=N / N=N=N / 线性 N → sp。修复脱离 H 的
   水(O+1H → sp)——此前误标 sp2 键参数差 25%。
5. **键规则补全**:N-sp2 键强 bstren(2)×1.04、CO ×0.90、F-F shift 0.22、3-环/醛 CH
   fxh(ctype)、重-重原子 fcn。
6. **扭转规则补全**:pi-sp3 中心键 f1=0.5(N-π 时 0.2!)、sp3-sp3 的 N-N/N-O/O-O 特例
   (60/90° 与 3-25 倍)、酰胺 C-N fij×1.3、α-C=O fij×1.3、饱和 N 惩罚需 piadr=0。
   硝基甲烷扭转从 5× 过强降到精确(FC 0.039 对齐)。
7. 新验证:Cp⁻(7 项)、硝基甲烷(7 项 + 硝基 flag)、NH₃ improper(cos 双阱)、
   脱离 H 的三片段水二聚体。
8. **"键检测阈值 ~1% 差"已澄清为误报**(第四轮):差异源于手算用了错误的 qa——
   实际掉键点由**两遍 qloop 的自洽**决定:第一遍(qa=0)掉键 → 脱离拓扑的
   片段 EEQ 使 O 电荷更负 → 第二遍半径变大、键被重新加回。实测两引擎在完全
   相同的网格区间掉键(水 O-H:s=1.360 保留 / 1.362 掉;甲烷 C-H:1.30 保留 /
   1.33 掉),s=1.360 处(重连拓扑)能量五项全部 ≤5e-7。已固化为
   bond_drop_point_vs_xtb 回归测试。

## 2025-09 review 修复 · 第二轮(qa / 键检测 / egbond_hb / eg2_rnr)

1. **拓扑电荷 qa 的片段约束**:第一遍 EEQ 原先只加总电荷约束(电荷可在片段间泄漏,
   H2CO···H2O 中 O 被拉到 −0.375);xtb 是逐片段约束(goodeckera 每片段一行)。
   修复后复合物 qa 与 xtb 打印值完全一致(H2CO···H2O:−0.327/+0.094/+0.117 |
   −0.530/+0.177),HB eg3 / XB / bATM 的 3-4% 残差全部消除,断言全部收紧到 2e-6。
2. **键检测完整移植**(替换 rcov×1.25 近似):rco = gfnffrab(cn=normcn) −(qa_i+qa_j)·
   rqshrink(0.23),再乘 fat 元素因子(H/O 1.02 等);r < rthr(1.25)·rco;nb 列表另加
   超配位过滤(nbf>4/6)。拉伸 O-H 到 1.3×(1.245 Å)键仍保留,能量 vs xtb 全程 5e-8
   (旧判据 1.15× 即断键,能量跳 0.07 Eh)。
3. **egbond_hb(X-H 键软化)**:给体为 N/O 且存在 N/O 氢键伙伴的 X-H 键,指数
   alp_eff = (1−(1−vbond_scale)·hb_cn(H))·alp,hb_cn 用 kn=27.5、rcov×1.78 的
   陡峭 erf 计数。水二聚体 bond 从 −0.543048 → −0.543578939(xtb −0.543578922535,
   差 1.7e-8)——二聚体全部能量项现已精确。梯度含 dE/dhb_cn 链(H 与受体 B 上的散射)。
4. **eg2_rnr 移植**(芳香 N 受体,孤对位置项):lp = B −(0.50−0.018·repz(B))·
   unit(Σ环邻向量),额外 A···lp···B out-of-line 阻尼(hblpcut=56),lp 位置的
   梯度链(gii 矩阵);xtb 原样保留其 dgb 丢弃的字面行为。water···pyridine hb
   −0.011372527 vs xtb −0.011372530(3e-9),FD 梯度 6e-5。
   至此 eg1/eg2new/eg2_rnr/eg3 四条 HB 路径全部精确。

## 2025-09 review 修复(第一轮)

逐行重读了 gxtb 源码(/tmp/gxtb-src)+ 用本地 xtb 二进制单点对照,修复:

1. **`create_dlogcn` 双重 bug**(梯度 CN 链):代数反演写错(乘应为除)+ 调用方传 raw cn
   而非 logCN → dlogCN/dCN 在 cn≈4 处小 30 倍。修复后水分子梯度 vs xtb 从
   0.099 → **0.000137 kcal/mol/Å**(范数 0.062773724 vs 0.062773545)。
2. **improper 能量形式反了**:原 cos²ω 以平面为极大;xtb 是 fc·(1−cosφ)(平面=0),
   sat-N 双阱是 (cosφ−cos80°)²;邻原子须按距离排序;fc 补羰基×38/卤×10/硝基×10f2。
   甲醛平面/锥形化能量与 xtb 到 1e-9,梯度含 domegadrPBC 移植。
3. **HB 重写为忠实移植**(eg1/eg2new/eg3,含列表构建 gfnff_hbset):水二聚体 hb
   −0.003319534 与 xtb 精确一致;eg3(羰基受体)含 tors_hb/bend_hb 乘积链。
   原实现把键合 H 同时计入 eg1+eg3(过强 2.2 倍)且用几何电荷而非 qa。
4. **bATM**:三体枚举改为无序 1-4 对×两端邻居;**zb3atm 表 off-by-one**(−i 应为
   −(i+1)·∛batmscal,苯差 43%!);梯度组装除以长度(Fortran invsr2ij=1/r 非平方)。
   苯 batm −0.003628897 vs xtb −0.003628901480,丁烷 4 位一致。
5. **HMO 奇数电子**:na=⌈nel/2⌉/nb=⌊nel/2⌋ 分别 smear(occu 语义),总数守恒;
   双基整占据;ε(HOMO)>0.4 eV → Nel−1 重试(硝基/自由基)。carbene C 跳过 +1。
6. **能量/梯度统一入口**:hb/xb/batm/improper 每项一个函数同时算 E 和 ∇,
   energy().total() ≡ energy_and_gradient().total()(原先差 4.54 kcal/mol,
   生产 ForceField 静默丢掉全部 HB/XB/bATM)。
7. clippy 0 警告(原 79);删除被取代的 adjusted_hbbas/aci、is_amide_h、xhaci_coh。

已知残差(记录在测试注释):HB/XB 对 xtb 差 3-4%,源自拓扑电荷 qa 的简化链
(羰基 O −0.375 vs −0.327 等),非公式错误;eg2_rnr(芳香 N 受体孤对项)未移植,
回落 eg2new。

## 当前验证状态(vs `xtb --gfnff --verbose`)

**能量:全部项 vs xtb 机器精度(≤1e-7 Eh)—— HB 四路径、XB、bATM、egbond_hb 软化键、
带电体系(Cp⁻/硝基甲烷)、NH₃ improper、脱离 H 拓扑(17 个分子级测试)**
**梯度:水 vs xtb 差 0.000137 kcal/mol/Å;全部验证分子解析 vs FD ≤ 1.4e-5;
energy().total() ≡ energy_and_gradient().total() < 1e-12;O-H 拉伸扫描(1.0–1.3×)
逐点 vs xtb ≤ 5e-7;掉键点(O-H 1.360/1.362,C-H 1.30/1.33)与 xtb 逐网格一致**

| 项 | 能量(水/甲烷) | 能量(8 双原子) | 能量(乙烷/甲醇扭转) | 梯度(修复后) |
|---|---|---|---|---|
| bond | 1e-8 ✅ | — | — | ✅(vs FD) |
| angle | 精确 ✅ | — | — | ✅(含阻尼导数) |
| es (EEQ) | 1e-7 ✅ | — | — | ✅(驻点定理:不含 dq/dr) |
| rep | 3e-8 ✅ | — | — | ✅ |
| disp | 1e-8 ✅ | 1e-8 ✅ | — | ✅(含 dc6/dcn 链) |
| torsion | — | — | 1e-9 / 1.4% | ✅(dphidrPBC 直译) |

梯度验证：vs `xtb --gfnff --grad` 的梯度文件(注意文件有两块，
第二块才是真梯度，与打印的 |dE/dxyz| 一致)最大差 0.000137 kcal/mol/Å，
范数 0.062773724 vs 0.062773545 Eh/Bohr(相对差 3e-6)。
此前的 0.099 kcal/mol/Å 残差 **不是** r0(cn) 链近似,而是 `create_dlogcn`
的双重 bug(代数反演写错 + 调用方传的是 raw cn 而非 logCN,dlogCN/dCN
在 cn≈4 处小了 30 倍);修复后水分子梯度与 xtb 吻合到机器精度量级。

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
2. ~~解析梯度 + ForceField trait~~ ✅(vs xtb ≤0.0002 kcal/mol/Å;GfnffForceField 可直接优化)
3. ~~π 体系 Hückel~~ ✅(苯 piBO=0.666 与 xtb 一致;bond/es/disp 精确;rep 差 3% 待环检测)
4. ~~环检测(getring36):环扭转 FC(fr3-fr6)、角 φ0 环规则、ringf、1-3/1-4 H...H 排除~~ ✅(第四轮补完环扭转;苯 rep 精确)
5. 离平面 improper(sp2 中心,omega 势)——依赖 π 检测
6. HB/XB 项 + bonded ATM(凝聚相/生物分子重要)
7. 金属/eta、PBC;eg1 未键合 H 的非零能量场景(已移植+FD 验证,未找到 xtb 非零单点)
8. 接入 WebMM 的 L-BFGS(optimize() 目前限 MMFFForceField,需泛化或适配)
9. 金属、PBC(可无限期推后)

## 文件

- `src/gfnff/mod.rs` — 参数加载 + 拓扑 + EEQ + 能量(全部已实现项)
- `scripts/extract_gfnff_params.py` — 参数提取(从 xtb 源码,重跑需源码快照)
- `data/gfnff_params.json` — 元素表×103 + 生成器常数 + D4 参考系数据
- 测试:`cargo test --release --lib gfnff`(水分解验证)
