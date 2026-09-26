# Plan: MMFF 优化 opt1 差距审计(3× vs RDKit)→ 目标 v1.2.8

## 背景

嵌入与管线已达 RDKit 持平;剩余差距集中在纯 MMFF94s 单分子优化
(opt1:aspirin 3.4×、ibuprofen 3.0×)。v1.2.5 已把单点 E+G 提速
5×(解析梯度),需重新审计剩余 3× 的构成。

## 任务

1. **分解审计**(插桩,OPT_ITERS 门控):
   - opt1 时间 = FF 构建 + L-BFGS 循环(能量求值次数 / 梯度求值
     次数 / 迭代数 / 线搜索拒绝率)
   - 剖析 compute_energy_and_gradient_into 内部:键合项
     (X_energy + X_gradient 分离调用 = 几何量算两遍)、非键循环、
     分配
2. **与 RDKit 求值次数对照**(其总时间 / 单点成本估计)
3. **按数据实施 1–2 项**(候选:键合项 E+G 融合——bond 的 r、
   angle 的 θ 等几何量只算一次;线搜索参数;优化器内分配)
4. **门禁**:cargo test、benchmark 230/230(能量逐位不变——只动
   梯度/调用结构)、ensemble 6/6、clippy 0、fmt;原生 + wasm
   交错 A/B
5. **发布**:1.2.7→1.2.8;CODE_STATUS/PLAN;commit+tag;冒烟
   必须命中引擎输出行

## 验收(实施后实测记录)

- `cargo test` 280/280 全绿(+1 oop 解析-FD 测试);benchmark 230/230
  与基线逐字节一致(能量算式逐位未动);ensemble 6/6;clippy 0;
  fmt;wasm(node 冒烟 1.2.8,引擎输出行确认);API 零变化
- **审计发现**:vdW 组合参数(R* 的 exp、ε 的 2 sqrt)与静电
  qq·scale 每对每次求值重算(构建期常量!);vdW 与静电各自独立
  计算距离(2× sqrt/pair);oop 仍是中心 FD(24 求值/项)
- **实施**:
  1. 对列表携带预计算 (r*, ε, qq·scale)(构建期用同一函数→
     逐位一致);非键循环距离单算 + vdw/静电力合并单系数
  2. MMFF oop 解析梯度(归一化版 asin 链式,与 ETKDG 版不同:
     此处 χ=asin(clamp(û·n̂)) 带符号;|s|=1 饱和区零梯度);
     修正过程中抓出自身初版 bug(atom2/atom3 的 dn 必须分离)
- **严格交错 A/B(3 轮,load ~8–15)**:
  - E+G:aspirin 20.7–24.1→8.8–12.8 μs(**~2.2×**);ibuprofen
    31.5–32.3→14.7–15.0 μs(**~2.1×**)
  - opt:aspirin 5.3→2.8 ms(**1.9×**);ibuprofen 14.7–15.2→7.6
    ms(**1.95×**——与 RDKit 原生 7.1–7.5 ms 持平)
  - E-only:ibuprofen 11.7→4.9 μs(2.4×)
- **wasm vs RDKit 原生(同窗口)**:opt1 ibuprofen 21.95→9.82 ms
  (3.0×→**1.38×**);aspirin 7.15→3.23(3.4×→1.52×);ethanol
  1.08× 近持平;**pipe30 ibuprofen 14.48 vs 14.78——wasm 反超**
