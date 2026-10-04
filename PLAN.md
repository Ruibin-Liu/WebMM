# Plan: 平台入口链接 + SA score Rust 移植 + scaffold-hop spike(main 线;全部完成)

## 1. 入口链接(小)✅

- Workbench 页导航加 Platform 链接;README 功能清单加平台页条目
- (Demo/Playground 导航不动——教育线与平台线分离)

## 2. SA score 移植(大,验证先行)

Ertl SA score = RDKit contrib sascorer.py + fpscores 数据表。移植
三件:①fpscores 表定位/导出(homebrew/pinned RDKit contrib)②
Rust 侧展开型 Morgan(radius 2, counts)环境 ID **逐位对齐 RDKit**
(先对拍后定论——bit 不齐则查找表静默错分,不许)③打分公式按
sascorer.py 原文移植 + Python 金标对拍(我们语料 50+ 分子)。
失败判据(诚实 no-go→换路):Rust 位对不齐且不可修。
落点:analog 表加 SA 列(与 RotB/ΔMW 并列,不做组合分)。

## 3. scaffold-hop spike(限时,go/no-go)

单芳环交换(benzene↔pyridine↔pyrimidine↔thiophene),≤2 取代位,
取代模式按环距映射(ortho/meta/para);父本 RGD 拆解→新骨架重装
→shape+color(投影位点)对父本打分。输出:若干真实例 + go/no-go
结论(质量 vs 成本)写入 CODE_STATUS;不做产品化 UI。

## 验收

- 链接:Workbench 导航可见平台页,m0 回归
- SA:Python 金标对拍通过则 wasm 导出 + analog 列 + 测试;
  不齐则数据存档 + 诚实 no-go
- spike:例证数字 + go/no-go 记录
- 七套件 + 平台 Node 全绿


## 实施与验收(本轮三件全落)

1. **入口链接 ✅**:Workbench 导航加 Platform 链接;README Features
   加平台页条目(含探索器/SA/aza-hop 摘要)
2. **SA score 移植 ✅(v1.9.0→v1.10.0)**:src/sascore.rs 全新模块
   ——Ertl SA = 展开型 Morgan(半径 2,稀疏计数)**逐位复刻**(gboost
   32 位 hash_combine/向量-对子哈希/BondType 枚举权重 AROMATIC=12/
   键集掩码持久去重/层种子)+ 705 292 项片段表懒加载(app/fpscores.bin
   5.6MB,sa_load_table_wasm)+ 融合环系统芳香感知(Hückel 4n+2 于
   共享原子环系,外环双键计 sp2;萘/吲哚/嘌呤需要系统级)+ RDKit
   隐氢价模型(原始 Kekulé 键级上累计,芳香过价钳制 +0.1 取整);
   **金标 98 分子位集逐位全等(72 golden + 26 probe),wasm 端与
   Python 参考最大差 1.32(唯一 = 类固醇的立体/桥头罚项近似,文档化
   ——位集相同故片段分严格相等)**
   - 排雷链(全存档):AROMATIC=13→12(枚举位);总度含全氢
     (getTotalDegree≠重原子度);环检测 parent 跳过; surgeries
     在含氢 sdf3d 上做→无氢 molblock 重建;C→N 换元素在 type-4
     输入上失败→Kekulé 基底;隐氢按原始键级(五元 [nH] 芳香 N:
     1.5 键累计算 H 错,单+单=2→1H 对)
   - SA 列进探索器表(与 RotB/ΔMW 并列);sa_score_wasm 直读
     molblock 原始键级(引擎 parse_sdf 会芳香化键型——契约破坏)
3. **scaffold-hop spike ✅(go)**:aza-scan——对当前分子每个含氢
   芳环碳做 C→N 元素行手术(31-34 列),canonical 去重,喂同一
   B→C→D 漏斗。paracetamol → 2 个吡啶基异构体,T 83.4%/82.8%,
   flex 87.6%/87.0%,SA 2.03/2.29。**结论 go**:嫁接机器泛化到
   骨架交换成立;正式化(取代位映射、更多骨架)另立项

## 验收数字

- m0-m6 七套件 37/10/11/10/32/44/**26** 全绿;平台 Node 34/34;
  cargo 323+4(golden 位级);clippy 0;v1.10.0 双重建;零 page error
