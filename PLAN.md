# Plan: Top-3 —— 原子映射贯通 + 点选约束(工作台交互设计)

## 范围
1. **原子账本(atom ledger)**:确立并实证不变量 —— Ketcher molblock →
   RDKit → WebMM V2000 解析(文件序)→ H 追加在重原子后 → 3Dmol
   serial(1 基)= SDF 序。每次拾取把 3Dmol 报告的元素与 SDF 行符号线
   核对,不符即拒绝并明示(Astra 红线:功能正常但指错原子比没有更危险)。
2. **拾取 UI**(工作台):⌖约束 按钮开关拾取;3Dmol setClickable →
   pickAtom;2/3/4 原子 → 距离/角度/二面角;实时几何读数 + 目标值输入
   (默认当前几何);再点取消;Esc 清选;橙色高亮。
3. **约束集**:chips 列表(可单个 ✕ 移除);冻结所选(freeze 数组);
   restraintJson 组装 {freeze/distance/angle/dihedral}(0 基序,镜像
   lib.rs set_restraints 文档);随「优化」生效,优化往返保留(可迭代
   目标值);molEpoch 变更(新分子)即清——挂在 process() 的递增点,
   无需等嵌入。
4. **worker 通道**:conf.worker 'optimize' 消息加可选 restraints 字段 →
   opts.set_restraints(向后兼容:不传行为不变)。

## 不做(记录)
Ketcher 侧桥接(selection/getSmarts 需重建 31MB bundle 扩展
WebMMKetcher wrapper)—— 下一轮;锚定叠合/核心固定构象生成 —— 依赖
Ketcher 选区,随后。

## 验收
- 新 CDP 套件 m9_constraints 10 断言:账本拒错(元素不符拒拾取)/
  实时几何/约束 chip/**端到端金断言**(d(1,3)=3.2Å 约束使该原子对
  实际位移 0.19Å 向目标)/冻结原子逐位不动/chips 跨优化存活/新分子
  清陈旧约束/零 page error。
- run.sh 全量 22 套件。
