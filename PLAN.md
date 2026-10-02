# Plan: 多构象 Shape 检索(Search 标签,v1.4.0 起文档化的单构象限制收口)

## 背景

Shape (3D) 检索目前每条目单构象——柔性分子(如本库 9-RotB 分子)排序
取决于运气构象;药效团筛选已多构象化,shape 是最后一块。基础设施齐
备:generate_optimized_conformers_wasm 逐条目 N 构象(药效团实测 55
库@10 构象 2.5s)、引擎侧 prepared-shape 缓存按 SDF 哈希(逐构象天然
命中)、两段式(screen 代理 ~7ms/对 → 全量 44ms/对)控制成本。

## 语义

- **条目得分 = max over conformers**(combo 排序、shape T 阈值、
  Pharm 过滤全部作用于最佳构象);行点击回载该构象的对齐姿态
- 表格增 Conf 列(shape 模式显示,1 基)——与药效团面板同款
- **N=1 = 原路径逐位不动**(单构象 ensureEntry3D、阈值决定两段式);
  N≥2 = 系综路径:生成 N 构象(seed 42+idx、MMFF94s、iter 100,药效团
  先例)→ phase 1 screen 对齐全部构象取条目最佳代理分 → top-50 条目
  × 每条目 top-3 构象(phase 1 排序)全量 shape+color 重打分 → 胜者
  构象做 pharmMatch 与回载姿态
- **自匹配注入**:查询即库成员时,该条目构象 0 的 SDF 替换为查询
  构象 → 对齐恒等 → 自匹配精确 1.0(m5 既有断言语义保持)
- 位点缓存从条目级(`e.sites`)改为构象级(胜者惰性感知)
- 查询侧维持单构象(ROCS 惯例,文档化)

## 实施(app/index.html 纯 JS,引擎零改动)

1. 搜索行增 `confs` 数字输入(1-50,默认 10,`shapeConfsWrap` 仅
   shape 模式显示;onSearchModeChange 同步)
2. `ensureEntryConfs(e, n)` 帮助函数:批量生成+逐构象切片+缓存
   `e.confs`;查询成员注入
3. runShapeSearch 分支:N=1 原代码路径;N≥2 新系综路径(分块异步、
   阶段状态行含构象数、两段式注释更新)
4. renderSearchResults:shape 分支增 Conf 单元格;表头 `confCol`
   随 shape 模式显隐
5. 性能护栏:55 库@10 构象实测目标 ≤20s(生成 2.5s + screen 550×
   ~7ms + 全量 150×44ms);超标则降默认或降 C

## 验收

- m5 增补:①confs=1 旧断言逐项不回归(自匹配 100%、三列、≥60%
  过滤 {aspirin, salicylic});②confs=10 系综:自匹配精确 100%、
  aspirin 居首、Conf 列在位、状态行含构象数、分数全在 [0,1+]、
  salicylic 在结果中;③行点击回载对齐构象(Conf 号一致)
- m2/m4 抽查回归;390px;零 page error;node --check
- 性能实测记录(55 库@10、@5)

## 验收结果(实施后)

- m5 **32/32**(26 旧 + 6 新):系综默认下自匹配经注入精确 100% 居首
  (Conf=1)、Conf 列显示且胜者构象正确、状态行文档化 "10 confs/entry
  (best over conformers) · two-phase (screen 55×10 → full 50×3)"、
  top-50 候选上限;confs=1 旧路径逐项不回归(自匹配居首、无 confs
  注记、Conf 列隐藏——修掉 renderSearchResults 覆盖 sync 的 bug);
  ≥60% 过滤系综下仍 salicylic 留/glucose 出;行点击回载对齐构象
- 质量实测:ibuprofen 查询 naproxen 59.7%→**71.9%**;warfarin 查询
  diclofenac 22.1%→**72.0%**——单构象对柔性分子的系统性低估被收口
- 性能:55 库@10 构象 14.6s(生成 2.5s+screen 550×~7ms+全量 150×
  44ms)、@5 构象 9.3s、@1 旧路径 2.2s;>300 条大库沿用两段式阈值
  语义不变(N≥2 恒两段式)
- 回归:m0 37/37、m1 10/10、m2 11/11、m3 10/10、m4 32/32;390px
  shape 模式无溢出、confs 输入可见;零 page error;7 内联脚本
  node --check
