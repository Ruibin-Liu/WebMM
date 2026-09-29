# Plan: 大库形状检索两段式筛选(screen 初筛 + Top-N 全质量重排;v1.6.0)

## 成本实测(决策依据)

ibuprofen(最慢用例,wasm):26 起点 526ms;10 帧起点 245ms;
**max_iter 5 vs 200 ≈ 233 vs 251ms——L-BFGS 迭代近乎免费,成本在每起点的
全量包含-排斥重打分(~24ms/起点)**。故初筛模式 = 跳过 IE 重打分、按
成对代理重叠直接排序(预计 ~15-20ms/条 = 全质量 1/25-1/30)。

## 设计

### T1 引擎(v1.5.0 → 1.6.0,additive)

- `AlignOptions.screen: bool`(opts_json `"screen": true`):
  - 仅帧起点(10 个,无随机)、无 polish
  - **不做任何 IE/颜色重打分**:按 optimize_pose 返回的代理重叠 O 选优
  - 返回 tanimoto=0(未计算,文档化),surrogate_overlap = 排序键
- Rust 单测:screen 模式 surrogate > 0、tanimoto==0、速度不回归

### T2 页面两段式

- 库 > 300 条自动分层(状态行明示两阶段):
  - **Phase 1**:`shape_align_wasm(sdf, sdf, '{"screen":true}')` 全库,
    按代理重叠取 **top 50**(代理是绝对重叠、偏大分子——初筛偏保守,
    文档化);进度分块
  - **Phase 2**:top 50 走现有全质量联合 shape+color 路径(16 随机起点
    + polish + 颜色),结果排序/表列不变
- 库 ≤ 300:现有单段路径不变

### T3 验证

1. **召回率测量**(Playwright,demo 库 55 条强制两段):全质量单段
   top-20(combo)vs 两段 top-20 → 重合 ≥ 18/20;差异条目如实记录
2. 性能:screen 相 ≥ 全质量相每条耗时比 ≤ 1/10(实测断言)
3. 端到端回归:小库(≤300)路径行为不变;390px;零 page error
4. 门禁:292 测试 + 新增、clippy(1.98,`-- -D warnings` 强制全查)、fmt

## 边界

- top-N 固定 50(不做 UI 设置);颜色不参与初筛(代理是纯形状——文档化);
  screen 阈值 300 条固定
