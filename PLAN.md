# Plan: Workbench 三标签按模式切换 section 显隐 + 输入框按模式隔离(纯 app 端)

## 背景

用户指出:switchMode 只切了 tab 高亮、batchBar/searchBar、singleActions 与
placeholder,**结果 section 与输入框内容不随模式走**:

1. **#output(2D Structure/Properties/3D Structure/Representations)在
   Batch/Search 模式下仍然挂着**——单分子结果与批量/检索语义无关;
   反过来 #batchPanel/#searchPanel 跑出结果后切去别的标签也依然可见。
2. **textarea#input 三模式共享一份内容**:Single 里贴的 SMILES 切到
   Batch/Search 还在框里,而那里语义是批量清单/库文本(会被
   parseBatchInputs/parseLibraryText 当成多行输入解析)。

既有的行点击回载(loadBatchRow/loadSearchRow/scaffold 行 →
switchMode('single') → 写 input → process)依赖"回 Single",因此
输入框不能一切就清空——**按模式各自记住自己的内容**,切走不串、
切回还在;section 同理按模式记住显隐。

## 修复(app/index.html,纯 JS/标记,零引擎/零 API)

1. switchMode 重构:
   - 切出前:把当前 textarea 值存入 `modeInputs{single,batch,search}`,
     当前模式主结果区(#output/#batchPanel/#searchPanel 之一)的显隐存入
     `modeShown{...}`
   - 切入后:恢复目标模式的输入框内容;`#output` 仅 single 且
     modeShown.single 时显示;`#batchPanel` 仅 batch 且 modeShown.batch;
     `#searchPanel` 仅 search 且 modeShown.search;rgd/scaffold/pharm
     面板维持既有 search 条件逻辑
   - 切标签时 clearError()(错误提示不跨模式残留)
2. 三处会使 searchPanel 失效的入口同步清 modeShown.search:
   loadSearchLibrary / loadDemoLibrary / clearSearchLibrary
3. drop 拖放按模式分流:仅 single 模式自动 process();batch/search 只填
   文本框(search 额外提示去点 Load library)——现存量行为是在批量/检索
   模式下把多行库文本当单分子 process 报错
4. 不动:rgd/scaffold/pharm 既有显隐条件、restoreSearchLibrary、
   行点击回载链路(switchMode('single') 先于写 input,缓冲在下次切走时
   刷新)、URL hash/History/JSME 写入路径(均在 single 模式可达)

## 验收(Playwright,localhost:8901;desktop 1440 + 390px)

1. single 处理 ibuprofen → output 可见;切 batch:output 隐、输入框空、
   batchBar 显;贴 3 条 SMILES 跑 batch → batchPanel 显
2. 切 search:batchPanel 隐、output 隐、输入框空;装 demo 库 + 跑一次
   similarity → searchPanel 显;点行回 single:输入=该行分子、output 显
3. 回 batch:批量输入文本恢复、batchPanel 恢复;回 single:ibuprofen
   恢复、output 显且不重复 process(input===lastProcessedInput)
4. clearSearchLibrary 后切走再切回:searchPanel 不复活
5. 390px 无横向溢出;全程零 page error;node --check 全部内联脚本;
   m0/m4 CDP 冒烟;m3 导航回归;cargo test 295(引擎零改动 sanity)

## 验收结果(实施后)

- Playwright 33/33 通过(localhost:8901,桌面 1440 + 390px):
  ① single 处理后切 batch/search → #output 隐藏、输入框不串;
  ② batch 3 条跑完切走再切回 → 输入文本与 batchPanel 都恢复;
  ③ search 装 demo 库检索 → 点行回 single 正常,searchPanel 不残留在
  single;回 search 面板恢复;
  ④ clearSearchLibrary 后切走切回 searchPanel 不复活;
  ⑤ search 模式 drop 只填框 + 提示,不再误触发单分子 process 报错;
  ⑥ 切标签 clearError;⑦ 390px 三模式 scrollW==clientW==390;
  ⑧ 零 page error;截图目检 batch 页无单分子 section、single 恢复完整
- 回归:m4_batch 7/7;m3 9/10(axe nested-interactive 存量,文献在案);
  m0 33/37——4 个 "props identical" 失败经 git stash A/B 证实存量:
  参照站 molecule-clipboard 无 QED 行使整表比对不等,与本次改动无关
- node --check 7 个内联脚本全过;cargo test 295 全绿(引擎零改动 sanity)
