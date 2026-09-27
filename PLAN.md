# Plan: Workbench file:// 打开时的裸 TypeError 修复(RDKit 未加载防御 + 可操作提示)

## 背景

用户更正复现路径:并非 file://,而是在 Workbench 里点 History 的一条
记录时看到 `TypeError: Cannot read properties of null (reading
'get_mol')`。

真实根因(代码审查 + 竞态复现确认):History 按钮与弹窗在 RDKit
初始化完成前就可用(showHistory 直接读 localStorage,不经 RDKit);
vendored RDKit_minimal.wasm 有 7.3MB,http/gh-pages 首访下载需数秒——
该窗口内点 history 条目(loadFromHistory → process() →
rdkitModule.get_mol)即裸抛 TypeError。file:// 只是让该窗口变成
『永远』(上一轮已修),慢网络下的竞态是独立且更普遍的入口。
history 弹窗内的 Export CSV(exportHistoryToCSV 直接调 get_mol)
同样可达;其余调用点均在已守卫的 process()/runBatch() 下游或
disabled 按钮之后。

## 任务(app/index.html,纯 JS 防御层,上一轮基础上补齐)

1. rdkitInitFailed 状态位:区分『仍在加载』(竞态窗口,提示
   several MB / try again in a moment)与『加载失败』(file:// 归因 +
   `python3 -m http.server` 指引,或保留原始错误);
2. exportHistoryToCSV 顶部加 rdkitAvailable() 守卫(history 弹窗
   内的最后一个未守卫入口);
3. initRDKit 成功后清空 #error(去掉残留的 still-loading 提示)。

## 验收(实施后实测记录)

- Playwright 竞态复现(route 将 wasm 延迟 8s + localStorage 预置
  history):加载窗口内点 history 条目——无 TypeError,提示为
  still loading / try again in a moment;同窗口 Export CSV——
  无 TypeError;wasm 到位后一切正常、#error 自动清空、全程零
  page error(7/7 + 复核);
- file:// 场景保持:失败归因 + http.server 指引,无 TypeError;
- http 正常加载:行为零变化;
- cargo test 281/281、clippy 0、fmt(零 Rust 改动,门禁例行)。
