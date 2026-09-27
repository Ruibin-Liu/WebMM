# Plan: Demo 页文案/单位标注补全 + 三页 a11y 与 Open Graph 元数据

## 背景

第二轮评审(文本提取版)与首轮(截图版)交叉后确认的新增有效项,本轮收口。
已核实的事实基础:
- "~1.4× RDKit speed" 徽章语义过时且有歧义:CODE_STATUS v1.2.8 时该数字
  是"wasm 比 RDKit 慢 1.38×(opt1 ibuprofen)";v1.2.10 后实测为持平至
  反超(aspirin 1.03×/ibuprofen 反超 1.28×/ethanol 持平,30 构象管线
  14.48 vs 14.78 ms 反超)。README 无 1.4× 出处,仅徽章一处。
- Metad CV 参数无单位:引擎侧确认 Hill Height = kcal/mol、
  Hill Width = 弧度(二面角 CV)/Å(距离 CV)(src/lib.rs:1149-1150)、
  CV Atoms 为 0 基索引(src/metad/mod.rs 测试注释)、Deposit Every 单位
  为 MD 步。MD Options 组已有单位(dt (fs)/Temp (K)/Friction (1/ps)),
  仅 Metad CV 组缺失。
- "Advanced options" summary 与内层 `<label>MD Options</label>` 语义重复。
- 三页 nav active 链接只有 class 高亮,无 aria-current="page";Demo 播放
  按钮(◀/►/▶)有 title 无 aria-label;三页均无 Open Graph 标签
  (meta description 三页齐备,og:description 直接复用)。
- 副标题三句 68 词,与徽章/步骤卡重复 marketing 信息。

不做(记录理由):⧉→📋 emoji 替换——CODE_STATUS 记录设计决策曾刻意
从操作按钮去 emoji,且按钮已有 title 提示与"copied"成功反馈;
"Load your own" 位置/拖拽上传、图表 tab 化——布局级改动,另行立项;
播放控件换图标集——超出本轮。

## 任务

1. **文案(site/index.html)**:
   - 副标题精简为单句:"MMFF94s energies match RDKit to <0.01 kcal/mol
     on 230/230 validation molecules — computed entirely in your
     browser."(五步叙事/Python 复现留在步骤卡,隐私主张留在徽章+页脚)
   - 徽章 2 `~1.4× RDKit speed` → `≈1× RDKit opt speed`,title 补基准
     说明(WASM vs RDKit 原生,单分子 MMFF94s 优化,v1.2.10 基准:
     aspirin 1.03×/ibuprofen 1.28×/ethanol 持平)——修正为与实测一致
   - 徽章 4 `+MD & metadynamics` → `7 energy terms decomposed`(数字型
     徽章,步骤 2 证据框页内可验证;MD/MetaD 能力由步骤 4/5 卡呈现)
2. **高级选项(site/index.html)**:
   - `CV Atoms` → `CV Atoms (0-based)`,title 注明格式(逗号分隔,
     dihedral 4 原子/distance 2 原子)
   - `Hill Height` → `Hill Height (kcal/mol)`
   - `Hill Width` → `Hill Width (rad/Å)`,title 注明随 CV 类型取弧度或 Å
   - `Deposit Every` → `Deposit Every (steps)`
   - 删除与 summary 重复的内层 `<label>MD Options</label>`
     (保留 "Metadynamics CV" 分组头)
3. **a11y + meta(3 页)**:
   - nav 当前页链接补 `aria-current="page"`(site/index Demo、
     site/playground Playground、app/index Workbench)
   - Demo 播放三按钮补 aria-label("Previous frame"/"Play or pause"/
     "Next frame"),保留 title
   - 三页 `<head>` 补 og:title/og:description/og:type/og:url
     (description 复用各自 meta description;og:type=website)

不改动:引擎/WASM/任何 Rust 代码;布局与样式;JS 逻辑。

## 验收(实施后实测记录)

- 文案:页面渲染文本含 "≈1×"/"(kcal/mol)"/"(rad/Å)"/"(0-based)"/
  "(steps)"/"7 energy terms",不含 "~1.4×"/孤立 "MD Options" label;
  tagline 为单句;
- aria:三页 `a[aria-current="page"]` 恰 1 个且为当前页;播放按钮
  accessible name 非空;
- og:三页 og:title/og:description/og:type/og:url 齐备;
- CDP m0–m4 全绿(8901 服务 + /tmp/caff24.sdf);
- 桌面/移动布局测量不变(390px scrollW==clientWidth,1440px nav 高 54);
- `cargo test` 281/281、clippy 0(零 Rust 改动,门禁例行)。
