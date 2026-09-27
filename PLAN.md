# Plan: 三页全功能真实可用性审查(v1.3.2 后综合回归)

## 背景

连续多轮功能落地后,用户要求对**所有功能的真实可用性**做一轮完整审查
——不是代码审查,而是模拟真实用户路径的端到端行为验证(点每个按钮、
走每条链路、验证每个声称的输出)。

## 审查矩阵(Playwright headless,全程捕获 page error/console error)

**Demo(site/index.html)**:
1. 10 个 preset 逐一:加载 → 奇偶 MATCH(顺带全量验证冻结参考);
2. 五步走查:Embed 3D → 奇偶 → Optimize → MD(5000 步,轨迹/播放器/
   坐标图)→ MetaD(3000 步裸跑);
3. butane 实验(60k):evidence/坐标/标注 + Use frame → Optimize 落
   已知盆地;
4. 自定义 SDF:改 textarea → 奇偶走无参考分支;
5. 视图重置、390px(裸载 + 实验后);
6. 全程零 page error。

**Playground(site/playground.html)**:
1. 加载 → live MD 自跑,overlay E/T 随时间变化;
2. 拖拽原子(pointer 事件)、温度拉高后 T 上升、力发光开关;
3. 坐标图、视图重置、preset 切换、分享链接反馈;
4. 全程零 page error。

**Workbench(app/index.html)**:
1. CCO → 描述符/规则/2D 图;URL 参数加载;
2. History:入史/弹窗/恢复/删除;
3. JSME 弹窗打开/关闭;
4. Embed 3D → 3D 面板/能量表;MMFF94s Optimize 收敛;GFN-FF 小分子;
5. Conformers(乙醇 N=5):列表/选中/导出 SDF;
6. Batch(3 SMILES):表格/CSV 导出;
7. 2D 导出四件套(SVG/PNG/MOL/SDF 下载事件)、3D 导出;
8. 视图重置;全程零 page error。

## 通过标准

- 上述每一项按页面自身语义判定成功(状态栏 success、evidence 内容、
  数据行、下载事件、数值断言);
- 零 page error(设计内 404 探测类 console 警告不算);
- 发现的缺陷:小缺陷当场修复并复验;大缺陷如实记录并单独立项。

## 验收(实施后实测记录)

Playwright 三页全流程审计 **49/49 通过**:

- **Demo 17/17**:10 preset 奇偶全 MATCH(冻结参考全量复核)、五步
  走查(Embed/奇偶/Optimize/MD 5000 步含播放器与坐标图/裸 MetaD)、
  butane 实验端到端(gauche 帧→优化 −4.2938 逐位)、自定义 SDF 无参考
  分支、视图重置、390px、零 page error;
- **Playground 11/11**:live MD 自跑(overlay E/T 实时变化)、拖拽原子
  (能量响应且无爆炸重置)、温度拉高 T→1129K、力发光开关、坐标图、
  视图重置、preset 切换、零 page error;
- **Workbench 21/21**:CCO 描述符/InChIKey、history 入史/恢复/删除、
  JSME 弹窗、Embed 3D + 能量表、MMFF94s 优化(E=−1.34 converged)、
  GFN-FF 优化(converged)、构象系综(N=5→prune 2,图表+轴+导出)、
  batch 3 行 + CSV 导出、2D 三件套导出、3D SDF/XYZ 导出、视图重置、
  零 page error。

**发现并修复 1 项**:Playground 分享链接在剪贴板不可用(权限/非安全
上下文)时仅报 'Copy failed' 无降级——补 window.prompt 降级(镜像
Demo 页先例),实测降级触发且携带 ?mol= 链接、零错误。

**审查中排除的假阳性**(审计脚本问题,非产品缺陷):window.traj 访问
方式(let 不挂 window)、构象面板选择器猜错(实为 #confPanel/#confChart)。

`cargo test` 281/281、clippy 0、fmt 干净(引擎零改动)。
