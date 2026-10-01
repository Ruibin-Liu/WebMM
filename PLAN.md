# Plan: 收掉两项存量测试失败(m0 props 比对交集化 + m3 axe nested-interactive)

## 背景

两项"文献在案存量"已多轮 git stash A/B 证实非回归,本轮收口:

1. **m0 4×"props identical"**:m0 把 Workbench 与参照站 MC(molecule-
   clipboard,同 vendored RDKit)整表逐行比对;LBDD 里程碑给 app 属性表
   加了 QED(mean/max)两行,MC 没有 → 整表长度不等恒败。共同字段
   (MW/cLogP/TPSA/HBD/HBA/RotB)实际全一致。
2. **m3 axe nested-interactive(severe)**:#viewer3d 标 role="img",视图
   重置按钮嵌在其内 = 交互元素嵌进纯展示角色(a11e 树被截断)。
   v1.3.2 加框内重置按钮时带入。

## 修复

1. **m0(tests/cdp/m0_core.test.js)**:propsEqual 从"等长逐行"改为
   "按键交集"——MC 侧每行必须在 app 侧存在且值相等(保留 HBA 版本
   漂移豁免);app 多出的行(QED)合法。detail 提示同步改用键查。
2. **m3(app/index.html)**:新增 .viewer-wrap(position:relative 包裹层)
   承载重置按钮为 #viewer3d 的**兄弟节点**(视觉位置不变:按钮锚到
   wrap 右下角);#viewer3d 保留 role="img"。按钮移出后容器重写
   不再误删按钮 → clear3DViewer/initViewer 两处"先取引用再 append"
   的补丁代码删除(简化)。site 两页无 role="img" 不涉。

## 验收

- m0 37/37;m3 10/10(axe clean)
- 3D 流程回归:embed→按钮可见且点击无错;m1/m2 全绿;390px 无溢出;
  零 page error
- 更新 CODE_STATUS(存量清单清零)

## 验收结果(实施后)

- **m0 37/37(历史首次全绿)**:propsEqual 改键交集后共同字段全过
  (HBA 漂移豁免保留、detail 提示改键查)
- **m3 10/10(axe clean)**:按钮移入 .viewer-wrap 后无 nested-interactive;
  按钮锚定 wrap 右下实测在 viewer 盒内、点击无错、容器重写后存活
  (clear3DViewer/initViewer 的补丁代码随之简化)
- 回归:m1 10/10、m2 11/11、m4 32/32、m5 26/26——**m0-m5 六套件
  首次全部通过**;390px 带 3D 无溢出;零 page error
