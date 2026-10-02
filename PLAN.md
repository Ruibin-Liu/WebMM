# Plan: 第三方归属常规操作(MIT 保持,补齐声明)

## 背景
用户决定维持 MIT;按常规操作补齐 LGPL 衍生模块(xtb/GFN-FF)的归属
声明与 vendor 组件的上游许可文本(模块头 src/gfnff/mod.rs 已有出处
声明,不动)。

## 实施
1. LICENSE 追加 Third-party components and attribution 节(gfnff→LGPL-3
   + xtb 出处与方法引用;RDKit/3Dmol/JSME→BSD-3 指向文本;MMFF 参数
   与 LBDD 数据→RDKit 提取物 + Halgren 出处)
2. app/vendor/ 补 RDKit-LICENSE.txt、3Dmol-LICENSE.txt(BSD-3 标准文本
   + upstream 链接;JSME 已有 license.txt 不动)
3. README License 节加一句 LGPL 模块与第三方声明指引

## 验收
文本就位、无代码改动、零功能影响。

## 验收结果(实施后)
LICENSE/README/vendor 两文件就位;纯文档,零代码 diff。
