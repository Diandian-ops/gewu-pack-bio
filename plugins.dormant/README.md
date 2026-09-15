# plugins.dormant/ — 封存插件区（BLOCK-86 收口后测试纪律定型）

本目录存放**暂不适配、不参与装载**的 14 个分析插件（cellchat / celltypist /
deg-standardizer / expression-qc / gene-correlation / go-kegg / gsea /
pca-explorer / sample-metadata-validator / scanpy-advanced / scvi-tools /
seurat-standard / volcano-plot / wgcna）。

## 为什么封存

用户裁决（2026-09-15，原话「选择1 开始这个」）：15 个插件逐个适配测试浪费
时间，先封存大部分活跃度，**保留 `deseq2` 为唯一活跃参考实现**——它持续
验证平台契约（R 运行时 + tool-job 管道 + 产物链 + BLOCK-86 装卸编排），
其余 14 个是冻结的内容资产。软件版本成型后（时点由用户裁决）再激活适配。

## 加载语义

插件加载器只扫描 `plugins/` 目录；本目录是普通目录，**零装载、零目录投影、
零测试义务**。核心仓一行不改（插件化纪律）。

## 激活适配 checklist（单个插件）

1. `git mv plugins.dormant/<id> plugins/<id>`；
2. 把 `plugins.dormant/installed.json` 中该条目移回 `plugins/installed.json`；
3. 重启 dev 栈（boot 装载）；
4. 跑该插件的真实工具作业（demo 数据）+ 检查产物链；
5. 在本 README 勾掉该插件，pack.json version bump。

## 全量激活里程碑

「软件版本成型」时点由用户裁决。届时按上表逐个激活适配（或整批），
并恢复 `plugins/manifest.json` 与实际装载面的一致性核对。
