# scfastq-qc：面向单细胞长读长 FASTQ 文库的结构感知质控与过滤工具

[作者姓名和单位待补充]  
通讯作者：[通讯作者邮箱待补充]  
稿件类型：Bioinformatics 软件 / 方法论文

## 摘要

### 研究动机

单细胞长读长测序文库保留了分子层面的结构信息，而这些信息很难仅通过通用 FASTQ 质控工具进行检查。在这类实验中，一条 read 是否可用于后续分析，不仅取决于读长和碱基质量，还取决于预期的分子特征是否存在于正确的位置、顺序和方向中，例如末端接头、细胞条形码相关序列、引物或 poly(A) 片段。现有通用长读长质控报告能够总结测序产量和质量分布，但不能直接评估每条 read 是否符合用户定义的单细胞文库结构。

### 结果

我们开发了 scfastq-qc，这是一个基于 Python 的命令行工具，用于单细胞长读长 FASTQ 文件的结构感知质控和 read 过滤。该方法接受 JSON 配置文件，用于描述预期 anchor motif 及其顺序；随后在每条 read 及其反向互补序列中搜索这些 anchor，分配文库结构类别，并将结构信息与读长、Phred 质量值、末端偏移、N 碱基比例和非法碱基检查结合起来。scfastq-qc 可生成自包含的交互式 HTML 报告、机器可读的 JSON 汇总、可选 CSV 导出，以及可选的 passed/failed FASTQ 文件。在仓库附带的两个 10,000 条 read 示例数据集中，scfastq-qc 分别在 56.26% 和 54.78% 的 reads 中识别到完整预期结构，并检测到低质量、anchor 缺失、内部接头、截短和 concatemer-like 等失败模式，同时揭示了两个示例之间明显的方向性差异。这些输出表明，结构感知质控可在下游单细胞长读长分析前辅助诊断文库构建和测序问题。

### 可用性与实现

scfastq-qc 基于 Python 3.10+ 实现，并为固定 anchor 搜索提供可选 Rust 加速器。当前源代码位于随附仓库中；投稿前仍需补充公开仓库 URL、软件许可证、归档 DOI 和软件分发信息。

## 引言

质量控制是测序数据分析中的关键检查点，因为未被发现的技术伪影可能会传递到比对、定量、isoform 发现、细胞分配以及下游生物学解释中。对于单细胞 RNA-seq 及相关实验而言，当质控报告能够反映实验文库结构，而不仅仅是通用产量和质量指标时，其价值通常更高。以 QCatch 为代表的 Bioinformatics 软件论文模式强调：质控工具应将复杂处理结果整合为清晰、可复现、且针对特定数据类型和工作流的报告。

单细胞长读长测序带来了额外的质控问题。长读长 reads 可能包含全长或接近全长的分子信息，但其可解释性取决于能否恢复文库的分子结构。一条可用 read 可能需要按特定顺序包含 5-prime 接头、条形码相邻序列、转录本主体和 3-prime poly(A) 相关信号。相反，即使一条 read 的碱基质量可以接受，如果其末端 anchor 缺失、anchor 出现在内部、motif 顺序异常，或 read 的最佳解释是反向互补方向，它仍然可能无法用于可靠分析。这些失败模式不能仅靠总体读长、N50 或每碱基质量汇总充分捕获。

FastQC 等通用工具以及 NanoPlot 等面向长读长的工具，对于广义测序诊断非常有用，包括质量值和读长分布分析。然而，这些工具并不旨在评估用户定义的单细胞长读长文库“分子语法”。类似地，某些文库特异性预处理工具可以对特定协议的 reads 进行修剪或定向，但其决策结果并不总是以交互式、样本级质控报告的形式呈现，也不一定能直接暴露结构失败模式及其过滤后果。

为解决这一问题，我们开发了 scfastq-qc，一个轻量、可配置的结构感知 FASTQ 质控命令行工具。其核心设计是将实验知识与程序逻辑分离：用户在 JSON 配置文件中定义 anchor 名称、匹配规则、预期顺序和质控阈值，而软件对每条 read 应用一致的流式分析流程。该工具旨在帮助用户回答几个实际问题：预期 motif 是否可检测；它们是否以正确顺序出现；文库是否存在截短或内部接头特征；reads 是否混合了不同方向；以及哪些 reads 应保留用于下游分析。

## 方法

### 数据集

仓库包含两个示例 FASTQ 文件，每个文件含 10,000 条 reads，并配有一个 JSON 配置文件。该配置定义了两个 anchor：一个固定的 5-prime 接头序列 `CGACATGGCTACGATCCGACTT`，允许 1 个 mismatch；以及一个正则表达式 poly(A) anchor，`A{8,}`。预期结构为 5-prime 接头位于 poly(A) anchor 之前。示例阈值将长读长最小长度设为 200 bp，将最小 read Q-score 设为 13。这些示例在本文中用于展示报告内容和方法行为，而不是作为生物学性能 benchmark。

### 设计与实现

scfastq-qc 以 Python 包形式分发，并提供命令行入口 `scfastq-qc`。命令行支持单文件运行和批处理运行。在单文件模式中，用户提供 FASTQ 或 FASTQ.GZ 文件、JSON 配置、输出目录，以及可选的 CSV 导出和 passed/failed FASTQ 输出参数。在批处理模式中，该工具可接受输入目录、文件列表，或配置文件内嵌的 samples 列表，并可并行处理多个样本。软件以流式方式读取 FASTQ 记录，因此不需要将整个输入文件加载到内存中。

配置模型包含四个主要部分：样本元数据、anchor 定义、预期结构和质控阈值。Anchor 可以是固定序列，也可以是正则表达式。固定 anchor 支持基于 bounded mismatch 参数的近似匹配；正则 anchor 支持同聚物或简并序列模式等 motif。阈值控制 read 级别决策，包括最小读长、最小 read Q-score、最大允许末端 anchor 偏移、高 N 比例，以及样本级报告中的 warn/fail 判定阈值。

为提高速度，scfastq-qc 包含一个可选 Rust 共享库，用于固定 anchor 搜索。当共享库可用时，Python 层通过 ctypes 加载该库。如果共享库不存在或无法构建，程序会回退到纯 Python 实现。回退实现对 0 mismatch anchor 使用精确搜索；对于近似固定 anchor，则采用 segment-seeding 策略：将 motif 分割为 `max_mismatches + 1` 个片段，通过精确片段匹配生成候选位置，并仅接受 bounded Hamming distance 不超过配置 mismatch 限制的候选窗口。正则 anchor 使用 Python 正则表达式进行评估。

### Read 级结构分类

对于每条 FASTQ 记录，scfastq-qc 会在原始序列和反向互补序列中搜索所有配置的 anchor。Anchor hits 按位置排序，并汇总为检测到的 anchor 名称有序列表。随后，工具分别对正向和反向互补解释进行评分，评分依据包括顺序有效性、结构标签优先级、检测到的 anchor 数量、hit 数量和总 mismatch 负担。得分更高的方向会被保留为该 read 的最佳结构解释。

分类器会分配一组可解释的结构标签。若 read 中所有预期 anchor 均按预期顺序出现，则标记为 `full_structure`。若没有检测到 anchor hit，则标记为 `no_anchor_detected`。对于部分结构 reads，工具会尽可能根据缺失的末端 anchor 进行标记；对于重复或内部 motif，则分配如 `duplicated_anchor`、`internal_5p_anchor`、`internal_3p_anchor` 或 `concatemer_candidate` 等标签。这一标签体系旨在暴露对文库构建评估或下游 read 过滤具有可操作意义的失败模式。

### Read 级 QC bucket 与样本级 verdict

结构分类会与常规 FASTQ 检查相结合，生成主要的 read 级 QC bucket。判定逻辑会检查空 reads、最小长度、read Q-score、高 N 比例、非法碱基、缺失 anchor、异常 anchor 顺序、内部接头证据、concatemer-like 模式、末端截短和密集 anchor hits。工具通过优先级顺序将每条 read 分配到一个主要 bucket；只有同时满足基本质量和结构条件的 reads 才会被标记为 `pass`。

在样本级别，scfastq-qc 汇总总 reads 数、总碱基数、平均和中位读长、N50、最小和最大读长、中位 read Q-score、long high-quality ratio、reversed-read ratio、截短比例、高 N 比例、非法碱基比例、anchor 检出比例、anchor mismatch 负担、结构计数、方向计数和 QC bucket 计数。配置文件中的 warn/fail 阈值会应用于关键指标，包括 long high-quality ratio、correct anchor-order ratio、no-anchor ratio、reversed-read ratio 和 high-N read ratio，最终生成整体 pass、warn 或 fail 状态。

### 报告与导出

主要输出是一个自包含 HTML 报告。报告包含运行配置、基于阈值的 QC verdict、概览卡片、读长与 Q-score 分布、anchor 检出和 mismatch 汇总、anchor 位置轨迹、结构分类图、QC bucket 图以及 read 结构热图。该报告不依赖外部网络资源，便于本地查看和归档。工具还会写出 `summary.json` 以供程序化复用。根据用户请求，scfastq-qc 还可导出 summary CSV、read 级 CSV、anchor-hit CSV、`passed.fastq`、`failed.fastq`，以及按非完整结构类别分组的 FASTQ 文件。

## 结果

### 示例运行汇总

我们将 scfastq-qc 应用于仓库中包含的两个 10,000 条 read 示例 FASTQ 文件。两个样本均在批处理模式下成功完成，批处理报告汇总了每个样本的总 reads 数、long high-quality ratio、correct-anchor-order ratio 和整体 QC 状态。第一个示例包含 8.65 Mb 序列，中位读长为 745 bp，N50 为 926 bp，中位 read Q-score 为 15.58。第二个示例包含 7.48 Mb 序列，中位读长为 685 bp，N50 为 769 bp，中位 read Q-score 为 15.07。

**表 1. 示例数据集汇总**

| 指标 | real_10k_1 | real_10k_2 |
| --- | --- | --- |
| 总 reads 数 | 10,000 | 10,000 |
| 总碱基数 | 8,650,110 | 7,475,020 |
| 中位读长 | 745 bp | 685 bp |
| N50 | 926 bp | 769 bp |
| 中位 read Q-score | 15.58 | 15.07 |
| Long high-quality ratio | 86.06% | 82.23% |
| Full-structure reads | 56.26% | 54.78% |
| Reversed-read ratio | 2.85% | 49.17% |
| 整体 verdict | warn | warn |

该示例也说明了结构感知报告的必要性。两个样本均通过了配置中的长读长质量标准，long high-quality ratio 均超过 80%。然而，由于 anchor 顺序正确的 reads 比例低于配置的 warning 阈值，两个样本均被标记为整体 warn。仅凭读长和质量分布，很难推断出这一结构层面的差异。

### Anchor 检出与结构失败模式

在两个示例中，5-prime 接头的检出比例分别为 73.24% 和 73.80%，poly(A) anchor 的检出比例分别为 95.53% 和 95.98%。固定 5-prime 接头的 mean best-hit mismatch burden 在两个示例中均约为 0.125，说明在允许 1 个 mismatch 的规则下，大多数检测结果接近精确匹配。poly(A) anchor 的 multi-hit ratio 更高，分别为 19.02% 和 21.20%，这反映了同聚物样序列模式的预期歧义，也说明 motif multiplicity 需要与检出比例一起解释。

结构标签将 reads 分为可解释的群体。在 `real_10k_1` 中，5,626 条 reads 被分类为 `full_structure`，2,110 条为 `missing_5p_anchor`，1,817 条为 `internal_3p_anchor`，354 条为 `missing_3p_anchor`，93 条为 `no_anchor_detected`。在 `real_10k_2` 中，5,478 条 reads 被分类为 `full_structure`，2,049 条为 `missing_5p_anchor`，1,928 条为 `internal_3p_anchor`，292 条为 `missing_3p_anchor`，108 条为 `concatemer_candidate`，36 条为 `internal_5p_anchor`，4 条为 `anchor_order_invalid`，105 条为 `no_anchor_detected`。这些标签为区分 anchor 缺失、内部接头样信号和 concatemer-like 模式提供了直接词汇。

### 方向性感知 QC

两个示例在方向组成上差异明显。在 `real_10k_1` 中，只有 2.85% 的 reads 最适合用反向互补方向解释，且几乎所有 full-structure reads 均为正向。在 `real_10k_2` 中，49.17% 的 reads 最适合解释为反向，其中 full-structure reads 中有 2,753 / 5,478 条为反向。该样本级差异被 reversed-read verdict 捕获：第一个示例保持 pass，第二个示例为 warn。当测序或预处理流程可能产生混合链方向表示时，这类方向性感知汇总具有实用价值。

### 过滤输出

启用 passed 和 failed FASTQ 导出后，scfastq-qc 为 `real_10k_1` 写出 5,095 条 passed reads 和 4,905 条 failed reads，为 `real_10k_2` 写出 4,573 条 passed reads 和 5,427 条 failed reads。第一个示例中主要的 non-pass buckets 为 `internal_adapter`、`missing_5p` 和 `low_read_q`；第二个示例中主要为 `low_read_q`、`internal_adapter` 和 `missing_5p`。由于程序同时写出机器可读汇总，用户可以在最终下游过滤前调节阈值或检查单条 read 级别判定。

## 讨论

scfastq-qc 通过将文库结构作为一类核心 QC 对象，解决了单细胞长读长工作流中的一个实际缺口。它不是仅将 FASTQ 质控视作读长和碱基质量问题，而是进一步判断 reads 是否携带预期分子特征，以及这些特征是否出现在预期方向和顺序中。对于那些下游解释依赖于从噪声长读长 reads 中恢复末端接头、条形码邻近片段、转录本主体和 poly(A) 相关序列特征的实验，该设计尤其有用。

该方法有意保持高度可配置。用户可以根据不同文库设计调整 anchor 定义和预期顺序，而无需修改源代码。固定 anchor 支持近似匹配，这对存在替换错误的长读长数据非常重要；正则 anchor 则允许以简单形式表示 poly(A) tract 等 motif。方向性感知分类、末端偏移检查、内部 anchor 标签和批处理报告的组合，为样本质量提供了简洁但细致的视图。

当前实现仍有一些限制，需要在正式投稿前进一步处理。第一，示例数据展示了功能，但并不构成覆盖不同协议、化学体系、测序平台或生物样本的广泛 benchmark。第二，当前固定 anchor 的近似匹配使用 segment seeding 后的 bounded Hamming distance，并不像 edit-distance aligner 那样建模插入和删除。第三，最佳方向评分是启发式的，用户应针对具有特殊或重复 motif 结构的协议验证其行为。第四，投稿前还应补充软件测试、运行时间 benchmark、公开安装说明、许可证信息和归档 DOI。

未来可从多个方向扩展 scfastq-qc。协议模板可降低常见单细胞长读长实验的配置负担。Edit-distance 匹配或 minimizer-based motif 搜索可提升 indel 容忍度。Barcode-aware 模块可将结构级 QC 与细胞级汇总连接起来。最后，将报告输出集成到工作流系统中，可使结构感知 QC 更容易作为可复现测序流程中的常规步骤运行。

## 可用性与实现

该软件基于 Python 3.10+ 实现，并提供可选 Rust 加速器。当前代码位于随附源代码仓库中。正式投稿前，应补充公开 GitHub URL、永久归档 DOI、软件许可证、软件版本、安装渠道、文档 URL，以及在干净环境中验证过的示例命令。

## 数据可用性

本文草稿所用示例 FASTQ 文件和配置位于本地仓库 `examples2/` 下。如果这些数据不适合公开发布，应在投稿前替换为公开数据集，或将合适的演示数据集存入稳定的数据仓库。

## 作者贡献

[待根据最终作者列表补充。建议 CRediT 角色包括：Conceptualization、Methodology、Software、Validation、Visualization、Writing - original draft、Writing - review and editing。]

## 基金

[待补充。]

## 利益冲突

作者声明无竞争性利益。[投稿前请确认。]

## 参考文献

Andrews, S. (2010) FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics. [请核对引用格式和 URL。]

De Coster, W. et al. (2018) NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. [请核对最终文献信息。]

Gao, Y., He, D. and Patro, R. (2026) QCatch: a framework for quality control assessment and analysis of single-cell sequencing data. Bioinformatics, 42(5), btag184.

Oxford Nanopore Technologies. Pychopper: identify, orient and trim full-length cDNA reads. [如果作为直接比较工具使用，请添加版本化引用或仓库归档。]

[请根据最终生物学应用场景补充协议特异性和单细胞长读长测序参考文献。]

## 作者修订清单

正式投稿前，请替换所有占位符，补充公开软件 URL 和许可证，提供版本化 release 元数据，在干净环境中完成安装测试，加入代表性公开数据集，分别 benchmark Rust 加速器启用和未启用时的运行时间，与相关通用 QC 或协议特异性预处理工具进行比较，并从 HTML 报告中整理最终论文图件。
