# metaTP
metaTP: a pipeline for analyzing meta-transcriptome.metaTP is a pipeline that integrated bioinformatics tools for analyzing metatranscriptomic data comprehensively.It includes quality control, non-coding RNA removal, transcript expression quantification, differential gene expression analysis, functional annotation, co-expression network analysis.<br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14997338.svg)](https://zenodo.org/records/14997338)
## 引用
### 版本控制
本研究代码已通过 Zenodo 存档，DOI 为 [10.5281/zenodo.14997338](https://zenodo.org/records/14997338)。  
**稿件使用版本**: v1.0.0([永久链接](https://zenodo.org/records/14997338))。

## Prerequisites
Download metaTP project
```Python
git clone https://github.com/nanbei45/metaTP.git
```
Configure the environment through yaml files
```Python
cd metaTP
conda env create -f metaTP.yaml
```
Install snakemake
```Python
conda activate metaTP
pip install snakemake
```
用户需要将samplesheet.csv文件中的数据替换为自己的数据(只需提供sra_id号，管道会自动进行下)），样本元数据表格式如下：
```Bash
# samplesheet.csv
sample,sra
sample1,SRR22411016
sample2,SRR22411017
...
```
sample: 样本唯一名称
sra: SRA 编号（流程自动下载）

修改 config/config.yaml 以自定义分析参数：
```Bash
#输入输出配置
samplesheet: "samplesheet.csv"   # 样本表路径
output_dir: "results"            # 输出根目录
fastq_dir: "{output_dir}/fastq"  # 原始 FASTQ 文件目录
qc_dir: "{output_dir}/qc"        # 质量控制结果
quant_dir: "{output_dir}/quant"  # 定量结果
deg_dir: "{output_dir}/deg"      # 差异表达结果
emapper_dir: "{output_dir}/emapper"  # 功能注释结果
#分析参数
threads: 24       # 并行线程数
pvalue: 0.05      # 差异分析显著性阈值
fold_change: 1    # 差异倍数阈值
#数据库路径
eggnog_db: "/path/to/eggnog-mapper_database"      # eggNOG-mapper 数据库目录
dmnd_db: "/path/to/eggnog_proteins.dmnd"          # DIAMOND 数据库路径
```
注意：需提前下载并配置 eggNOG-mapper 数据库[参考文档](https://github.com/eggnogdb/eggnog-mapper/wiki)。

## Execute
Dry run: Use --dry-run to see what tasks Snakemake will perform without actually running them:
```Bash
snakemake --dry-run #Does not execute anything, just shows what the process will do
```
<img src="https://github.com/nanbei45/metaTP/blob/master/img/1.png" width="200px"> <br>
Run the following command to get all the analysis results
```Bash
snakemake --cores 4  #The maximum number of CPU cores/jobs to use for parallelization.
```
本流程可通过 Snakemake 的集群模式 适配 PBS/Torque、SGE 等集群系统，支持动态资源分配和作业调度。<br>
根据您的集群类型，修改config/cluster.yaml文件中的cluster字段
```yaml
cluster: "sbatch --mem {resources.mem}G --cpus-per-task {threads} --time {resources.time} --output slurm_logs/%j.out --error slurm_logs/%j.err"
default-resources:
  mem: 16
  time: "04:00:00"
  threads: 4
```
关键参数说明<br>
• {resources.mem}：任务内存需求（单位：GB）。<br>
• {resources.time}：任务最大运行时间（格式：HH:MM:SS）。<br>
• {threads}：任务请求的 CPU 核心数。<br>
• --output/--error 或 -o/-e：作业日志输出目录（需提前创建，如 slurm_logs）。<br>
 提交作业到集群<br>
 ```Bash
# 使用 --profile 指定集群配置，--jobs 控制最大并发任务数
snakemake --profile config/cluster --jobs 100
```
You can view the complete workflow diagram by running the following command.
```Python
snakemake --dag | dot -Tpng > dag.png
```
<img src="https://github.com/nanbei45/metaTP/blob/master/img/dag.png" width="200px"> <br>
<br>
At the same time, we can perform the tasks we want in the overall process.
#### 1. Download the sra sequence according to the ACC number
The metaTP pipeline integrates data download options using the SRA toolkit, run the following command to implement this function.
```Bash
snakemake prefetch_sra2fastq
```
#### 2. Sequence quality test
Assess the quality of FASTQ files using FastQC,run the following command to implement this function.<br>
```Python
snakemake QC_test
```
#### 3. Sequence quality control, rmrRNA contig cds
（1）Quality control (QC_control):<br>
• Trim the FASTQ file using Trimmomatic, including removing adapters, low-quality sequences, and short sequences.<br>
• Distinguish between single-end sequencing data (SE) and paired-end sequencing data (PE), and call different trimming commands for each.<br>
（2）rRNA removal (rmrRNA):<br>
• Use Bowtie2 to remove rRNA sequences and retain non-rRNA reads.<br>
（3）Transcript assembly (to_contigs):<br>
• Use megahit to assemble non-rRNA reads to generate transcript groups.<br>
（4）Coding region prediction (contigs_change_id):<br>
• Use TransDecoder to predict the coding region (CDS) in the transcript, and modify the ID of the assembly result.<br>
（5）De-redundancy (rmdup):<br>
• Use the seqkit tool to deduplicate the predicted coding region sequence, and finally generate a file without redundant sequences<br>
Run the following command to complete the above process.<br>
```Python
snakemake QC_rmrRNA_contigs_cds
```
#### 4. transcript_index
This is a function for creating transcript indexes, which mainly relies on the salmon tool. Its function is to generate index files for transcript sequences (FASTA format) in RNA-Seq data for subsequent quantitative expression analysis.The output folder will contain the Salmon index files<br>
Run the following command to complete the function.<br>
```Python
snakemake transcript_index
```

#### 5. gene_expression_quant
This process is to filter out low-expression genes by checking the expression values ​​of each group of samples in the gene expression data table, eliminating genes with low expression levels in the samples, and retaining genes with a certain expression level.Run the following command to implement this function<br>
```Python
snakemake gene_expression_quant
```

```bash
     library("tidyverse")
     table_filter<- read.csv("transcript_abundance_quantification_table_filter.csv",header=T,row.names=1)
     S_id <- read.csv("Streptophyta_id.csv",header=F,row.names = 1)
     rownames(S_id)

     table_filter2<- table_filter[-which(rownames(table_filter) %in% rownames(S_id)),]
     #table_filter2<- table_filter %>% filter(rownames(table_filter) %in% rownames(S_id))
     nrow(table_filter2)
     nrow(table_filter)
     nrow(S_id)
     write.csv(data.frame(GeneID=rownames(table_filter2),table_filter2),"transcript_abundance_quantification_table_filter2.csv",row.names=F)
```

#### 6. DEG_analysis

```Python
snakemake DEG_analysis
snakemake up_regulated_gene
snakemake down_regulated_gen
```
(1)snakemake DEG_analysis<br>
a.Call the R script DEG_analysis.r to perform differential gene analysis:the output is the differentially expressed genes result file differential_genes.csv: contains gene ID, logFC, P value and gene expression change status (such as upregulation Up, downregulation Down)<br>
b.Filter up-regulated or down-regulated genes from the differential_genes.csv file and generate the gene ID file differential_genes_id.txt.<br>
c.Use the seqkit tool to extract the DNA sequences of differentially expressed genes from the input gene sequence file (FASTA format).<br>
d.Translate the extracted DNA sequence into a protein sequence，The output is differential_gene_sequence.fastadifferential_gene_sequence.pep, which represent the DNA sequence of the differential gene and the translated protein sequence of the differential gene, respectively.<br>
Run the following command to automatically execute the above functions：
```Bash
snakemake DEG_analysis
```

(2)snakemake up_regulated_gene<br>
a.Screening down-regulated genes: Extracting related gene IDs based on expression change status (Down).<br>
b.Extract gene sequence: Extract DNA sequence from sequence file based on gene ID.The output is differential_genes_id_down.txt: a list of IDs of downregulated genes.The output is differential_gene_sequence_down.fasta: the DNA sequence of the downregulated gene.<br>
c.Translate protein sequence: The DNA sequence of the downregulated gene is translated into protein sequence.The output is differential_gene_sequence_down.pep: the protein sequence of the downregulated gene.<br>
Run the following commands to perform downstream analyses such as gene function annotation, enrichment analysis, and structure prediction.<br>
```Bash
snakemake up_regulated_gene
```

(3)snakemake down_regulated_gen<br>
This process mainly involves extracting the sequences (DNA and protein) of up-regulated genes from the results of differentially expressed gene analysis for subsequent analysis, which is similar to the above steps.<br>
Run the following commands to perform downstream analyses such as gene function annotation, enrichment analysis, and structure prediction.<br>
```Bash
snakemake down_regulated_gen
```

#### 7. emapper.py
The integrated eggNOG -mapper provides several key features including: 1) de novo gene prediction based on raw alignments, 2) integrated pairwise homology prediction, and 3) rapid protein domain detection.<br>
```Python
snakemake emapper
```

## Analyze
We use 8 ontologies and rhizosphere metatranscriptome samples from Mendes et al. as examples for case analysis. The following are the specific steps and input and output of the analysis:
#### 1. Downstream analysis.R
Run the script directly using the R language command:
```Python
 Rscript downstream_analysis.R
```
(1)Dimensionality Reduction Analysis<br>
input: <br>
	•transcript_abundance_quantification_table_filter2.csv:Transcript abundance table containing expression values ​​for each sample.<br>
output:<br>
	•The data visualization results after dimensionality reduction, principal coordinate analysis (PCoA) plot, showed significant separation between rhizosphere samples and bulk samples at different doses, indicating significant transcriptional differences in different soil environments.<br>
(2)Venn analysis<br>
input:<br>
	•venn_setdiff_rhi.emapper.annotations.csv：Gene function data annotated by eggNOG, including gene ID and its functional classification information.<br>
output:<br>
	•venn_setdiff_bulk_go.pdf：Venn diagram showing the shared and unique genes between the rhizosphere and rhizosphere samples<br>
(3)Functional Annotation <br>
input：<br>
	•final_table_sequence.emapper.annotations.csv：Functional annotation table generated from the eggNOG annotation tool, including gene ID, functional classification, GO annotation, KEGG pathway, etc.<br>
	•cog_funclass.tab：The COG functional classification table defines the meaning of each COG classification.<br>
output：<br>
	•cog.pdf：Histogram of COG functional categories showing the number of genes included in each functional category.<br>
	•Intermediate result table of functional annotation（such as gene2cog），Mapping of genes to COG functions is documented.<br>
 (4)Functional Enrichment Analysis<br>
 input：<br>
	•differential_genes_id_up.txt：A list of differentially expressed gene IDs.<br>
	•go_term.list：GO glossary, defining the function and classification of each GO term.<br>
	•final_table_sequence.emapper.annotations.csv：GO annotation data of genes.<br>
	•ko00001.json：Pathway annotation information downloaded from KEGG (for KO enrichment analysis).<br>
output：<br>
	•go_rich.txt：The GO functional enrichment analysis result table records the significant GO classifications and the corresponding gene numbers.<br>
	•go_rich_bar.pdf and go_rich_dot.pdf：Bar graph and dot plot of GO functional enrichment.<br>
	•Ko_rich.txt：The KEGG functional enrichment result table records the significant KEGG pathways and corresponding genes.<br>
	•ko_rich_bar.pdf and ko_rich_dot.pdf：Bar graph and dot plot of KEGG functional enrichment.<br>


#### 2. Gene co-expression network
This flow describes the construction process of the gene co-expression network.<br>
（1）Correlation analysis: remove classification information and keep ID (1,1) as the name of the sample and species. Samples are used as rows and species are used as columns. Use the Python script correlation_analysis.py to analyze the correlation.<br>
（2）RMT screening: filter out significant correlation data, extract important structural information, and output an RMT_result.csv file.<br>
（3）Adjacency matrix construction: use the significant correlation matrix to construct the adjacency matrix adjacency_matrx.csv of the gene co-expression network.<br>
（4）Meta network analysis: analyze the network connections between different species at different levels and generate a Meta network.<br>
（5）Subnetwork analysis: analyze the network structure of different samples and groups, extract topological structure features, and generate network information corresponding to each sample.<br>
（6）Subnetwork analysis of each group: analyze the subnetwork structure corresponding to each group (such as different altitudes, different species, etc.) and draw the network of each group.<br>

```Python
  all codes in Gene_co-expression_network.zip
```



