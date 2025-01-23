# metaTP
metaTP: a pipeline for analyzing meta-transcriptome.metaTP is a pipeline that integrated bioinformatics tools for analyzing metatranscriptomic data comprehensively.It includes quality control, non-coding RNA removal, transcript expression quantification, differential gene expression analysis, functional annotation, co-expression network analysis.
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
## Execute
Dry run: Use --dry-run to see what tasks Snakemake will perform without actually running them:
```Python 
snakemake --dry-run
```
Dry run: Use --dry-run to see what tasks Snakemake will perform without actually running them:
```Python
snakemake -j 4 
```
#### 1. Download the sra sequence according to the ACC number
```Python
snakemake prefetch_sra2fastq
```
#### 2. Sequence quality test
```Python
snakemake QC_test
```
#### 3. Sequence quality control, rmrRNA contig cds
```Python
snakemake QC_rmrRNA_contigs_cds
```
#### 4. transcript_index
```Python
snakemake transcript_index
```
#### 5. gene_expression_quant
```Python
snakemake gene_expression_quant
```
#### 6. DEG_analysis
```Python
snakemake DEG_analysis
snakemake up_regulated_gene
snakemake down_regulated_gen
```
#### 7. emapper.py
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
(1)Dimensionality Reduction Analysis
input: 
	•	transcript_abundance_quantification_table_filter2.csv:Transcript abundance table containing expression values ​​for each sample.br>
output:
	•	降维后的数据可视化结果，主坐标分析（PCoA）图，显示根际样品和不同剂量的本体样品之间存在显著分离，表明在不同土壤环境中存在显著的转录差异。br>
(2)Venn analysis
input:
	•	venn_setdiff_rhi.emapper.annotations.csv：通过eggNOG注释的基因功能数据，包括基因ID及其功能分类信息。br>
output:
	•	venn_setdiff_bulk_go.pdf：展示本体样本和根际样本共享和特有基因的维恩图。
(3)Functional Annotation br>
input：
	•	final_table_sequence.emapper.annotations.csv：从eggNOG注释工具生成的功能注释表，包含基因ID、功能分类、GO注释、KEGG路径等。br>
	•	cog_funclass.tab：COG功能分类表，定义了各COG分类的含义。br>
output：
	•	cog.pdf：COG功能分类的柱状图，显示每个功能类别中包含的基因数量。br>
	•	功能注释的中间结果表（如gene2cog），记录了基因到COG功能的映射。br>
 (4)Functional Enrichment Analysis
 input：
	•	differential_genes_id_up.txt：差异表达基因的ID列表。br>
	•	go_term.list：GO术语表，定义每个GO术语的功能和分类。br>
	•	final_table_sequence.emapper.annotations.csv：基因的GO注释数据。br>
	•	ko00001.json：从KEGG下载的路径注释信息（用于KO富集分析）。br>
output：
	•	go_rich.txt：GO功能富集分析结果表，记录了显著GO分类及对应基因数。br>
	•	go_rich_bar.pdf 和 go_rich_dot.pdf：GO功能富集的柱状图和点图。br>
	•	Ko_rich.txt：KEGG功能富集结果表，记录了显著的KEGG路径及对应基因。br>
	•	ko_rich_bar.pdf 和 ko_rich_dot.pdf：KEGG功能富集的柱状图和点图。br>


#### 1. Gene co-expression network
```Python
  all codes in Gene_co-expression_network.zip
```
