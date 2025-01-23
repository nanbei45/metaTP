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


#### 1. Gene co-expression network
```Python
  all codes in Gene_co-expression_network.zip
```
