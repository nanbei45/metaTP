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
The metaTP pipeline integrates data download options using the SRA toolkit.
```Bash
snakemake prefetch_sra2fastq
```
input: SRR_Acc_List.txt<br>
output: ./test_sra_data/fastq<br>
#### 2. Sequence quality test
Assess the quality of FASTQ files using FastQC.<br>
```Python
snakemake QC_test
```
input: ./test_sra_data/fastq<br>
output: ./test_sra_data/QC_before_result<br>
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
input: ./test_sra_data/fastq<br>
output: ./test_sra_data/QC_control;./test_sra_data/rmrRNA;./test_sra_data/magahit;<br>
#### 4. transcript_index
Coding sequences obtained from Transdecoder were used to build the index, create a decoy-aware transcriptome file for Salmon transcript quantification, and then create a transcriptome index. Gene expression levels were normalized to transcript length and library size (TPM).
```Python
snakemake transcript_index
```
input: ./test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta<br>
output: ./test_sra_data/transcripts_index<br>

#### 5. gene_expression_quant
```Python
snakemake gene_expression_quant
```
input: ./test_sra_data/rmrRNA、./test_sra_data/transcripts_index<br>
output: ./test_sra_data/transcripts_quant/transcript_abundance_quantification_table.csv<br>
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
input:transcript_abundance_quantification_table_filter.csv<br>
output:transcript_abundance_quantification_table_filter2.csv<br>
#### 6. DEG_analysis
```Python
snakemake DEG_analysis
snakemake up_regulated_gene
snakemake down_regulated_gen
```
(1)snakemake DEG_analysis
input:./test_sra_data/transcripts_quant/transcript_abundance_quantification_table_filter.csv<br>
output:./test_sra_data/DEG_result0.05<br>
(2)snakemake up_regulated_gene
input:./test_sra_data/DEG_result0.05<br>
output:/home/mne/metaTP/test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta <br>
(3)snakemake down_regulated_gen
input:./test_sra_data/DEG_result0.05<br>
output:/home/mne/metaTP/test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta<br>
#### 7. emapper.py
The integrated eggNOG -mapper provides several key features including: 1) de novo gene prediction based on raw alignments, 2) integrated pairwise homology prediction, and 3) rapid protein domain detection.
```Python
snakemake emapper
```
input:./test_sra_data/DEG_result0.05/differential_gene_sequence_up.fasta<br>
output:./test_sra_data/differential_gene_sequence_up<br>
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
