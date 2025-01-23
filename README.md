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
