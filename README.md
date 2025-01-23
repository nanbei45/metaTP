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
```Python
snakemake -j 4  #Use the -j parameter to set the number of tasks to execute in parallel. For example, to use 4 parallel tasks:
snakemake --dry-run #Dry run: Use --dry-run to see what tasks Snakemake will perform without actually running them:
```
### 1. Download the sra sequence according to the ACC number
```Python

```
### 2. Sequence quality test
```Python

```
### 3. Sequence quality control, rmrRNA contig cds
```Python

```
### 4. transcript_index
```Python

```
### 5. gene_expression_quant
```Python

```
### 6. DEG_analysis
```Python

```
### 7. emapper.py
```Python

```
## Analyze
