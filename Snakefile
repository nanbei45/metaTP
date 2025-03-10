import pandas as pd

# Loading a Configuration File
configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samplesheet"])
SAMPLES = samples_df["sample"].tolist()
SRAS = samples_df["sra"].tolist()

# Defining dynamic output paths
output_dir = config["output_dir"]
fastq_dir = config["fastq_dir"]
qc_dir = config["qc_dir"]
quant_dir = config["quant_dir"]
deg_dir = config["deg_dir"]
emapper_dir = config["emapper_dir"]

#Global rules: define all target outputs
rule all:
    input:
        expand(f"{deg_dir}/{{sample}}_deg_results", sample=SAMPLES),
        expand(f"{emapper_dir}/{{sample}}_emapper", sample=SAMPLES)

# Rule 1: Download SRA data and convert to FASTQ
rule prefetch_sra2fastq:
    output:
        fastq = f"{fastq_dir}/{{sra}}.fastq"
    params:
        sra_id = lambda wildcards: samples_df[samples_df["sra"] == wildcards.sra]["sra"].iloc[0]
    shell:
        """
        prefetch {params.sra_id} && \
        fastq-dump {params.sra_id} \
            --outdir {fastq_dir} \
            --gzip
        """

# Rule 2: Quality Control
rule QC_test:
    input:
        fastq = f"{fastq_dir}/{{sra}}.fastq.gz"
    output:
        qc_report = directory(f"{qc_dir}/{{sra}}_qc")
    shell:
        "python scripts/2.QC_test.py -i {input.fastq} -o {output.qc_report}"

# Rule 3: Remove rRNA and assemble contigs
rule QC_rmrRNA_contigs_cds:
    input:
        fastq = f"{qc_dir}/{{sra}}_qc/clean.fastq.gz"
    output:
        rmrna = directory(f"{qc_dir}/{{sra}}_rmrRNA"),
        contigs = directory(f"{qc_dir}/{{sra}}_contigs")
    shell:
        "python scripts/3.QC_rmrRNA_contigs_cds.py -i {input.fastq} -o {output.rmrna} {output.contigs}"

# Rule 4: Transcript Indexing
rule transcript_index:
    input:
        contigs = f"{qc_dir}/{{sra}}_contigs/contigs.fasta"
    output:
        index = directory(f"{quant_dir}/{{sra}}_index")
    shell:
        "salmon index -t {input.contigs} -i {output.index}"

# Rule 5: Quantification of gene expression
rule gene_expression_quant:
    input:
        fastq = f"{qc_dir}/{{sra}}_rmrRNA/clean.fastq.gz",
        index = f"{quant_dir}/{{sra}}_index"
    output:
        quant = directory(f"{quant_dir}/{{sra}}_quant")
    params:
        threads = config["threads"]
    shell:
        "salmon quant -i {input.index} -l A -r {input.fastq} -o {output.quant} --threads {params.threads}"

# Rule 6: Differential Expression Analysis
rule DEG_analysis:
    input:
        quant_dirs = expand(f"{quant_dir}/{{sra}}_quant", sra=SRAS),
        groups = lambda wildcards: samples_df[samples_df["sample"] == wildcards.sample]["group"].iloc[0]
    output:
        deg_results = directory(f"{deg_dir}/{{sample}}_deg_results")
    params:
        pvalue = config["pvalue"],
        fold_change = config["fold_change"]
    shell:
        "python scripts/6.DEG_analysis.py -i {input.quant_dirs} -g {input.groups} -o {output.deg_results} -p {params.pvalue} -f {params.fold_change}"

# Rule 7: Functional Annotation (emapper)
rule emapper:
    input:
        deg_genes = f"{deg_dir}/{{sample}}_deg_results/upregulated_genes.fasta"
    output:
        emapper_out = directory(f"{emapper_dir}/{{sample}}_emapper")
    params:
        eggnog_db = config["eggnog_db"],
        threads = config["threads"]
    shell:
        "emapper.py -i {input.deg_genes} -o {output.emapper_out} --data_dir {params.eggnog_db} --cpu {params.threads} -m diamond"
