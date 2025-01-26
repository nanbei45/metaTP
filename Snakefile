rule all:
    input:
        "test_sra_data2/differential_gene_sequence_up",
        "test_sra_data2/up_regulated_gene",
        "test_sra_data2/down_regulated_gene"

# 规则 1: 下载并转换 SRA 到 FASTQ
rule prefetch_sra2fastq:
    input:
        sra_list="SRR_Acc_List.txt"
    output:
        fastq_list=directory("test_sra_data2/fastq")
    shell:
        "python 1.prefetch_sra2fastq.py -i {input} -o {output}"

# 规则 2: 质量控制
rule QC_test:
    input:
        fastq_list=rules.prefetch_sra2fastq.output.fastq_list # 引用规则1的输出
    output:
        qc_result=directory("test_sra_data2/QC_before_result")

    shell:
        "python 2.QC_test.py -i {input} -o {output}"

# 规则 3: 运行 QC_rmrRNA_contigs_cds.py
rule QC_rmrRNA_contigs_cds:
    input:
        fastq_dir=rules.prefetch_sra2fastq.output.fastq_list
    output:
        QC_control=directory("test_sra_data2/QC_control"),
        rmrRNA = directory("test_sra_data2/rmrRNA"),
        megahit = directory("test_sra_data2/megahit")

    shell:
        "python 3.QC_rmrRNA_contigs_cds.py -i {input.fastq_dir} -o {output.QC_control} {output.rmrRNA} {output.megahit}"

# 规则 4: 运行 transcript_index.py
rule transcript_index:
    input:
        fasta_file=rules.QC_rmrRNA_contigs_cds.output.megahit
    output:
        index_dir=directory("test_sra_data2/transcripts_index"),
    shell:
        "python 4.transcript_index.py -i {input.fasta_file} -o {output.index_dir}"

# 规则 5: 运行 gene_expression_quant.py
rule gene_expression_quant:
    input:
        fastq_dir=rules.QC_rmrRNA_contigs_cds.output.rmrRNA,
        index_dir=rules.transcript_index.output.index_dir
    output:
        quant_file= directory("test_sra_data2/transcripts_quant")
    params:
        p=24
    shell:
        "python 5.gene_expression_quant.py -i {input.fastq_dir} - index {input.index_dir} -p {params.p} -o {output.quant_file}"

# 规则 6: 运行 DEG_analysis.py
rule DEG_analysis:
    input:
        quant_file=rules.gene_expression_quant.output,
        sample_group = "sample_group.txt",
        bulk = "rhizosphere/bulk",
        megahit = rules.QC_rmrRNA_contigs_cds.output.megahit
    output:
        deg_result= directory("test_sra_data2/DEG_result0.05")
    params:
        pvalue=0.05, fold_change=1
    shell:
        "python 6.DEG_analysis.py -i {input.quant_file} -g {input.sample_group} -n {input.bulk} -s {input.megahit} -o {output.deg_result} -p {params.pvalue} -f {params.fold_change}"

# 规则 7: 运行 up_regulated_gene.py
rule up_regulated_gene:
    input:
        deg_result=rules.DEG_analysis.output.deg_result,
        megahit= rules.QC_rmrRNA_contigs_cds.output.megahit
    output:
        up_genes=directory("test_sra_data2/up_regulated_gene")
    shell:
        "python 6.up_regulated_gene.py -i {input.deg_result} -s {input.megahit} -o {output.up_genes}"

# 规则 8: 运行 down_regulated_gene.py
rule down_regulated_gene:
    input:
        deg_result=rules.DEG_analysis.output.deg_result,
        megahit= rules.QC_rmrRNA_contigs_cds.output.megahit
    output:
        down_genes=directory("test_sra_data2/down_regulated_gene")
    shell:
        "python 6.down_regulated_gene.py -i {input.deg_result} -s {input.megahit} -o {output.down_genes}"

# 规则 9: 运行 emapper.py
rule emapper:
    input:
        up_gnene = rules.DEG_analysis.output.deg_result
    output:
        emapper_output=directory("test_sra_data2/differential_gene_sequence_up")
    params:
        cpu=20,
        eggnog_db="/home/mne/metaTP/eggnog-mapper_database",
        dmnd_db="/home/mne/metaTP/eggnog-mapper_database/eggnog_proteins.dmnd"
    shell:
        "emapper.py -m diamond -i {input.up_genes} --itype CDS --translate --cpu {params.cpu} --data_dir {params.eggnog_db} --dmnd_db {params.dmnd_db} -o {output.emapper_output} --block_size 0.4 --override"



