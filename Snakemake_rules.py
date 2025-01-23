# 设置输出目录
output_dir = "./test_sra_data"

# 规则 1: 运行 prefetch_sra2fastq.py
rule prefetch_sra2fastq:
    input:
        sra_list="SRR_Acc_List.txt"
    output:
        fastq_dir=output_dir + "/fastq"
    shell:
        "python 1.prefetch_sra2fastq.py -i {input.sra_list} -o {output.fastq_dir}"

# 规则 2: 运行 QC_test.py
rule QC_test:
    input:
        fastq_dir=output_dir + "/fastq"
    output:
        qc_result=output_dir + "/QC_before_result"
    shell:
        "python 2.QC_test.py -i {input.fastq_dir} -o {output.qc_result}"

# 规则 3: 运行 QC_rmrRNA_contigs_cds.py
rule QC_rmrRNA_contigs_cds:
    input:
        fastq_dir=output_dir + "/fastq"
    output:
        result_dir=output_dir
    shell:
        "python 3.QC_rmrRNA_contigs_cds.py -i {input.fastq_dir} -o {output.result_dir}"

# 规则 4: 运行 transcript_index.py
rule transcript_index:
    input:
        fasta_file=output_dir + "/megahit/all_longest_orfs_cds_rmdup_id.fasta"
    output:
        index_dir=output_dir
    shell:
        "python 4.transcript_index.py -i {input.fasta_file} -o {output.index_dir}"

# 规则 5: 运行 gene_expression_quant.py
rule gene_expression_quant:
    input:
        fastq_dir=output_dir + "/rmrRNA",
        index_dir=output_dir + "/transcripts_index"
    output:
        quant_file=output_dir + "/transcripts_quant/transcript_abundance_quantification_table_filter.csv"
    params:
        p=24
    shell:
        "python 5.gene_expression_quant.py -i {input.fastq_dir} -index {input.index_dir} -p {params.p} -o {output.quant_file}"

# 规则 6: 运行 DEG_analysis.py
rule DEG_analysis:
    input:
        quant_file=output_dir + "/transcripts_quant/transcript_abundance_quantification_table_filter.csv",
        sample_group="./sample_group.csv"
    output:
        deg_result=output_dir + "/DEG_result0.05"
    params:
        pvalue=0.05, fold_change=1
    shell:
        "python 6.DEG_analysis.py -i {input.quant_file} -g {input.sample_group} -n rhizosphere/bulk -s {output.deg_result} -o {output.deg_result} -p {params.pvalue} -f {params.fold_change}"

# 规则 7: 运行 up_regulated_gene.py
rule up_regulated_gene:
    input:
        deg_result=output_dir + "/DEG_result0.05"
    output:
        up_genes=output_dir + "/DEG_result0.05/differential_gene_sequence_up.fasta"
    shell:
        "python 6.up_regulated_gene.py -i {input.deg_result} -s {output.up_genes}"

# 规则 8: 运行 down_regulated_gene.py
rule down_regulated_gene:
    input:
        deg_result=output_dir + "/DEG_result0.05"
    output:
        down_genes=output_dir + "/DEG_result0.05/differential_gene_sequence_down.fasta"
    shell:
        "python 6.down_regulated_gene.py -i {input.deg_result} -s {output.down_genes}"

# 规则 9: 运行 emapper.py
rule emapper:
    input:
        up_genes=output_dir + "/DEG_result0.05/differential_gene_sequence_up.fasta"
    output:
        emapper_output=output_dir + "/emapper_result"
    params:
        cpu=20,
        eggnog_db="/home/mne/metaTP/eggnog-mapper_database",
        dmnd_db="/home/mne/metaTP/eggnog-mapper_database/eggnog_proteins.dmnd"
    shell:
        "emapper.py -m diamond -i {input.up_genes} --itype CDS --translate --cpu {params.cpu} --data_dir {params.eggnog_db} --dmnd_db {params.dmnd_db} --output_dir {output.emapper_output} -o differential_gene_sequence_up --block_size 0.4 --override"