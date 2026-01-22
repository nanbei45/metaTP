# <metaTP>
   
# metaTP

metaTP: a pipeline for analyzing meta-transcriptome.metaTP is a pipeline that integrated bioinformatics tools for analyzing metatranscriptomic data comprehensively.It includes quality control, non-coding RNA removal, transcript expression quantification, differential gene expression analysis, functional annotation, co-expression network analysis. 

## **Prerequisites**
### **Install conda software**
### **Downloaded the conda environment in our OBS** 

### **This is [OBS](https://e-share.obs-website.cn-north-1.myhuaweicloud.com?token=AoEwdFl40jDPD8Kz5t/mbW0buwwa1OGzg1o18CZSpn+B7FoikdQ+uXQlXfu7v8RnT4JnRsEOHvDXFbuJ468nCwZ154OvWGbP2iMZWd26IgnN74cKjlF1xiBCGStN/o/g1jiY5fKU2tv9Gv0JWZ9CHTqMBMZyQ953BsGzXZPkBAnCkdxKb3j7Mj29wcH87geAyJRIhGJh5/w26qITd1LtdNxxXbn31x97Bo0bmIfd5Zl3CwwDUATfBc6t5zZcQ5zpUCuxw7yDWEviAURWlVFCuIV/++d0r106aFmf+JrnXh9zRwVbDY+1z0ARAz0t72WgZTu+jyt67mG0nIVAb/3Brunmpc8D+zXlWc4ed4+Z9TVnHYc6fWMZOJLLly5pwMXLqiJi+fVqQb0alqNjFW+j4HksD5JCxBmGTIDziJRa217n+0V2oqhIBGril2cGe+MYEJr3MfDdADDJ2BGr5gPYg0aMdLVHN/GkidLMBzY0AZCyEtrmEii57ptiHHqGV4t8/Z5OIBUIStn4FzhUlzOclA==) link**
### **password:**
### **123456**

### **Codes also can be download from [my google drive](https://drive.google.com/drive/folders/1OMdrzCSU_e81C3xLBAPkKDp0LlM4OjTc?usp=sharing)**

      get -r metaTP.tar.tgz
      tar zxvf  metaTP.tar.tgz -C ./   
      conda create -n metaTP
      tar -xzvf metaTP.tar.gz -C /home/xxx/anaconda3/envs/metaTP
      conda activate metaTP
   
### **The codes can be downloaded from our OBS**
      
      get -r metaTP_pipeline.tar.tgz
      tar -xzvf metaTP_pipeline.tar.tgz
      cd metaTP_pipeline

### **1. Download the sra sequence according to the ACC number**

     python 1.prefetch_sra2fastq.py -i SRR_Acc_List.txt -o test_sra_data  

       # -i (SRR_Acc_List.txt)  
       # -o (output_dir)   
       # output_dir：test_sra_data； test_sra_data/fastq  
 
### **2. Sequence quality test**

      python 2.QC_test.py -i ./test_sra_data/fastq -o ./test_sra_data/QC_before_result  
      
      # -i (fastq files)  
      # -o (output_dir)  
      # output_dir: test_sra_data/QC_before_result  

### **3. Sequence quality control, rmrRNA contig cds**

       python 3.QC_rmrRNA_contigs_cds.py -i ./test_sra_data/fastq -o ./test_sra_data  
   
      #Required_dir: trimmomatic_adapters, rRNA_databases 
      #output_dir：test_sra_data/QC_control;test_sra_data/rmrRNA;test_sra_data/magahit  
      # all_longest_orfs_cds_rmdup.fasta  
      #-i (input folder contain fastq )  
      #-threads (thread setting,default=24)  
      #-phred （you can choose "-phred33 or -phred64",default="-phred33"）  
      #-ILLUMINACLIP （paths for adapters,default="trimmomatic_adapters/TruSeq3"）  
      #-LEADING （Remove front-end low-quality sequences or N bases,",type =int,default=5）  
      #-SLIDINGWINDOW (SLIDINGWINDOW:4:15, scan sequences using a 4-base width sliding window and trim when the average mass per base is below 15",
                      default="4:15")  
      #-MINLEN (MINLEN:25, Remove the sequences that less than 25 bases",type =int,
                      default=25)  
      # -m minimum protein length (default: 100)  
   
      #-o (the output folder)  
	
	'error 
	/xxxx/xxxx/bin/megahit_core: /xxx/metaTP/bin/../lib/libstdc++.so.6: version `CXXABI_1.3.11' not found (required by /xxx/anaconda3/envs/metaTP2.0/bin/megahit_core)  
    /xxxxx/xxxxx/anaconda3/envs/metaTP/bin/megahit_core: /xxxx/metaTP/bin/../lib/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /xxxx/anaconda3/envs/metaTP2.0/bin/megahit_core)  
    /xxxxx/xxxxx/anaconda3/envs/metaTP/bin/megahit_core: /xxx/envs/metaTP/bin/../lib/libstdc++.so.6: version `GLIBCXX_3.4.22' not found (required by /xxxx/anaconda3/envs/metaTP2.0/bin/megahit_core)  
    /xxxxx/xxxx/anaconda3/envs/metaTP/bin/megahit_core: /xxxx/anaconda3/envs/metaTP/bin/../lib/libstdc++.so.6: version `CXXABI_1.3.13' not found (required by /xxxxx/anaconda3/envs/metaTP2.0/bin/megahit_core)  

	conda install -c anaconda libstdcxx-ng  
	delete anaconda3/envs/metaTP/lib/libstdc++.so， Keep one version 'libstdc++.so.6.xxx.xxx'  

### **4. transcript_index**

	python 4.transcript_index.py -i ./test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta -o ./test_sra_data    
	
        # output_dir：transcripts_index    
	# -i (input folder contain fastq )    
	# -k (kmer setting,default=31)   
	# -p（threads,default=24）  
	# -o (the output folder)   
	# salmon index -t transcripts.fa -i transcripts_index --type quasi -k 31  
	# salmon index -t ./TransDecoder/all_final_cds_rmdup.fa -i ./salmon_index/transcripts_index  
	# salmon index -t ./TransDecoder/all_final_cds_rmdup.fa -i ./TransDecoder/salmon_index/transcripts_index -p 24  
	# salmon index -t ./all_final.contigs_cds_rmdup_id.fa -i ./TransDecoder/salmon_index/transcripts_index  -p 24 -k 31  
	# wget https://anaconda.org/bioconda/salmon/1.1.0/download/linux-64/salmon-1.1.0-hf69c8f4_0.tar.bz2  
	# conda install --use-local salmon-1.1.0-hf69c8f4_0.tar.bz2  

### **5. gene_expression_quant**

	python 5.gene_expression_quant.py -i ./test_sra_data/rmrRNA -index ./test_sra_data/transcripts_index -p 24 -o ./test_sra_data  
   
	#output_dir：./test_sra_data/transcripts_quant/transcript_abundance_quantification_table.csv  
	#-i (input folder containing fastq)   
	#-index (transcripts_index obtained in the previous step)  
	#-p (threads,default=24)  
	#-o (the output folder)  
   
	# fastq后缀文件   
	# 对双端测序数据reads表达量的估计    
	# salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant    
	# 对单端测序数据reads表达量的估计    
	# salmon quant -i transcripts_index -l <LIBTYPE> -r reads.fq -o transcripts_quant    
	### 命令quant均适用于这两个index(quasi-mapping or SMEM-based)，此外，Salmon能够自动检测到使用的是哪种index，从而采用与之匹配的比对方法。    
	### 注意：参数-l必须指定在参数-1，-2和-r的前面。   
	### 生成一个目录，内含文件quant.sf   

	可选步骤过滤低表达的基因初步过滤低表达基因与保存counts数据   
 	由于transcript_abundance_quantification_table.csv数据量比较大，所以选择过滤掉样本序列数（设定一定阈值进行过滤）   
  
        # seqkit grep --pattern-file filter_id.txt ./test_sra_data2/megahit/all_longest_orfs_cds_rmdup_id.fasta > ./test_sra_data2/final_table.fasta   

### **5.1 filter_gene.R**

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

### **6. DEG_analysis**

	python 6.DEG_analysis.py -i ./test_sra_data/transcripts_quant/transcript_abundance_quantification_table_filter.csv -g ./sample_group.csv -n rhizosphere/bulk -s /home/mne/metaTP/test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta -o ./test_sra_data/DEG_result0.05 -p 0.05 -f 1  
	python 6.up_regulated_gene.py -i ./test_sra_data/DEG_result0.05 -s /home/mne/metaTP/test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta  
	python 6.down_regulated_gene.py -i ./test_sra_data/DEG_result0.05 -s /home/mne/metaTP/test_sra_data/megahit/all_longest_orfs_cds_rmdup_id.fasta      
   
	# Required file: Rscript DEG_analysis.r;   
	# -i transcript_abundance_quantification_table.csv # Gene expression matrix based on TPM from Step 5；  
	# -g Required file: sample_group.csv(The first column is the sample name,column name is "sample"; the second column is the group name; colum name is "group")   
	# -n grouping information; The name must be the same with sample_group.csv  
	# -s all_final.contigs_cds_rmdup_id.fa; the cds sequences from Step 3 #absolute path  
	# -o output directory  
	# -p cut_off_pvalue  
	# -f cut_off_logFC  
 
### **7. emapper.py**

	emapper.py -m diamond -i ./test_sra_data/DEG_result0.05/differential_gene_sequence_up.fasta --itype CDS --translate --cpu 20 --data_dir /home/mne/metaTP/eggnog-mapper_database --dmnd_db /home/mne/metaTP/eggnog-mapper_database/eggnog_proteins.dmnd --output_dir test_sra_data -o differential_gene_sequence_up --block_size 0.4 --override

	# /home/wys/anaconda3/envs/metaTP/bin/download_eggnog_data.py
	# wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
	# wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
	# wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz

	# 如果内存比较小可以设置小一些的block_size (--block_size 0.4),默认值是block_size 为2
	# wget https://anaconda.org/bioconda/eggnog-mapper/2.1.12/download/noarch/eggnog-mapper-2.1.12-pyhdfd78af_0.tar.bz2
	# conda install --use-local eggnog-mapper-2.1.12-pyhdfd78af_0.tar.bz2

### **8. Downstream analysis.R**

      Dimensionality reduction analysis
      Venn analysis
      Functinal annotation
      Functional enrichment analysis

### **9. Gene co-expression network**

      all codes in Gene_co-expression_network.zip
 
## **If you use metaTP, please cite:** 
    (1) metaTP: a pipeline for analyzing meta-transcriptome
    (2) Yansu Wang, Murong Zhou, Quan Zou, Lei Xu. Machine learning for phytopathology: from the molecular scale towards the network scale. Briefings in Bioinformatics. 2021, Doi: 10.1093/bib/bbab037
   
