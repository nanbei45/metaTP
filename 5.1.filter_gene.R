#### 初步过滤低表达基因 ####
##过滤掉每组中至少一半的样品大于1的行

counts <- read.csv("./test_sra_data/transcripts_quant/transcript_abundance_quantification_table.csv",header=T,row.names=1)

cutoff = 0.5
num =1 
group1 <- data.frame(counts[,1:4][which(apply(counts[,1:4], 1, function(x){length(which (x>=1))/length(x)}) > cutoff),])
group2 <- data.frame(counts[,5:8][which(apply(counts[,5:8], 1, function(x){length(which (x>=1))/length(x)}) > cutoff),])

gene_filter1 <- counts[union(rownames(group1),rownames(group2)),]

#gene_filter <- data.frame(counts[which(apply(counts, 1, function(x){length(which (x>=1))/length(x)}) >= 0.5),])

write.csv(data.frame(GeneID=rownames(gene_filter1),gene_filter1),"./test_sra_data/transcripts_quant/transcript_abundance_quantification_table_filter.csv",row.names=F)