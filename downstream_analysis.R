###############################################################
##################过滤掉Streptophyta##########################
############ df[,-which(names(df)%in%c("a","b")]##############
##############################################################
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


#########################################################
       ##################################
#######Dimensionality reduction analysis#################
       ###################################

library(vegan)
library(ape)
library(ggplot2)
library(dplyr)

gene_dist <- vegdist(t(table_filter2), "bray")
##分组信息，第一列样本名称，第二列
group<- read.csv("sample_group.csv",header = T,row.names = 1)

res <- pcoa(gene_dist)
biplot(res)
res2 <- res[c('values','vectors')]
vectors<- res2$vectors
x_value<- res2$values[1,2]
#> x_value
#[1] 0.6043274
y_value<- res2$values[2,2]
#> y_value
#[1] 0.1536484

###合并vectors 与group

vectors_df <- data.frame(sampleID=rownames(vectors),vectors)
group_df <- data.frame(sampleID=rownames(group),group)
vectors_group<-left_join(vectors_df, group_df, by = "sampleID")

#write.csv(res$values,"effectome_unifrac_pcoa_values.csv")
#write.csv(res$vectors,"effectome_unifrac_pcoa_vectors.csv")
#pcoa<-read.csv("secretome_unifrac_pcoa_vectors.csv",header=T,row.names = 1)

p<- ggplot(vectors_group, aes(Axis.1,Axis.2)) + 
  geom_point(aes(colour = factor(group)), size = 5) + 
  xlab("PCoA1(60.43%)") + 
  ylab("PCoA2(15.36%)") +
  geom_vline(xintercept=0, linetype="dotted",size=1.2)+
  geom_hline(yintercept=0, linetype="dotted",size=1.2)+
  theme_set(theme_bw())+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",linewidth=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=20),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=20))

ggsave(p, file="gene_pcoa.pdf", width=9, height=7,device=cairo_pdf)

#组间差异分析
library(vegan)
group1<- vectors_group$group
anosim(t(table_filter2),group1,permutations = 999, distance = "bray") 

#########################################################
######################################
###########   Venn   #################
#####################################
# The script was downloaded from https://pan.baidu.com/s/10h8VqkFE_1Z_60GxsPSX0w#list/path=%2F
#venn
library (VennDiagram)
#table_filter2_t <- read.delim("transcript_abundance_quantification_table_filter2.txt", header = T, row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
group <- read.csv("sample_group.csv", header = T, stringsAsFactors = F)
#View(group)
group_id <- unique(group[,2])
names(group) <- c('sample', 'group')
#View(group_id)
#View(group)
result <- NULL
len = 0
for (i in group_id) {
  group_i <- subset(group, group == i)
  table_filter2_i <- table_filter2[group_i$sample]
  table_filter2_i <- subset(table_filter2_i, rowSums(table_filter2_i) > 0)
  table_filter2_select <- c(i, rownames(table_filter2_i))
  if (length(table_filter2_select) > len) len <- length(table_filter2_select)
  length(table_filter2_select) = nrow(table_filter2)
  result <- cbind(result, table_filter2_select)
}
result <- result[1:len, ]
write.table(result, "venn_group.txt", sep = '\t', col.names = F,row.names = F, na = '', quote = F)
venn_dat <- read.delim('venn_group.txt', header = T, sep = '\t', stringsAsFactors = F, check.names = F)
venn_list <- list(venn_dat[,1], venn_dat[,2])
names(venn_list) <- colnames(venn_dat)
write.csv(venn_list,"venn_list.csv")
#作图
library(venn)
pdf("venn_7.pdf",width=9,height=7)
#png('venn_7.pdf', width = 1500, height = 1500, res = 200, units = 'px')
venn(venn_list,zcolor='style')
dev.off()

##############################################
##########  functinal annotation   ###############
##############################################
###特异基因功能#############
#加载dplyr包
library(dplyr)
#直接利用dplyr包里面的intersect函数对数据框取交集
result2=intersect(venn_list$rhizosphere,venn_list$bulk)
#保存交集结果
write.table(file="venn_intersect.txt",result2,quote=F,row.names = F,sep="\t")
##补集 setdiff(a,b) a与b的补集：setdiff(b,a)
result3=setdiff(venn_list$rhizosphere,venn_list$bulk)
write.table(file="venn_setdiff_rhi.txt",result3,quote=F,row.names = F,sep="\t")
result4=setdiff(venn_list$bulk,venn_list$rhizosphere)
write.table(file="venn_setdiff_bulk.txt",result4,quote=F,row.names = F,sep="\t")

####根据差异基因提取功能注释信息################

emapper.annotations_up <- read.csv("differential_gene_sequence_up.emapper.annotations.csv",header = T)

###上调基因###

table_filter_up<- table_filter2[which(row.names(table_filter2) %in% emapper.annotations_up["X.query"]),]

emapper.annotations_up_df <- data.frame(emapper.annotations_up)
table_filter2_df <- data.frame(geneID =rownames(table_filter2),table_filter2)
table_filter_up<-left_join(table_filter2_df, emapper.annotations_up_df, by =c("geneID" = "X.query"))
table_filter_up<-table_filter_up[complete.cases(table_filter_up),]

View(table_filter_up)

write.csv(table_filter_up,"table_filter_up.csv")

merge_otu <- function(data, formula, e_species=NULL, FUN){
  
  require(dplyr)
  
  df <- as.data.frame(data)
  
  if(!is.null(e_species)){
    
    if(is.character( data[,1])){
      rownames(df) <- df[, 1]
      df <- df[,-1]
      df[is.na(df)] <- 0
      df <- df %>% select( -e_species )
      result <- aggregate(formula , data = df, FUN)
      return(result)
    }
    
    if(! is.character(data[,-1])){
      df[is.na(df)] <- 0
      df <- df %>% select( -e_species )
      result <- aggregate(formula , data = df, FUN)
      return(result)
    }
  }
  
  if(is.null(e_species)){
    
    if(is.character( data[,1])){
      rownames(df) <- df[, 1]
      df <- df[,-1]
      df[is.na(df)] <- 0
      result <- aggregate(formula , data = df, FUN)
      return(result)
    }
    
    if(! is.character(data[,-1])){
      df[is.na(df)] <- 0
      result <- aggregate(formula , data = df, FUN)
      return(result)
    }
  }
}

#table_filter_rhi_sum<- aggregate(table_filter_rhi,by=list(table_filter_rhi$Description),sum,na.rm=TRUE)
rownames(table_filter_up) <- table_filter_up[,1]
table_filter_up <- table_filter_up[,-1]
table_filter_up <- table_filter_up[,c(1,2,3,4,5,6,7,8,15)]
table_filter_up_sum <- merge_otu(data =table_filter_up, formula =. ~Description, e_species = NULL, sum)
rownames(table_filter_up_sum) <- table_filter_up_sum[,1]
table_filter_up_sum <- table_filter_up_sum[,-1]
table_filter_up_sum <- table_filter_up_sum[-1,]
table_filter_up_sum <- table_filter_up_sum[which(rowSums(table_filter_up_sum) >= 50), ]   
write.csv(table_filter_up_sum,"table_filter_up_sum.csv")
table_filter_up_sum <- read.csv("table_filter_up_sum.csv",header = T,row.names = 1)

library(pheatmap)
library(RColorBrewer)
data <- as.data.frame(table_filter_up_sum)
#x  <- as.matrix(dat_heat)
x1<-scale(data,center=FALSE,scale=T)

par(no.readonly=TRUE)
par(mgp=c(1.6,0.6,0),mar=c(3,3,12,1))
pheatmap(x1,color = colorRampPalette(c("blue", "yellow", "red"))(50),
         border_color = "black",
         cluster_row = FALSE, cluster_col = FALSE,
         scale="none",cellwidth=20,cellheight=20,
         breaks =NA,
         fontsize_row=12,
         fontsize_col=12,
         angle_col="315",
         filename = "DEG_up_heatmap.pdf")

###################################################
#####################################################
###对于无参数据可以自己搭建数据库功能富集分析#########
####################################################
#####################################################

##根据emapper.annotations提取注释结果######
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationForge")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("wiggleplotr")

library(tidyverse)
library(stringr)
library(AnnotationForge)
library(formattable)
library(clusterProfiler)
library("wiggleplotr")
library("dplyr")

emapper <- read.delim("final_table_sequence.emapper.annotations.csv",sep=",") %>%
  dplyr::select(GID=X.query,Gene_Symbol=Preferred_name, 
                GO=GOs,KO=KEGG_ko,Pathway =KEGG_Pathway, 
                OG=eggNOG_OGs,COG = COG_category,Gene_Name =seed_ortholog)

#####提取eggnog信息######
#<<<<<<<<<<<<<<<<<<<<获得gene2eggnog<<<<<<<<<<<<<<<<<<<<<<<<
gene_info <- dplyr::select(emapper,GID,Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))

#<<<<<<<<<<<<<<<<<<<<获得gene2go<<<<<<<<<<<<<<<<<<<<<<<<

gene2go <- dplyr::select(emapper,GID,GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  dplyr::filter(GO!="",!is.na(GO))%>%   #这是只提取有GO注释信息的行，判断的标准时GO信息不是NA，这也就是为什么前面必须将“-”替换为NA，不替换就无法进行有效过滤。
  mutate(EVIDENCE = 'A')     #硬生生加了1列EVIDENCE，全部赋值A,凑数的。
dim(gene2go)    #查看数据维度。
gene2go<- gene2go[-which(gene2go$GO=="-"),]

#<<<<<<<<<<<<<<<<<<<<获得gene2cog<<<<<<<<<<<<<<<<<<<<<<<<

cog_info <- read_delim("cog_funclass.tab",
                       "\t", escape_double = FALSE, trim_ws = TRUE)
cogs <- emapper %>%
  dplyr::select(GID, COG = `COG`) %>%
  na.omit()

gene2cog = data.frame(GID = character(),
                      COG = character())

df_temp <- list()
for (row in 1:nrow(cogs)) {
  the_gid <- cogs[row,"X.query"][[1]]
  the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
  
  df_temp[[row]] <- tibble(GID = rep(the_gid, length(the_cogs)),
                               COG = the_cogs)
}
gene2cog <- bind_rows(df_temp)

gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")

gene2cog<- gene2cog[-which(gene2cog$COG=="-"),]
gene2cog<- gene2cog[!is.na(gene2cog$COG_Name),]

#<<<<<<<<<<<<<<<<<<<<获得gene2ko<<<<<<<<<<<<<<<<<<<<<<<<
gene2ko<- dplyr::select(emapper,GID,KO) %>%
  separate_rows(KO, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(KO))
dim(gene2ko)

#提取GID与Pathway信息，这里只有2列信息

gene2pathway<- dplyr::select(emapper,GID,Pathway) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  dplyr::filter(!is.na(Pathway))
dim(gene2pathway)

##########################################
#至此获得gene_eggnog,gene2go,gene2ko,gene2cog
####制作自己的Orgdb######################
AnnotationForge::makeOrgPackage(gene_info=gene_info,
                                go=gene2go,
                                ko=gene2ko,
                                pathway=gene2pathway,
                                maintainer='test <test@hormail.com>',#格式最好不要变
                                author='test <test@hormail.com>',
                                version="0.1" ,
                                outputDir=".", 
                                tax_id="0",
                                genus="M",
                                species="y",
                                goTable="go")
##重启.rs.restartR()
# R CMD build org.My.eg.db
dir.create("R_library")
install.packages("org.My.eg.db", repos = NULL, lib = "R_library",type = "source")
library(org.My.eg.db,lib.loc = "./R_library")
columns(org.My.eg.db)

##########################################################
	
gene_list <- read.delim("venn_setdiff_bulk.txt",header=T,stringsAsFactors = FALSE) 

names(gene_list) <- c('gene_id') 
gene_select <- gene_list$gene_id  
go_bp <- groupGO(gene     = gene_select,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "BP",
                 level    = 2,
                 readable = FALSE)

go_bp <- as.data.frame(go_bp)
go_bp$GO_Class <- "Biological Process"

go_cc <- groupGO(gene     = gene_select,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "CC",
                 level    = 2,
                 readable = FALSE)

go_cc <- as.data.frame(go_cc)
go_cc$GO_Class <- "Cellular Component"

go_mf <- groupGO(gene     = gene_select,
                 OrgDb    = org.My.eg.db,
                 keyType  = "GID",
                 ont      = "MF",
                 level    = 2,
                 readable = FALSE)
go_mf <- as.data.frame(go_mf)
go_mf$GO_Class <- "Molecular Function"

go_all <- rbind(go_bp, go_cc, go_mf)

p <- ggplot(go_all) + 
  geom_histogram(aes(x = Description, 
                     y = Count,
                     fill = GO_Class),
                 stat = "identity") + facet_wrap(~GO_Class, scales = "free_x") + 
  labs(title = "GO function classification", y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("venn_setdiff_bulk_go.pdf", p, width = 20, height = 7)


#################################################################				   

#gene2cog$COG_Name = paste("(", gene2cog$COG, ")", gene2cog$COG_Name, sep = " ")

emapper_up <- read.delim("venn_setdiff_rhi.emapper.annotations.csv",sep=",") %>%
  dplyr::select(GID=X.query,Gene_Symbol=Preferred_name, 
                GO=GOs,KO=KEGG_ko,Pathway =KEGG_Pathway, 
                OG=eggNOG_OGs,COG = COG_category,Gene_Name =seed_ortholog)

cog_info <- read_delim("cog_funclass.tab",
                       "\t", escape_double = FALSE, trim_ws = TRUE)
cogs <- emapper_up %>%
  dplyr::select(GID, COG = `COG`) %>%
  na.omit()

gene2cog = data.frame(GID = character(),
                      COG = character())

df_temp <- list()
for (row in 1:nrow(cogs)) {
  the_gid <- cogs[row,"X.query"][[1]]
  the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
  
  df_temp[[row]] <- tibble(GID = rep(the_gid, length(the_cogs)),
                           COG = the_cogs)
}
gene2cog <- bind_rows(df_temp)

gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")

gene2cog<- gene2cog[-which(gene2cog$COG=="-"),]
gene2cog<- gene2cog[!is.na(gene2cog$COG_Name),]

p <- ggplot(data = gene2cog) + 
  geom_bar(aes(x = COG, 
               fill = COG_Name)) +
  labs(title = "COG/KOG Function Classification ", 
       x = "",
       y = "Number of genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size=unit(1,"line"),
        legend.text = element_text(size = 7.5)) +
  guides(fill=guide_legend(ncol=1))
ggsave("cog.pdf", p, width = 16, height = 7)

##########################################################
################ GO 富集分析##############################
##########################################################

library(tidyverse)
library(clusterProfiler)
go_class <- read.delim("go_term.list", header=FALSE, stringsAsFactors =FALSE) 
names(go_class) <- c('GO_ID','Description','Ontology') 

egg <- read.csv("final_table_sequence.emapper.annotations.csv",header=T, stringsAsFactors =FALSE) 
egg<- egg[,c("X.query","GOs")]
names(egg) <- c('gene_id','GO_ID') 

go_anno<- egg %>%
  mutate(GO_ID = strsplit(GO_ID, ",")) %>% 
  unnest(GO_ID)

gene_list <- read.delim("differential_genes_id_up.txt",header=T,stringsAsFactors = FALSE) 

names(gene_list) <- c('gene_id') 
gene_select <- gene_list$gene_id  

#Go富集分析 
#cat("         >> go_rich analysis starting...", fill=TRUE)
#cat("", fill = TRUE)
go_rich <- enricher(gene = gene_select,
                    TERM2GENE = go_anno[c('GO_ID','gene_id')],
                    TERM2NAME = go_class[c('GO_ID','Description')],
                    pvalueCutoff = 1,
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 1) 
write.table(go_rich,paste('go_rich.txt',sep=""), sep = '\t', row.names = FALSE, quote = FALSE)
p1<- barplot(go_rich,drop=T,showCategory = 30) 
ggsave(p1, file=paste('go_rich_bar.pdf',sep =""), width=8, height=6)
p2<- dotplot(go_rich, showCategory=30,font.size = 10) + ggtitle("dotplot for ORA")
ggsave(p2, file=paste('go_rich_dot.pdf',sep =""), width=8, height=6)
#cat("  

###################################################
#####################KO富集#######################

###pathway2name

##################STEP4-2: 得到pathway2name, ko2pathway##########################
# 需要下载 json文件(这是是经常更新的)
# https://www.genome.jp/kegg-bin/get_htext?ko00001
# 代码来自：http://www.genek.tv/course/225/task/4861/show
if(F){
  library(jsonlite)
  library(purrr)
  library(RCurl)
  
  update_kegg <- function(json = "ko00001.json") {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    
    kegg <- fromJSON(json)
    
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
    
    save(pathway2name, ko2pathway, file = "kegg_info.RData")
  }
  
  update_kegg(json = "ko00001.json")
  
}
pathway2name

pathway2gene <- AnnotationDbi::select(org.My.eg.db, 
                                      keys = keys(org.My.eg.db), 
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(GID,KO)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

egg <- read.csv("final_table_sequence.emapper.annotations.csv",header=T, stringsAsFactors =FALSE) 
egg<- egg[,c("X.query","KEGG_Pathway")]
names(egg) <- c('gene_id','KO_ID') 

ko_anno<- egg %>%
  mutate(KO_ID = strsplit(KO_ID, ",")) %>% 
  unnest(KO_ID)

gene_list <- read.delim("differential_genes_id_up.txt",header=T,stringsAsFactors = FALSE) 

names(gene_list) <- c('gene_id') 
gene_select <- gene_list$gene_id 

ko_rich <- enricher(gene = gene_select,
                    TERM2GENE = ko_anno[c('KO_ID','gene_id')],
                    TERM2NAME = pathway2name[c("Pathway","Name")],
                    pvalueCutoff = 1,
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 1) 
write.table(ko_rich,paste('Ko_rich.txt',sep=""), sep = '\t', row.names = FALSE, quote = FALSE)
p1<- barplot(ko_rich,drop=T,showCategory = 30) 
ggsave(p1, file=paste('ko_rich_bar.pdf',sep =""), width=8, height=6)
p2<- dotplot(ko_rich, showCategory=30,font.size = 10) + ggtitle("dotplot for ORA")
ggsave(p2, file=paste('ko_rich_dot.pdf',sep =""), width=8, height=6)
#cat(" 

###################################
####plot_bar#######################
###################################

library(ggplot2)
databar1<- read.table("clipboard",header = T,sep ="\t")
p1<- ggplot(data=databar1,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Node_counts") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar2<- read.table("clipboard",header = T,sep ="\t")
p2<- ggplot(data=databar2,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Edge_counts") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar3<- read.table("clipboard",header = T,sep ="\t")
p3<- ggplot(data=databar3,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Average_degree") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar4<- read.table("clipboard",header = T,sep ="\t")
p4<- ggplot(data=databar4,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Average_path_length") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar5<- read.table("clipboard",header = T,sep ="\t")
p5<- ggplot(data=databar5,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Clustering_coefficient") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar6<- read.table("clipboard",header = T,sep ="\t")
p6<- ggplot(data=databar6,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Cluster_number") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar7<- read.table("clipboard",header = T,sep ="\t")
p7<- ggplot(data=databar7,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Modularity") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar8<- read.table("clipboard",header = T,sep ="\t")
p8<- ggplot(data=databar8,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Diameter") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar9<- read.table("clipboard",header = T,sep ="\t")
p9<- ggplot(data=databar9,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Degree_assortativity") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))

databar10<- read.table("clipboard",header = T,sep ="\t")
p10<- ggplot(data=databar10,mapping=aes(x =group, y =value,fill=group))+
  geom_bar(stat="identity",width=0.6,position=position_dodge(0.90))+ 
  #geom_hline(aes(yintercept=0.930),colour="#990000", linetype="dashed")+
  #geom_text(aes(label = y, vjust = -0.8, hjust = 0.5), show.legend = TRUE)+
  labs(x = "", y = "Density") +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))+
  theme_set(theme_bw())+
  theme(legend.position = 'none')+
  theme(axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'), 
        axis.text = element_text(face = "bold",colour="black",size=15),
        axis.text.x = element_text(vjust=0.8,hjust=0.7,angle=20,size=15),
        axis.title.x = element_text(face = "bold",colour="black",size=15),
        axis.title.y = element_text(face = "bold",colour="black",size=15))
library(ggpubr)
library(grid)
library(gridExtra)
p0<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 3)
ggsave(p0, file="gene_netwotk.pdf", width=9, height=7,device=cairo_pdf)

