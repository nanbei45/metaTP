#!/usr/bin/env Rscript
# R code for DEG_analysis
# This function will be imported into Python program
# author: Wang,Yansu 
# ///
# standalone usage:
# R --slave --no-restore --args --source_dir sourceDir -i gene_expression.csv -g group_file \ 
# -n experiment/control -p cut_off_pvalue -f cut_off_logFC -o DEG_analysis_result < sourceDir/DEG_analysis.r
 
# print help string:
# R --slave --args -h --source_dir sourceDir < sourceDir/DEG_analysis.r
#

# load R packages and suppress the Warnings and Messages 
suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("rstatix")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggrepel")))

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i","--inputfile"), type="character",
        help="The gene expression matrix; csv-based file format is necessary [required]."),
    make_option(c("-g", "--group_informaton"), type="character",
        help="The group information; csv-based file format is necessary [required]."),
    make_option(c("-n", "--group_name"), type="character",
         help="Experimental group/control group;which is the same as that in the groupfile;
		      such as Rhizosphere/Bulk[required]."),
    make_option(c("-p", "--cut_off_pvalue"), type="double",default = 0.001,
        help="the value of cut_off_pvalue;"),
	make_option(c("-f", "--cut_off_logFC"), type="double",default = 1,
        help="the value of cut_off_logFC;"),
    make_option(c("-o", "--output_dir"), type="character",
       help="The output foldersaving differential gene list,volcano plots[required]"),
    make_option(c('-v','--verbose'), action='store_true', default=FALSE,
        help="Print information about execution [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), 
    args=commandArgs(trailing=TRUE))

if (!dir.exists(opts$output_dir)){
dir.create(opts$output_dir)
} else {
    print("Dir already exists!")
}
# loading data 
gene_file <- read.csv(opts$inputfile, header=TRUE)
# loading group information
group_file <- read.csv(opts$group_informaton, header=TRUE)
# grouping and sorting
df = gene_file %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  
  left_join(group_file,by=c("sample" = "sample")) # merge with group data 
# calculating FC
group1 <- strsplit(opts$group_name, "/", fixed= T)[[1]][1]
group2 <- strsplit(opts$group_name, "/", fixed= T)[[1]][2]
dfFC1 = df %>%
  group_by(GeneID,group)%>%
  summarise(mean = mean(value,na.rm=T)) %>% 
  pivot_wider(names_from = group,values_from = mean)%>%
  select(any_of(c(group1,group2)))

dfFC1= as.data.frame(dfFC1)
#dfFC1[,3]<-dfFC1[,3]+0.00001
dfFC1[,2][dfFC1[,2] == '0']<- 0.000000001
dfFC1[,3][dfFC1[,3] == '0']<- 1
dfFC1 <- mutate(dfFC1, FC = dfFC1[,2]/dfFC1[,3])
dfFC1 <- mutate(dfFC1, logFC = log2(dfFC1[,2]/dfFC1[,3]))
dfFC <- dfFC1[,c("GeneID","FC","logFC")]
##转换数据类型
dfFC <- tibble::as_tibble(dfFC)
# Calculating p-value
dfP = df %>%
  group_by(GeneID) %>%
  t_test(value ~ group,var.equal=T)
# FDR correction for P-values
dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  #the algorithm chooses BH
#####################################################################
dfdata = dfP %>%
  left_join(dfFC1) %>%
  left_join(dfFC) %>%
  left_join(dfP_FDR)

#print(R.utils::getAbsolutePath(opts$output_dir))
write.csv(dfdata,paste(opts$output_dir,"/","differential_genes.csv", sep=""))
cat("         >> Differential gene analysis finished...", fill=TRUE)
cat("", fill = TRUE)
###volcano plot#####
cat("         >> Volcano plot starting...", fill=TRUE)
cat("", fill = TRUE)

library(ggplot2)
dataset<- read.csv(paste(opts$output_dir,"/","differential_genes.csv", sep=""),header = T,row.names = 1)
cut_off_pvalue = opts$cut_off_pvalue  
cut_off_logFC = opts$cut_off_logFC
dataset$change = ifelse(dataset$FDR < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC,
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
write.csv(dataset,paste(opts$output_dir,"/","differential_genes.csv", sep=""))
		
p <- ggplot(
  dataset, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","grey48","red"))+
  scale_x_continuous(limits = c(-2.5,10))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))
ggsave(p, file=paste(opts$output_dir,"/",'differential_genes_volcano.pdf',sep =""), width=8, height=6)
while (!is.null (dev.list ())) dev.off ()

p2 <- ggplot(
  dataset, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=1, size=4) +
  scale_color_manual(values=c("blue","grey48","red"))+
  scale_x_continuous(limits = c(-2.5,10))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  geom_label_repel(aes(label =GeneID),data=dataset, nudge_y = 2,  alpha = 0.6  )+
  theme_set(theme_bw())+
  theme(legend.position="right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "white", color="black",size=2),
        panel.grid = element_blank(),
        text=element_text(face = "bold",colour="black",size=15),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(face = "bold",colour="black",size=15))
ggsave(p2, file=paste(opts$output_dir,"/",'differential_gene_volcano_label.pdf',sep =""), width=8, height=6)
while (!is.null (dev.list ())) dev.off ()

cat("         >> Volcano plot finished...", fill=TRUE)
cat("", fill = TRUE)
