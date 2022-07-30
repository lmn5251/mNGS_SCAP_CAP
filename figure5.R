##SCAP survived vs nonsurvived
setwd("D:\\SCAP")
pheno_SCAP=readxl::read_xlsx('pheno_SCAP-1.xlsx') 
severe_pneumonia=read.table("dds_normalized_counts.txt",header=F,row.names=1)
a=match(as.matrix(pheno_SCAP[,7]),as.matrix(severe_pneumonia[1,]))
normalize_SCAP=severe_pneumonia[,a]
colnames(normalize_SCAP)=normalize_SCAP[1,]
normalize_SCAP=as.matrix(normalize_SCAP[-1,])
normalize_SCAP=apply(normalize_SCAP,1,as.numeric)
immune_gene_1=readxl::read_xlsx("immune_gene.express.xlsx",col_names=F)
b=match(as.matrix(pheno_SCAP[,7]),as.matrix(immune_gene_1[1,]))
immune_gene=as.matrix(immune_gene_1[,b])
rownames(immune_gene)=as.matrix(immune_gene_1[,1])
colnames(immune_gene)=immune_gene[1,]
immune_gene=as.matrix(immune_gene[-1,])
immune_gene=apply(immune_gene,1,as.numeric)
condition=as.matrix(pheno_SCAP[which(!is.na(a)),60])
condition=gsub("1","Balf_death",condition)
condition=gsub("2","Sputum_death",condition)
condition=gsub("3","Balf_survive",condition)
condition=gsub("4","Sputum_survive",condition)
pheno_survive=condition
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
#fig 5a
dgelist <- DGEList(counts = t(normalize_SCAP),group=condition) 
#（2）
keep <- rowSums(cpm(dgelist) > 1 ) >= 1
dgelist <- dgelist[keep, keep.lib.sizes = FALSE]
#（3）
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#
design <- model.matrix(~0+factor(pheno_survive))
#
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#
fit <- glmQLFit(dge, design, robust = TRUE)
colnames(design)=c("Balf_death","Balf_survive","Sputum_death","Sputum_survive")
fit1= glmQLFTest(fit, contrast=c(1,-1,0,0))
fit2= glmQLFTest(fit, contrast=c(0,0,1,-1))
lrt <- topTags(fit1, n=Inf,adjust.method="BH")
lrt1 <- topTags(fit2, n=Inf,adjust.method="BH")
DEG=as.data.frame(lrt)
DEG1=as.data.frame(lrt1)
length(which(DEG$FDR < 0.1 & abs(DEG$logFC)>1))
length(which(DEG1$FDR < 0.1 & abs(DEG1$logFC)>1))
a=DEG[which(DEG$FDR < 0.1 &  abs(DEG$logFC)>1),]
b=DEG1[which(DEG1$FDR < 0.1 &  abs(DEG1$logFC)>1),]
k=intersect(rownames(a),rownames(b)) #
write.table(lrt, 'control_treat.glmLRT.txt', sep = '\t', col.names = NA, quote = FALSE)

#GSEA
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
msigdb_1 <- read.gmt("c2.cp.reactome.v7.5.1.symbols.gmt")
msigdb_2 <-read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")  
msigdb_3 <- read.gmt("h.all.v7.5.1.symbols.gmt")
kegmt <- rbind(msigdb_1,msigdb_2)
kegmt1 <- rbind(kegmt,msigdb_3)
geneList = a$logFC
names(geneList) = rownames(a)
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
set.seed(1)
library(enrichplot)
egmt<-GSEA(geneList,TERM2GENE = kegmt1,  pvalueCutoff = 0.05,
           pAdjustMethod = "BH") 
egmt_result_df <- as.data.frame(egmt)
save(egmt,egmt_result_df,file = "GSEA_deg_h.all.rda")
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_balf_death_survive_gsea.png", width = 7, height = 6,units = "in",dpi = 800,bg="transparent")
gseaplot2(egmt, geneSetID = c(1:1), pvalue_table = F,color="#7876B1FF")
dev.off()
#fig 5b
u1=DEG[which(rownames(DEG)=="LTF"),]
u2=DEG[which(rownames(DEG)=="PTPRB"),]
u3=DEG[which(rownames(DEG)=="SLC2A5"),]
u4=DEG[which(rownames(DEG)=="CAMP"),]
u5=DEG[which(rownames(DEG)=="ORM2"),]
u6=DEG[which(rownames(DEG)=="RNASE3"),]
u7=DEG[which(rownames(DEG)=="MS4A3"),]
u8=DEG[which(rownames(DEG)=="MPO"),]
u9=DEG[which(rownames(DEG)=="ORM1"),]
u10=DEG[which(rownames(DEG)=="RETN"),]
u11=DEG[which(rownames(DEG)=="DEFA4"),]
o1=rbind(u1,u2);o2=rbind(o1,u3);o3=rbind(o2,u4);o4=rbind(o3,u5);o5=rbind(o4,u6);
o6=rbind(o5,u7);o7=rbind(o6,u8);o8=rbind(o7,u9);o9=rbind(o8,u10);o10=rbind(o9,u11)
write.table(o10,file="o10.txt",sep="\t")
setwd("D:\\SCAP")
o10=read.table("o10.txt",sep="\t",header=T)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "NEU_SCAP_balf_FC.png", width = 5, height = 7,units = "in",dpi = 800,bg="transparent")
p=ggplot(o10, aes(x=gene, y=logFC, fill=sample)) +
  geom_bar(stat="identity", color="black", width=.6,position = position_dodge(width=0.6)) 
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
#fig 5d
#SCAP SOFA-decreased vs SOFA-nondecreased
pheno_SCAP=readxl::read_xlsx('pheno_SCAP-1.xlsx') 
a=match(as.matrix(pheno_SCAP[,4]),as.matrix(RNA.SDSMRN.species[1,])) 
condition=as.matrix(pheno_SCAP[which(!is.na(a)),65])
condition=gsub("0","SOFA_bad",condition)
condition=gsub("2","SOFA_good",condition)
pheno_survive=condition
library(edgeR)
#
dgelist <- DGEList(counts = t(normalize_SCAP),group=condition) 
#
keep <- rowSums(cpm(dgelist) > 1 ) >= 1
dgelist <- dgelist[keep, keep.lib.sizes = FALSE]
#
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#
design <- model.matrix(~0+factor(pheno_survive))
#（1）
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
colnames(design)=c("SOFA_bad","SOFA_good")
fit1= glmQLFTest(fit, contrast=c(-1,1))
lrt <- topTags(fit1, n=Inf,adjust.method="BH")
DEG=as.data.frame(lrt)
length(which(DEG$FDR < 0.05 & abs(DEG$logFC)>1))
a=DEG[which(DEG$FDR < 0.05 &  abs(DEG$logFC)>1),]
# 
DEG$change = as.factor(ifelse(DEG$FDR < 0.05 & abs(DEG$logFC) > 1,
                               ifelse(DEG$logFC > 1 ,'Up','Down'),'Stable'))
table(DEG1$change)
this_tile <- paste0('Cutoff for logFC is ',1,
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='Up',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='Down',]))
g <- ggplot(data=DEG, 
            aes(x=logFC, 
                y=-log10(FDR),color=change)) + 
  geom_point(alpha=0.4, size=1.75) + 
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") + 
  ggtitle( this_tile ) + 
  theme_bw(base_size=20) +  
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','black','red'))+ ## 与分组一一对应
  theme(axis.text.x=element_text(vjust=1,size=18,face = "bold"))+
  theme(axis.text.y=element_text(vjust=1,size=18,face = "bold"))+
  theme(axis.title.x =element_text(size=18,face = "bold"), axis.title.y=element_text(size=18,face = "bold"))
g

#GSEA
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
msigdb_1 <- read.gmt("c2.cp.reactome.v7.5.1.symbols.gmt")
msigdb_2 <-read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")  
msigdb_3 <- read.gmt("h.all.v7.5.1.symbols.gmt")
kegmt <- rbind(msigdb_1,msigdb_2)
kegmt1 <- rbind(kegmt,msigdb_3)
j=match(colnames(immune_gene),rownames(DEG))  
j1=DEG[j,]
geneList =j1$logFC
names(geneList) = rownames(j1)
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
set.seed(1)
library(enrichplot)
egmt<-GSEA(geneList,TERM2GENE = kegmt1,  pvalueCutoff = 0.05,
           pAdjustMethod = "BH") 
egmt_result_df <- as.data.frame(egmt)
save(egmt,egmt_result_df,file = "GSEA_deg_h.all.rda")
gseaplot2(egmt, geneSetID = c(1:2),color=c("#BC3C29FF","#0072B5FF"),pvalue_table = F)
gseaplot2(egmt, geneSetID = c(1:5), subplots = 1:3,pvalue_table = T)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_SOFA_gsea.png", width = 14.5, height = 7,units = "in",dpi = 800,bg="transparent")
dotplot(egmt,split=".sign",showCategory = 10,font.size = 12, title = "",
        label_format = 30,)+facet_grid(~.sign)+
  #edit legends
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE),
    #reverse size order (higher diameter on top) 
    size = guide_legend(reverse = TRUE))+
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()
#fig 5e

setwd("D:\\脓毒症\\结果整理\\SCAP") 
t=read.table("IL10_45gene.txt",header=T)
t1=match(t[,1],colnames(normalize_SCAP))
t2=normalize_SCAP[,t1]
write.table(t2,file="t2.txt",sep="\t")
t2=read.table("t2.txt",sep="\t",header=T,row.names=1)
annotation=read.table("annotation_IL10.txt",sep="\t",header=T,row.names=1)
t2=log2(t2)
#rownames(four_gene_1)=colnames(four_gene)
library(RColorBrewer)
q1=t2[,which(colnames(t2)=="CXCL10")];
q2=t2[,which(colnames(t2)=="CCL20")];
q3=t2[,which(colnames(t2)=="CXCL1")];
q4=t2[,which(colnames(t2)=="CCL19")];
q5=cbind(q1,q2);q6=cbind(q5,q3);q7=cbind(q6,q4)
colnames(q7)=c("CXCL10","CCL20","CXCL1","CCL19")
rownames(q7)=rownames(t2)
pheatmap(q7,cluster_rows=F, cluster_cols = F,annotation_row = annotation,
         color = colorRampPalette(colors = c("lightblue","white","red"))(900),
         display_numbers =F,fontsize = 16,face="bold",border_color = NA)+
  theme_bw(base_size=18)

t3=match(t[,1],rownames(DEG))
t4=DEG[t3,]
write.table(t4,file="t4.txt",sep="\t")
t4=read.table("t4.txt",header=T,sep="\t")
setwd("D:\\figure\\v3_transparent") 
Cairo::CairoPNG(filename = "SCAP_IL10_45_FC.png", width = 10, height = 7,units = "in",dpi = 800,bg="transparent")
p=ggplot(t4, aes(x=gene, y=logFC)) +
  geom_bar(stat="identity", color="black",width=0.6,fill="#E18727FF") 
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
#fig 5f
###ISG DEGs
setwd("D:\\SCAP")
ISG=read.table('ISG_gene.txt',header=F,sep="\t",row.names=1) 
c=match(as.matrix(pheno_SCAP[,7]),as.matrix(ISG[1,])) 
colnames(ISG)=ISG[1,]
ISG=ISG[-1,]
ISG_gene=ISG[,c]
condition=as.matrix(pheno_SCAP[which(!is.na(c)),65])
condition=gsub("0","SOFA_bad",condition)
condition=gsub("2","SOFA_good",condition)
pheno_survive=condition
library(edgeR)
ISG_gene=apply(ISG_gene,2,as.numeric)
#
dgelist <- DGEList(counts = t(log2(ISG_gene)),group=condition) 
#
keep <- rowSums(cpm(dgelist) > 1 ) >= 1
dgelist <- dgelist[keep, keep.lib.sizes = FALSE]
#
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#
design <- model.matrix(~0+factor(pheno_survive)) 
#
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#
fit <- glmQLFit(dge, design, robust = TRUE)
colnames(design)=c("SOFA_bad","SOFA_good")

fit1= glmQLFTest(fit, contrast=c(-1,1))
lrt <- topTags(fit1, n=Inf,adjust.method="BH")
DEG=as.data.frame(lrt)
length(which(DEG$FDR < 0.05 & abs(DEG$logFC)>1))
a=DEG[which(DEG$FDR < 0.05 &  abs(DEG$logFC)>1),]
w1=ISG_gene[,which(colnames(ISG_gene)=="GZMB")]
w2=ISG_gene[,which(colnames(ISG_gene)=="IL15RA")]
w3=ISG_gene[,which(colnames(ISG_gene)=="CTCFL")]
w4=ISG_gene[,which(colnames(ISG_gene)=="IFNLR1")]
w5=ISG_gene[,which(colnames(ISG_gene)=="CXCL11")]
w6=ISG_gene[,which(colnames(ISG_gene)=="CCL8")]
w7=ISG_gene[,which(colnames(ISG_gene)=="HESX1")]
w8=ISG_gene[,which(colnames(ISG_gene)=="ANXA2R")]
w9=ISG_gene[,which(colnames(ISG_gene)=="VEGFC")]
w10=ISG_gene[,which(colnames(ISG_gene)=="CXCL9")]
w11=ISG_gene[,which(colnames(ISG_gene)=="CXCL10")]
r1=cbind(w1,w2);r2=cbind(r1,w3);r3=cbind(r2,w4);r4=cbind(r3,w5);r5=cbind(r4,w6);
r6=cbind(r5,w7);r7=cbind(r6,w8);r8=cbind(r7,w9);r9=cbind(r8,w10);r10=cbind(r9,w11)
colnames(r10)=c("GZMB","IL15RA","CTCFL","IFNLR1","CXCL11","CCL8","HESX1","ANXA2R","VEGFC","CXCL9","CXCL10")
write.table(r10,file="r10.txt",sep="\t")
setwd("D:\\SCAP") 
r10=read.table("r10.txt",sep="\t",header=T,row.names=1)
r10=log2(r10)
annotation_r10=read.table("annotation_r10.txt",row.names = 1,header=T)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_SOFA_ISG_pheatmap.png", width = 7, height = 7,units = "in",dpi = 800,bg="transparent")
pheatmap(r10,annotation_row =data.frame(annotation_r10),cluster_rows=F, cluster_cols = F,
         color = colorRampPalette(colors = c("lightblue","white","red"))(1500),
         display_numbers =F,fontsize = 16,face="bold",border_color = NA,show_rownames = F)+
  theme_bw(base_size=18)
dev.off()
##FC of ISG DEGs
setwd("D:\\SCAP") 
e=readxl::read_xlsx("SOFA_ISG_DEG.xlsx")
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_SOFA_ISG_FC.png", width = 4, height = 7,units = "in",dpi = 800,bg="transparent")
p=ggplot(e, aes(x=gene, y=logFC)) +
  geom_bar(stat="identity", color="black",width=0.6,fill="#E18727FF") 
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
#fig 5c
#CIBERSORT
#SCAP survived vs nonsurvived
cibersort=read.table("CIBERSORT-Results.txt",header=T,sep="\t")
a=match(as.matrix(pheno_SCAP[,7]),as.matrix(cibersort[,1])) 
cibersort1=cibersort[na.omit(a),]
condition=as.matrix(pheno_SCAP[which(!is.na(a)),64])
condition=as.matrix(condition)
cibersort_1=cibersort1[which(condition=="0"),] #survived
cibersort_2=cibersort1[which(condition=="1"),] #nonsurvived
cibersort_survive=rbind(cibersort_1,cibersort_2)
cibersort_1=cibersort_survive[,1:23]
pheno_condition_1=as.matrix(condition[which(condition=="0"),])
pheno_condition_2=as.matrix(condition[which(condition=="1"),])
pheno_condition=rbind(pheno_condition_1,pheno_condition_2)
pheno_condition=as.data.frame(pheno_condition)
pheno_condition$group = as.factor(ifelse(as.numeric(pheno_condition$V1) ==1,1,0))
pheno_condition_1=pheno_condition$group
cibersort_2=cbind(cibersort_1,pheno_condition_1)
rownames(cibersort_2)=cibersort_2[,1]
cibersort_2=cibersort_2[,-1]
library(ggpubr)
library(RColorBrewer)
dat <- cibersort_2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample,-pheno_condition_1)
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "cibersort.png", width = 10, height = 7,units = "in",dpi = 800,bg="transparent")
p2=ggplot(dat,aes(Cell_type,Proportion,fill = pheno_condition_1 )) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw(base_size = 18) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,hjust = 1,face="bold",size=14))+
  theme(axis.text.y = element_text(face="bold",size=14))+
  theme(axis.title.y = element_text(face="bold",size=14))+
  theme(axis.title.x = element_text(face="bold",size=14))+
  scale_fill_npg()+ 
  stat_compare_means(aes(group = pheno_condition_1,label = ..p.signif..),method = "kruskal.test")
dev.off()
p2
#fig 5g
##
cibersort=read.table("CIBERSORT-Results.txt",header=T,sep="\t")
a=match(as.matrix(pheno_SCAP[,7]),as.matrix(cibersort[,1])) 
cibersort1=cibersort[na.omit(a),]
condition=as.matrix(pheno_SCAP[which(!is.na(a)),68])
condition=as.matrix(condition)
cibersort_1=cibersort1[which(condition=="0"),] #0 nondecreased
cibersort_2=cibersort1[which(condition=="2"),] #2 decreased
cibersort_survive=rbind(cibersort_1,cibersort_2)
cibersort_1=cibersort_survive[,1:23]
pheno_condition_1=as.matrix(condition[which(condition=="0"),])
pheno_condition_2=as.matrix(condition[which(condition=="2"),])
pheno_condition=rbind(pheno_condition_1,pheno_condition_2)
pheno_condition=as.data.frame(pheno_condition)
pheno_condition$group = as.factor(ifelse(as.numeric(pheno_condition$V1) <2,1,0))
pheno_condition_1=pheno_condition$group
cibersort_2=cbind(cibersort_1,pheno_condition_1)
rownames(cibersort_2)=cibersort_2[,1]
cibersort_2=cibersort_2[,-1]
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
dat <- cibersort_2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample,-pheno_condition_1)
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "cibersort-SOFA.png", width = 10, height = 7,units = "in",dpi = 800,bg="transparent")
p2=ggplot(dat,aes(Cell_type,Proportion,fill = pheno_condition_1 )) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw(base_size = 18) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,hjust = 1,face="bold",size=14))+
  theme(axis.text.y = element_text(face="bold",size=14))+
  theme(axis.title.y = element_text(face="bold",size=14))+
  theme(axis.title.x = element_text(face="bold",size=14))+
  scale_fill_npg()+ 
  stat_compare_means(aes(group = pheno_condition_1,label = ..p.signif..),method = "kruskal.test")
dev.off()
p2
