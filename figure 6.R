#fig 6a
setwd("D:\\SCAP") 
pheno_all=readxl::read_xlsx('pheno-all.xlsx')
diversity=read.table("diversity_all.txt",sep="\t",header=T,row.names = 1)
gene6=read.table("6-gene.txt",header=T,sep="\t",row.names=1)
m_SCAP_spearman= rcorr(as.matrix(gene6[,1:6]),as.matrix(diversity[,1:6]),type="spearman") 
write.table(m_SCAP_spearman$r,file="relation_gene_diversity_r_spearman.txt",sep="\t")
write.table(m_SCAP_spearman$P,file="relation_gene_diversity_p_spearman.txt",sep="\t")
library(corrplot)#
library(Hmisc)
library(RColorBrewer)
library(pheatmap)
setwd("D:\\SCAP")
corr_all=read.table("relation_gene_diversity_r_spearman.txt",header=T,sep="\t",row.names = 1)
corr_all_p=read.table("relation_gene_diversity_p_spearman.txt",header=T,sep="\t",row.names = 1)
corr_all_1=t(apply(corr_all,2,as.numeric))
corr_all_p=t(apply(corr_all_p,2,as.numeric))
colnames(corr_all_1)=rownames(corr_all)
corr_all_2=t(corr_all_1)
corr_all_p=t(corr_all_p)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "corelation_gene_micro_pheatmap1.png", width = 7, height = 8,units = "in",res = 800,bg="transparent")
pheatmap(corr_all_2[1:6,],
         scale = "none",cluster_rows=F, cluster_cols = F, #4DBBD5FF
         color = colorRampPalette(colors = c("white","#E64B35FF"))(1000),
         display_numbers =T,
         fontsize =20,face="bold",angle=90)+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))
dev.off()
#fig 6b 
setwd("D:\\correlation")
q1=which(colnames(all_count_1)=="IL15RA")
IL15RA=all_count_1[,q1]
q2=which(colnames(all_count_1)=="CTCFL")
CTCFL=all_count_1[,q2]
q3=which(colnames(all_count_1)=="IFNLR1")
IFNLR1=all_count_1[,q3]
q4=which(colnames(all_count_1)=="CXCL11")
CXCL11=all_count_1[,q4]
q5=which(colnames(all_count_1)=="CCL8")
CCL8=all_count_1[,q5]
q6=which(colnames(all_count_1)=="HESX1")
HESX1=all_count_1[,q6]
q7=which(colnames(all_count_1)=="ANXA2R")
ANXA2R=all_count_1[,q7]
q8=which(colnames(all_count_1)=="VEGFC")
VEGFC=all_count_1[,q8]
q9=which(colnames(all_count_1)=="CXCL9")
CXCL9=all_count_1[,q9]
q10=which(colnames(all_count_1)=="CXCL10")
CXCL10=all_count_1[,q10]
q11=which(colnames(all_count_1)=="CCL19")
CCL19=all_count_1[,q11]
q12=which(colnames(all_count_1)=="CCL20")
CCL20=all_count_1[,q12]
q13=which(colnames(all_count_1)=="CXCL1")
CXCL1=all_count_1[,q13]

IFN1=cbind(IL15RA,CTCFL);IFN2=cbind(IFN1,IFNLR1);IFN3=cbind(IFN2,CXCL11);
IFN4=cbind(IFN3,CCL8);IFN5=cbind(IFN4,HESX1);IFN6=cbind(IFN5,ANXA2R);IFN7=cbind(IFN6,VEGFC);IFN8=cbind(IFN7,CXCL9);
IFN9=cbind(IFN8,CXCL10);IFN10=cbind(IFN9,CCL19);IFN11=cbind(IFN10,CCL20);IFN12=cbind(IFN11,CXCL1);

result=cbind(gene_6,IFN12)
result1=cbind(result,gene10)
m_SCAP_spearman= rcorr(as.matrix(result1[,1:30]),as.matrix(diversity[,1:6]),type="spearman") 
write.table(m_SCAP_spearman$r,file="relation_IFNgene_diversity_r_spearman.txt",sep="\t")
write.table(m_SCAP_spearman$P,file="relation_IFNgene_diversity_p_spearman.txt",sep="\t")
library(corrplot)#
library(Hmisc)
library(RColorBrewer)
library(pheatmap)
setwd("D:\\correlation")
corr_all=read.table("relation_IFNgene_diversity_r_spearman.txt",header=T,sep="\t",row.names = 1)
corr_all_p=read.table("relation_IFNgene_diversity_p_spearman.txt",header=T,sep="\t",row.names = 1)
corr_all_1=t(apply(corr_all,2,as.numeric))
colnames(corr_all_1)=rownames(corr_all)
corr_all_2=t(corr_all_1)
corr_all_p=as.matrix(corr_all_p)
#corr_all_p=t(corr_all_p)
setwd("D:\\v3_transparent")
Cairo::CairoPNG(filename = "corelation_IFNgene_micro_pheatmap1.png", width = 7, height = 8,units = "in",res = 800,bg="transparent")
pheatmap(corr_all_2[1:13,],
         scale = "none",cluster_rows=F, cluster_cols = F, #4DBBD5FF
         color = colorRampPalette(colors = c("#4DBBD5FF","white","#E64B35FF"))(1000),
         display_numbers = T,
         fontsize = 16,face="bold")+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15,angle=45),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))
dev.off()
#fig 6c

species=readxl::read_xlsx("SCAP_CAP_top_species.xlsx")
k2=match(as.matrix(species),colnames(all_DNA))
all_DNA_1=all_DNA[,k2]
gene6=read.table("6-gene.txt",sep="\t",header=T)
species=rownames(u1)
k2=match(as.matrix(species),colnames(all_DNA))
all_DNA_1=all_DNA[,k2]
k4=rcorr(as.matrix(all_DNA_1),as.matrix(gene6),type="spearman")
library(pheatmap)
relation_clinical_gene_r=read.table("relation_clinical_gene_r_spearman.txt",sep="\t",header=T,row.names=1)
relation_clinical_gene_p=read.table("relation_clinical_gene_p_spearman.txt",sep="\t",header=T,row.names=1)

setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "DEG6_gene_abundance_species_correlation.png", width = 7, height =8,units = "in",dpi = 800,bg="transparent")
pheatmap(relation_clinical_gene_r[1:12,],
         scale = "none",cluster_rows=F, cluster_cols = F,
         color = colorRampPalette(colors = c("#4DBBD5FF","white","#E64B35FF"))(50),
         display_numbers = matrix(ifelse(relation_clinical_gene_p<0.05,ifelse(relation_clinical_gene_p<0.01,ifelse(relation_clinical_gene_p<0.001,"***", '**'),'*'), ' '), nrow(relation_clinical_gene_p)),
         fontsize = 16,fontface="italic",)+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))
dev.off()
#fig 6d
setwd("D:\\correlation")
w1=c("LTF","PTPRB","SLC2A5","CAMP","ORM2","RNASE3","MS4A3","MPO","ORM1","RETN","DEFA4")
w2=c("Acinetobacter_baumannii","Staphylococcus_aureus","Streptococcus_pneumoniae",
     "Rothia_mucilaginosa","Haemophilus_influenzae","Streptococcus_salivarius","Prevotella_oris",
     "Prevotella_pallens","Streptococcus_mitis","Gemella_haemolysans")
r1=match(w1,colnames(normalize_SCAP))
r2=normalize_SCAP[,r1]
t1=match(w2,colnames(SCAP_DNA))
t2=SCAP_DNA[,t1]
library(Hmisc);library(corrplot)
w4=rcorr(as.matrix(r2),as.matrix(t2),type="spearman")
write.table(w4$r,file="corelation_gene_species_balf_NEU_r_spearman.txt",sep="\t")
write.table(w4$P,file="corelation_gene_species_balf_NEU_p_spearman.txt",sep="\t")
setwd("D:\\SCAP")
corr_r=read.table("corelation_gene_species_balf_NEU_r_spearman.txt",sep="\t",header=T,row.names = 1)
corr_p=read.table("corelation_gene_species_balf_NEU_p_spearman.txt",sep="\t",header=T,row.names = 1)
library(corrplot)
col <- colorRampPalette(c("red", "white", "cyan"))(500)
corr_r=as.matrix(corr_r)
library(pheatmap)
corr_r=t(corr_r)
corr_p=t(corr_p)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "NEU_gene_abundance_species_correlation.png", width = 7, height =8,units = "in",dpi = 800,bg="transparent")
pheatmap(corr_r[1:10,],
         scale = "none",cluster_rows=F, cluster_cols = F,
         color = colorRampPalette(colors = c("#4DBBD5FF","white","#E64B35FF"))(50),
         display_numbers = matrix(ifelse(corr_p<0.05,ifelse(corr_p<0.01,ifelse(corr_p<0.001,"***", '**'),'*'), ' '), nrow(corr_p)),#数值>2时显示星号, ≤2时显示·,
         fontsize = 18,face="bold",fontface="italic")+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))
dev.off()
