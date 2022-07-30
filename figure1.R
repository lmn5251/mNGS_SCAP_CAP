setwd("D:\\community_diversity\\SCAP_CAP") 
install.packages("vegan")
library(vegan)  
install.packages("ggsignif")
library(ggsignif) 
pheno_SCAP=readxl::read_xlsx("pheno_SCAP-1.xlsx") 
pheno_CAP=readxl::read_xlsx("pheno_CAP-1.xlsx") 
DNA.SDSMRN.species=read.table("DNA.SDSMRN.179.filter.species.txt",header=F,row.names=1)
c=match(as.matrix(pheno_CAP[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
d=match(as.matrix(pheno_SCAP[,3]),as.matrix(DNA.SDSMRN.species[1,]))
d=as.numeric(na.omit(d))
pheno_SCAP_DNA=DNA.SDSMRN.species[,d]
pheno_CAP_DNA=DNA.SDSMRN.species[,na.omit(c)]
SCAP_group_DNA=as.matrix(pheno_SCAP[which(!is.na(d)),1])
CAP_group_DNA=as.matrix(pheno_CAP[which(!is.na(c)),1])
##############################################################################################################
#SCAP vs CAP 
pheno_all=readxl::read_xlsx('pheno-all.xlsx') 
a=match(as.matrix(pheno_all[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
DNA=DNA.SDSMRN.species[,na.omit(a)]
condition=as.matrix(pheno_all[which(!is.na(a)),62])
condition=gsub("1","SCAP_BALF",condition)
condition=gsub("2","SCAP_Sputum",condition)
condition=gsub("3","CAP_BALF",condition)
condition=gsub("4","CAP_Sputum",condition)
DNA=as.data.frame(t(DNA)) 
rownames(DNA)=DNA[,1]
DNA=DNA[,-1]
DNA_1=as.data.frame(apply(DNA,2,as.numeric))
alpha_diversity <- function(x){
  library(vegan) 
  shannon = diversity(x,index="shannon") #shannon index
  simpson = diversity(x, index = "simpson") #simpson index
  Chao1 = estimateR(x)[2,] #Chao1 index
  ACE = estimateR(x)[4,] #ACE index
}
diversity=alpha_diversity(DNA_1)
data_ggplot=data.frame(diversity[,1]) )##1: shannon,2:simpson;3:Chao1; 4:ACE
data_ggplot_DNA=data_ggplot
data_ggplot=cbind(condition,data_ggplot_DNA)
colnames(data_ggplot)=c("Group","data_norm_shannon")
compared_list = list(c("SCAP_BALF", "CAP_BALF"),c("SCAP_Sputum","CAP_Sputum")) 
library(ggsci)
#ggplot2
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "SCAP_CAP_shannon.png", width = 6.3, height = 7,units = "in",dpi = 800,bg="transparent")
alpha_boxplot=ggplot(data_ggplot, aes(x=Group, y=data_norm_shannon, fill=Group))+
  geom_boxplot()+
  labs(x="Group", y="ACE Index",face="bold",size=18)+
  theme(plot.title=element_text(hjust=0), legend.title=element_blank())+
  theme_bw(base_size=18)+
  scale_fill_npg()
alpha_boxplot  + geom_signif(comparisons = compared_list, test = wilcox.test, step_increase = 0.1, map_signif_level=T,size=0.8,textsize=6)+
  theme(axis.text.x = element_text(face="bold", size=14,angle=45,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=20))+
  geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=1)    +
  scale_fill_npg()
dev.off()
data_ggplot
##########################################################################################################
##bray-curtis distance 
pheno_all=readxl::read_xlsx('pheno-all.xlsx') 
a=match(as.matrix(pheno_all[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
DNA=DNA.SDSMRN.species[,na.omit(a)]
condition=as.matrix(pheno_all[which(!is.na(a)),61])
condition=gsub("0","SCAP",condition)
condition=gsub("2","CAP",condition)
DNA=as.data.frame(t(DNA)) 
rownames(DNA)=DNA[,1]
DNA=DNA[,-1]
DNA_1=as.data.frame(apply(DNA,2,as.numeric))
bray_dis <- vegdist(DNA_1, method = 'bray')
bray_dis=as.matrix(bray_dis)
condition=data.frame(condition)
data1=cbind(bray_dis,condition)
dune.adonis <- adonis2(bray_dis~as.matrix(condition),data=data1, permutations = 999) #PERANOVA
dune.adonis
write.table(bray_dis, file = "dist_matrix_SCAP_CAP.txt", row.names = T, quote = F, sep="\t")
distance=read.table("dist_matrix_SCAP_CAP.txt",sep="\t",header=T)
m=distance[which(condition=="SCAP_Sputum"),which(condition=="SCAP_Sputum")]
m[!upper.tri(m, diag = TRUE)]=0
m=as.numeric(as.matrix(m))
write.table(m,file="m.txt",sep="\t")
m1=distance[which(condition=="SCAP_Balf"),which(condition=="SCAP_Balf")]
m1[!upper.tri(m1, diag = TRUE)]=0
m1=as.numeric(as.matrix(m1))
write.table(m1,file="m1.txt",sep="\t")
m2=distance[which(condition=="CAP_Sputum"),which(condition=="CAP_Sputum")]
m2[!upper.tri(m2, diag = TRUE)]=0
m2=as.numeric(as.matrix(m2))
write.table(m2,file="m2.txt",sep="\t")
m3=distance[which(condition=="CAP_Balf"),which(condition=="CAP_Balf")]
m3[!upper.tri(m3, diag = TRUE)]=0
m3=as.numeric(as.matrix(m3))
write.table(m3,file="m3.txt",sep="\t")
distance=readxl::read_xlsx("distance.xlsx")
compared_list=list(c("SCAP_Sputum","CAP_Sputum"),c("SCAP_BALF","CAP_BALF"))
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_CAP_bray_distance.png", width = 7, height = 7,units = "in",dpi = 800,bg="transparent")
boxplot1=ggplot(distance, aes(x=sample, y=distance, fill=sample))+
  geom_boxplot()+
  labs(x="sample", y="Bray-Curtis distance",face="bold",size=18)+
  theme(plot.title=element_text(hjust=0), legend.title=element_blank())+
  theme_bw(base_size=18)
boxplot1  + geom_signif(comparisons = compared_list, test = wilcox.test, step_increase = 0.1, map_signif_level=T,size=0.8,textsize=6)+
  theme(axis.text.x = element_text(face="bold", size=14,angle=45,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=20))+
  geom_jitter(aes(fill=sample),width =0.2,shape = 21,size=1)    +
  scale_fill_nejm()
dev.off()
####################################################################################################################################
#PCoA analysis of SCAP vs CAP
library(ape)
DNA_1=as.data.frame(apply(DNA,2,as.numeric))
bray_dis <- vegdist(DNA_1, method = 'bray')
bray_dis=as.matrix(bray_dis)
df.pcoa<-pcoa(bray_dis,correction = "cailliez")
df.pcoa$vectors
df.pcoa$values
condition=as.matrix(pheno_all[which(!is.na(a)),62])
condition=gsub("1","SCAP_BALF",condition)
condition=gsub("2","SCAP_Sputum",condition)
condition=gsub("3","CAP_BALF",condition)
condition=gsub("4","CAP_Sputum",condition)
bray_dis=as.data.frame(bray_dis);condition=as.data.frame(condition)
dune.adonis <- adonis2(bray_dis~condition,data=condition, permutations = 999) #PERANOVA
dune.adonis
summary(dune.adonis)
library(usethis)
library(devtools)
devtools::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))
df.plot<-data.frame(df.pcoa$vectors)
x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "SCAP_CAP_PCoA.png", width = 9, height = 7,units = "in",dpi = 800,bg="transparent")
p=ggplot(data=round(df.plot,4),aes(x=round(Axis.1,4),y=round(Axis.2,4), color= as.matrix(condition)))+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))
p+stat_ellipse(data=round(df.plot,4),
               geom = "polygon",
               aes(fill=as.matrix(condition)),
               alpha=0.3)+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))+
  scale_color_nejm()+
  scale_fill_nejm()
dev.off()
#############################################################################################################################  
#relative abundances of top15 species comparing SCAP with CAP in BALF and sputum samples
mean_data=matrix(0,10,2871) 
m=as.matrix(apply(DNA_1,1,sum))
DNA_2=matrix(0,177,2871)
for (i in 1:177) {
  for (j in 1:2871) {
  DNA_2[i,j]=DNA_1[i,j]/m[i,]
  } 
}
DNA_2=as.matrix(DNA_2) 
colnames(DNA_2)=colnames(DNA_1) 
write.table(DNA_2,file="DNA_all_relative_abundance.txt",sep="\t")
colnames(mean_data)=colnames(DNA_2)
DNA_3=DNA_2
for (i in 1:2871) {
  mean_data[1,i]=mean((DNA_2[which(condition=="SCAP_Sputum"),i]))
  mean_data[2,i]=mean((DNA_2[which(condition=="SCAP_Balf"),i]))
  mean_data[3,i]=mean((DNA_2[which(condition=="CAP_Sputum"),i]))
  mean_data[4,i]=mean((DNA_2[which(condition=="CAP_Balf"),i]))
  mean_data[5,i]=sd((DNA_2[which(condition=="SCAP_Sputum"),i]))
  mean_data[6,i]=sd((DNA_2[which(condition=="SCAP_Balf"),i]))
  mean_data[7,i]=sd((DNA_2[which(condition=="CAP_Sputum"),i]))
  mean_data[8,i]=sd((DNA_2[which(condition=="CAP_Balf"),i]))
  mean_data[9,i]=mean_data[1,i]-mean_data[3,i] #SCAP-sputum  -   CAP-sputum
  mean_data[10,i]=mean_data[2,i]-mean_data[4,i]  #SCAP-balf  -  CAP-balf
}
p_data=matrix(0,2871,2) 
for (i in 1:2871) {
  p_data[i,1]=wilcox.test(as.matrix(DNA_2[which(condition=="SP_Sputum"),i]),as.matrix(DNA_2[which(condition=="NSP_Sputum"),i]))$p.value
  p_data[i,2]=wilcox.test(as.matrix(DNA_2[which(condition=="SP_Balf"),i]),as.matrix(DNA_2[which(condition=="NSP_Balf"),i]))$p.value
}
write.table(mean_data,file="mean_data.txt",sep="\t")
write.table(p_data,file="p_data.txt",sep="\t")
#top15 species of relative abundances in sputum samples
top15=readxl::read_xlsx("SCAP_CAP_sputum_top15.xlsx")
top15$type=factor(top15$type,levels = c("Rothia_mucilaginosa",
                                        "Acinetobacter_baumannii",
                                        "Pseudomonas_aeruginosa",
                                        "Streptococcus_parasanguinis",
                                        "Weissella_viridescens",
                                        "Candida_albicans",
                                        "Prevotella_melaninogenica",
                                        "Human_alphaherpesvirus_1",
                                        "Clavispora_lusitaniae",
                                        "Klebsiella_pneumoniae",
                                        "Elizabethkingia_anophelis",
                                        "Streptococcus_salivarius",
                                        "Moraxella_catarrhalis",
                                        "Staphylococcus_epidermidis",
                                        "Enterococcus_faecalis"
)) 
library(ggpubr) 
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "SCAP_CAP_sputum_top15.png", width = 9, height = 6,units = "in",dpi = 800,bg="transparent")
p=ggplot(top15, aes(x=type, y=percentage, fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  #geom_signif(comparisons = list(c("SCAP","CAP")),test = t.test, step_increase = 0.0, map_signif_level=T,  y_position=0.5)
  stat_compare_means(comparisons =list(c("SCAP","CAP")))
p+theme_bw(base_size=18)+scale_fill_npg()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
###top15 species of relative abundances in BALF samples
top15_balf=readxl::read_xlsx("SCAP_CAP_balf_top15-1.xlsx")
top15_balf$type=factor(top15_balf$type,levels = c("Acinetobacter_baumannii",
                                                  "Haemophilus_influenzae",
                                                  "Scedosporium_apiospermum",
                                                  "Klebsiella_pneumoniae",
                                                  "Rothia_mucilaginosa",
                                                  "Veillonella_parvula",
                                                  "Corynebacterium_striatum",
                                                  "Stenotrophomonas_maltophilia",
                                                  "Campylobacter_mucosalis",
                                                  "Prevotella_melaninogenica",
                                                  "Abiotrophia_defectiva",
                                                  "Human_betaherpesvirus_5",
                                                  "Neisseria_subflava",
                                                  "Streptococcus_pneumoniae",
                                                  "Human_alphaherpesvirus_1"
)) 
library(ggpubr)
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "SCAP_CAP_balf_top15.png", width = 9, height = 6,units = "in",dpi = 800,bg="transparent")
p=ggplot(top15_balf, aes(x=type, y=percentage, fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  #geom_signif(comparisons = list(c("SCAP","CAP")),test = t.test, step_increase = 0.0, map_signif_level=T,  y_position=0.5)
  stat_compare_means(comparisons =list(c("SCAP","CAP")))
p+theme_bw(base_size=18)+scale_fill_npg()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
 dev.off()
# geom_errorbar(aes(ymin=percentage-percentage_sd, ymax=percent +percent_sd),position=position_dodge(.6), width=.2)
################################################################################################
#the relationship between 12 species and clinical indicators    SCAP vs CAP
mds <- cmdscale(bray_dis, k = 2, eig = TRUE)
mds_point <- data.frame(mds$points)   # 
colnames(mds_point) <- c('X1','X2') 
eig <- mds$eig
library(vegan)
condition=as.matrix(pheno_all[which(!is.na(a)),61])
condition=gsub("0","SCAP",condition)
condition=gsub("2","CAP",condition)
condition=as.data.frame(condition)
group <- as.factor(condition$condition)
#fit <- envfit(mds, immune_gene,permutations = 999)
m=as.matrix(apply(DNA_1,1,sum))
all_DNA=matrix(0,177,2871)
for (i in 1:177) {
  for (j in 1:2871) {
    all_DNA[i,j]=DNA_1[i,j]/m[i,]
  } 
}
colnames(all_DNA)=colnames(DNA_1) 
all_DNA_1=apply(all_DNA,2,sum)
all_DNA_1=as.matrix(all_DNA_1)
k1=all_DNA[,which(all_DNA_1[,1]>0.1)]
setwd("D:\\\correlation")
pheno_all=readxl::read_xlsx('pheno-all-1.xlsx')
a=match(as.matrix(pheno_all[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
q1=pheno_all[!is.na(a),37:55]#
q2=pheno_all[!is.na(a),22]#APACHEII
q3=pheno_all[!is.na(a),35]#SOFA
q4=pheno_all[!is.na(a),64]#No.Comorbidity

q7=cbind(q1,q2);q8=cbind(q7,q3);q9=cbind(q8,q4)
fit <- envfit(mds,k1 ,permutations = 999)
fit_val <- scores(fit,display = c("vectors"))
fit_val <- fit_val*vegan::ordiArrowMul(fit_val, fill = 1.5)
u=fit$vectors
u1=u[["arrows"]][which(u[["pvals"]]<0.05),]
color <- c(brewer.pal(4,"Set1"))
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "SCAP_CAP_PCoA_1.png", width = 9, height = 7,units = "in",dpi = 800,bg="transparent")
ggplot(mds_point, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 4, alpha = .6) +
  #  scale_fill_manual(values=color) +
  theme_bw(base_size=16)+
  scale_color_manual(values = color) +
  geom_segment(data=data.frame(u1), 
               aes(x=0,y=0,xend=Dim1, yend=Dim2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color='black',alpha=1)  + 
  geom_label_repel(data=data.frame(u1), aes(Dim1, Dim2, label=rownames(u1)),
                   color='black',alpha=1,
                   segment.color = 'grey35',
                   point.padding = unit(0.1,"lines")) +
  labs(x = paste("PCoA 1 (", format(100*eig[1]/sum(eig), digits = 4), "%)",sep = ""), 
       y = paste("PCoA 2 (", format(100*eig[2]/sum(eig), digits = 4), "%)",sep = "")) 
dev.off()
species=readxl::read_xlsx("SCAP_CAP_top_species.xlsx")
k2=match(as.matrix(species),colnames(all_DNA))
all_DNA_1=all_DNA[,k2]
k3=rcorr(as.matrix(all_DNA_1),as.matrix(q9[,1:22]))
k4=rcorr(as.matrix(all_DNA_1),as.matrix(q9[,1:22]),type="spearman")
write.table(k3$r,file="relation_clinical_species_r.txt",sep="\t")
write.table(k3$P,file="relation_clinical_species_p.txt",sep="\t")
write.table(k4$r,file="relation_clinical_species_r_spearman.txt",sep="\t")
write.table(k4$P,file="relation_clinical_species_p_spearman.txt",sep="\t")
library(pheatmap)
relation_clinical_species_r=read.table("relation_clinical_species_r.txt",sep="\t",header=T,row.names=1)
#relation_clinical_species_r_spearman=read.table("relation_clinical_species_r_spearman.txt",sep="\t",header=T,row.names=1)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "correlation_SCAP_CAP.png", width = 15, height = 7,units = "in",dpi = 800,bg="transparent")
relation_clinical_species_r_1=t(relation_clinical_species_r)
pheatmap(relation_clinical_species_r[1:12,],
         scale = "none",cluster_rows=F, cluster_cols = F,
         color = colorRampPalette(colors = c("#4DBBD5FF","white","#E64B35FF"))(50),
         display_numbers =F,fontsize = 15,face="bold",angle=45)+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))
dev.off()
