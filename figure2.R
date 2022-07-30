setwd("D:\\community_diversity\\SCAP_CAP") 
install.packages("vegan")
library(vegan)  
install.packages("ggsignif")
library(ggsignif) 
pheno_SCAP=readxl::read_xlsx("pheno_SCAP-1.xlsx") 
d=match(as.matrix(pheno_SCAP[,3]),as.matrix(DNA.SDSMRN.species[1,]))
d=as.numeric(na.omit(d))
pheno_SCAP_DNA=DNA.SDSMRN.species[,d]
SCAP_group_DNA=as.matrix(pheno_SCAP[which(!is.na(d)),1])
##########################################################################################################
#metaDNA
pheno_SCAP=readxl::read_xlsx('pheno_SCAP-1.xlsx') 
a=match(as.matrix(pheno_SCAP[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
DNA=DNA.SDSMRN.species[,na.omit(a)]
condition=as.matrix(pheno_SCAP[which(!is.na(a)),61])
condition=gsub("1","BALF_nonsurvived",condition)
condition=gsub("2","Sputum_nonsurvived",condition)
condition=gsub("3","BALF_survived",condition)
condition=gsub("4","Sputum_survived",condition)
DNA=as.data.frame(t(DNA)) 
rownames(DNA)=DNA[,1]
DNA=DNA[,-1]
DNA_1=as.data.frame(apply(DNA,2,as.numeric))
alpha_diversity <- function(x){
  library(vegan) 
  shannon = diversity(x,index="shannon") #shannon
  simpson = diversity(x, index = "simpson") #simpson
  Chao1 = estimateR(x)[2,] #Chao1 index
  ACE = estimateR(x)[4,] #ACE index
  return(a)
}
diversity=alpha_diversity(DNA_1)
data_ggplot=data.frame(diversity[,1])##1: shannon,2:simpson;3:Chao1; 4:ACE
data_ggplot_DNA=data_ggplot
data_ggplot=cbind(condition,data_ggplot_DNA)
colnames(data_ggplot)=c("Group","data_norm_shannon")
compared_list = list(c("BALF_nonsurvived", "BALF_survived"),c("Sputum_nonsurvived","Sputum_survived")) 
library(ggsci)
library(ggplot2)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_survived_nonsurvived_shannon.png", width = 7.5, height = 8,units = "in",dpi = 800,bg="transparent")
alpha_boxplot=ggplot(data_ggplot, aes(x=Group, y=data_norm_shannon, fill=Group))+
  geom_boxplot()+
  labs(x="Group", y="Shannon Index",face="bold",size=18)+
  theme(plot.title=element_text(hjust=0), legend.title=element_blank())+
  theme_bw(base_size=18)
alpha_boxplot  + 
  geom_signif(comparisons = compared_list, test = wilcox.test, step_increase = 0.1, map_signif_level=T,size=0.8,textsize=6)+
  theme(axis.text.x = element_text(face="bold", size=14,angle=45,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=20))+
  geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=1)    +
  scale_fill_npg()
dev.off()
data_ggplot
write.table(diversity,file="diversity_SCAP.txt",sep="\t")
##########################################################################################################################################################
#the relative abundances of each species when compaaring SCAP nonsurvived with SCAP survived patients
m=as.matrix(apply(DNA_1,1,sum))
SCAP_DNA=matrix(0,136,2871)
for (i in 1:136) {
  for (j in 1:2871) {
    SCAP_DNA[i,j]=DNA_1[i,j]/m[i,]
  } 
}
colnames(SCAP_DNA)=colnames(DNA_1) 
condition=as.matrix(pheno_SCAP[which(!is.na(a)),61])
condition=gsub("1","BALF_nonsurvived",condition)
condition=gsub("2","Sputum_nonsurvived",condition)
condition=gsub("3","BALF_survived",condition)
condition=gsub("4","Sputum_survived",condition)
mean_data=matrix(0,10,2871) 
colnames(mean_data)=colnames(SCAP_DNA)
for (i in 1:2871) {
  mean_data[1,i]=mean((SCAP_DNA[which(condition=="BALF_nonsurvived"),i]))
  mean_data[2,i]=mean((SCAP_DNA[which(condition=="BALF_survived"),i]))
  mean_data[3,i]=mean((SCAP_DNA[which(condition=="Sputum_nonsurvived"),i]))
  mean_data[4,i]=mean((SCAP_DNA[which(condition=="Sputum_survived"),i]))
  mean_data[5,i]=sd((SCAP_DNA[which(condition=="BALF_nonsurvived"),i]))
 mean_data[6,i]=sd((SCAP_DNA[which(condition=="BALF_survived"),i]))
  mean_data[7,i]=sd((SCAP_DNA[which(condition=="Sputum_nonsurvived"),i]))
  mean_data[8,i]=sd((SCAP_DNA[which(condition=="Sputum_survived"),i]))
  mean_data[9,i]=mean_data[1,i]-mean_data[2,i] #BALF_nonsurvived -  BALF_survived
  mean_data[10,i]=mean_data[3,i]-mean_data[4,i] #Sputum_nonsurvived  -   Sputum_survived
}
p_data=matrix(0,2871,2) 
for (i in 1:2871) {
  p_data[i,1]=wilcox.test(as.matrix(SCAP_DNA[which(condition=="BALF_nonsurvived"),i]),as.matrix(SCAP_DNA[which(condition=="BALF_survived"),i]))$p.value
  p_data[i,2]=wilcox.test(as.matrix(SCAP_DNA[which(condition=="Sputum_nonsurvived"),i]),as.matrix(SCAP_DNA[which(condition=="Sputum_survived"),i]))$p.value
} 
write.table(mean_data,file="mean_data.txt",sep="\t")
write.table(p_data,file="p_data.txt",sep="\t")
#############################################################################################
#BALF top10 species that the relative abundances are significanltly when comparing SCAP nonsurvived with survived ones 
setwd("D:\\re-shannon-index\\SCAP")
top10=readxl::read_xlsx("balf_survive_death_top10_pvalue.xlsx")
top10$type=factor(top10$type,levels = c("Acinetobacter_baumannii",
                                        "Staphylococcus_aureus",
                                        "Streptococcus_pneumoniae",
                                        "Rothia_mucilaginosa",
                                        "Haemophilus_influenzae",
                                        "Streptococcus_salivarius",
                                        "Prevotella_oris",
                                        "Prevotella_pallens",
                                        "Streptococcus_mitis",
                                        "Gemella_haemolysans"
)) 
library(ggpubr)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_top10_balf_significant.png", width = 9, height = 6,units = "in",dpi = 800,bg="transparent")
p=ggplot(top10, aes(x=type, y=percentage, fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  #geom_signif(comparisons = list(c("SCAP","CAP")),test = t.test, step_increase = 0.0, map_signif_level=T,  y_position=0.5)
  stat_compare_means(comparisons =list(c("BALF-nonsurvived","BALF-survived")))
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
############################################################################################################################
#Sputum top10 species that the relative abundances are significanltly when comparing SCAP nonsurvived with survived ones 
setwd("D:\\SCAP\\important_results")
SCAP_sputum=readxl::read_xlsx("sputum_survive_death_top10_pvalue.xlsx")
SCAP_sputum$type=factor(SCAP_sputum$type,levels = c("Moraxella_catarrhalis",
                                                    "Corynebacterium_striatum",
                                                    "Pseudomonas_fluorescens",
                                                    "Elizabethkingia_meningoseptica",
                                                    "Parvimonas_micra",
                                                    "Porphyromonas_endodontalis",
                                                    "Streptococcus_constellatus",
                                                    "Streptococcus_anginosus",
                                                    "Streptococcus_mutans",
                                                    "Streptococcus_salivarius"))
library(ggpubr)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "SCAP_top10_sputum_significant.png", width = 9, height = 6,units = "in",dpi = 800,bg="transparent")
p=ggplot(SCAP_sputum, aes(x=type, y=percentage, fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  stat_compare_means(comparisons =list(c("BALF-nonsurvived","BALF-survived")))
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()