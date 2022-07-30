setwd("D:\\SCAP") 
install.packages("vegan")
library(vegan)  
install.packages("ggsignif")
library(ggsignif) 
DNA.SDSMRN.species=read.table("DNA.SDSMRN.179.filter.species.txt",header=F,row.names=1)
RNA.SDSMRN.species=read.table("RNA.SDSMRN.179.filter.species.txt",header=F,row.names=1)
a=match(as.matrix(pheno_SCAP[,4]),as.matrix(RNA.SDSMRN.species[1,]))
b=match(as.matrix(pheno_CAP[,4]),as.matrix(RNA.SDSMRN.species[1,]))
c=match(as.matrix(pheno_CAP[,3]),as.matrix(DNA.SDSMRN.species[1,]))
d=match(as.matrix(pheno_SCAP[,3]),as.matrix(DNA.SDSMRN.species[1,]))
pheno_SCAP_DNA=DNA.SDSMRN.species[,d]
pheno_CAP_DNA=DNA.SDSMRN.species[,na.omit(c)]
SCAP_group_RNA=as.matrix(pheno_SCAP[which(!is.na(a)),1])
CAP_group_RNA=as.matrix(pheno_CAP[which(!is.na(b)),1])
SCAP_group_DNA=as.matrix(pheno_SCAP[which(!is.na(d)),1])
CAP_group_DNA=as.matrix(pheno_CAP[which(!is.na(c)),1])
##########################################################################################################
# SCAP  vs CAP
pheno_all=readxl::read_xlsx('pheno-all.xlsx')
all_count=read.table("dds_normalized_counts.txt",header=F,row.names=1)
a=match(as.matrix(pheno_all[,7]),as.matrix(all_count[1,]))
all_count=all_count[,a]
colnames(all_count)=all_count[1,]
all_count=as.matrix(all_count[-1,])
all_count=apply(all_count,1,as.numeric)

immune_gene_1=readxl::read_xlsx("immune_gene.express.xlsx",col_names=F)
b=match(as.matrix(pheno_all[,7]),as.matrix(immune_gene_1[1,]))
immune_gene=as.matrix(immune_gene_1[,b])
rownames(immune_gene)=as.matrix(immune_gene_1[,1])
colnames(immune_gene)=immune_gene[1,]
immune_gene=as.matrix(immune_gene[-1,])
immune_gene=apply(immune_gene,1,as.numeric)

condition=as.matrix(pheno_all[which(!is.na(a)),62])
condition=gsub("1","SCAP_balf",condition)
condition=gsub("2","SCAP_sputum",condition)
condition=gsub("3","CAP_balf",condition)
condition=gsub("4","CAP_sputum",condition)
pheno_survive=condition
library(edgeR)
#（1）
dgelist <- DGEList(counts = t(all_count),group=condition) 
#（2）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, keep.lib.sizes = FALSE]
#（3）
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#DEGs analysis
design <- model.matrix(~0+factor(pheno_survive))
#（1）
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）
fit <- glmQLFit(dge, design, robust = TRUE) 
colnames(design)=c("CAP_balf","CAP_sputum","SCAP_balf","SCAP_sputum")
fit1= glmQLFTest(fit, contrast=c(-1,0,1,0))
fit2= glmQLFTest(fit, contrast=c(0,-1,0,1))
lrt <- topTags(fit1, n=Inf,adjust.method="BH")
lrt1 <- topTags(fit2, n=Inf,adjust.method="BH")
DEG=as.data.frame(lrt) 
DEG1=as.data.frame(lrt1)
length(which(DEG$FDR < 0.05 & abs(DEG$logFC)>1))
length(which(DEG1$FDR < 0.05 & abs(DEG1$logFC)>2))
a=DEG[which(DEG$FDR < 0.05 &  abs(DEG$logFC)>1),]
b=DEG1[which(DEG1$FDR < 0.05 &  abs(DEG1$logFC)>2),]
k=intersect(rownames(a),rownames(b)) #intersected DEGs 
write.table(a,file="a_SCAP_CAP_balf.txt",sep="\t")
write.table(b,file="b_SCAP_CAP_SP.txt",sep="\t")
library(clusterProfiler)
#fig 4b and 4c
ego_BP <- enrichGO(gene          = row.names(b),   
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",  #
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
setwd("D:\\\figure\\v3_transparent")
png(filename = "SCAP_CAP_GO_sputum.png", width = 13, height = 10,units = "in",res = 800,bg="transparent")
enrichplot::cnetplot(ego_BP,circular=FALSE,colorEdge = TRUE,showCategory =5)
dev.off()
write.csv(ego_ALL,'enrichGO_all.csv')
View(ego_BP)
ggplot(data=ego_BP, aes(x=Description,y=Count)) + 
  geom_bar(stat="identity", width=0.8,fill='salmon1') + 
  coord_flip() +  xlab("GO term") + ylab("Num of Genes") + 
  theme_bw()+
  theme_classic(base_size = 20)+
  theme(axis.text.x=element_text(vjust=1,size=20,face = "bold"))+
  theme(axis.text.y=element_text(vjust=1,size=20,face = "bold"))+
  theme(axis.title.x =element_text(size=20,face = "bold"), axis.title.y=element_text(size=20,face = "bold"))
#DEGs FC of SCAP vs CAP fig4a
p1=which(rownames(a)=="OPRPN");o1=a[p1,] 
p2=which(rownames(a)=="HTN1");o2=a[p2,] 
p3=which(rownames(a)=="HTN3");o3=a[p3,] 
p4=which(rownames(a)=="LCN1");o4=a[p4,] 
p5=which(rownames(a)=="FDCSP");o5=a[p5,] 
p6=which(rownames(a)=="STATH");o6=a[p6,]
y1=which(rownames(b)=="OPRPN");j1=b[y1,]
y2=which(rownames(b)=="HTN1");j2=b[y2,]
y3=which(rownames(b)=="HTN3");j3=b[y3,]
y4=which(rownames(b)=="LCN1");j4=b[y4,]
y5=which(rownames(b)=="FDCSP");j5=b[y5,]
y6=which(rownames(b)=="STATH");j6=b[y6,]
e1=rbind(o1,o2);e2=rbind(e1,o3);e3=rbind(e2,o4);e4=rbind(e3,o5);e5=rbind(e4,o6);e6=rbind(e5,o7)
w1=rbind(j1,j2);w2=rbind(w1,j3);w3=rbind(w2,j4);w4=rbind(w3,j5);w5=rbind(w4,j6);w6=rbind(w5,j7)
write.table(e6,file="e6.txt",sep="\t")
write.table(w6,file="w6.txt",sep="\t")
setwd("D:\\SCAP") 
e6=read.table("e6.txt",sep="\t",header=T)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "DEG_SCAP_CAP_FC.png", width = 5, height = 7,units = "in",dpi = 800,bg="transparent")
p=ggplot(e6, aes(x=gene, y=logFC, fill=sample)) +
  geom_bar(stat="identity", color="black", width=.6,position = position_dodge(width=0.6)) 
p+theme_bw(base_size=18)+scale_fill_aaas()+
  theme(axis.text.x = element_text(face="bold", size=14,angle=60,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))
dev.off()
annotation_6_gene=read.table("annotation_6_gene.txt",sep="\t",header=T,row.names=1)
gene_6=read.table("6-gene.txt",header=T,sep="\t",row.names=1)
gene_6=log2(gene_6)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "DEG_SCAP_CAP_pheatmap.png", width = 6.0, height = 7.5,units = "in",res = 800,bg="transparent")
pheatmap(gene_6,annotation_row =data.frame(annotation_6_gene),cluster_rows=F, cluster_cols = F,
         color = colorRampPalette(colors = c("lightblue","white","red"))(150),
         display_numbers =F,fontsize = 12.1,face="bold",border_color = NA,show_rownames = F)+
  theme_bw()+ 
  theme(axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold",  size=15),
        axis.title.x=element_text(size=15,face="bold"),axis.title.y=element_text(size=15,face="bold"))+
  scale_fill_aaas()  
dev.off()  