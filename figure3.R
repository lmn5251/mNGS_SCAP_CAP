setwd("D:\\re-shannon-index\\SCAP\\model")
m = ifelse(as.matrix(SCAP_DNA)>0 ,'1','0')
m=as.matrix(m)
write.table(m,file="m.txt",sep="\t")
m=read.table("m.txt",header=T,sep="\t")
pheno_SCAP=readxl::read_xlsx("pheno_SCAP-1.xlsx") 
diversity=read.table("diversity_SCAP.txt",sep="\t")
p1=m[,which(colnames(m)=="Acinetobacter_baumannii")]
p2=m[,which(colnames(m)=="Staphylococcus_aureus")]
p3=m[,which(colnames(m)=="Streptococcus_pneumoniae")]
p4=m[,which(colnames(m)=="Rothia_mucilaginosa")]
p5=m[,which(colnames(m)=="Haemophilus_influenzae")]
p6=m[,which(colnames(m)=="Streptococcus_salivarius")]
p7=m[,which(colnames(m)=="Prevotella_oris")]
p8=m[,which(colnames(m)=="Prevotella_pallens")]
p9=m[,which(colnames(m)=="Streptococcus_mitis")]
p10=m[,which(colnames(m)=="Gemella_haemolysans")]
p11=m[,which(colnames(m)=="Moraxella_catarrhalis")]
p12=m[,which(colnames(m)=="Corynebacterium_striatum")]
p13=m[,which(colnames(m)=="Pseudomonas_fluorescens")]
p14=m[,which(colnames(m)=="Elizabethkingia_meningoseptica")]
p15=m[,which(colnames(m)=="Parvimonas_micra")]
p16=m[,which(colnames(m)=="Porphyromonas_endodontalis")]
p17=m[,which(colnames(m)=="Streptococcus_constellatus")]
p18=m[,which(colnames(m)=="Streptococcus_anginosus")]
p19=m[,which(colnames(m)=="Streptococcus_mutans")]
p20=m[,which(colnames(m)=="Streptococcus_salivarius")]
p21=pheno_SCAP[,17]
p22=cbind(p1,p2);p23=cbind(p22,p3);p24=cbind(p23,p4);p25=cbind(p24,p5);p26=cbind(p25,p6);
p27=cbind(p26,p7);p28=cbind(p27,p8);p29=cbind(p28,p9);p30=cbind(p29,p10);p31=cbind(p30,p11);
p32=cbind(p31,p12);p33=cbind(p32,p13);p34=cbind(p33,p14);p35=cbind(p34,p15);p36=cbind(p35,p16);
p37=cbind(p36,p17);p38=cbind(p37,p18);p39=cbind(p38,p19);p40=cbind(p39,p20);p41=cbind(p40,p21)
library(glmnet) 
library(foreign)
DNA.SDSMRN.species=read.table("DNA.SDSMRN.179.filter.species.txt",header=F,row.names=1)
a=match(as.matrix(pheno_SCAP[,3]),as.matrix(DNA.SDSMRN.species[1,])) 
library(caret) 
set.seed(123456) #set seed
pheno_SCAP_2=readxl::read_xlsx("pheno_SCAP-2.xlsx") 
q1=pheno_SCAP_2[!is.na(a),38:56]#clinical indicators
q2=pheno_SCAP_2[!is.na(a),22]#APACHEII
q3=pheno_SCAP_2[!is.na(a),36:37]#clinical indicators
q4=diversity
q5=cbind(p41)
q6=pheno_SCAP_2[!is.na(a),64]#outcome
q7=cbind(q1,q2);q8=cbind(q7,q3);q9=cbind(q8,q4);q10=cbind(q9,q5);q11=cbind(q10,q6)
f1 = glmnet(data.matrix(q11[,1:48]),data.matrix(q11[,49]), family="binomial", nlambda=1000, alpha=1) 
#set trainning and validation sets
sam<- createDataPartition(q11[,49], p = .7,list = FALSE) 
head(sam)
k=data.matrix(q11[,1:48])
k1=data.matrix(q11[,49])
train <- data.matrix(k[sam,])
test <- data.matrix(k[-sam,])
train_meta <-data.matrix( k1[sam,])
test_meta <- data.matrix(k1[-sam,])
x=train
y = train_meta
cv_fit <- cv.glmnet(x=x, y=y,alpha = 1,family="binomial", nlambda=1000)
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "lambda.png", width = 7, height = 7,units = "in",dpi = 800,bg="transparent")
plot(cv_fit)
dev.off()
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
lasso.prob <- predict(cv_fit, newx=test, s=c(cv_fit$lambda.min),type="response")
lasso.prob1 <- predict(cv_fit, newx=test, s=c(cv_fit$lambda.1se),type="response")

choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
choose_gene_1se
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_min

testdata=as.data.frame(cbind(test,test_meta))
traindata=as.data.frame(cbind(train,train_meta))
mod3=glm(train_meta~Hb+WBC+NEU+AST+Cr+PCT+PT+APTT+use_of_immunosuppressor+APACHE+SOFA1+SOFA_delta +simpson+p3+p11+p12+p15+p16, family = binomial(),data=traindata)
summary(mod3)
exp(confint(mod3))  #CI
exp(coef(mod3))   ##OR
mod=glm(train_meta~APTT+SOFA_delta+simpson+p3+use_of_immunosuppressor, family = binomial(),data=traindata)
summary(mod)
exp(confint(mod))  #CI
exp(coef(mod))   ##OR
coef(mod)
#univariate analysis
mod4=glm(train_meta~Hb, family = binomial(),data=traindata)
summary(mod4)
exp(confint(mod4))  #CI
exp(coef(mod4))   
#logit probability
setwd("D:\\figure\\v3_transparent")
ggsave(filename = "logit-prob.png", width =5, height = 7,units = "in",dpi = 800,bg="transparent")
logit.prob1 <- predict(mod, newx=q11[,1:48],,type="response")
re2=cbind(q11[,49] ,logit.prob1)
head(re2)
re2=as.data.frame(re2)
colnames(re2)=c("event","probability")
re2$event=as.factor(re2$event)
library(ggpubr) 
compared_list=c("0","1")
p2= ggboxplot(re2, x = "event", y = "probability",
              color = "event", palette = "jco",
              add = "jitter")+ stat_compare_means()+
  theme(axis.text.x = element_text(face="bold", size=14,hjust=1),
        axis.text.y = element_text(face="bold",  size=20),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=20))+
  scale_color_aaas()
p2
dev.off()

j1=cbind(train[,16],train[,22])
j2=cbind(j1,train[,25]);
j3=cbind(j2,train[,30])
j4=cbind(j3,train[,19])
train_1=as.data.frame(j4)
colnames(train_1)=c("APTT","SOFA_delta","simpson","p3","use_of_immunosuppressor")
pre <- predict(mod,type='response',newdata=train_1) 
library(pROC)
modelroc1 <- roc(train_meta,pre) 
k1=cbind(test[,16],test[,22])
k2=cbind(k1,test[,25]);k3=cbind(k2,test[,30]);k4=cbind(k3,test[,19])
test_1=as.data.frame(k4)
colnames(test_1)=c("APTT","SOFA_delta","simpson","p3","use_of_immunosuppressor")
pre <- predict(mod,type='response',newdata=test_1) 
modelroc <- roc(test_meta,pre)
#plot
setwd("D:\\figure\\v3_transparent")
Cairo::CairoPNG(filename = "AUC.png", width = 7, height = 7,units = "in",dpi = 800,bg="transparent")
plot(modelroc, print.auc=F, auc.polygon=F,
     col='#E18727FF', xlab= "1-Specificity",ylab="Sensitivity",
     lwd = 2, cex.main=1.3, cex.lab=1.7, cex.axis=1.6, font=1.5,legacy.axes=T)
plot(modelroc1, print.auc=F, auc.polygon=F,
     col='#20854EFF', xlab= "1-Specificity",ylab="Sensitivity",
     lwd = 2, cex.main=1.3, cex.lab=1.7, cex.axis=1.6, font=1.5,legacy.axes=T,add=T)
pal_nejm("default")(8)
dev.off()
