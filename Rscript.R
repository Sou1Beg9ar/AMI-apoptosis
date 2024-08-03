#############数据集预处理##################

GSE48060=read.table('GSE48060_series_matrix.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

GSE66360=read.table('GSE66360_series_matrix.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

GSE97320=read.table('GSE97320_series_matrix.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

GPL=read.table('GPL570-55999.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

colnames(GSE48060)[1]='ID'
colnames(GSE66360)[1]='ID'  ##改为ID
colnames(GSE97320)[1]='ID'

GSE48060[,1]=gsub('["]', '', GSE48060[,1])
GSE66360[,1]=gsub('["]', '', GSE66360[,1])     ##去掉“” 
GSE97320[,1]=gsub('["]', '', GSE97320[,1])

GPL_1=GPL[,c(1,11)]       ##只取ID和Genesymbol

GSE48060=merge(GSE48060,GPL_1,by='ID')
GSE66360=merge(GSE66360,GPL_1,by='ID')  #根据ID的名称把基因名对号入座，合并到最后一列
GSE97320=merge(GSE97320,GPL_1,by='ID')


GSE48060=GSE48060[!duplicated(GSE48060[,54]),]   #去掉Genesymbol里的重复值
GSE66360=GSE66360[!duplicated(GSE66360[,101]),]
GSE97320=GSE97320[!duplicated(GSE97320[,8]),]

row.names(GSE48060)=GSE48060[,54]
row.names(GSE66360)=GSE66360[,101]  #修改行名
row.names(GSE97320)=GSE97320[,8]

GSE48060=GSE48060[,-54]
GSE66360=GSE66360[,-101]  #去掉最后一列的Genesymbol
GSE97320=GSE97320[,-8]

GSE48060=GSE48060[,-1]
GSE66360=GSE66360[,-1]   #去掉没用的ID
GSE97320=GSE97320[,-1]

write.csv(GSE48060, file = 'GSE48060.csv')
write.csv(GSE66360, file = 'GSE66360.csv')  #输出
write.csv(GSE97320, file = 'GSE97320.csv')

#引用包
###########去批次效应#################

BiocManager::install('sva')


library(limma)
library(sva)
outFile="merge.txt"       #输出文件
setwd("D:\\生信工程\\李梦瑶\\原始数据\\训练集")      #设置工作目录

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
  if(file==outFile){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()  #做数据表
batchType=c()  #做向量

for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  rt=normalizeBetweenArrays(rt)
  
  #数据合并
  if(i==1){
    allTab=rt[row.names(rt)==interGenes,]
  }else{
    allTab=cbind(allTab, rt[row.names(rt)==interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}



#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)


##########差异基因做热图##############
####导入差异基因

DEG_1.5=read.table('DEG1.5-select.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

select_DEG=DEG_1.5[1:50,]   ##根据P值升序排列

merge_matrix=read.table('merge.txt', sep = '\t',quote = "",fill = T,comment.char = "!",header = T)
row.names(merge_matrix)=merge_matrix[,1]   ##导入merge后修改行名
merge_matrix=merge_matrix[,-1]

##在merge中挑选出差异基因和其表达量
choose_matrix=merge_matrix[which(row.names(merge_matrix)%in%select_DEG$Tag),]
write.csv(choose_matrix,file = 'heatmap_matrix.csv')

#################

GSE775=read.table('GSE775_series_matrix.txt',sep = '\t',quote = "",fill = T,comment.char = "!",header = T)

rownames(GSE775)=GSE775[,1]
GSE775=GSE775[,-1]
GSE775=log2(GSE775+1)

write.csv(GSE775,file = 'GSE775matrix.csv')


BiocManager::install('sva')


library(limma)
library(sva)
outFile="merge.txt"       #输出文件
setwd("D:\\生信工程\\李梦瑶\\原始数据\\训练集")      #设置工作目录

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
  if(file==outFile){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()  #做数据表
batchType=c()  #做向量

for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  rt=normalizeBetweenArrays(rt)
  
  #数据合并
  if(i==1){
    allTab=rt[row.names(rt)==interGenes,]
  }else{
    allTab=cbind(allTab, rt[row.names(rt)==interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}



#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)


write.csv(GSE775,file = "GSE775.csv")

install.packages("glmnet")


set.seed(123)
library(glmnet)                   #???ð?
inputFile="diffGeneExp.txt"       #?????ļ?
setwd("D:\\生信工程\\缺血再灌焦亡\\6.机器学习筛选\\GSE160516")       #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#????ģ??
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#????ɸѡ??????????
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

pdf(file="lambda.pdf",width=6,height=5.5)

plot(fit, xvar = "lambda", label = TRUE)
dev.off()

pdf(file="dev.pdf",width=6,height=5.5)

plot(fit, xvar = "dev", label = TRUE)
dev.off()

install.packages("e1071")
install.packages("kernlab")
install.packages("caret")


#???ð?
library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="diffGeneExp.txt"        #?????ļ?
setwd("D:\\生信工程\\缺血再灌焦亡\\6.机器学习筛选\\GSE160516")      #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

#SVM-RFE????
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#????ͼ??
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#??ע??????֤??????С?ĵ?
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#????ѡ???Ļ???
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

install.packages("randomForest")


#引用???
library(randomForest)
set.seed(123)

inputFile="diffGeneExp.txt"       #输入文件
setwd("D:\\生信工程\\缺血再灌焦亡\\6.机器学习筛选\\GSE160516")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

#随机森林???
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的???
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要???
importance=importance(x=rf2)

#绘制基因的重要性图
pdf(file="geneImportance.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

#挑选疾病特征基???
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>=2])      #挑选重要性评分大???2的基???
#rfGenes=names(rfGenes[1:5])           #挑选重要性评分最高的5个基???
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

#输出重要基因的表达量
sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)

install.packages("rms")
install.packages("rmda")

remove.packages("Hmisc")
remove.packages("rmda")

install.packages("Hmisc")
install.packages("rmda")


#引用包
library(rms)
library(rmda)

inputFile="train_data.txt" #输入文件
inputFile="test_data.txt"

setwd("D:\\生信工程\\李梦瑶\\17.临床决策曲线")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
write.csv(data,file = 'cluster.csv')
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type ~ CDKN1A+BCL10+PMAIP1+IL1B+GNA15+CD14, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.1,0.3,0.5,0.7,0.9),
              lp=F, funlabel="Risk of Disease")
#输出列线图
pdf("Nom.pdf", width=9, height=6)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=6, height=6)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

#绘制决策曲线
rt$Type=ifelse(rt$Type=="CON", 0, 1)
dc=decision_curve(Type ~ CDKN1A+BCL10+PMAIP1+IL1B+GNA15+CD14, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
#输出DCA图形
pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
                    curve.names="Crucial genes",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

#绘制临床影响曲线
pdf(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("red", "blue"))
dev.off()

