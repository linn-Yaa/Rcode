#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("kernelshap")
#install.packages("pROC")
#install.packages("shapviz")
#install.packages("xgboost")
#install.packages("klaR")


#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
#names(getModelInfo())

set.seed(12345)      #设置种子

bestMethod="Stepglm[both]+NaiveBayes"     #机器学习方法
shapMethod="nb"                           #SHAP分析方法
inputFile="merge.normalize.txt"        #表达数据文件
geneFile="model.genes.txt"              #模型基因的文件


#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取模型基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
geneRT=geneRT[geneRT$algorithm==bestMethod,]
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

#获取样品分组信息(对照组和实验组)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data=as.data.frame(data)
data$Type=group

#机器学习的模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
data$Type=ifelse(data$Type=="Control", 0, 2)
data$Type <- as.factor(data$Type)
model=train(Type ~., data = data, method = shapMethod, trControl=control)

#计算SHAP值
custom_predict <- function(object, newdata) {
  # 返回第二类的概率（通常是正类的概率）
  predict(object, newdata, type = "prob")[, 1]
}
fit <- permshap(model, data[,-ncol(data)], 
                pred_fun = custom_predict,
                bg_X = data[,-ncol(data)][sample(1:nrow(data), 200), ])
fit=permshap(model, data[,-ncol(data)])          #基因比较少的时候运行这命令
#fit=kernelshap(model, data[,-ncol(data)])     #基因比较多的时候运行这命令
shp <- shapviz(fit, X_pred = data[,-ncol(data)], X=data[,-ncol(data)], interactions=T)
#根据贡献程度对基因进行排序
important=sort(colMeans(abs(shp$S)), decreasing=T)
showVars=names(important)
write.table(important, file="important.genes1.txt", sep="\t", quote=F, col.names=F)

#绘制柱状图
pdf(file="barplot1.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE)+theme_bw()
dev.off()

#绘制蜂群图
pdf(file="bee1.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE)+theme_bw()
dev.off()

#散点图(依赖图)
pdf(file="dependence1.pdf", width=9, height=6)
sv_dependence(shp, v = showVars)+theme_bw()
dev.off()

#瀑布图
pdf(file="waterfall1.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 10)
dev.off()

#单样品力图
pdf(file="force1.pdf", width=9, height=5)
sv_force(shp, row_id = 10)
dev.off()



