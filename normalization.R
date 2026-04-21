library(GEOquery)
library(stringr)
library(dplyr)
library(glmnet)
library(tidyverse)
library(tinyarray)
library(hgu133plus2.db)
library(tidyr)
library(sva)
library(patchwork)
library(ggplot2)
library(illuminaHumanv4.db)
library(limma)

##数据集1
geoID="GSE14520"
eSet1 <- getGEO("GSE14520", 
                destdir = '.', 
                getGPL = F)
eSet1 = eSet1[[1]]
exp1 <- exprs(eSet1)
# exp1 = log2(exp1+1)
# exp1=limma::normalizeBetweenArrays(exp1)
if(any(is.na(exp1))) {
  exp1 <- exp1[complete.cases(exp1), ]}
exp1 <- scale(t(exp1))
exp1 <- t(exp1)
head(exp1[, 1:5])
pd1 <- pData(eSet1)
p1 = identical(rownames(pd1),colnames(exp1));p1
gpl_number1 <- eSet1@annotation
ids1 <- toTable(hgu133plus2SYMBOL)
ids1 = ids1[!duplicated(ids1$symbol),]
exp1 <- as.data.frame(exp1)
exp1 <- mutate(exp1,probe_id=rownames(exp1))
exp1 <- inner_join(exp1,ids1,by="probe_id")
rownames(exp1) <- exp1$symbol
exp1 <- exp1[, -(446:447) ]

pd1$characteristics_ch1[pd1$characteristics_ch1 == 'Tissue: Liver Tumor Tissue'] <- '1'
pd1$characteristics_ch1[pd1$characteristics_ch1 == 'Tissue: Liver Non-Tumor Tissue'] <- '0'
pd1$characteristics_ch1[pd1$characteristics_ch1 == 'tissue: Liver Tumor Tissue'] <- '1'
pd1$characteristics_ch1[pd1$characteristics_ch1 == 'tissue: Liver Non-Tumor Tissue'] <- '0'
pd1 <- subset(pd1, select = characteristics_ch1)
pd1 <- mutate(pd1,sample=rownames(pd1))
con_pd1 <- as.matrix(pd1[pd1$characteristics_ch1 == 0, -1])
tumor_pd1 <- as.matrix(pd1[pd1$characteristics_ch1 == 1, -1])
conName1=gsub("^ | $", "", as.vector(con_pd1[,1]))
tumorName2=gsub("^ | $", "", as.vector(tumor_pd1[,1]))
conData=exp1[,conName1]
tumorData=exp1[,tumorName2]
data1=cbind(conData,tumorData)
conNum=ncol(conData)
tumorNum=ncol(tumorData)
Type=c(rep("Control",conNum),rep("Tumor",tumorNum))
outData1=rbind(id=paste0(colnames(data1),"_",Type),data1)
write.table(outData1,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

##数据集2
geoID="GSE36376"
eSet2 <- getGEO("GSE36376", 
                destdir = '.', 
                getGPL = F)
eSet2 = eSet2[[1]]
exp2 <- exprs(eSet2)
if(any(is.na(exp2))) {
  exp2 <- exp2[complete.cases(exp2), ]}
exp2 <- scale(t(exp2))
exp2 <- t(exp2)
head(exp2[, 1:5])
# exp2 = log2(exp2+1)
# exp2=limma::normalizeBetweenArrays(exp2)
pd2 <- pData(eSet2)
p2 = identical(colnames(exp2),rownames(pd2),);p2
gpl_number2 <- eSet2@annotation
find_anno(gpl_number2)
ids2 <- toTable(illuminaHumanv4SYMBOL)
ids2 = ids2[!duplicated(ids2$symbol),]
exp2 <- as.data.frame(exp2)
exp2 <- mutate(exp2,probe_id=rownames(exp2))
exp2 <- inner_join(exp2,ids2,by="probe_id")
rownames(exp2) <- exp2$symbol
exp2 <- exp2[, -(434:435) ]

# sum(is.na(exp2))
# which(is.na(exp2), arr.ind = TRUE)
# na_ratio <- mean(is.na(exp2))
# na_ratio
# expr_imputed <- exp2
# expr_imputed <- expr_imputed[, -(434:435) ]
# expr_matrix_numeric <- as.matrix(expr_imputed)
# mode(expr_matrix_numeric) <- "numeric"  # 确保数据类型为数值型
# class(expr_matrix_numeric)# 应该返回 "matrix" 或 "array"
# typeof(expr_matrix_numeric) # 应该返回 "double"
# dim(expr_matrix_numeric)    # 检查维度是否正确
# for(i in 1:nrow(expr_matrix_numeric)) {
#   na_index <- is.na(expr_matrix_numeric[i, ])
#   if(any(na_index)) {
#     expr_matrix_numeric[i, na_index] <- mean(expr_matrix_numeric[i, !na_index], na.rm = TRUE)
#   }
# }
# typeof(expr_matrix_numeric)
# sum(is.na(expr_matrix_numeric))
# exp2 <- expr_matrix_numeric

pd2$characteristics_ch1[pd2$characteristics_ch1 == 'tissue: liver tumor'] <- '1'
pd2$characteristics_ch1[pd2$characteristics_ch1 == 'tissue: adjacent non-tumor liver'] <- '0'
pd2 <- subset(pd2, select = characteristics_ch1)
pd2 <- mutate(pd2,sample=rownames(pd2))
con_pd2 <- as.matrix(pd2[pd2$characteristics_ch1 == 0, -1])
tumor_pd2 <- as.matrix(pd2[pd2$characteristics_ch1 == 1, -1])
conName2=gsub("^ | $", "", as.vector(con_pd2[,1]))
tumorName2=gsub("^ | $", "", as.vector(tumor_pd2[,1]))
conData=exp2[,conName2]
tumorData=exp2[,tumorName2]
data2=cbind(conData,tumorData)
conNum=ncol(conData)
tumorNum=ncol(tumorData)
Type=c(rep("Control",conNum),rep("Tumor",tumorNum))
outData2=rbind(id=paste0(colnames(data2),"_",Type),data2)
write.table(outData2,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

##数据集3
geoID="GSE76427"
eSet3 <- getGEO("GSE76427", 
                destdir = '.', 
                getGPL = F)
eSet3 = eSet3[[1]]
exp3 <- exprs(eSet3)
# exp3 = log2(exp3+1)
# exp3=limma::normalizeBetweenArrays(exp3)
if(any(is.na(exp3))) {
  exp3 <- exp3[complete.cases(exp3), ]}
exp3 <- scale(t(exp3))
exp3 <- t(exp3)
head(exp3[, 1:5])
pd3 <- pData(eSet3)
p3 = identical(rownames(pd3),colnames(exp3));p3
gpl_number3 <- eSet3@annotation
ids3 <- toTable(illuminaHumanv4SYMBOL)
ids3 = ids3[!duplicated(ids3$symbol),]
exp3 <- as.data.frame(exp3)
exp3 <- mutate(exp3,probe_id=rownames(exp3))
exp3 <- inner_join(exp3,ids3,by="probe_id")
rownames(exp3) <- exp3$symbol
exp3 <- exp3[, -(168:169) ]

pd3$characteristics_ch1.2[pd3$characteristics_ch1.2 == 'tissue: primary hepatocellular carcinoma tumor'] <- '1'
pd3$characteristics_ch1.2[pd3$characteristics_ch1.2 == 'tissue: adjacent non-tumor liver tissue'] <- '0'
pd3 <- subset(pd3, select = characteristics_ch1.2)
pd3 <- mutate(pd3,sample=rownames(pd3))
con_pd3 <- as.matrix(pd3[pd3$characteristics_ch1.2 == 0, -1])
tumor_pd3 <- as.matrix(pd3[pd3$characteristics_ch1.2 == 1, -1])
conName3=gsub("^ | $", "", as.vector(con_pd3[,1]))
tumorName3=gsub("^ | $", "", as.vector(tumor_pd3[,1]))
conData=exp3[,conName3]
tumorData=exp3[,tumorName3]
data3=cbind(conData,tumorData)
conNum=ncol(conData)
tumorNum=ncol(tumorData)
Type=c(rep("Control",conNum),rep("Tumor",tumorNum))
outData=rbind(id=paste0(colnames(data3),"_",Type),data3)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

##数据集4
geoID="GSE10143"
eSet4 <- getGEO("GSE10143", 
                destdir = '.', 
                getGPL = F)
eSet4 = eSet4[[1]]
exp4 <- exprs(eSet4)
if(any(is.na(exp4))) {
  exp4 <- exp4[complete.cases(exp4), ]}
exp4 <- scale(t(exp4))
exp4 <- t(exp4)
head(exp4[, 1:5])
# exp4 = log2(exp4+1)
# exp4=limma::normalizeBetweenArrays(exp4)
pd4 <- pData(eSet4)
gpl_number4 <- eSet4@annotation
a = getGEO(gpl_number4,destdir = ".")
b = a@dataTable@table
colnames(b)
ids4 = b[,c("ID","Symbol")]
colnames(ids4) = c("probe_id","symbol")
ids4 = ids4[ids4$symbol!="" & !str_detect(ids4$symbol,"///"),]
ids4 = ids4[!duplicated(ids4$symbol),]
exp4 <- as.data.frame(exp4)
exp4 <- mutate(exp4,probe_id=rownames(exp4))
exp4 <- inner_join(exp4,ids4,by="probe_id")
rownames(exp4) <- exp4$symbol
exp4 <- exp4[, -(388:389) ]

pd4$source_name_ch1[pd4$source_name_ch1 == 'hepatocellular carcinoma'] <- '1'
pd4$source_name_ch1[pd4$source_name_ch1 == 'hepatitis/cirrhotic liver'] <- '0'
pd4 <- subset(pd4, select = source_name_ch1)
pd4 <- mutate(pd4,sample=rownames(pd4))
con_pd4 <- as.matrix(pd4[pd4$source_name_ch1 == 0, -1])
tumor_pd4 <- as.matrix(pd4[pd4$source_name_ch1 == 1, -1])
conName4=gsub("^ | $", "", as.vector(con_pd4[,1]))
tumorName4=gsub("^ | $", "", as.vector(tumor_pd4[,1]))
conData=exp4[,conName4]
tumorData=exp4[,tumorName4]
data4=cbind(conData,tumorData)
conNum=ncol(conData)
tumorNum=ncol(tumorData)
Type=c(rep("Control",conNum),rep("Tumor",tumorNum))
outData=rbind(id=paste0(colnames(data4),"_",Type),data4)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

###GSE98617
geoID="GSE98617"
eSet5 <- getGEO("GSE98617", 
                destdir = '.', 
                getGPL = F)
eSet5 = eSet5[[1]]
exp5 <- exprs(eSet5)
if(any(is.na(exp5))) {
  exp5 <- exp5[complete.cases(exp5), ]}
exp5 <- scale(t(exp5))
exp5 <- t(exp5)
head(exp5[, 1:5])
# exp5 = log2(exp5+1)
# exp4=limma::normalizeBetweenArrays(exp4)
pd5 <- pData(eSet5)
gpl_number5 <- eSet5@annotation
a = getGEO(gpl_number5,destdir = ".")
b = a@dataTable@table
colnames(b)
ids5 = b[,c("ID","ILMN_Gene")]
colnames(ids5) = c("probe_id","symbol")
ids5 = ids5[ids5$symbol!="" & !str_detect(ids5$symbol,"///"),]
ids5 = ids5[!duplicated(ids5$symbol),]
exp5 <- as.data.frame(exp5)
exp5 <- mutate(exp5,probe_id=rownames(exp5))
exp5 <- inner_join(exp5,ids5,by="probe_id")
rownames(exp5) <- exp5$symbol
exp5 <- exp5[, -(50:51) ]

pd5$source_name_ch1[pd5$source_name_ch1 == 'HCC tissue from cirrotic liver'] <- '1'
pd5$source_name_ch1[pd5$source_name_ch1 == 'Non-tumor cirrhotic liver'] <- '0'
pd5 <- subset(pd5, select = source_name_ch1)
pd5 <- mutate(pd5,sample=rownames(pd5))
con_pd5 <- as.matrix(pd5[pd5$source_name_ch1 == 0, -1])
tumor_pd5 <- as.matrix(pd5[pd5$source_name_ch1 == 1, -1])
conName5=gsub("^ | $", "", as.vector(con_pd5[,1]))
tumorName5=gsub("^ | $", "", as.vector(tumor_pd5[,1]))
conData=exp5[,conName5]
tumorData=exp5[,tumorName5]
data5=cbind(conData,tumorData)
conNum=ncol(conData)
tumorNum=ncol(tumorData)
Type=c(rep("Control",conNum),rep("Tumor",tumorNum))
outData=rbind(id=paste0(colnames(data5),"_",Type),data5)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)


###数据整合
library(limma)
library(sva)
files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

for(file in files){
  if(file=="merge.preNorm.txt"){next}
  if(file=="merge.normalize.txt"){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile=="merge.preNorm.txt"){next}
  if(inputFile=="merge.normalize.txt"){next}
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=data
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  
  #数据合并
  if(i==1){
    allTab=rt[interGenes,]
  }else{
    allTab=cbind(allTab, rt[interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#输出合并后的表达数据
outTab=rbind(geneNames=colnames(allTab), allTab)

#对合并后数据进行批次矫正，输出批次矫正后的表达数据
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)

BiocManager::install("hthgu133a.db")
# 数据整合 --------------------------------------------------------------------
##数据集1
eSet1 <- getGEO("GSE14520", 
                destdir = '.', 
                getGPL = F)
eSet1 = eSet1[[1]]
exp1 <- exprs(eSet1)
exp1 = log2(exp1+1)
exp1=limma::normalizeBetweenArrays(exp1)
pd1 <- pData(eSet1)
p1 = identical(rownames(pd1),colnames(exp1));p1
gpl_number1 <- eSet1@annotation
find_anno(gpl_number1)
library(hthgu133a.db)
ids1 <- toTable(hthgu133aSYMBOL)
ids1 = ids1[!duplicated(ids1$symbol),]
exp1 <- as.data.frame(exp1)
exp1 <- mutate(exp1,probe_id=rownames(exp1))
exp1 <- inner_join(exp1,ids1,by="probe_id")
rownames(exp1) <- exp1$symbol

Group1=ifelse(str_detect(pd1$characteristics_ch1,"Tissue: Liver Tumor Tissue"),
              "patient",
              "control")
Group1 = factor(Group1,
                levels = c("control","patient"))
Group1


##数据集2
eSet2 <- getGEO("GSE36376", 
                destdir = '.', 
                getGPL = F)
eSet2 = eSet2[[1]]
exp2 <- exprs(eSet2)
keep <- rowSums(exp2 > 5) >= ncol(exp2) * 0.1
exp2 <- exp2[keep, ]
exp2 = log2(exp2+1)
exp2=limma::normalizeBetweenArrays(exp2)
pd2 <- pData(eSet2)
p2 = identical(colnames(exp2),rownames(pd2),);p2
gpl_number2 <- eSet2@annotation
find_anno(gpl_number2)
ids2 <- toTable(illuminaHumanv4SYMBOL)
ids2 = ids2[!duplicated(ids2$symbol),]
exp2 <- as.data.frame(exp2)
exp2 <- mutate(exp2,probe_id=rownames(exp2))
exp2 <- inner_join(exp2,ids2,by="probe_id")
rownames(exp2) <- exp2$symbol

Group2=ifelse(str_detect(pd2$characteristics_ch1,"tissue: liver tumor"),
              "patient",
              "control")
Group2 = factor(Group2,
                levels = c("control","patient"))
Group2



# sum(is.na(exp2))
# which(is.na(exp2), arr.ind = TRUE)
# na_ratio <- mean(is.na(exp2))
# na_ratio
# expr_imputed <- exp2
# expr_imputed <- expr_imputed[, -(434:435) ]
# expr_matrix_numeric <- as.matrix(expr_imputed)
# mode(expr_matrix_numeric) <- "numeric"  # 确保数据类型为数值型
# class(expr_matrix_numeric)# 应该返回 "matrix" 或 "array"
# typeof(expr_matrix_numeric) # 应该返回 "double"
# dim(expr_matrix_numeric)    # 检查维度是否正确
# for(i in 1:nrow(expr_matrix_numeric)) {
#   na_index <- is.na(expr_matrix_numeric[i, ])
#   if(any(na_index)) {
#     expr_matrix_numeric[i, na_index] <- mean(expr_matrix_numeric[i, !na_index], na.rm = TRUE)
#   }
# }
# typeof(expr_matrix_numeric)
# sum(is.na(expr_matrix_numeric))
# exp3 <- expr_matrix_numeric
# 合并数据集 -------------------------------------------------------------------
##表达量
table(rownames(exp1) %in% rownames(exp2))
common_genes <- intersect(rownames(exp1), rownames(exp2))
exp1 <- exp1[common_genes, ]
exp2 <- exp2[common_genes, ]
table(rownames(exp1) %in% rownames(exp2))
exp1 <- exp1[, -(446:447) ]
exp2 <- exp2[, -(434:435) ]
exp <- cbind(exp1,exp2)
#boxplot(exp)
Group <- c(Group1,Group2)
Group <- factor(Group,levels = c("control","patient"))

batch <- c(rep("GSE14520",445),rep("GSE36376",433))
mod <- model.matrix(~Group)
exp_B <- ComBat(dat = exp,batch = batch,mod = mod,par.prior = TRUE,ref.batch = "GSE36376")
#boxplot(exp_B)

##临床数据
common_cols <- intersect(colnames(pd1), colnames(pd2))
cpd <- rbind(pd1[, common_cols], pd2[, common_cols])
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'Tissue: Liver Tumor Tissue'] <- '1'
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'Tissue: Liver Non-Tumor Tissue'] <- '0'
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'tissue: liver tumor'] <- '1'
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'tissue: adjacent non-tumor liver'] <- '0'
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'tissue: Liver Tumor Tissue'] <- '1'
cpd$characteristics_ch1[cpd$characteristics_ch1 == 'tissue: Liver Non-Tumor Tissue'] <- '0'
pd <- subset(cpd, select = characteristics_ch1)
colnames(pd)[1] <- "status" 
pd <- mutate(pd,sample=rownames(pd))

# PCA ---------------------------------------------------------------------
pca1 <- prcomp(t(exp), scale.=TRUE)
plot_data1 <- data.frame(PC1=pca1$x[,1], PC2=pca1$x[,2], Batch=as.factor(batch))
ggplot(plot_data1, aes(PC1, PC2, color=Batch)) + geom_point() +
  stat_ellipse(level = 0.9)

pca2 <- prcomp(t(exp_B), scale.=TRUE)
plot_data2 <- data.frame(PC1=pca2$x[,1], PC2=pca2$x[,2], Batch=as.factor(batch))
ggplot(plot_data2, aes(PC1, PC2, color=Batch)) + geom_point() +
  stat_ellipse(level = 0.9)

# 热图 ----------------------------------------------------------------------
cg=names(tail(sort(apply(exp_B,1,sd)),100))
n=exp[cg,]

# 直接画热图，对比不鲜明
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col
)

# 按行标准化
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3.5,3.5,length.out = 100)
) 

# 查看annotation_col的结构
str(annotation_col)
head(annotation_col)

# 按分组排序
if ("Group" %in% colnames(annotation_col)) {
  # 按Group列排序
  sorted_order <- order(annotation_col$Group)
} else if ("Type" %in% colnames(annotation_col)) {
  # 按Type列排序
  sorted_order <- order(annotation_col$Type)
} else {
  # 使用第一列排序
  sorted_order <- order(annotation_col[, 1])
}

# 按排序顺序重新排列数据矩阵和注释
n_sorted <- n[, sorted_order]
annotation_col_sorted <- annotation_col[sorted_order, , drop = FALSE]

# 绘制热图
pheatmap(n_sorted,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = annotation_col_sorted,
         scale = "row",
         breaks = seq(-3.5, 3.5, length.out = 100),
         cluster_cols = FALSE  # 禁止列聚类，保持手动排序
)



# 差异分析 --------------------------------------------------------------------
library(limma)
design=model.matrix(~Group)
fit=lmFit(exp_B,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
deg <- mutate(deg,probe_id=rownames(deg))

logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
#filter_up <- subset(deg, P.Value < P.Value_t & logFC > logFC_t) #过滤上调基因
#filter_down <- subset(deg, P.Value < P.Value_t & logFC < -logFC_t) #过滤下调基因
table(deg$change)
write.csv(deg, file = "deg_results.csv")

p <- ggplot(data = deg, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=1, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

# 富集分析 --------------------------------------------------------------------
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)

gene <- deg %>%
  filter(abs(logFC) > logFC_t & P.Value < P.Value_t) %>%
  rownames(.)
genelist <- bitr(gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

#GO分析
ego <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
write.table(ego_res,file="GO.txt",sep="\t",quote=F,row.names = F)

ego_CC <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "CC",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
dotplot(ego_CC)

ego_BP <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "BP",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
dotplot(ego_BP)

ego_MF <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "MF",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
dotplot(ego_MF)


#KEGG分析
kk <- enrichKEGG(gene = genelist$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.05)#可调整
kk_res <- kk@result
write.table(kk_res,file="KEGG.txt",sep="\t",quote=F,row.names = F)
dotplot(kk, showCategory = 20, orderBy = "GeneRatio", label_format=50, color = "pvalue")

# GSEA
gene_list <- deg$logFC
names(gene_list) <- deg$probe_id
gene_list <- sort(gene_list, decreasing = TRUE)
gsea_result <- gseGO(
  geneList     = gene_list,           # 排序的基因列表
  OrgDb        = org.Hs.eg.db,        # 物种数据库
  keyType      = "SYMBOL",            # 明确指定基因ID类型
  ont          = "ALL",                # 本体：BP, MF, CC
  minGSSize    = 10,                  # 最小基因集大小
  maxGSSize    = 500,                 # 最大基因集大小
  pvalueCutoff = 0.05,                # p值阈值
  pAdjustMethod = "BH",               # 多重检验校正方法
  verbose      = FALSE                # 不显示详细过程
) 
sig_results <- gsea_result@result[gsea_result@result$p.adjust < 0.05, ]
sig_results <- sig_results[order(sig_results$NES, decreasing = TRUE), ]
write.csv(sig_results, "gsea_significant_pathways.csv", row.names = FALSE)
# 绘制多个通路的GSEA图
gseaplot2(gsea_result, 
          geneSetID = 1:3,            # 选择前3个通路
          pvalue_table = FALSE)        # 显示p值表格

gseaplot2(gsea_result, geneSetID = 1:10) 



# 绘制网络图（Cnetplot）
cnetplot(gsea_result, showCategory = 5,
         categorySize = "pvalue", 
         foldChange = gene_list) +
  ggtitle("Gene-Concept Network")


















# ####TCGA
# load("TCGA_data.RData")
# data5 <- mutate(combined_data,sample=rownames(combined_data))
# interGenes <- read.table("interGenes.txt", header = FALSE, sep = "")
# interGenes <- interGenes[, 1]  
# data5 <- data5[, interGenes]
# data5$status <- factor(data5$status,
#                          levels = c("0", "1"),  # 原始水平
#                          labels = c("Control", "Tumor"))  # 有效的R变量名
# data5 <- mutate(data5,id=rownames(data5))
# data5$id <- paste(data5[, 162], data5[, 161], sep = "_")  # 使用"_"作为分隔符
# data5 <- data5[, -161 ]
# data5 <- as.data.frame(t(data5))
# data5 <- data5[c(nrow(data5), 1:(nrow(data5)-1)), ]
# write.table(data5,file="TCGA1.normalize.txt",sep="\t",quote=F,col.names=F)

# 准备验证集和测试集 ---------------------------------------------------------------
geneFile="interGenes.txt"      #基因列表文件
#获取目录下所有"normalize.txt"结尾的文件
files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

#读取所有表达数据文件中的基因信息，保存到geneList
for(file in files){
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#获取表达数据的交集基因
interGenes=Reduce(intersect, geneList)
deg_results <- subset(deg, !change %in% c("stable"))
interGenes <- intersect(interGenes, deg_results$probe_id)
write.table(interGenes,file="interGenes.txt",sep="\t",quote=F,row.names = F)

#数据合并
allTab=data.frame()
batchType=c()
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
  rt=data
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  
  #数据合并
  if(i==1){
    allTab=rt[interGenes,]
  }else{
    allTab=cbind(allTab, rt[interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#提取交集基因的表达量
svaTab=allTab
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
geneTab=svaTab[intersect(row.names(svaTab), as.vector(geneRT[,1])),]
geneTab=t(geneTab)

#提取训练组和测试组的数据
train=grepl("^merge", rownames(geneTab), ignore.case=T)
trainExp=geneTab[train,,drop=F]
testExp=geneTab[!train,,drop=F]
rownames(trainExp)=gsub("merge_", "Train.", rownames(trainExp))
trainType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(trainExp))
testType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(testExp))
trainType=ifelse(trainType=="Control", 0, 1)
testType=ifelse(testType=="Control", 0, 1)
trainExp=cbind(trainExp, Type=trainType)
testExp=cbind(testExp, Type=testType)

#输出训练组和测试组的数据
trainOut=rbind(id=colnames(trainExp), trainExp)
write.table(trainOut, file="data.train.txt", sep="\t", quote=F, col.names=F)
testOut=rbind(id=colnames(testExp), testExp)
write.table(testOut, file="data.test.txt", sep="\t", quote=F, col.names=F)





