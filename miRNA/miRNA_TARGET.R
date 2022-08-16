##options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")) 
#options(download.file.method = 'libcurl')
#options(url.method='libcurl')
#BiocManager::install("miRNAtap",ask = F,update = F)
#BiocManager::install("topGO",ask = F,update = F)
#BiocManager::install("miRNAtap.db",ask = F,update = F)
rm(list = ls())

library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)
mir = 'miR-125a-5p'#这里只需要输入miRNA名字就可以了。其他的无脑跑代码，最终将得到一个miRNA名字的CSV文件
predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)

head(predictions)
predictions<-as.data.frame(predictions)
predictions$"ENTREZID"<-rownames(predictions)

gene<-rownames(predictions)
gene.df2 <- bitr(gene, fromType = "ENTREZID",
                toType = c("SYMBOL"),
                OrgDb = org.Hs.eg.db)

c=dplyr::inner_join(predictions,gene.df2,by ="ENTREZID")
write.csv(c,file = paste(mir,".csv",sep=""))




fef
ffs



















options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/"))
options("repos" = c(CRAN="https://mirrors.aliyun.com/CRAN/"))
options(download.file.method = 'libcurl')
options(url.method='libcurl')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multiMiR",ask = F,update = F)
library(multiMiR)
db.ver = multimir_dbInfoVersions()
