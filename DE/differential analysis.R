
library(edgeR)
group<-c(rep(1,6),rep(2,63))  #1表示对照，2表示实验组                      
y<-DGEList(counts = df,group=group)
keep<-filterByExpr(y)
y<-y[keep, ,keep.lib.sizes=FALSE]
y<-calcNormFactors(y)
y<-estimateDisp(y)
et<-exactTest(y)
et<-topTags(et,n=20000)
et<-as.data.frame(et)
et<-cbind(rownames(et),et)
colnames(et)<-c("gene_id","log2foldchange","log2CPM","pvalue","FDR")

et[which(et$pvalue<0.05 & et$log2foldchange>2),"change"]<-"up"
et[which(et$pvalue<0.05 & et$log2foldchange< -2),"change"]<-"down"
et[which(et$pvalue<0.05 & abs(et$log2foldchange)< 2),"change"]<-"not"
et[which(et$pvalue>0.05 ),"change"]<-"not"
et$logP<--log10(et$pvalue)
etsig<-et[which(et$pvalue<0.05 & abs(et$log2foldchange)>2),]

write.csv(etsig,file = "FD=2mRNA.csv")
write.csv(et,file = "et.csv")

