
for(gene in list[1:480]){#特别要注意这里，是1：lengthha还是直接gene in colnames(rt)
  
  x<-as.numeric(dna1[,gene]) 
  
  y<-as.numeric(rna1[,gene])
  
  r<-round(cor(x,y,method="pearson"),digits = 2)
  tiff(file=paste(gene,".cor.tiff",sep=""),
       width = 14,            #图片的宽度
       height =14,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
  plot(x,y, main=paste(gene,"(r=",r,")",sep=""),xlab = "DNA_methylation",ylab = "Gene_expression",cex.axis=1.5,cex.lab=1.5,cex=1,col = "black",pch=16)
  abline(lm(y~x),col = "black",lwd=2)
  dev.off()
  print(gene)
}
































png(filename="xianguan.png",width=680*3,height=480*3,res=72*3)
dna2<-dna1[,list]
rna2<-rna1[,list]
list<-cor_data_sig$symbol
opar<-par(no.readonly = TRUE)
par(mfrow=c(3,4))
par(mai=c(0.3,0.3,0.3,0.2))
par(col=c("blue","red"))
for(i in 1:length(list)){
x<-dna2[,i]
y<-rna2[,i]
r<-round(cor(x,y,method="pearson"),digits = 2)
plot(x,y, main=paste(list[i],"(r=",r,")",sep=""),xlab = "",ylab = "",
        cex=1)
abline(lm(y~x))
}
par(opar)
dev.off()


x1<-dna2[,1]
y1<-rna2[,1]
plot(x1,y1,main = "xxxx",xlab = "",ylab = "")
abline(lm(y1~x1)) #lm
r<-cor(x,y,method="pearson")