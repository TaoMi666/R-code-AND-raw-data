for(gene in rt$M2相关与差异.预后基因交集){#特别要注意这里，是1：lengthha还是直接gene in colnames(rt)
  x<-trt[,gene]
  y<-TME$`Macrophages M2`
  r<-round(cor(x,y,method="pearson"),digits = 2)
  tiff(file=paste(gene,".cor.tiff",sep=""),
       width = 14,            #图片的宽度
       height =14,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
  plot(x,y, main=paste(gene,"(r=",r,")",sep=""),xlab = gene,ylab = "Macrophages M2",cex=1,col = "black",pch=20)
  abline(lm(y~x),col = "red")
  dev.off()
  print(gene)
}



