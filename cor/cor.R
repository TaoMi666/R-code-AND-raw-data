for(gene in rt$M2��������.Ԥ����򽻼�){#�ر�Ҫע�������1��lengthha����ֱ��gene in colnames(rt)
  x<-trt[,gene]
  y<-TME$`Macrophages M2`
  r<-round(cor(x,y,method="pearson"),digits = 2)
  tiff(file=paste(gene,".cor.tiff",sep=""),
       width = 14,            #ͼƬ�Ŀ���
       height =14,            #ͼƬ�ĸ߶�
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
  plot(x,y, main=paste(gene,"(r=",r,")",sep=""),xlab = gene,ylab = "Macrophages M2",cex=1,col = "black",pch=20)
  abline(lm(y~x),col = "red")
  dev.off()
  print(gene)
}


