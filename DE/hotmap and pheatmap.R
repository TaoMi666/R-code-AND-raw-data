 library("ggpubr")
library(ggthemes)
 et$logP<- -log10(et$FDR)
 et[which(et$pvalue>0),"change"]<-"not"
 et[which(et$log2foldchange>1& et$FDR<0.05),"change"]<-"up"
 et[which(et$log2foldchange< -1 & et$FDR<0.05),"change"]<-"down"

 #加标签，我们去p前10上调和下调的基因为例
et$label<-""

et<-et[order(abs(et$pvalue)),]
upgene<-head(et$gene_id[which(et$`change`=="up")],10)
downgene<-head(et$gene_id[which(et$`change`=="down")],10)
top10<-c(as.character(upgene),as.character(downgene))
et$label[match(top10,et$gene_id)]<-top10#
#正式画火山图
png("火山label1.png",width=800*2,height=800*2,res=100*3)
ggscatter(et,x="log2foldchange",y="logP",color = "up.down",palette=c("blue","gray","red"),size = 1,label =et$label,font.label = 6,repel = T)
#ggscatter(et,x="log2foldchange",y="logP",color = "up-down",palette=c("green","gray","red"),size = 1 )+theme_base() 
dev.off() 
library(ggrepel) #标签用
library(ggplot2)
ggplot(et, aes(log2foldchange, logp1))+geom_point(aes(color = change))+ labs(tittle = "volcanoplot", x = expression(log2FoldChange), y = expression(-log10(pvalue)))+scale_color_manual(values =c("navy","black", "red"))+ geom_hline(yintercept = -log10(0.05),linetype="dotted")+geom_vline(xintercept = c(-1,1),linetype="dotted")+theme(axis.title.x =element_text(size=18), axis.title.y=element_text(size=18),legend.title =element_text(size=18),legend.text =element_text(size=14),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"),,panel.background=element_rect(color="black"))

#p<-ggplot(et, aes(log2foldchange, logP))+geom_point(aes(color = change))+ labs(tittle = "volcanoplot", x = expression(log2FoldChange), y = expression(-log10(pvalue)))+scale_color_manual(values =c("navy","black", "red"))+ geom_hline(yintercept = -log10(0.05),linetype="dotted")+geom_vline(xintercept = c(-1,1),linetype="dotted")+theme(axis.title.x =element_text(size=18), axis.title.y=element_text(size=18),legend.title =element_text(size=18),legend.text =element_text(size=14),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"),panel.background=element_rect(fill="white",color="grey50")，panel.background = element_rect(fill = "lightblue", colour = "red", size=3)) 
p<-ggplot(et, aes(log2foldchange, logP))+geom_point(aes(color = change))+ labs(tittle = "volcanoplot", x = expression(log2FoldChange), y = expression(-log10(pvalue)))+
        scale_color_manual(values =c("green","black", "red"))+ 
        geom_hline(yintercept = -log10(0.05),linetype="dotted")+geom_vline(xintercept = c(-2,2),linetype="dotted")+theme(axis.title.x =element_text(size=18), axis.title.y=element_text(size=18),legend.title =element_text(size=18),legend.text =element_text(size=14),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"),panel.grid.major =element_line(size=1,linetype =3,color="grey70"),,panel.background=element_rect(color="black"))+
        geom_text_repel(
                data = et[et$pvalue<0.05&abs(et$log2foldchange)>1,],
                aes(label = label),
                size = 3,
                color = "black",
                segment.color = "black",force = 7, segment.size=0.1,show.legend = FALSE )

ggsave(p,file = "火山-甲基化需要1.png", width = 6, height =5.5, type = "cairo", dpi = 700) 
    
pheatmap(hot50,kmeans_k = NA,breaks = NA,gaps_row = NULL,border = F,
         border_color = "grey" ,labels_row = NULL,annotation_col=annotation_col,
         cluster_cols = FALSE,show_rownames = T,show_colnames = F,
         treeheight_col=13,treeheight_row =10, 
         color = colorRampPalette(c("springgreen4", "white", "darkred"))(50))





p<-pheatmap(n6,annotation_col=annotation_col,annotation_row =annotation_row ,cluster_row = FALSE,show_rownames = T,show_colnames = F,treeheight_col=10,treeheight_row =10,cutree_cols=2)
 ggsave(p,file = "hotnew2.png", width = 5, height =6, type = "cairo", dpi = 800) 
 xpheatmap(mird1)
 hot50<-hot[1:50,]
 