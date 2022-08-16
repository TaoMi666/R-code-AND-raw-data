library(tinyarray)
library(tidyverse)

setwd("D:\\shenxin3\\GSE87212_mRNA\\immune")
source("cibersort.R")
tpms<-read.csv("rt_TPM.csv")
tpms[1:3,1:3]
exp2 = as.data.frame(tpms)
exp2 = rownames_to_column(exp2)
write.table(exp2,file = "exp.txt",row.names = F,quote = F,sep = "\t")

TME.results = CIBERSORT("LM22.txt", 
                          "exp.txt" , 
                          perm = 1000, 
                          QN = T)
save(TME.results,file = "cibersort.Rdata")

re <- TME.results[,-(23:25)]
library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
Group = c(rep("normal",6),rep("tumor",124))
an = data.frame(group = Group,
                row.names = colnames(exp))
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))




library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

p<-ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

ggsave(p,file = "hotnew.png", width = 10, height =5.5, type = "cairo", dpi = 800) 

#展示免疫细胞之间的比较。
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

p<-ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
ggsave(p,file = "箱式图.png", width = 10.1, height =5.5, type = "cairo", dpi = 800) 

#肿瘤与正常之间的比较
dat$Group = ifelse(as.numeric(str_sub(dat$Sample,18,19))<10,"tumor","normal")
library(ggpubr)
p<-ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave(p,file = "T_N.png", width = 10.1, height =5.5, type = "cairo", dpi = 600) 
