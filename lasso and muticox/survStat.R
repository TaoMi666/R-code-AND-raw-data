setwd("C:\\Users\\lexb4\\Desktop\\tcgaLncRNA\\19.survStat")
rt$vital_status<-as.character(rt$vital_status)
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$vital_status)
color[color==1]="red"
color[color==0]="green"
tiff(file="survStat.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=300)
plot(rt$days_to_last_follow_up,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topright", 
       c("Dead", "Alive"),
       pch=19,
       col=c("red","green"))
abline(v=lowLength,lty=2)
dev.off()

rt$vital_status<-factor(rt$vital_status,
                       levels = c('1','2'),
                       labels = c("Alive","Dead"))
library(ggplot2)
p<-ggplot(rt, aes(x = `riskScore`, y = `days_to_last_follow_up`))+geom_point(aes(color=`vital_status`),shape=16,size=4)+theme_bw() + theme(panel.grid=element_blank())
ggsave(p,file = "Éú´æ×ªÌ¬.tiff", width = 6, height =4, type = "cairo", dpi = 600) 

