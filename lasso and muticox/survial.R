
#install.packages("survival")

setwd("C:\\Users\\lexb4\\Desktop\\tcgaLncRNA\\14.survival")   #工作目录（需修改）
library(survival)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(days_to_last_follow_up, vital_status) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
#pValue=round(pValue,3)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ risk, data = rt)
summary(fit)    #查看五年生存率
tiff(file="survival.tiff",
       width = 14,            #图片的宽度
       height =14,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()


