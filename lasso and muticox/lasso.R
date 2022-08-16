#��Ҫ�����ļ���һ���Ǳ��������һ��genelist���ƣ�����ҵı��������c,����listΪgene

c$days_to_last_follow_up<-c$days_to_last_follow_up/365
c1<-c[,-c(1,2)]
c1<-log2(c1+1)
c2<-c[,c(1,2)]
c<-cbind(c1,c2)


afterRSF<-read.csv("afterRSF_�Ѿ��ߵ�������.csv")
library(caret)
library("glmnet")
library("survival")
set.seed(1)
sam<- createDataPartition(c$vital_status, p = .7,list = FALSE)
head(sam)
train <- c[sam,]
test <- c[-sam,]

train<-afterRSF
x<-train[,3:ncol(train)]
x<-as.matrix(x)#x����Ϊ����

#�鿴����һЩ�ٴ������и����
#prop.table(table(c$status))
#prop.table(table(c$status)) 


y <- data.matrix(Surv(train$days_to_last_follow_up,train$vital_status==1))


cv.fit <- cv.glmnet(x, y, type.measure="C", family="cox",nlambda = 100)

plot(cv.fit)

png(filename="lasso1.png",width=680*3,height=480*3,res=150*3)
plot(cv.fit)
dev.off()


fit <- glmnet(x, y, type.measure="C",family = "cox",nlambda = 100)

png(filename="lasso2.png",width=680*3,height=480*3,res=150*3)
plot(fit)
dev.off()

Coefficients <- coef(fit, s = cv.fit$lambda.min)

Active.Index <- which(Coefficients != 0)

Active.Coefficients <- Coefficients[Active.Index]


Active.Index

Active.Coefficients

row.names(Coefficients)[Active.Index]

gene<-row.names(Coefficients)[Active.Index]
afterlasso<-cbind(train[,c(1,2)],train[,gene])
save(afterlasso,file = "afterlasso.Rdata")
write.csv(afterlasso,file = "afterlasso.csv")
