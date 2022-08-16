library(randomForestSRC)

surv.rf <- rfsrc(Surv(days_to_last_follow_up, vital_status) ~ ., data = rt, ntree = 100,block.size=1, importance = TRUE)

#rf.model <- rfsrc(Surv(days_to_last_follow_up, vital_status) ~ ., data = rt, ntree = 100, block.size=1, importance = TRUE)




err.rate <-surv.rf$err.rate

#err.rate <-rf.model$err.rate

rf.model[["importance"]]

plot(rf.model)

oo <- subsample(rf.model, verbose = FALSE)


ntree=100
err.rate <-surv.rf$err.rate
par(bty = "o", mgp = c(1.5,.33,0),mar = c(3,4,1,2),las = 1,tcl = -.25)
plot(1:ntree,err.rate,xlab = "Number of Trees",ylab = "",type = "l",las = 1,cex = 1.5)
mtext("ErrorRate",side = 2,line = 2.5,las = 3)


par(bty = "o", mgp = c(1.5,.33,0),mar = c(3,4,1,2),las = 1,tcl = -.25)
plot(1:ntree,err.rate,xlab = "Number of Trees",ylab = "",type = "l",las = 1,cex = 1.5)
mtext("ErrorRate",side = 2,line = 2.5,las = 3)


raw.imp <- surv.rf$importance

names(raw.imp) <- gsub("_","-",names(raw.imp))
#
normalize <- function(x){return((x-min(x))/(max(x)-min(x)))}

rel.imp <- normalize(raw.imp)# 输出重要性矩阵

imp.res <- data.frame(gene = names(raw.imp), raw.importance = raw.imp, rel.importance = rel.imp, stringsAsFactors = F)

imp.cutoff <- 0.3

rel.imp.sel <-rel.imp[rel.imp > imp.cutoff]

rel.imp.sel <-sort(rel.imp.sel)

rel.imp.sel

xrange <- range(pretty(range(rel.imp.sel)))
yrange <- c(1,length(rel.imp.sel))
par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,6,1,2), las = 1, tcl = -.25)
plot(NULL,NULL, xlim = xrange, ylim = yrange, xlab = "Variable Relative Importance", ylab = "", yaxt = "n", las = 1)
axis(side = 2,at = 1:length(rel.imp.sel),names(rel.imp.sel))
for (i in 1:length(rel.imp.sel)) { lines(c(xrange[1],rel.imp.sel[i]),c(i,i), lwd = 2.5, col = "steelblue")}

imgene<-as.data.frame(rel.imp.sel)
imgene$name<-rownames(imgene)
imgene<-imgene[order(imgene$rel.imp.sel,decreasing = TRUE),]
write.csv(imgene,file = "基因重要性排序从高到底.csv")

afterRSF<-cbind(rt[,c(1,2)],rt[,imgene$name])

write.csv(afterRSF,file = "afterRSF_已经高到底排序.csv")

save(afterRSF,file = "afterRSF_已经高到底排序.Rdata")
