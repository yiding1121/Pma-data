library(gplots)
library(RColorBrewer)
y = read.csv("~/Desktop/Pma_data/DEseq/Pma_transporters_fpkm.csv", sep=",")
rnames<-y[,1]
z<-data.matrix(y[,2:ncol(y)])
distance = dist(z, method = "maximum")
cluster = hclust(distance, method = "ward.D2")
#sidecolours <- colorRampPalette(c("white", "black"))(2)[y$Secreted]
mycolor=colorRampPalette(rev(brewer.pal(9,'RdYlBu')))
colors=c(seq(-5, 0, length=88),seq(0.1,5, length=88),seq(5.1,10, length=91))
pdf(file='Pma_transporters.pdf', height=15, width=7)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(10,28),trace=c('none'),col=mycolor,cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression",breaks=colors)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(10,28),trace=c('none'),cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression")

dev.off()

y = read.csv("~/Desktop/Pma_data/Analysis_by_HXG/DEseq/Pma_SP_fpkm.csv", sep=",")
rnames<-y[,1]
z<-data.matrix(y[,2:ncol(y)])
distance = dist(z, method = "maximum")
cluster = hclust(distance, method = "ward.D2")
#sidecolours <- colorRampPalette(c("white", "black"))(2)[y$Secreted]
mycolor=colorRampPalette(rev(brewer.pal(9,'RdYlBu')))
colors=c(seq(-12, -4, length=50),seq(-3.9,4, length=200),seq(4.1,12, length=150))
pdf(file='Pma_transporters.pdf', height=15, width=7)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(10,28),trace=c('none'),col=mycolor,cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression",breaks=colors)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(10,28),trace=c('none'),cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression")

y = read.csv("~/Desktop/Pma_data/Analysis_by_HXG/DEseq/Pma_CAN_fpkms.csv", sep=",")
rnames<-y[,1]
z<-data.matrix(y[,2:ncol(y)])
distance = dist(z, method = "maximum")
cluster = hclust(distance, method = "ward.D2")
#sidecolours <- colorRampPalette(c("white", "black"))(2)[y$Secreted]
mycolor=colorRampPalette(rev(brewer.pal(9,'RdYlBu')))
colors=c(seq(-10, 0, length=20),seq(0.1,4, length=40),seq(4.1,10, length=60))
pdf(file='Pma_transporters.pdf', height=15, width=7)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(15,35),trace=c('none'),col=mycolor,cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression",breaks=colors)
heatmap.2(z,scale='none',dendrogram="row",Colv=NA,density.info = c('none'),margins=c(10,28),trace=c('none'),cexRow=0.001,cexCol=1,symbreaks=TRUE,symkey=FALSE,keysize=1, key.title=NA, key.par=list(mar=c(3.5,0.5,3,0)), key.xlab = "log2 expression")
