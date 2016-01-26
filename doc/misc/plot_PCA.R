trainset = read.delim("trainset.txt")
trainset = subset(trainset, complete.cases(trainset))
pca = prcomp(trainset[,-1], scale=TRUE)

simCTCF = as.matrix(read.delim("simCTCF.txt")[,-9])
simH3K36 = as.matrix(read.delim("simH3K36me3.txt")[,-9])


pdf("PCA.pdf", height=6, width=11)
par(mfrow=c(1,2))
plot(pca$x, col=2-trainset$OK, pch=19, cex=.7)
points((simCTCF %*% pca$rotation), col=3, pch=19, cex=.7)
points((simH3K36 %*% pca$rotation), col=4, pch=19, cex=.7)
plot(pca$x[,3:2], col=2-trainset$OK, pch=19, cex=.7)
points((simCTCF %*% pca$rotation)[,3:2], col=3, pch=19, cex=.7)
points((simH3K36 %*% pca$rotation)[,3:2], col=4, pch=19, cex=.7)
dev.off()
