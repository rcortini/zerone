trainset = read.delim("trainset.txt")
trainset = subset(trainset, complete.cases(trainset))
pca = prcomp(trainset[,-1], scale=TRUE)
# Change orientation of the first PC.
pca$x[,1] = -pca$x[,1]
pca$rotation[,1] = -pca$rotation[,1]

#simCTCF = as.matrix(read.delim("simCTCF.txt")[,-9])
#simH3K36 = as.matrix(read.delim("simH3K36me3.txt")[,-9])

pdf("PCA.pdf", height=4.5, width=8)
par(mfrow=c(1,2))
plot(pca$x, col=2-trainset$OK, pch=19, cex=.7,
     yaxp=c(-6,6,4), xaxp=c(-2,6,4), xlab="", ylab="")
title(xlab="First Principal Component",
      ylab="Second Principal Component", line=2.2)
legend(x="topleft", inset=0, pch=19, pt.cex=.7,
      col=c(1,2), legend=c("Positive", "Negative"), bty="n")
#points((simCTCF %*% pca$rotation), col=3, pch=19, cex=.7)
#points((simH3K36 %*% pca$rotation), col=4, pch=19, cex=.7)
plot(pca$x[,3:2], col=2-trainset$OK, pch=19, cex=.7,
     yaxp=c(-6,6,4), xaxp=c(-2,6,4), xlab="", ylab="")
title(xlab="Third Principal Component",
      ylab="Second Principal Component", line=2.2)
#points((simCTCF %*% pca$rotation)[,3:2], col=3, pch=19, cex=.7)
#points((simH3K36 %*% pca$rotation)[,3:2], col=4, pch=19, cex=.7)
dev.off()
