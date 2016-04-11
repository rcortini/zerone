J1 = subset(read.table("jarid.out.gz", comm="#"), V1 == "chr1")
# Set many values to NA in order to reduce the size of the figure.

n = nrow(J1)
J1$V6[(1:n %% 10 != 0) & (J1$V6 < 11)] = NA
J1$V7[(1:n %% 10 != 0) & (J1$V7 < 11)] = NA

pdf("jarid.pdf")
par(mfrow=c(2,1))
par(mar=c(2.9,3.5,1,0))

plot(J1$V6[1:200000], type='h', col="#00000080",
     ylim=c(0,120), xaxt="n", yaxt="n", bty="n", ylab="")
axis(side=2, at=c(0,40,80,120), cex.axis=1.2)
axis(side=1, at=(0:6)*1e5/3, labels=(0:6)*10, cex.axis=1.2, padj=-.5)
title(ylab="Read count", line=2.4, cex.lab=1.3)
title(xlab="Position on chr1 (Mb)", line=2, cex.lab=1.3)
text(x=-7000, y=115, "JARID1A (rep1)", pos=4, cex=1.3)
par(mar=c(3,3.5,0,0))
plot(J1$V7[1:200000], type='h', col="#00000080",
     ylim=c(0,120), bty="n", yaxt="n", xaxt="n", cex.axis=0.8,
     ylab="", xlab="")
title(ylab="Read count", line=2.4, cex.lab=1.3)
title(xlab="Position on chr1 (Mb)", line=2, cex.lab=1.3)
axis(side=1, at=(0:6)*1e5/3, labels=(0:6)*10, cex.axis=1.2, padj=-.5)
axis(side=2, at=c(0,40,80,120), cex.axis=1.2)
text(x=-7000, y=115, "JARID1A (rep2)", pos=4, cex=1.3)

dev.off()
