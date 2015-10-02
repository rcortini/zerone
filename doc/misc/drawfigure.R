# Output of the fit.
# NB: alpha/size = 1.983, prob = 0.421
# ZINB: alpha/size = 3.842, pi = 0.880, prob = 0.554
dat = read.table("wgEncodeSydhTfbsK562InputRawData.out")

y = table(dat$V3)
ypois = table(rpois(lambda=mean(dat$V3), nrow(dat)))
ynb = table(rnbinom(size=1.983, prob=0.421, nrow(dat)))
yzinb = table(c(rep(0, (1-0.880)*nrow(dat)),
          rnbinom(size=3.842, prob=0.554, 0.880*nrow(dat))))

pdf("figure.pdf", width=6, height=6)
plot(0:20, y[1:21] / sum(y), type="n", ylim=c(0,.25),
    xlab="Read count (300 bp windows)", ylab="Frequency")
rect(xleft=0, xright=20, ybottom=0, ytop=.25, border=NA, col="grey96")
segments(x0=rep(0,4), y0=c(.05,.1,.15,.2),
    x1=rep(20,4), lwd=2, col="white")
lines(0:20, y[1:21] / sum(y), type="h", lwd=3)
lines(0:20+.15, ypois[1:21] / sum(y), type="h", lwd=3, col="grey75")
lines(0:20+.3, ynb[1:21] / sum(y), type="h", lwd=3, col="grey50")
lines(0:20+.45, yzinb[1:21] / sum(y), type="h", lwd=3, col="grey25")
legend(legend=c("Observed", "Poisson", "NB", "ZINB"),
    x="topright", inset=0.04, lwd=3,
    col=c("black", "grey75", "grey50", "grey25"), bg="white")
dev.off()
