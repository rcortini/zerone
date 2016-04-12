data.dir = "/home/pcusco/src/C/zerone/docker/"
output.dir = "/home/pcusco/Zerone/figures/"
setwd(output.dir)

require(GenomicRanges)

# Read profile in 300 bp windows.
CTCF = read.table("zerone_windows/CTCF.txt.gz")
H3K36me3 = read.table("zerone_windows/H3K36me3.txt.gz")
H3K4me3 = read.table("zerone_windows/H3K4me3.txt.gz")
Pol2 = read.table("zerone_windows/Pol2.txt.gz")

# Compute cumulative sums.
xCTCF = cumsum(sort(CTCF$V4, decreasing=TRUE))
xCTCF = xCTCF[seq(from=1,to=length(xCTCF),by=1000)]
nCTCF = xCTCF[length(xCTCF)]

xH3K36me3 = cumsum(sort(H3K36me3$V4, decreasing=TRUE))
xH3K36me3 = xH3K36me3[seq(from=1,to=length(xH3K36me3),by=1000)]
nH3K36me3 = xH3K36me3[length(xH3K36me3)]

xH3K4me3 = cumsum(sort(H3K4me3$V4, decreasing=TRUE))
xH3K4me3 = xH3K4me3[seq(from=1,to=length(xH3K4me3),by=1000)]
nH3K4me3 = xH3K4me3[length(xH3K4me3)]

xPol2 = cumsum(sort(Pol2$V4, decreasing=TRUE))
xPol2 = xPol2[seq(from=1,to=length(xPol2),by=1000)]
nPol2 = xPol2[length(xPol2)]

# Create GRanges objects for overlap.
gCTCF = GRanges(Rle(CTCF$V1), IRanges(start=CTCF$V2, end=CTCF$V3), count=CTCF$V4)
gH3K36me3 = GRanges(Rle(H3K36me3$V1), IRanges(start=H3K36me3$V2, end=H3K36me3$V3), count=H3K36me3$V4)
gH3K4me3 = GRanges(Rle(H3K4me3$V1), IRanges(start=H3K4me3$V2, end=H3K4me3$V3), count=H3K4me3$V4)
gPol2 = GRanges(Rle(Pol2$V1), IRanges(start=Pol2$V2, end=Pol2$V3), count=Pol2$V4)

################################################################################################

COL=c("#798698", "#bb9966", "#843e24", "black")

pdf("~/curves.pdf", height=4.5, width=17, useDingbats = FALSE)
layout(matrix(1:4, ncol=4))
par(mar = c(3.5, 3, 2, 1), cex = 1)


# CTCF figure. #####################################
# Read discretization data for all the tools.
setwd(data.dir)
BayesPeak1 <- read.table('BayesPeak/Ctcf1.bed', stringsAsFactors=FALSE)
BayesPeak2 <- read.table('BayesPeak/Ctcf2.bed', stringsAsFactors=FALSE)
JAMM <- read.table('JAMM/Ctcf/peaks/filtered.peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
JAMM <- JAMM[nchar(JAMM$V1) <= 5, ]
#MACS1 <- read.table('MACS/Ctcf1_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
#MACS2 <- read.table('MACS/Ctcf2_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
MACS1 <- read.table('MACS/Ctcf1.bed', stringsAsFactors=FALSE)[1:3]
MACS2 <- read.table('MACS/Ctcf2.bed', stringsAsFactors=FALSE)[1:3]
MACS1 <- MACS1[nchar(MACS1$V1) <= 5, ]
MACS2 <- MACS2[nchar(MACS2$V1) <= 5, ]
Zerone <- read.table('zerone/Ctcf.bed', stringsAsFactors=FALSE)

# Create GRange objects. Note that only the midpoint of the peak is
# considered so that each peak maps to a single 300 bp window.
gBayesPeak1 <- GRanges(seqnames = BayesPeak1$V1,
                       ranges = IRanges(start = (BayesPeak1$V2+BayesPeak1$V3)/2, width=1))
gBayesPeak2 <- GRanges(seqnames = BayesPeak2$V1,
                       ranges = IRanges(start = (BayesPeak2$V2+BayesPeak2$V3)/2, width=1))
gJAMM <- GRanges(seqnames = JAMM$V1,
                 ranges = IRanges(start = (JAMM$V2+JAMM$V3)/2, width=1))
gMACS1 <- GRanges(seqnames = MACS1$V1,
                  ranges = IRanges(start = (MACS1$V2+MACS1$V3)/2, width=1))
gMACS2 <- GRanges(seqnames = MACS2$V1,
                  ranges = IRanges(start = (MACS2$V2+MACS2$V3)/2, width=1))
gZerone <- GRanges(seqnames = Zerone$V1,
                   ranges = IRanges(start = Zerone$V2, end = Zerone$V3))
ovBP1 = countOverlaps(gCTCF, gBayesPeak1) > 0
ovBP2 = countOverlaps(gCTCF, gBayesPeak2) > 0
ovJ = countOverlaps(gCTCF, gJAMM) > 0
ovM1 = countOverlaps(gCTCF, gMACS1) > 0
ovM2 = countOverlaps(gCTCF, gMACS2) > 0
ovZ = countOverlaps(gCTCF, gZerone) > 0

setwd(output.dir)
plot(xCTCF / nCTCF, type='n', ylab="", xlab="", xaxt="n", cex.axis=1.2)
axis(side=1, at=c(0,length(xCTCF)), labels=c("First","Last"), cex.axis=1.2)
title(xlab="Ordered genomic windows", line=2.2, cex.lab=1.2)
title(main="CTCF ChIP-seq", line=.7)
polygon(x=c(1:length(xCTCF),length(xCTCF)),
        y=c(xCTCF / nCTCF,0), col="grey95", border=NA)
lines(xCTCF / nCTCF)
legend(x="topright", inset=0.02, col=COL, pch=19, pt.cex=1.3,
   legend=c("BayesPeak", "JAMM", "MACS", "Zerone"), box.col="grey50", bg="white")
points(sum(ovBP1)/1000, sum(gCTCF$count[ovBP1]) / nCTCF, pch=19, col=COL[1], cex=1.4)
points(sum(ovJ)/1000, sum(gCTCF$count[ovJ]) / nCTCF, pch=19, col=COL[2], cex=1.4)
points(sum(ovM1)/1000, sum(gCTCF$count[ovM1]) / nCTCF, pch=19, col=COL[3], cex=1.4)
points(sum(ovBP2)/1000, sum(gCTCF$count[ovBP2]) / nCTCF, pch=19, col=COL[1], cex=1.4)
#points(sum(ovM2)/1000, sum(gCTCF$count[ovM2]) / nCTCF, pch=19, col=COL[3], cex=1.4)
points(sum(ovZ)/1000, sum(gCTCF$count[ovZ]) / nCTCF, cex=1.3, pch=19, col=COL[4])
#points(sum(ovZ)/1000, sum(gCTCF$count[ovZ]) / nCTCF, cex=1.8, col=COL[4])

# Pol2 figure. #############################
# Read discretization data for all the tools.
setwd(data.dir)
BayesPeak1 <- read.table('BayesPeak/Pol21.bed', stringsAsFactors=FALSE)
BayesPeak2 <- read.table('BayesPeak/Pol22.bed', stringsAsFactors=FALSE)
#JAMM <- read.table('JAMM/Pol2/peaks/filtered.peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
#JAMM <- JAMM[nchar(JAMM$V1) <= 5, ]
JAMM <- read.table('benchmark/idrCode/JAMM/Pol2-overlapped-peaks.txt')
JAMM <- subset(JAMM, IDR < 0.05)[, 1:3]
colnames(JAMM) <- c('V1', 'V2', 'V3')
#MACS1 <- read.table('MACS/Pol21_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
#MACS2 <- read.table('MACS/Pol22_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
MACS1 <- read.table('MACS/Pol21.bed', stringsAsFactors=FALSE)[1:3]
MACS2 <- read.table('MACS/Pol22.bed', stringsAsFactors=FALSE)[1:3]
MACS1 <- MACS1[nchar(MACS1$V1) <= 5, ]
MACS2 <- MACS2[nchar(MACS2$V1) <= 5, ]
Zerone <- read.table('zerone/Pol2.bed', stringsAsFactors=FALSE)

# Create GRange objects. Note that only the midpoint of the peak is
# considered so that each peak maps to a single 300 bp window.
gBayesPeak1 <- GRanges(seqnames = BayesPeak1$V1,
                       ranges = IRanges(start = (BayesPeak1$V2+BayesPeak1$V3)/2, width=1))
gBayesPeak2 <- GRanges(seqnames = BayesPeak2$V1,
                       ranges = IRanges(start = (BayesPeak2$V2+BayesPeak2$V3)/2, width=1))
gJAMM <- GRanges(seqnames = JAMM$V1,
                 ranges = IRanges(start = (JAMM$V2+JAMM$V3)/2, width=1))
gMACS1 <- GRanges(seqnames = MACS1$V1,
                  ranges = IRanges(start = (MACS1$V2+MACS1$V3)/2, width=1))
gMACS2 <- GRanges(seqnames = MACS2$V1,
                  ranges = IRanges(start = (MACS2$V2+MACS2$V3)/2, width=1))
gZerone <- GRanges(seqnames = Zerone$V1,
                   ranges = IRanges(start = Zerone$V2, end = Zerone$V3))
ovBP1 = countOverlaps(gPol2, gBayesPeak1) > 0
ovBP2 = countOverlaps(gPol2, gBayesPeak2) > 0
ovJ = countOverlaps(gPol2, gJAMM) > 0
ovM1 = countOverlaps(gPol2, gMACS1) > 0
ovM2 = countOverlaps(gPol2, gMACS2) > 0
ovZ = countOverlaps(gPol2, gZerone) > 0

setwd(output.dir)
plot(xPol2 / nPol2, type='n', ylab="", xlab="", cex.axis=1.2, xaxt="n")
axis(side=1, at=c(0,length(xCTCF)), labels=c("First","Last"), cex.axis=1.2)
title(xlab="Ordered genomic windows", line=2.2, cex.lab=1.2)
title(main="Pol2 ChIP-seq", line=.7)
polygon(x=c(1:length(xPol2),length(xPol2)),
        y=c(xPol2 / nPol2,0), col="grey95", border=NA)
lines(xPol2 / nPol2)
points(sum(ovBP1)/1000, sum(gPol2$count[ovBP1]) / nPol2, pch=19, col=COL[1], cex=1.4)
points(sum(ovBP2)/1000, sum(gPol2$count[ovBP2]) / nPol2, pch=19, col=COL[1], cex=1.4)
points(sum(ovJ)/1000, sum(gPol2$count[ovJ]) / nPol2, pch=19, col=COL[2], cex=1.4)
points(sum(ovM1)/1000, sum(gPol2$count[ovM1]) / nPol2, pch=19, col=COL[3], cex=1.4)
points(sum(ovM2)/1000, sum(gPol2$count[ovM2]) / nPol2, pch=19, col=COL[3], cex=1.4)
points(sum(ovZ)/1000, sum(gPol2$count[ovZ]) / nPol2, cex=1.3, pch=19, col=COL[4])

# H3K4me3 ######################################
setwd(data.dir)
BayesPeak1 <- read.table('BayesPeak/H3k4me31.bed', stringsAsFactors=FALSE)
BayesPeak2 <- read.table('BayesPeak/H3k4me32.bed', stringsAsFactors=FALSE)
#JAMM <- read.table('JAMM/H3k4me3/peaks/filtered.peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
#JAMM <- JAMM[nchar(JAMM$V1) <= 5, ]
JAMM <- read.table('benchmark/idrCode/JAMM/H3k4me3-overlapped-peaks.txt')
JAMM <- subset(JAMM, IDR < 0.05)[, 1:3]
colnames(JAMM) <- c('V1', 'V2', 'V3')
#MACS1 <- read.table('MACS/H3k4me31_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
#MACS2 <- read.table('MACS/H3k4me32_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
MACS1 <- read.table('MACS/H3k4me31.bed', stringsAsFactors=FALSE)[1:3]
MACS2 <- read.table('MACS/H3k4me32.bed', stringsAsFactors=FALSE)[1:3]
MACS1 <- MACS1[nchar(MACS1$V1) <= 5, ]
MACS2 <- MACS2[nchar(MACS2$V1) <= 5, ]
Zerone <- read.table('zerone/H3k4me3.bed', stringsAsFactors=FALSE)

# Create GRange objects. Note that only the midpoint of the peak is
# considered so that each peak maps to a single 300 bp window.
gBayesPeak1 <- GRanges(seqnames = BayesPeak1$V1,
                       ranges = IRanges(start = (BayesPeak1$V2+BayesPeak1$V3)/2, width=1))
gBayesPeak2 <- GRanges(seqnames = BayesPeak2$V1,
                       ranges = IRanges(start = (BayesPeak2$V2+BayesPeak2$V3)/2, width=1))
gJAMM <- GRanges(seqnames = JAMM$V1,
                 ranges = IRanges(start = (JAMM$V2+JAMM$V3)/2, width=1))
gMACS1 <- GRanges(seqnames = MACS1$V1,
                  ranges = IRanges(start = (MACS1$V2+MACS1$V3)/2, width=1))
gMACS2 <- GRanges(seqnames = MACS2$V1,
                  ranges = IRanges(start = (MACS2$V2+MACS2$V3)/2, width=1))
gZerone <- GRanges(seqnames = Zerone$V1,
                   ranges = IRanges(start = Zerone$V2, end = Zerone$V3))
ovBP1 = countOverlaps(gH3K4me3, gBayesPeak1) > 0
ovBP2 = countOverlaps(gH3K4me3, gBayesPeak2) > 0
ovJ = countOverlaps(gH3K4me3, gJAMM) > 0
ovM1 = countOverlaps(gH3K4me3, gMACS1) > 0
ovM2 = countOverlaps(gH3K4me3, gMACS2) > 0
ovZ = countOverlaps(gH3K4me3, gZerone) > 0

setwd(output.dir)
plot(xH3K4me3 / nH3K4me3, type='n', ylab="", xlab="", xaxt="n", cex.axis=1.2)
axis(side=1, at=c(0,length(xH3K4me3)), labels=c("First","Last"), cex.axis=1.2)
title(xlab="Ordered genomic windows", line=2.2, cex.lab=1.2)
title(main="H3K4me3 ChIP-seq", line=.7)
polygon(x=c(1:length(xH3K4me3),length(xH3K4me3)),
        y=c(xH3K4me3 / nH3K4me3,0), col="grey95", border=NA)
lines(xH3K4me3 / nH3K4me3)
points(sum(ovBP1)/1000, sum(gH3K4me3$count[ovBP1]) / nH3K4me3, pch=19, col=COL[1], cex=1.4)
points(sum(ovBP2)/1000, sum(gH3K4me3$count[ovBP2]) / nH3K4me3, pch=19, col=COL[1], cex=1.4)
points(sum(ovJ)/1000, sum(gH3K4me3$count[ovJ]) / nH3K4me3, pch=19, col=COL[2], cex=1.4)
points(sum(ovM1)/1000, sum(gH3K4me3$count[ovM1]) / nH3K4me3, pch=19, col=COL[3], cex=1.4)
points(sum(ovM2)/1000, sum(gH3K4me3$count[ovM2]) / nH3K4me3, pch=19, col=COL[3], cex=1.4)
points(sum(ovZ)/1000, sum(gH3K4me3$count[ovZ]) / nH3K4me3, cex=1.3, pch=19, col=COL[4])

# H3K36me3 figure. ############################################
# Read discretization data for all the tools.
setwd(data.dir)
BayesPeak1 <- read.table('BayesPeak/H3k36me31.bed', stringsAsFactors=FALSE)
BayesPeak2 <- read.table('BayesPeak/H3k36me32.bed', stringsAsFactors=FALSE)
JAMM <- read.table('JAMM/H3k36me3/peaks/filtered.peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
JAMM <- JAMM[nchar(JAMM$V1) <= 5, ]
#MACS1 <- read.table('MACS/H3k36me31_peaks.broadPeak', stringsAsFactors=FALSE)[1:3]
#MACS2 <- read.table('MACS/H3k36me32_peaks.broadPeak', stringsAsFactors=FALSE)[1:3]
MACS1 <- read.table('MACS/H3k36me31.bed', stringsAsFactors=FALSE)[1:3]
MACS2 <- read.table('MACS/H3k36me32.bed', stringsAsFactors=FALSE)[1:3]
MACS1 <- MACS1[nchar(MACS1$V1) <= 5, ]
MACS2 <- MACS2[nchar(MACS2$V1) <= 5, ]
Zerone <- read.table('zerone/H3k36me3.bed', stringsAsFactors=FALSE)

# Create GRange objects. Note that only the midpoint of the peak is
# considered so that each peak maps to a single 300 bp window.
gBayesPeak1 <- GRanges(seqnames = BayesPeak1$V1,
   ranges = IRanges(start = (BayesPeak1$V2+BayesPeak1$V3)/2, width=1))
gBayesPeak2 <- GRanges(seqnames = BayesPeak2$V1,
   ranges = IRanges(start = (BayesPeak2$V2+BayesPeak2$V3)/2, width=1))
gJAMM <- GRanges(seqnames = JAMM$V1,
   ranges = IRanges(start = (JAMM$V2+JAMM$V3)/2, width=1))
gMACS1 <- GRanges(seqnames = MACS1$V1,
   ranges = IRanges(start = (MACS1$V2+MACS1$V3)/2, width=1))
gMACS2 <- GRanges(seqnames = MACS2$V1,
   ranges = IRanges(start = (MACS2$V2+MACS2$V3)/2, width=1))
gZerone <- GRanges(seqnames = Zerone$V1,
   ranges = IRanges(start = Zerone$V2, end = Zerone$V3))
ovBP1 = countOverlaps(gH3K36me3, gBayesPeak1) > 0
ovBP2 = countOverlaps(gH3K36me3, gBayesPeak2) > 0
ovJ = countOverlaps(gH3K36me3, gJAMM) > 0
ovM1 = countOverlaps(gH3K36me3, gMACS1) > 0
ovM2 = countOverlaps(gH3K36me3, gMACS2) > 0
ovZ = countOverlaps(gH3K36me3, gZerone) > 0

setwd(output.dir)
plot(xH3K36me3 / nH3K36me3, type='n', ylab="", xlab="", xaxt="n", cex.axis=1.2)
axis(side=1, at=c(0,length(xH3K36me3)), labels=c("First","Last"), cex.axis=1.2)
title(xlab="Ordered genomic windows", line=2.2, cex.lab=1.2)
title(main="H3K36me3 ChIP-seq", line=.7)
polygon(x=c(1:length(xH3K36me3),length(xH3K36me3)),
        y=c(xH3K36me3 / nH3K36me3,0), col="grey95", border=NA)
lines(xH3K36me3 / nH3K36me3)
points(sum(ovBP1)/1000, sum(gH3K36me3$count[ovBP1]) / nH3K36me3, pch=19, col=COL[1], cex=1.4)
points(sum(ovBP2)/1000, sum(gH3K36me3$count[ovBP2]) / nH3K36me3, pch=19, col=COL[1], cex=1.4)
points(sum(ovJ)/1000, sum(gH3K36me3$count[ovJ]) / nH3K36me3, pch=19, col=COL[2], cex=1.4)
points(sum(ovM1)/1000, sum(gH3K36me3$count[ovM1]) / nH3K36me3, pch=19, col=COL[2], cex=1.4)
points(sum(ovM2)/1000, sum(gH3K36me3$count[ovM2]) / nH3K36me3, pch=19, col=COL[3], cex=1.4)
points(sum(ovZ)/1000, sum(gH3K36me3$count[ovZ]) / nH3K36me3, cex=1.3, pch=19, col=COL[4])

###################################################
dev.off()
