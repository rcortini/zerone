require(GenomicRanges)

# Read profile in 300 bp windows.
txn = read.table("wg2.txt")

# Compute cumulative sums.
xtxn = cumsum(sort(txn$V4, decreasing=TRUE))
xtxn = xtxn[seq(from=1,to=length(xtxn),by=1000)]
ntxn = xtxn[length(xtxn)]

# Create GRanges objects for overlap.
gtxn = GRanges(Rle(txn$V1), IRanges(start=txn$V2, end=txn$V3), count=txn$V4)

pdf("pareto_front_txn.pdf")

# Read discretization data for all the tools.
BayesPeak1 <- read.table('H3K36me31_output.bed', stringsAsFactors=FALSE)
BayesPeak2 <- read.table('H3K36me32_output.bed', stringsAsFactors=FALSE)
JAMM <- read.table('filtered.peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
JAMM <- JAMM[nchar(JAMM$V1) <= 5, ]
MACS1 <- read.table('H3K36me3_MACS_Rep1_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
MACS2 <- read.table('H3K36me3_MACS_Rep2_peaks.narrowPeak', stringsAsFactors=FALSE)[1:3]
MACS1 <- MACS1[nchar(MACS1$V1) <= 5, ]
MACS2 <- MACS2[nchar(MACS2$V1) <= 5, ]
Zerone <- subset(read.table('H3K36.txt.gz', stringsAsFactors=FALSE), V4 == 2)

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
   ranges = IRanges(start = (Zerone$V2+Zerone$V3)/2, width=1))

ovBP1 = countOverlaps(gtxn, gBayesPeak1) > 0
ovBP2 = countOverlaps(gtxn, gBayesPeak2) > 0
ovJ = countOverlaps(gtxn, gJAMM) > 0
ovM1 = countOverlaps(gtxn, gMACS1) > 0
ovM2 = countOverlaps(gtxn, gMACS2) > 0
ovZ = countOverlaps(gtxn, gZerone) > 0

n = nrow(txn)
xbar = mean(txn$V4)

SS = nrow(txn) * var(txn$V4)
print (sum(sum(ovBP1)*(mean(txn$V4[ovBP1]) - xbar)^2 +
    sum(!ovBP1)*(mean(txn$V4[!ovBP1]) - xbar)^2) / SS)
print (sum(sum(ovBP2)*(mean(txn$V4[ovBP2]) - xbar)^2 +
    sum(!ovBP2)*(mean(txn$V4[!ovBP2]) - xbar)^2) / SS)
print (sum(sum(ovJ)*(mean(txn$V4[ovJ]) - xbar)^2 +
    sum(!ovJ)*(mean(txn$V4[!ovJ]) - xbar)^2) / SS)
print (sum(sum(ovM1)*(mean(txn$V4[ovM1]) - xbar)^2 +
    sum(!ovM1)*(mean(txn$V4[!ovM1]) - xbar)^2) / SS)
print (sum(sum(ovM2)*(mean(txn$V4[ovM2]) - xbar)^2 +
    sum(!ovM2)*(mean(txn$V4[!ovM2]) - xbar)^2) / SS)
print (sum(sum(ovZ)*(mean(txn$V4[ovZ]) - xbar)^2 +
    sum(!ovZ)*(mean(txn$V4[!ovZ]) - xbar)^2) / SS)

plot(xtxn / ntxn, type='n', ylab="", xlab="Ordered genomic windows",
     main="H3K36me3 ChIP-seq", xaxt="n", cex.lab=2, cex.main=2, cex.axis=1.2)
     title(ylab="Cumulative fraction of reads (RNA)", line=2.2, cex.lab=2)
axis(side=1, at=c(0,length(xtxn)), labels=c(0,1), cex=2)
polygon(x=c(1:length(xtxn),length(xtxn)),
        y=c(xtxn / ntxn,0), col="grey85", border=NA)
lines(xtxn / ntxn)
points(sum(ovBP1)/1000, sum(gtxn$count[ovBP1]) / ntxn, pch=19, cex=2)
points(sum(ovBP2)/1000, sum(gtxn$count[ovBP2]) / ntxn, pch=19, cex=2)
text(sum(ovBP2)/1000, sum(gtxn$count[ovBP2]) / ntxn,
     labels="BayesPeak", pos=4, offset=1, cex=2)
points(sum(ovJ)/1000, sum(gtxn$count[ovJ]) / ntxn, pch=19, cex=2)
text(sum(ovJ)/1000, sum(gtxn$count[ovJ]) / ntxn,
     labels="JAMM", pos=4, offset=1, cex=2)
points(sum(ovM1)/1000, sum(gtxn$count[ovM1]) / ntxn, pch=19, cex=2)
points(sum(ovM2)/1000, sum(gtxn$count[ovM2]) / ntxn, pch=19, cex=2)
text(sum(ovM2)/1000, sum(gtxn$count[ovM2]) / ntxn,
     labels="MACS", pos=4, offset=1, cex=2)
points(sum(ovZ)/1000, sum(gtxn$count[ovZ]) / ntxn, pch=19, cex=1.6)
points(sum(ovZ)/1000, sum(gtxn$count[ovZ]) / ntxn, cex=3.5)
text(sum(ovZ)/1000, sum(gtxn$count[ovZ]) / ntxn,
     labels="Zerone", pos=4, offset=1, cex=2)

dev.off()
