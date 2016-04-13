setwd('/')

total.tss <- unique(read.table('benchmark/hg19_TSS.bed'))
gro.plus  <- read.table('benchmark/GSM1480325_K562_GROseq_plus.bed')
gro.minus <- read.table('benchmark/GSM1480325_K562_GROseq_minus.bed')

# total.jamm.idr <- read.table('benchmark/idrCode/JAMM/Pol2-overlapped-peaks.txt')
# total.macs.idr <- read.table('benchmark/idrCode/MACS_q0.05/Pol2/IDR_Pol2-overlapped-peaks.txt')
# pass.idr.jamm <- subset(total.jamm.idr, IDR < 0.05)
# pass.idr.macs <- subset(total.macs.idr, IDR < 0.05)
#
# total.bayespeak1 <- read.table('BayesPeak/Pol21.bed')
# total.bayespeak2 <- read.table('BayesPeak/Pol22.bed')
# total.jamm       <- read.table('JAMM/Pol2/peaks/filtered.peaks.narrowPeak')[1:3]
# total.macs1      <- read.table('MACS/Pol21.bed')[1:3]
# total.macs2      <- read.table('MACS/Pol22.bed')[1:3]
# total.zerone     <- read.table('zerone/Pol2.bed')
# total.jamm.idr1 <- pass.idr.jamm[, 1:3]
# total.jamm.idr2 <- pass.idr.jamm[, 5:7]
# total.macs.idr1 <- pass.idr.macs[, 1:3]
# total.macs.idr2 <- pass.idr.macs[, 5:7]

total.jamm.idr <- read.table('benchmark/idrCode/JAMM/H3k4me3-overlapped-peaks.txt')
total.macs.idr <- read.table('benchmark/idrCode/MACS_q0.05/H3k4me3/IDR_H3k4me3-overlapped-peaks.txt')
pass.idr.jamm <- subset(total.jamm.idr, IDR < 0.05)
pass.idr.macs <- subset(total.macs.idr, IDR < 0.05)

total.bayespeak1 <- read.table('BayesPeak/H3k4me31.bed')
total.bayespeak2 <- read.table('BayesPeak/H3k4me32.bed')
total.jamm       <- read.table('JAMM/H3k4me3/peaks/filtered.peaks.narrowPeak')[1:3]
total.macs1      <- read.table('MACS/H3k4me31.bed')[1:3]
total.macs2      <- read.table('MACS/H3k4me32.bed')[1:3]
total.zerone     <- read.table('zerone/H3k4me3.bed')
total.jamm.idr1 <- pass.idr.jamm[, 1:3]
total.jamm.idr2 <- pass.idr.jamm[, 5:7]
total.macs.idr1 <- pass.idr.macs[, 1:3]
total.macs.idr2 <- pass.idr.macs[, 5:7]

names(total.tss)        <- c('chrom', 'start', 'end', 'strand')
names(gro.plus)         <- c('chrom', 'start', 'end', 'tags')
names(gro.minus)        <- c('chrom', 'start', 'end', 'tags')
names(total.bayespeak1) <- c('chrom', 'start', 'end')
names(total.bayespeak2) <- c('chrom', 'start', 'end')
names(total.jamm)       <- c('chrom', 'start', 'end')
names(total.jamm.idr1)  <- c('chrom', 'start', 'end')
names(total.jamm.idr2)  <- c('chrom', 'start', 'end')
names(total.macs1)      <- c('chrom', 'start', 'end')
names(total.macs2)      <- c('chrom', 'start', 'end')
names(total.macs.idr1)  <- c('chrom', 'start', 'end')
names(total.macs.idr2)  <- c('chrom', 'start', 'end')
names(total.zerone)     <- c('chrom', 'start', 'end')

library(GenomicRanges)
gr.t.tss     <- GRanges(total.tss$chrom, IRanges(total.tss$start, total.tss$end), strand = total.tss$strand)
gr.gro.plus  <- GRanges(gro.plus$chrom,  IRanges(gro.plus$start,  gro.plus$end),  strand = '+', tags =  gro.plus$tags)
gr.gro.minus <- GRanges(gro.minus$chrom, IRanges(gro.minus$start, gro.minus$end), strand = '-', tags = -gro.minus$tags)
gr.t.gro <- c(gr.gro.plus, gr.gro.minus)
gr.t.bayespeak1 <- GRanges(total.bayespeak1$chrom, IRanges(total.bayespeak1$start, total.bayespeak1$end))
gr.t.bayespeak2 <- GRanges(total.bayespeak2$chrom, IRanges(total.bayespeak2$start, total.bayespeak2$end))
gr.t.jamm       <- GRanges(total.jamm$chrom,       IRanges(total.jamm$start,       total.jamm$end))
gr.t.jamm.idr1  <- GRanges(total.jamm.idr1$chrom,  IRanges(total.jamm.idr1$start,  total.jamm.idr1$end))
gr.t.jamm.idr2  <- GRanges(total.jamm.idr2$chrom,  IRanges(total.jamm.idr2$start,  total.jamm.idr2$end))
gr.t.macs1      <- GRanges(total.macs1$chrom,      IRanges(total.macs1$start,      total.macs1$end))
gr.t.macs2      <- GRanges(total.macs2$chrom,      IRanges(total.macs2$start,      total.macs2$end))
gr.t.macs.idr1  <- GRanges(total.macs.idr1$chrom,  IRanges(total.macs.idr1$start,  total.macs.idr1$end))
gr.t.macs.idr2  <- GRanges(total.macs.idr2$chrom,  IRanges(total.macs.idr2$start,  total.macs.idr2$end))
gr.t.zerone     <- GRanges(total.zerone$chrom,     IRanges(total.zerone$start,     total.zerone$end))

ol.bayespeak1 <- findOverlaps(gr.t.tss, gr.t.bayespeak1)
ol.bayespeak2 <- findOverlaps(gr.t.tss, gr.t.bayespeak2)
ol.jamm       <- findOverlaps(gr.t.tss, gr.t.jamm)
ol.jamm.idr1  <- findOverlaps(gr.t.tss, gr.t.jamm.idr1)
ol.jamm.idr2  <- findOverlaps(gr.t.tss, gr.t.jamm.idr2)
ol.macs1      <- findOverlaps(gr.t.tss, gr.t.macs1)
ol.macs2      <- findOverlaps(gr.t.tss, gr.t.macs2)
ol.macs.idr1  <- findOverlaps(gr.t.tss, gr.t.macs.idr1)
ol.macs.idr2  <- findOverlaps(gr.t.tss, gr.t.macs.idr2)
ol.zerone     <- findOverlaps(gr.t.tss, gr.t.zerone)

ol.active.tss <- findOverlaps(resize(gr.t.tss, fix = 'center', width = 1000), gr.t.gro)
ol.active.tss <- tapply(gr.t.gro[subjectHits(ol.active.tss)]$tags, queryHits(ol.active.tss), sum) > 10
ol.active.tss <- as.vector(as.numeric(names(ol.active.tss)))[ol.active.tss]

aol.bayespeak1 <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.bayespeak1)
aol.bayespeak2 <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.bayespeak2)
aol.jamm       <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.jamm)
aol.jamm.idr1  <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.jamm.idr1)
aol.jamm.idr2  <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.jamm.idr2)
aol.macs1      <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.macs1)
aol.macs2      <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.macs2)
aol.macs.idr1  <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.macs.idr1)
aol.macs.idr2  <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.macs.idr2)
aol.zerone     <- findOverlaps(gr.t.tss[ol.active.tss], gr.t.zerone)

precision <- function(cov, total) length(unique(cov)) / length(total)
recall <- function(cov, tss) length(cov) / tss
f1 <- function(cov, total, tss) {
    p <- precision(cov, total)
    r <- recall(cov, tss)
    print(length(total))
    print(length(unique(cov)))
    print(p)
    print(r)
    print(2 * p * r / (p + r))
}

n.tss <- nrow(total.tss)
f1(subjectHits(ol.bayespeak1), gr.t.bayespeak1, n.tss)
f1(subjectHits(ol.bayespeak2), gr.t.bayespeak2, n.tss)
f1(subjectHits(ol.jamm),       gr.t.jamm,       n.tss)
f1(subjectHits(ol.jamm.idr1),  gr.t.jamm.idr1,  n.tss)
f1(subjectHits(ol.jamm.idr2),  gr.t.jamm.idr2,  n.tss)
f1(subjectHits(ol.macs1),      gr.t.macs1,      n.tss)
f1(subjectHits(ol.macs2),      gr.t.macs2,      n.tss)
f1(subjectHits(ol.macs.idr1),  gr.t.macs.idr1,  n.tss)
f1(subjectHits(ol.macs.idr2),  gr.t.macs.idr2,  n.tss)
f1(subjectHits(ol.zerone),     gr.t.zerone,     n.tss)

n.active.tss <- length(ol.active.tss)
f1(subjectHits(aol.bayespeak1), gr.t.bayespeak1, n.active.tss)
f1(subjectHits(aol.bayespeak2), gr.t.bayespeak2, n.active.tss)
f1(subjectHits(aol.jamm),       gr.t.jamm,       n.active.tss)
f1(subjectHits(aol.jamm.idr1),  gr.t.jamm.idr1,  n.active.tss)
f1(subjectHits(aol.jamm.idr2),  gr.t.jamm.idr2,  n.active.tss)
f1(subjectHits(aol.macs1),      gr.t.macs1,      n.active.tss)
f1(subjectHits(aol.macs2),      gr.t.macs2,      n.active.tss)
f1(subjectHits(aol.macs.idr1),  gr.t.macs.idr1,  n.active.tss)
f1(subjectHits(aol.macs.idr2),  gr.t.macs.idr2,  n.active.tss)
f1(subjectHits(aol.zerone),     gr.t.zerone,     n.active.tss)
