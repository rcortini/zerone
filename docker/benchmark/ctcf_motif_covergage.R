setwd('/')

total.motif <- read.table('benchmark/ctcf_motifs.bed')

total.macs.idr <- read.table('benchmark/idrCode/MACS_q0.05/Ctcf/IDR_Ctcf-overlapped-peaks.txt')
total.jamm.idr <- read.table('benchmark/idrCode/JAMM/Ctcf-overlapped-peaks.txt')
pass.idr.jamm <- subset(total.jamm.idr, IDR < 0.05)
pass.idr.macs <- subset(total.macs.idr, IDR < 0.05)

total.bayespeak1 <- read.table('BayesPeak/Ctcf1.bed')
total.bayespeak2 <- read.table('BayesPeak/Ctcf2.bed')
total.jamm       <- read.table('JAMM/Ctcf/peaks/filtered.peaks.narrowPeak')[1:3]
total.macs1      <- read.table('MACS/Ctcf1.bed')[1:3]
total.macs2      <- read.table('MACS/Ctcf2.bed')[1:3]
total.zerone     <- read.table('zerone/Ctcf.bed')
total.jamm.idr1 <- pass.idr.jamm[, 1:3]
total.jamm.idr2 <- pass.idr.jamm[, 5:7]
total.macs.idr1 <- pass.idr.macs[, 1:3]
total.macs.idr2 <- pass.idr.macs[, 5:7]

names(total.motif)      <- c('chrom', 'start', 'end')
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
gr.t.motif      <- GRanges(total.motif$chrom,      IRanges(total.motif$start,      total.motif$end))
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

ol.bayespeak1 <- findOverlaps(gr.t.motif, gr.t.bayespeak1)
ol.bayespeak2 <- findOverlaps(gr.t.motif, gr.t.bayespeak2)
ol.jamm       <- findOverlaps(gr.t.motif, gr.t.jamm)
ol.jamm.idr1  <- findOverlaps(gr.t.motif, gr.t.jamm.idr1)
ol.jamm.idr2  <- findOverlaps(gr.t.motif, gr.t.jamm.idr2)
ol.macs1      <- findOverlaps(gr.t.motif, gr.t.macs1)
ol.macs2      <- findOverlaps(gr.t.motif, gr.t.macs2)
ol.macs.idr1  <- findOverlaps(gr.t.motif, gr.t.macs.idr1)
ol.macs.idr2  <- findOverlaps(gr.t.motif, gr.t.macs.idr2)
ol.zerone     <- findOverlaps(gr.t.motif, gr.t.zerone)

precision <- function(cov, total) length(unique(cov)) / length(total)
recall <- function(cov) length(cov) / nrow(total.motif)
f1 <- function(cov, total) {
    p <- precision(cov, total)
    r <- recall(cov)
    print(length(total))
    print(length(unique(cov)))
    print(p)
    print(r)
    print(2 * p * r / (p + r))
}

f1(subjectHits(ol.bayespeak1), gr.t.bayespeak1)
f1(subjectHits(ol.bayespeak2), gr.t.bayespeak2)
f1(subjectHits(ol.jamm),       gr.t.jamm)
f1(subjectHits(ol.jamm.idr1),  gr.t.jamm.idr1)
f1(subjectHits(ol.jamm.idr2),  gr.t.jamm.idr2)
f1(subjectHits(ol.macs1),      gr.t.macs1)
f1(subjectHits(ol.macs2),      gr.t.macs2)
f1(subjectHits(ol.macs.idr1),  gr.t.macs.idr1)
f1(subjectHits(ol.macs.idr2),  gr.t.macs.idr2)
f1(subjectHits(ol.zerone),     gr.t.zerone)
