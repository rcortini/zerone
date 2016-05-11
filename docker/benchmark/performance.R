running.times <- matrix(  # BayesPeak
                        c(0,
                          mean(c(336218.775, 332847.403)),
                          0,
                          mean(c(338192.231, 336234.144)),
                          0,
                          mean(c(60 * 5413 + 26.109, 60 * 5416 + 24.009)),
                          0,
                          mean(c(335221.933, 335329.401)),
                          # JAMM
                          # 24 * 60 * 60 + (60 * (8792 + 15)),
                          60 * 250,
                          60 *  653 + 35.299,
                          24 * 60 * 60 + (60 * (8792 + 15)),
                          60 *  549 + 42.371,
                          60 * (8792 + 15),
                          60 * 1064 + 50.512,
                          0,
                          60 * 1127 +  7.212, # rerun with -r region?
                          # MACS
                          60 * 3 + 45.374,
                          mean(c(60 * 11 + 35.446, 60 * 10 + 25.626)),
                          60 * 4 + 47.285,
                          mean(c(60 * 17 + 32.532, 60 * 17 + 25.122)),
                          60 * 2 + 20.222,
                          mean(c(60 * 13 + 42.804, 60 * 13 +  8.881)),
                          0,
                          mean(c(60 * 11 + 31.013, 60 * 11 + 48.244)), # rerun with --broad?
                          # Zerone
                          0,
                          60 * 3 + 29.667,
                          0,
                          60 * 3 +  0.624,
                          0,
                          60 * 3 + 29.317,
                          0,
                          60 * 3 + 26.677),
                          nrow=8, ncol=4)

colnames(running.times) <- c('BayesPeak', 'JAMM', 'MACS', 'Zerone')
rownames(running.times) <- rep(c('CTCF', 'Pol2', 'H3K4me3', 'H3K36me3'), each=2)

running.times[c(1, 3, 5), c(2, 3)] <- running.times[c(1, 3, 5), c(2, 3)] + running.times[c(2, 4, 6), c(2, 3)]
log.running.times <- log2(running.times)
log.running.times[!is.finite(log.running.times)] <- 0
log.running.times[c(1, 3, 5), c(2, 3)] <- log.running.times[c(1, 3, 5), c(2, 3)] - log.running.times[c(2, 4, 6), c(2, 3)]

used.memory <- matrix(  # BayesPeak
                      c(mean(c(2038888, 2017828)),
                        mean(c(3119264, 2967020)),
                        mean(c(1885588, 1713604)),
                        mean(c(1981264, 1779136)),
                        # JAMM
                        1644248,
                        8931512,
                        10781632, # may need to be repeated.
                        1518724, # rerun with -r region?
                        # MACS
                        mean(c(613064, 571604)),
                        mean(c(816888, 787448)),
                        mean(c(632112, 570808)),
                        mean(c(582816, 645516)), # rerun with --broad?
                        # Zerone
                        663932,
                        654444,
                        775260,
                        775320),
                      nrow=4, ncol=4) / 1024^2

colnames(used.memory) <- c('BayesPeak', 'JAMM', 'MACS', 'Zerone')
rownames(used.memory) <- c('CTCF', 'Pol2', 'H3K4me3', 'H3K36me3')

setwd('/benchmark/')
pdf('performance.pdf', height=3)
layout(matrix(1:4, ncol=4, byrow=TRUE))

# Running times
par(mar=c(6, 5, 4, 0), cex=1)
barplot(log.running.times[2:1, ], col=c('grey', 'dimgrey'), main='CTCF', las=2, ylab=expression('log'[2]*'(Running time) (s)'))
# par(mar=c(5, 3, 4, 1))
# barplot(log.running.times[4:3, ], col=c('grey', 'dimgrey'), main='Pol2', las=2)
# par(mar=c(5, 3, 4, 1))
# barplot(log.running.times[6:5, ], col=c('grey', 'dimgrey'), main='H3K4me3', las=2)
par(mar=c(6, 3, 4, 2))
barplot(log.running.times[8:7, ], col=c('grey', 'dimgrey'), main='H3K36me3', las=2)

# Memory footprint
par(mar=c(6, 4, 4, 1), cex=1)
barplot(used.memory[1, ], main='CTCF', las=2, ylab='Memory footprint (GB)')
# par(mar=c(6, 3, 3, 1))
# barplot(used.memory[2, ], main=NA, las=2)
# par(mar=c(6, 3, 3, 1))
# barplot(used.memory[3, ], main=NA, las=2)
par(mar=c(6, 3, 4, 2))
barplot(used.memory[4, ], main='H3K36me3', las=2)

par(mar=c(5, 4, 4, 2))
layout(1)
dev.off()