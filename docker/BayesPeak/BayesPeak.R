library(BayesPeak)
args <- commandArgs(trailingOnly=TRUE)
raw.output <- bayespeak(args[1], args[2])
output <- summarize.peaks(raw.output)
write.table(as.data.frame(output)[, 1:3],
            #paste0(args[1], '_tmp'),
            'tmp',
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
