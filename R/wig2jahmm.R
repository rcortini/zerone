wig2jahmm <- function (directory) {
   fnames <- list.files(directory)
   input <- grep("^Input", fnames, ignore.case=TRUE, value=TRUE)
   targets <- grep("^Input", fnames, ignore.case=TRUE, value=TRUE, invert=TRUE)
   if (length(input) != 1) stop("You need exactly one Input profile.")
   if (length(targets) < 1) stop("You need at least one ChIP profile.")

   dfInput <- wig2df(input)
   dfAll <- data.frame(dfInput, lapply(X=targets, FUN=wig2df, onlyScore=TRUE))

   library(jahmm)
   return(jahmm(dfAll))
}

wig2df <- function(fname, onlyScore=FALSE) {
   input <- readLines(fname)
   input <- grep("^(#|browser|track)", input, value=TRUE, invert=TRUE)
   fileLength <- length(input)

   wholeDecLines <- grep("^fixedStep", input, value=TRUE)
   decLines <- matrix(unlist(strsplit(wholeDecLines, " ")),
                      nrow=length(decLinesInd), byrow=TRUE)
   chrom <- matrix(unlist(strsplit(decLines[, 2], "=")),
                   ncol=2, byrow=TRUE)[, 2]
   start <- as.numeric(matrix(unlist(strsplit(decLines[, 3], "=")),
                              ncol=2, byrow=TRUE)[, 2])
   step <- as.numeric(matrix(unlist(strsplit(decLines[, 4], "=")),
                             ncol=2, byrow=TRUE)[, 2])
   if (dim(decLines)[2] > 4) {
     span <- as.numeric(matrix(unlist(strsplit(decLines[, 5], "=")),
                               ncol=2, byrow=TRUE)[, 2])
     } else span <- FALSE

   input <- as.numeric(grep("^fixedStep", input, value=TRUE, invert=TRUE))
   chromRep <- c(diff(decLinesInd) - 1, fileLength - max(decLinesInd))
   chromCol <- rep(chrom, chromRep)

   if (onlyScore) return(data.frame(score=input))
   else return(data.frame(chrom=chromCol, score=input)) 
}
