wig2jahmm <- function (directory) {
   fnames <- list.files(directory, pattern="\\.wig(\\.gz)?$")
   inputs <- grep("(input|control)", fnames,
                 ignore.case=TRUE, value=TRUE)
   targets <- grep("(input|control)", fnames,
                   ignore.case=TRUE, value=TRUE, invert=TRUE)
   if (length(inputs) < 1) { stop("You need at least one Input profile.")
   } else { write(c("File(s) taken as Input:", inputs), stdout()) }
   if (length(targets) < 1) { stop("You need at least one ChIP profile.")
   } else { write(c("File(s) taken as Target:", targets), stdout()) }

   write("Reading files...", stdout())
   dfChrom <- wig2df(fnames[1], isChrom=TRUE)
   dfInputs <- apply(sapply(X=inputs, FUN=wig2df), 1, sum)
   dfTargets <- data.frame(lapply(X=targets, FUN=wig2df))
   names(dfTargets) <- seq(length(dfTargets))
   dfAll <- data.frame(chrom=dfChrom, input=dfInputs, dfTargets)

   write("Running jahmm...", stdout())
   return(jahmm(dfAll))
}

wig2df <- function(fname, isChrom=FALSE) {
   input <- readLines(fname)
   input <- grep("^(#|browser|track)", input, value=TRUE, invert=TRUE)
   fileLength <- length(input)

   decLinesInd <- grep("^fixedStep", input)
   wholeDecLines <- input[decLinesInd]
   decLines <- matrix(unlist(strsplit(wholeDecLines, " ")),
                      nrow=length(wholeDecLines), byrow=TRUE)

   if (isChrom) {
      chrom <- matrix(unlist(strsplit(decLines[, 2], "=")),
                      ncol=2, byrow=TRUE)[, 2]
      chromRep <- c(diff(decLinesInd) - 1, fileLength - max(decLinesInd))
      chromCol <- rep(chrom, chromRep)
      return(chromCol)
   } else {
      scores <- as.numeric(grep("^fixedStep", input, value=TRUE, invert=TRUE))
      return(scores)
   }
}
