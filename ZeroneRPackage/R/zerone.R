zerone <- function(y) {
   # The data frame must contain seqname, mock and profiles.
   stopifnot(is.data.frame(y))
   stopifnot(ncol(y) > 2)
   # Use run length encoding to get block lengths.
   x <- rle(as.character(y[,1]))
   names <- as.character(x$values)
   size <- as.integer(x$lengths)
   retval <- .Call(zerone_R_call, names, size, y)
   names(retval) <- c("Q", "a", "pi", "p", "phi", "pem", "path", "l")
   return(retval)
}
