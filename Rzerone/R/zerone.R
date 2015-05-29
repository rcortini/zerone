zerone <- function(y) {
   stopifnot(is.data.frame(y))
   stopifnot(ncol(y) > 2)
   retval <- .Call(zerone_R_call, y)
   names(retval) <- c("Q", "a", "pi", "p", "phi", "pem", "path", "l")
   return(retval)
}
