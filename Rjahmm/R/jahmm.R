jahmm <- function(m, y) {
   stopifnot(m > 0)
   stopifnot(is.data.frame(y))
   stopifnot(ncol(y) > 1)
   retval <- .Call(jahmm_R_call, m, y)
   names(retval) <- c("Q", "a", "pi", "p", "phi", "pem", "path", "l")
   return(retval)
}
