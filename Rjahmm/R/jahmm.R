jahmm <- function(m, y) {
   stopifnot(m > 0)
   stopifnot(is.data.frame(y))
   stopifnot(ncol(y) > 1)
   .Call(jahmm_R_call, m, y)
}
