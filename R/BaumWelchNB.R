BaumWelch.NB <- function (x, y, m=2, Q=NULL, alpha=1, gamma=c(1,2),
   beta=1, initialProb=NULL, maxiter=500, tol=1e-05, dig=3) {

# XXX Development XXX
stopifnot(m == 2)

###############################################
#              OPTION PROCESSING              #
###############################################

   if (is.vector(x))
      x <- list(x)
   if (is.list(x)) {
      if (!all(sapply(X = x, FUN = is.vector))) {
         stop ("non vector element found in list x")
      }
      else if (!all(sapply(X = x, FUN = is.numeric))) {
         stop ("x must be numeric")
      }
      else {
         concat.x <- unlist(x)
      }
   }
   else {
      stop ("x must be a vector or a list")
   }

   if (!is.null(y)) {
      if (is.vector(y))
         y <- list(y)
      if (is.list(y)) {
         if (!all(sapply(X = y, FUN = is.vector))) {
            stop ("non vector element found in list y")
         }
         else if (!all(sapply(X = y, FUN = is.numeric))) {
            stop ("y must be numeric")
         }
         else if (!identical(sapply(X = x, FUN = length),
            sapply(X = y, FUN = length))) {
            stop ("lists x and y must have same structure")
         }
         else {
            concat.y <- unlist(y)
         }
      }
      else {
         stop ("y must be a vector or a list")
      }
   }
   else {
      concat.y <- rep(0, length(concat.x))
   }

   if (!is.null(Q)) {
      if (!all.equal(dim(Q), c(m,m)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(0.1/(m-1), ncol = m, nrow = m)
      diag(Q) <- 0.9
   }

#   if (is.null(theta.x)) {
#      init.x <- 2*(1:m)*mean(concat.x, na.rm = TRUE) / (m+1)
#      init.y <- mean(concat.y, na.rm = TRUE)
#      tilde.theta.x <- init.x / (init.x + alpha)
#      tilde.theta.y <- init.y / (init.y + alpha)
#      theta.x <- tilde.theta.x * (1 - tilde.theta.y) /
#         (1 - tilde.theta.x*tilde.theta.y)
#      theta.y <- tilde.theta.y * (1 - tilde.theta.x) /
#         (1 - tilde.theta.x*tilde.theta.y)
#      theta.alpha <- 1 - theta.x - theta.y
#   }

   adjustInitialProb <- ifelse(is.null(initialProb), TRUE, FALSE)




###############################################
#            VARIABLE DEFINITIONS             #
###############################################


   n <- length(concat.x)
   blockSizes <- sapply(X = x, FUN = length)

   oldLogLikelihood <- -Inf

   phi <- emissionProb <- matrix(NA_real_, nrow = n, ncol = m)

   # Tabulation saves a lot of time at the M step.
   counts <- c(sum(concat.x + concat.y == 0), tabulate(bin = concat.x
      + concat.y))
   levels <- (0:(length(counts)-1))[counts > 0]
   counts <- counts[counts > 0]
   



###############################################
#           FUNCTION DEFINITIONS              #
###############################################


   initial.steady.state.probabilities <- function () {
   # Compute the steady-state initial probabilities.

      spectralDecomposition <- eigen(t(Q))

      if (is.complex(spectralDecomposition$values[1]))
         return(rep(1/m,m))
      if (spectralDecomposition$values[1] > 0.99)
         return( spectralDecomposition$vectors[,1] /
            sum(spectralDecomposition$vectors[,1]))

      return(rep(1/m,m))

   }

   E.step <- function () {
   # Performs the E-step of the modified Baum-Welch algorithm

      # Emission probabilities are computed up to constant terms.
      # The only term that depend on the state are
      #
      #        gamma_i^z / (1+gamma_i+1/beta)^(alpha+y+z)
      #
      # The other terms are merged into a multiplicative constant
      # which is normalized away.
      for (i in 1:m)
#         emissionProb[,i] <<- (theta.alpha[i]**alpha) * (theta.x[i]**concat.x) *
#            (theta.y[i]**concat.y)
         emissionProb[,i] <<- (gamma[i]^concat.x) /
                  (1+gamma[i]+1/beta)^(alpha+concat.y+concat.x)

      # Emission probabilities of NAs are set to 1 for every state.
      emissionProb[is.na(emissionProb)] <<- 1

      if (adjustInitialProb)
         initialProb <- initial.steady.state.probabilities()

      logLikelihood <<- 0
      cumulative.n <- 0
      transitions <- matrix(double(m * m), nrow = m)

      counter = 0

      for (n.i in blockSizes) {

         counter = counter+1

         forwardback <- .Fortran("fwdb", as.integer(m), as.integer(n.i), initialProb, 
            emissionProb[(cumulative.n + 1):(cumulative.n + n.i),], Q, double(n.i),
            matrix(double(n.i * m), nrow = n.i), matrix(double(m^2), nrow = m), double(1),
            PACKAGE = "HummingBee")

         logLikelihood <<- logLikelihood + forwardback[[9]]
         phi[(cumulative.n + 1):(cumulative.n + n.i),] <<- forwardback[[7]]
         transitions <- transitions + forwardback[[8]]

         cumulative.n <- cumulative.n + n.i

      }

      Q <<- transitions / rowSums(transitions)

   }




###############################################
#                 MAIN LOOP                   #
###############################################


   ybar = mean(concat.y)

   for (iter in 1:maxiter) {

      cat(paste("iteration:", iter, "\n"))

      E.step()

      # M-step.

      sumPhi <- colSums(phi, na.rm = TRUE)
      sumPhi.x <- colSums(phi*concat.x, na.rm = TRUE)
      sumPhi.y <- colSums(phi*concat.y, na.rm = TRUE)
      mean.x <- sumPhi.x / sumPhi
      mean.y <- sumPhi.y / sumPhi

      new.alpha <- alpha + 2*tol
      while(abs(new.alpha - alpha) > tol) {
         alpha <- new.alpha
         # TODO: tabulate to gain speed.
         f <- -digamma(alpha) - log(1+ybar/alpha) +
               mean(digamma(alpha+concat.x+concat.y)) -
               sum(sumPhi * log(1+mean.x/(alpha+mean.y)))/n
         df <- -trigamma(alpha) + (ybar/alpha)/(alpha+ybar) +
               mean(trigamma(alpha+concat.x+concat.y)) +
               sum(sumPhi*mean.x/((alpha+mean.y)
                  *(alpha+mean.y+mean.x)))/n
         new.alpha <- alpha - f/df
      }

      # Note: the rounding prevents the oscillation of the estimates.
      alpha <- round(new.alpha, dig)
      beta <- ybar / alpha
      gamma <- (1+1/beta)*mean.x/(alpha+mean.y)
      

#      beta.plus.1 <- sum((sumPhi.x + sumPhi.y + alpha*sumPhi) / (1 + zeta)) / (n*alpha)
#      beta.gamma <- beta.plus.1 * zeta
#      theta.x <- round(beta.gamma / (beta.plus.1 + beta.gamma), dig)
#      theta.y <- round((beta.plus.1-1) / (beta.plus.1 + beta.gamma), dig)
#      theta.alpha <- 1 - theta.x - theta.y

#      cat(paste(c(alpha, theta.x), "\n"))

      if (abs(logLikelihood - oldLogLikelihood) < tol)
         break

      oldLogLikelihood <- logLikelihood

   } # for (iter in 1:maxiter)

   cat("\n")

   if (adjustInitialProb) 
      initialProb <- initial.steady.state.probabilities()
      
   vPath <- Viterbi(Q, initialProb, emissionProb, blockSizes)

#   return(list(logL = logLikelihood, Q = Q, alpha = alpha, theta.x =
#      theta.x, theta.y = theta.y, vPath = vPath, iterations = iter))
   return(list(logL=logLikelihood, Q=Q, alpha=alpha, gamma=gamma,
      beta=beta, vPath = vPath, iterations = iter))
}
