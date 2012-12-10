BaumWelch.NB <- function (data, m=2, Q=NULL, alpha=NULL, gamma=c(1,2),
   beta=1, initialProb=NULL, maxiter=500, tol=1e-9, dig=5) {


###############################################
#              OPTION PROCESSING              #
###############################################


   n <- nrow(data)
   concat.x = data[,2]
   concat.y = data[,3]

   if (!is.null(Q)) {
      if (!all.equal(dim(Q), c(m,m)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(0.1/(m-1), ncol = m, nrow = m)
      diag(Q) <- 0.9
   }

   # Needed??
   concat.y[concat.y > 20*median(concat.y)] <- NA
   ybar <- mean(concat.y, na.rm=TRUE)

   if (is.null(alpha)) {
      tab <- tabulate(concat.y + 1L)
      u <- 0:(length(tab) - 1L)
      a <- 1
      da <- 2*tol
      while (abs(da) > tol) {
         a <- a + da
         f <- sum(tab*digamma(a+u))/n -digamma(a) - log(1+ybar/a)
         f. <- sum(tab*trigamma(a+u))/n -trigamma(a) + ybar/a/(a+ybar)
         da <- -f/f.
         for (ii in 1:20) {
            if (a + da > 0) break
            da <- da/2
         }
         if (a + da < 0) {
            a <- runif(1)
            da <- 2*tol
         }
      }
      alpha <- a
   }
   
   if (is.null(beta))
      beta <- mean(concat.y) / alpha

   adjustInitialProb <- ifelse(is.null(initialProb), TRUE, FALSE)




###############################################
#            VARIABLE DEFINITIONS             #
###############################################


   #blockSizes <- sapply(X = x, FUN = length)
   blockSizes = tapply(X=rep(1,n), INDEX=data[,1], FUN=sum)

   oldLogLikelihood <- -Inf

   phi <- emissionProb <- matrix(NA_real_, nrow = m, ncol = n)

   # Tabulation saves a lot of time.
   tab.x <- tabulate(concat.x + 1L)
   tab.x.y <- tabulate(concat.x + concat.y + 1L)
   u.x <- 0:(length(tab.x)-1L)
   u.x.y <- 0:(length(tab.x.y)-1L)



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
      #        gamma_i^z / (1+gamma_i+1/beta)^(alpha+y+z)
      # The other terms are merged into a multiplicative constant
      # which is normalized away.
      for (i in 1:m) {
         n.values <- gamma[i]^u.x 
         d.values <- (1+1/beta+gamma[i])^(alpha+u.x.y)
         emissionProb[i,] <<- n.values[concat.x + 1L] /
               d.values[concat.x+concat.y + 1L]
      }

      undef <- colSums(emissionProb) == 0
      emissionProb[,undef] <<- 1

      # Emission probabilities of NAs are set to 1 for every state.
      emissionProb[is.na(emissionProb)] <<- 1

      if (adjustInitialProb)
         initialProb <- initial.steady.state.probabilities()

      logLikelihood <<- 0
      cumulative.n <- 0
      transitions <- matrix(double(m*m), nrow = m)


      forwardback <- .C(
         "R_fwdb",
         as.integer(m),
         as.integer(length(blockSizes)),
         as.integer(blockSizes),
         Q,
         initialProb,
         emissionProb,
         double(n*m),
         double(m*m),
         double(1),
         as.integer(0),
         PACKAGE = "HummingBee"
      )

      logLikelihood <<- logLikelihood + forwardback[[9]]
      phi <<- matrix(forwardback[[7]], ncol=n)
      transitions <- transitions + forwardback[[8]]
      Q <<- transitions / rowSums(transitions)

   }




###############################################
#                 MAIN LOOP                   #
###############################################


   for (iter in 1:maxiter) {

      cat(paste("iteration:", iter, "\r"))

      E.step()

      # M-step.

      sumPhi <- rowSums(phi, na.rm = TRUE)
      sumPhi.x <- colSums(t(phi)*concat.x, na.rm = TRUE)
      sumPhi.y <- colSums(t(phi)*concat.y, na.rm = TRUE)
      mean.x <- sumPhi.x / sumPhi
      mean.y <- sumPhi.y / sumPhi

      dalpha <- 2*tol
      while(abs(dalpha) > tol) {
         alpha <- alpha + dalpha
         f <- -digamma(alpha) -log(1 + ybar/alpha) +
               sum(tab.x.y * digamma(alpha + u.x.y))/n -
               sum(sumPhi * log(1 + mean.x/(alpha+mean.y)))/n
         df <- -trigamma(alpha) + ybar/alpha/(alpha+ybar) +
               sum(tab.x.y * trigamma(alpha + u.x.y))/n +
               sum(sumPhi*mean.x/(alpha+mean.y) /
                  (alpha+mean.y+mean.x))/n
         dalpha = -f/df
         for (ii in 1:20) {
            if (alpha + dalpha > 0) break
            dalpha = dalpha/2
         }
         if (alpha + dalpha < 0) {
            alpha <- runif(1)
            dalpha <- 2*tol
         }
      }

      beta <- ybar / alpha
      gamma <- (1+1/beta)*mean.x/(alpha+mean.y)
      
      if (abs(logLikelihood - oldLogLikelihood) < tol)
         break

      oldLogLikelihood <- logLikelihood

   } # for (iter in 1:maxiter)

   cat("\n")

   if (adjustInitialProb) 
      initialProb <- initial.steady.state.probabilities()

   vPath <- Viterbi(Q, initialProb, t(emissionProb), blockSizes)

   # XXX DEBUG XXX
   print(c(alpha, beta, gamma, logLikelihood))

   return(list(logL=logLikelihood, Q=Q, alpha=alpha, gamma=gamma,
      beta=beta, vPath = vPath, iterations = iter))
}
