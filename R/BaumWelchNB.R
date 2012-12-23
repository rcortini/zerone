BaumWelch.NB <- function (data, Q, alpha=1, beta, gamma=c(1,2),
   maxiter=1000, tol=1e-5) {

###############################################
#              OPTION PROCESSING              #
###############################################

   # Number of states is hard coded.
   m = 2

   n <- nrow(data)
   x = data[,2]
   y = data[,3]

   if (!missing(Q)) {
      if (!all.equal(dim(Q), c(m,m)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(0.05/(m-1), ncol=m, nrow=m)
      diag(Q) <- 0.95
   }

   # Censor outliers in 'y' (they bear a huge weight on the estimation
   # of the 'alpha' parameter).
   y[y > quantile(y, .9999)] <- NA
   ybar <- mean(y, na.rm=TRUE)
  
   if (missing(beta)) beta <- ybar / alpha



###############################################
#            VARIABLE DEFINITIONS             #
###############################################


   # Tabulation saves a lot of time in the M-step.
   tab.x <- tabulate(x + 1L)
   tab.x.y <- tabulate(x + y + 1L)
   u.x <- 0:(length(tab.x)-1L)
   u.x.y <- 0:(length(tab.x.y)-1L)

   blockSizes <- tapply(X=rep(1,n), INDEX=data[,1], FUN=sum)


###############################################
#           FUNCTION DEFINITIONS              #
###############################################

   steady_state_probs <- function (Q) {
      eig <- eigen(t(Q))
      if (is.complex(eig$values[1])) return(rep(1/m,m))
      if (eig$values[1] > 0.99)
         return(eig$vectors[,1] / sum(eig$vectors[,1]))
      return(rep(1/m,m))
   }



###############################################
#                 MAIN LOOP                   #
###############################################


   for (iter in 1:maxiter) {

      cat(paste("iteration:", iter, "\r"))

      # E-step.
      initialProb <- steady_state_probs(Q)

      ccall1 <- .C("compute_pratio",
         # input #
         as.integer(n),
         as.integer(x),
         as.integer(y),
         # params #
         as.double(alpha),
         as.double(beta),
         as.double(gamma),
         # output #
         double(n),                    # Probability ratio.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE,
         PACKAGE = "HummingBee"
      )

      ccall2 <- .C("fwdb",
         # input #
         as.integer(length(blockSizes)),
         as.integer(blockSizes),
         # params #
         diag(Q),
         initialProb,
         # output #
         ccall1[[7]],                  # Forward alphas.
         double(n),                    # Phi for first state.
         double(2*2),                  # Sum of transitions.
         # extra '.C()' arguments #
         DUP = FALSE,
         NAOK = FALSE
      )

      phi <- cbind(ccall2[[6]], 1-ccall2[[6]])

      # M-step.
      Q <- matrix(ccall2[[7]], nrow = m)
      Q <- Q / rowSums(Q)

      sumPhi <- colSums(phi)
      mean.x <- colSums(phi*x, na.rm=TRUE) / sumPhi
      mean.y <- colSums(phi*y, na.rm=TRUE) / sumPhi

      oldparams <- c(alpha, beta, gamma)

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
         for (j in 1:20) {
            if (alpha + dalpha > 0) break
            dalpha = dalpha/2
         }
         if (alpha + dalpha < 0) {
            alpha <- runif(1)
            dalpha <- 2*tol
         }
      }

      alpha <- round(alpha, 6)
      beta <- round(ybar / alpha, 6)
      gamma <- round((1+1/beta)*mean.x/(alpha+mean.y), 6)
      
      if (all(abs(oldparams - c(alpha, beta, gamma)) < tol)) break

   } # for (iter in 1:maxiter)

   cat("\n")

   # Viterbi algorithm.
   ccall1 <- .C("compute_pratio",
      # input #
      as.integer(n),
      as.integer(x),
      as.integer(y),
      # params #
      as.double(alpha),
      as.double(beta),
      as.double(gamma),
      # output #
      double(n),                    # Probability ratio.
      # extra '.C()' arguments #
      NAOK = TRUE,
      DUP = FALSE,
      PACKAGE = "HummingBee"
   )

   initialProb <- steady_state_probs(Q)
   pem <- cbind(ccall1[[7]], rep(1,n))
   pem <- pem / rowSums(pem)
   log_pem <- log(t(pem))
   log_pem[log_pem < -320] <- -320

   vitC <- .C("block_viterbi",
      # input #
      as.integer(length(blockSizes)),
      as.integer(blockSizes),
      log(initialProb),
      log_pem,
      log(Q),
      # output #
      integer(n)
   )

   return(list(Q=Q, alpha=alpha, beta=beta, gamma=gamma,
      vPath=vitC[[6]], emissionProb=pem, iterations=iter))

}
