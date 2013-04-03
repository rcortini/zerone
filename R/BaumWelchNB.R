BaumWelch.NB <- function (data, Q, alpha=1, beta, gammas,
   maxiter=2000, tol=1e-5, verbose=TRUE) {

###############################################
#              OPTION PROCESSING              #
###############################################

   n <- nrow(data)
   y <- data[,2]
   z <- data[,3:ncol(data), drop=FALSE]
   r <- ncol(z)

   if (!missing(Q)) {
      if (!all.equal(dim(Q), c(3,3)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(c(
         0.95,  0.025, 0.025,
         0.025, 0.95,  0.025,
         0.025, 0.025, 0.95), nrow=3)
   }

   # Censor outliers in 'y' (they bear a huge weight on the estimation
   # of the 'alpha' parameter).
   #y[y > quantile(y, .9999, na.rm=TRUE)] <- NA
   ybar <- mean(y, na.rm=TRUE)
  
   if (missing(beta)) beta <- ybar / alpha
   if (missing(gammas)) {
      gammas <- matrix(NA, ncol=3, nrow=r)
      for (i in 1:r) gammas[i,] <- quantile(z[,i],
         prob=c(.3, .6, .9), na.rm=TRUE) / (alpha * beta)
   }



###############################################
#            VARIABLE DEFINITIONS             #
###############################################

   # Coerce 'z' for C.
   #yz_as_matrix <- t(as.matrix(data[,-1]))
   yz_as_matrix <- t(as.matrix(data.frame(y,z)))


   # Get complete cases (useful for M-step).
   complete <- complete.cases(data.frame(y,z))
   y0 <- y
   z0 <- t(as.matrix(z))
   y0[!complete] <- 0
   z0[,!complete] <- 0

   # Tabulation saves a lot of time in the M-step.
   tab.z <- tabulate(rowSums(z) + 1L)
   tab.z.y <- tabulate(y + rowSums(z) + 1L)
   u.z <- 0:(length(tab.z)-1L)
   u.z.y <- 0:(length(tab.z.y)-1L)

   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   id <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[id]


###############################################
#                 MAIN LOOP                   #
###############################################


   index <- as.integer(rep(-1L, n))

   for (iter in 1:maxiter) {

      if (verbose) cat(paste("iteration:", iter, "\r"))

      # E-step.

      C_call_1 <- .C("compute_pratio",
         # input #
         as.integer(n),
         as.integer(r),
         as.integer(yz_as_matrix),
         # params #
         as.double(alpha),
         as.double(beta),
         as.double(gammas),
         # index #
         index,
         # output #
         double(2*n),                    # Probability ratio.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE,
         PACKAGE = "HummingBee"
      )

      #initialProb <- steady_state_probs(Q)
      initialProb <- rep(1/3,3)

      if (any(is.na(C_call_1[[8]]))) return (list(loglik=NA))

      C_call_2 <- .C("fwdb",
         # input #
         as.integer(length(blocks)),
         as.integer(blocks),
         # params #
         Q,
         initialProb,
         # output #
         C_call_1[[8]],                # Forward alphas.
         matrix(double(2*n), nrow=2),  # Phi for first two states.
         double(3*3),                  # Sum of transitions.
         # extra '.C()' arguments #
         DUP = FALSE,
         NAOK = TRUE
      )

      phi <- cbind(t(C_call_2[[6]]), 1-colSums(C_call_2[[6]]))

      # M-step.
      Q <- matrix(C_call_2[[7]], nrow = 3)
      Q <- Q / rowSums(Q)

      # FIXME
      sumPhi <- colSums(phi[complete,])

      normPhi <- phi
      normPhi <- scale(normPhi, center=FALSE, scale=colSums(normPhi))
      mean.y <- colSums(normPhi*y0)
      mean.zj <- z0 %*% normPhi
      mean.z <- colSums(mean.zj)

      oldparams <- c(alpha, beta, gammas)

      dalpha <- 2*tol
      while(abs(dalpha) > tol) {
         alpha <- alpha + dalpha
         f <- -digamma(alpha) -log(1 + ybar/alpha) +
               sum(tab.z.y * digamma(alpha + u.z.y))/n -
               sum(sumPhi * log(1 + mean.z/(alpha+mean.y)))/n
         df <- -trigamma(alpha) + ybar/alpha/(alpha+ybar) +
               sum(tab.z.y * trigamma(alpha + u.z.y))/n +
               sum(sumPhi*mean.z/(alpha+mean.y) /
                  (alpha+mean.y+mean.z))/n
         dalpha = -f/df
         # In case of numerical instability, keep old alpha.
         if (is.na(alpha + dalpha) || (alpha + dalpha < 0)) {
            alpha <- oldparams[1]
            break
         }
      }

      alpha <- round(alpha, 4)
      beta <- round(ybar / alpha, 4)
      tTij <- t(mean.zj) / (alpha+mean.y+mean.z)
      gammas <- t(tTij*(1+1/beta)*(1+mean.z/(alpha+mean.y)))
      
      if (all(abs(oldparams - c(alpha, beta, gammas)) < tol)) break

   } # for (iter in 1:maxiter)

   if (verbose) cat("\n")

   loglik_C <- .C("compute_loglik",
      # input #
      as.integer(n),
      as.integer(r),
      as.integer(yz_as_matrix),
      # params #
      as.double(Q),
      as.double(alpha),
      as.double(beta),
      as.double(gammas),
      # index #
      index,
      as.integer(c(tab.z.y, -1)),
      # output #
      double(1),                    # Probability ratio.
      # extra '.C()' arguments #
      NAOK = TRUE,
      DUP = FALSE,
      PACKAGE = "HummingBee"
   )

#   # Viterbi algorithm.
#   C_call_1 <- .C("compute_pratio",
#      # input #
#      as.integer(n),
#      as.integer(r),
#      as.integer(yz_as_matrix),
#      # params #
#      as.double(alpha),
#      as.double(beta),
#      as.double(gammas),
#      # index #
#      index,
#      # output #
#      double(2*n),                    # Probability ratio.
#      # extra '.C()' arguments #
#      NAOK = TRUE,
#      DUP = FALSE,
#      PACKAGE = "HummingBee"
#   )
#
#   initialProb <- steady_state_probs(Q)
#   pem <- cbind(t(matrix(C_call_1[[8]], nrow=2)), rep(1,n))
#   pem[is.infinite(pem)] <- .Machine$double.xmax
#   pem <- pem / rowSums(pem)
#   log_p <- log(t(pem))
#   log_p[log_p < -320] <- -320
#   log_Q <- log(Q)
#   log_Q[log_Q < -320] <- -320
#
#   vitC <- .C("block_viterbi",
#      # input #
#      as.integer(m),
#      as.integer(length(blocks)),
#      as.integer(blocks),
#      log(initialProb),
#      log_p,
#      log_Q,
#      # output #
#      integer(n)
#   )
#
#
#   return(list(loglik=loglik_C[[10]] / nrow(HMMdata), Q=Q,
#      alpha=alpha, beta=beta, gamma=gammas, vPath=vitC[[7]],
#      emissionProb=pem, iterations=iter))
   return(list(loglik=loglik_C[[10]] / nrow(HMMdata), Q=Q,
      alpha=alpha, beta=beta, gammas=gammas, iterations=iter,
      blocks=blocks, index=index))

}
