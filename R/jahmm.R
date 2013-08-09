BaumWelch <- function(m, yz, theta, alpha, C1, C2, blocks, index,
      verbose=TRUE, maxiter=1000, tol=1e-6) {

   nas <- is.na(colSums(yz, na.rm=FALSE))

   r <- nrow(yz)
   n <- ncol(yz)

   n_na.rm <- sum(!nas)
   yz_na.rm <- yz[,!nas]

   index <- rep(-1L, n)
   index_na.rm <- rep(-1L, n_na.rm)

   Q <- matrix(.05/(m-1), nrow=m, ncol=m)
   diag(Q) <- ifelse(m > 1, .95, 1.0)

   p <- matrix(NA, ncol=m, nrow=r+1)
   q <- matrix(NA, ncol=m, nrow=r+1)
   ybar <- mean(yz[1,], na.rm=TRUE)
   # Fill in p row-wise.
   p[1,] <- ybar * C1
   q[1,] <- ybar * C2
   q[2,] <- p[2,] <- ybar
   if (m == 1) mult <- 1
   if (m == 2) mult <- c(1,2)
   if (m == 3) mult <- c(.5,1,2)
   for (i in 2:r) {
      zbar <- mean(yz[i,], na.rm=TRUE)
      q[i+1,] <- p[i+1,] <- zbar * mult
   }
   # Assert that values are well-defined.
   if (any(p < .Machine$double.eps) || any(q < .Machine$double.eps)) {
      stop('p or q undefined: check input to jahmm')
   }
   p <- scale(p, center=FALSE, scale=colSums(p))
   q <- scale(q, center=FALSE, scale=colSums(q))

   new.p <- p
   new.q <- q

   for (iter in 1:maxiter) {

      if (verbose) cat(paste("iteration:", iter, "\r"), file=stderr())
      oldparams <- c(p,q)

      # E-step.

      # Compute (unnormalized) emission probabilities.
      C_call_1 <- .C(
         mnmultinom_prob,
         # input #
         as.integer(m),
         as.integer(n),
         as.integer(r),
         as.integer(yz),
         # params #
         as.double(theta),
         as.double(alpha),
         as.double(p),
         as.double(q),
         # index #
         as.integer(index),
         # control #
         as.integer(0+4),  # No warning, linear space unless underflow.
         # output #
         double(m*n),      # Emission probabilities.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      initialProb <- steady_state_probs(Q)

      # Compute smoothing phis by the forward-backward algorithm.
      C_call_2 <- .C(
         block_fwdb,
         # input #
         as.integer(m),
         as.integer(length(blocks)),
         as.integer(blocks),
         # params #
         as.double(Q),
         as.double(initialProb),
         # output #
         C_call_1[[11]],  # Forward alphas.
         double(m*n),     # Smoothing phis.
         double(m*m),     # Sum of transitions.
         double(1),       # Log-likelihood.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      # M-step.
      Q <- matrix(C_call_2[[8]], nrow=m)
      Q <- Q / rowSums(Q)

      phi <- matrix(C_call_2[[7]], ncol=m, byrow=TRUE)[!nas,]

      for (j in 1:10) {

         old.p <- new.p
         old.q <- new.q

         C_call_1 <- .C(
            mnmultinom_prob,
            # input #
            as.integer(m),
            as.integer(n_na.rm),
            as.integer(r),
            as.integer(yz_na.rm),
            # params #
            as.double(theta),
            as.double(alpha),
            as.double(old.p),
            as.double(old.q),
            # index #
            as.integer(index_na.rm),
            # control #
            as.integer(2+4),     # No warning, ratio of mixture states.
            # output #
            double(m*n_na.rm),   # Probability ratios.
            # extra '.C()' arguments #
            NAOK = FALSE,
            DUP = FALSE
         )

         theta_i <- matrix(C_call_1[[11]], byrow=TRUE, ncol=m)
         zi1 <- scale(yz_na.rm %*% (phi * theta_i), center=FALSE,
                  scale=colSums(phi * theta_i))
         zi0 <- scale(yz_na.rm %*% (phi * (1-theta_i)), center=FALSE,
                  scale=colSums(phi * (1-theta_i)))

         # Add 'alpha' to first row.
         zi1[1,] <- alpha + zi1[1,]
         zi0[1,] <- alpha + zi0[1,]
         new.p <- scale(zi1, center=FALSE, scale=colSums(zi1))
         new.q <- scale(zi0, center=FALSE, scale=colSums(zi0))
         new.p[1,] <- new.p[1,] / (C1+1)
         new.q[1,] <- new.q[1,] / (C2+1)
         new.p <- rbind(C1*new.p[1,], new.p)
         new.q <- rbind(C2*new.q[1,], new.q)

         if (all(abs(c(new.p,new.q)-c(old.p,old.q))< 1e-3)) break

      }

      p <- new.p
      q <- new.q

      if (all(abs(oldparams - c(p,q)) < tol)) break

   }

   if (verbose) cat("\n", file=stderr())

   return (list(Q=Q, p=new.p, q=new.q, index=index, iter=iter))

}


jahmm <- function (data, PSO=TRUE, verbose=TRUE, threshold=0.07, ...) {

###############################################
#              OPTION PROCESSING              #
###############################################

   n <- nrow(data)

   # EM for mixture distribution of the baseline.
   baseline_params <- EM.mnb(data[,2])

   theta <- baseline_params[3]
   alpha <- baseline_params[2]
   C1 <- baseline_params[4]/(1-baseline_params[4])
   C2 <- baseline_params[5]/(1-baseline_params[5])

   # Coerce input data to matrix for C calls.
   yz <- t(as.matrix(data[,-1]))



###############################################
#                 MAIN LOOP                   #
###############################################


   # Keep only the profiles that pass quality control.
   #KLpq <- KL(yz, c(theta, alpha, C1, C2))
   #QC_scores <- KLpq$KLa - KLpq$KLb
   #if (all(QC_scores < QC)) return(list(QC=QC_scores>QC))
   #yz <- yz[c(TRUE, QC_scores > QC),]
   #p0 <- KLpq$p[c(TRUE, TRUE, QC_scores > QC)]
   #q0 <- KLpq$q[c(TRUE, TRUE, QC_scores > QC)]

   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   sorted <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[sorted]

   loglik <- rep(NA, 3)
   for (m in c(1,3)) {
      BW <- BaumWelch(m, yz, theta, alpha, C1, C2, blocks,
         verbose=verbose, ...)

      Q <- BW$Q
      p <- BW$p
      q <- BW$q

      index <- BW$index

      if (PSO && m > 1) {
         if (verbose) cat("running PSO\n", file=stderr())
         pso_call <- .C(
            pso,
            # input #
            as.integer(m),
            as.integer(n),
            as.integer(nrow(yz)),
            as.integer(yz),
            # fixed params #
            as.double(theta),
            as.double(alpha),
            # params #
            as.double(Q),
            as.double(p),
            as.double(q),
            # index #
            as.integer(index),
            # output #
            double(1),
            # extra '.C()' arguments #
            NAOK = TRUE,
            DUP = FALSE
         )
         Q <- array(pso_call[[7]], dim=dim(Q))
         p <- array(pso_call[[8]], dim=dim(p))
         q <- array(pso_call[[9]], dim=dim(q))
      }

      C_call_1 <- .C(
         mnmultinom_prob,
         # input #
         as.integer(m),
         as.integer(n),
         as.integer(nrow(yz)),
         as.integer(yz),
         # params #
         as.double(theta),
         as.double(alpha),
         as.double(p),
         as.double(q),
         # index #
         as.integer(index),
         # control #
         as.integer(1+4+8),   # No warning, log space, constant terms.
         # output #
         double(m*n),         # Emission probabilities.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )


      initialProb <- steady_state_probs(Q)

      if (m > 1) {
         vitC <- .C(
            block_viterbi,
            # input #
            as.integer(m),
            as.integer(length(blocks)),
            as.integer(blocks),
            as.double(log(Q)),
            as.double(log(initialProb)),
            as.double(C_call_1[[11]]),
            # control #
            as.integer(1),
            # output #
            integer(n),
            # extra '.C()' arguments #
            NAOK = TRUE,
            DUP = FALSE
         )

         vPath <- vitC[[8]]

      }

      C_call_2 <- .C(
         block_fwdb,
         # input #
         as.integer(m),
         as.integer(length(blocks)),
         as.integer(blocks),
         # params #
         as.double(Q),
         as.double(initialProb),
         # output #
         C_call_1[[11]],  # Forward alphas.
         double(m*n),     # Smoothing phis.
         double(m*m),     # Sum of transitions.
         double(1),       # Log-likelihood.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE
      )

      loglik[m] <- C_call_2[[9]]

   }

   # Quality control.
   score <- (loglik[1]-loglik[3]) / loglik[1]
   return(list(
      success = score > threshold,
      vPath = vPath,
      alpha = alpha,
      theta = theta,
      Q = Q,
      p = p,
      q = q,
      score = score,
      loglik = loglik[3]
   ))

}
