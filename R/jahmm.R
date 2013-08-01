BaumWelch <- function(m, yz, theta, alpha, Q, p, q,
   blocks, index, verbose=TRUE, maxiter=1000, tol=1e-6) {

   nas <- is.na(colSums(yz, na.rm=FALSE))

   r <- nrow(yz)
   n <- ncol(yz)

   n_na.rm <- sum(!nas)
   yz_na.rm <- yz[,!nas]

   index <- rep(-1L, n)
   index_na.rm <- rep(-1L, n_na.rm)

   C1 <- p[2]/p[1]
   C2 <- q[2]/q[1]

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
         as.integer(4),   # No warning, linear space unless underflow.
         # output #
         double(m*n),     # Emission probabilities.
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
      Q <- matrix(C_call_2[[8]], nrow = m)
      Q <- Q / rowSums(Q)

      phi <- matrix(C_call_2[[7]], ncol=m, byrow=TRUE)[!nas,]
      new.p <- p
      new.q <- q

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
            as.integer(6),       # No warning, ratio of mixture states.
            # output #
            double(m*n_na.rm),   # Emission probabilities.
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

   return (list(loglik=C_call_2[[9]]/n,
      Q=Q, p=new.p, q=new.q, index=index))

}


jahmm <- function (data, PSO=TRUE, verbose=TRUE, ...) {

###############################################
#              OPTION PROCESSING              #
###############################################

   m <- 3
   n <- nrow(data)

   Q <- matrix(.05/(m-1), nrow=m, ncol=m)
   diag(Q) <- ifelse(m > 1, .95, 1.0)

   # EM for mixture distribution of the baseline.
   baseline_params <- EM.mnb(data[,2])

   theta <- baseline_params[3]
   alpha <- baseline_params[2]
   C1 <- baseline_params[4]/(1-baseline_params[4])
   C2 <- baseline_params[5]/(1-baseline_params[5])

   p <- matrix(NA, ncol=m, nrow=ncol(data))
   q <- matrix(NA, ncol=m, nrow=ncol(data))
   med <- median(data[,2], na.rm=TRUE)
   # Fill in p row-wise.
   p[1,] <- med * C1
   q[1,] <- med * C2
   q[2,] <- p[2,] <- med
   p_ <- p[,1,drop=FALSE]
   q_ <- q[,1,drop=FALSE]
   for (i in 3:ncol(data)) {
      zbar <- mean(data[,i], na.rm=TRUE)
      q[i,] <- p[i,] <- zbar * c(.5, 1, 2)
      q_[i,] <- p_[i,] <- zbar
   }
   # Assert that values are well-defined.
   # Note that by definition if 'p' and 'q' are well defined, then
   # 'p_' and 'q_' also are, so there is no need to check them.
   if (any(p < .Machine$double.eps) || any(q < .Machine$double.eps)) {
      stop('p or q undefined: check input to jahmm')
   }
   p <- scale(p, center=FALSE, scale=colSums(p))
   q <- scale(q, center=FALSE, scale=colSums(q))
   p_ <- scale(p_, center=FALSE, scale=colSums(p_))
   q_ <- scale(q_, center=FALSE, scale=colSums(q_))


   # Coerce input data to matrix for C calls.
   yz <- t(as.matrix(data[,-1]))


###############################################
#                 MAIN LOOP                   #
###############################################


   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   sorted <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[sorted]

   BW <- BaumWelch(m, yz, theta, alpha, Q, p, q, blocks,
      verbose=verbose, ...)

   Q <- BW$Q
   p <- BW$p
   q <- BW$q

   loglik <- BW$loglik
   index <- BW$index

   if (PSO) {
      if (verbose) cat("running PSO\n", file=stderr())
      pso_call <- .C(
         pso,
         # input #
         as.integer(m),
         as.integer(n),
         as.integer(ncol(data)-1),
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
      loglik <- pso_call[[11]]
   }

   C_call_1 <- .C(
      mnmultinom_prob,
      # input #
      as.integer(m),
      as.integer(n),
      as.integer(ncol(data)-1),
      as.integer(yz),
      # params #
      as.double(theta),
      as.double(alpha),
      as.double(p),
      as.double(q),
      # index #
      as.integer(index),
      # control #
      as.integer(5),   # No warning, log space.
      # output #
      double(m*n),     # Emission probabilities.
      # extra '.C()' arguments #
      NAOK = TRUE,
      DUP = FALSE
   )

   initialProb <- steady_state_probs(Q)

   vitC <- .C(
      block_viterbi,
      # input #
      as.integer(m),
      as.integer(length(blocks)),
      as.integer(blocks),
      as.double(log(Q)),
      as.double(log(initialProb)),
      C_call_1[[11]],
      # control #
      as.integer(1),
      # output #
      integer(n),
      NAOK = TRUE,
      DUP = FALSE
   )

   H0 <- BaumWelch(1, yz, theta, alpha, matrix(1), p_, q_,
      blocks, index, verbose=FALSE)

   # For the sake of speed, log-likelihood is computed without
   # the constant terms (the ones that do not depend on the
   # parameters). For QC score the absolute value of the
   # log-likielihood is needed so we need to add the correction.
   i_unique <- unique(index+1)
   i_table <- table(index+1)
   u_yz <- yz[,i_unique]
   # Use the index of the time series to save a lot of time.
   n_na.rm <- sum(!is.na(colSums(yz, na.rm=FALSE)))
   lg_terms <- lgamma(alpha+colSums(u_yz))-colSums(lgamma(u_yz+1))
   eps <- sum(i_table*(lg_terms), na.rm=TRUE)/n_na.rm - lgamma(alpha)

   QC <- (H0$loglik - loglik) / (H0$loglik+eps)
   OK <- QC > 0.01

   return(list(OK=OK, QC=QC, vPath=vitC[[8]], loglik=loglik,
      Q=Q, alpha=alpha, theta=theta, p=p, q=q,
      blocks=blocks, index=index))

}
