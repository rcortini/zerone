yahmm <- function (data, PSO=TRUE,
   maxiter=500, tol=1e-5, verbose=TRUE) {

###############################################
#              OPTION PROCESSING              #
###############################################

   m <- 3
   n <- nrow(data)

   Q <- matrix(.05/(m-1), nrow=m, ncol=m)
   diag(Q) <- .95

   # EM for mixture distribution of the baseline.:
   baseline_params <- EM.mnmultinom(data[,2])

   theta <- baseline_params[2]
   alpha <- baseline_params[3]
   C1 <- baseline_params[4]/(1-baseline_params[4])
   C2 <- baseline_params[5]/(1-baseline_params[5])

   p <- matrix(NA, ncol=m, nrow=ncol(data))
   q <- matrix(NA, ncol=m, nrow=ncol(data))
   medy <- median(data[,2], na.rm=TRUE)
   # Fill in p row-wise.
   p[1,] <- medy * C1
   q[1,] <- medy * C2
   q[2,] <- p[2,] <- medy
   for (i in 3:ncol(data)) {
      if (m == 2) {
         q[i,] <- p[i,] <- quantile(data[,i],
            prob=c(.5, .9), na.rm=TRUE)
      }
      else if (m == 3) {
         q[i,] <- p[i,] <- quantile(data[,i],
            prob=c(.3, .6, .9), na.rm=TRUE)
      }
   }
   p <- scale(p, center=FALSE, scale=colSums(p))
   q <- scale(q, center=FALSE, scale=colSums(q))


   # Coerce input data to matrix for C calls.
   yz <- t(as.matrix(data[,-1]))
   nas <- is.na(colSums(yz, na.rm=FALSE))
   n_na.rm <- sum(!nas)
   yz_na.rm <- yz[,!nas]


###############################################
#                 MAIN LOOP                   #
###############################################


   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   sorted <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[sorted]

   index <- as.integer(rep(-1L, n))
   index_na.rm <- as.integer(rep(-1L, n_na.rm))

   # Particle Swarm Optimization.
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

      # PSO can scramble the order of the states.
      permute <- order(p[1,], decreasing=TRUE)
      p <- p[,permute]
      q <- q[,permute]
   }

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
         as.integer(0),
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
         DUP = FALSE,
         NAOK = TRUE
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
            as.integer(ncol(data)-1),
            as.integer(yz_na.rm),
            # params #
            as.double(theta),
            as.double(alpha),
            as.double(old.p),
            as.double(old.q),
            # index #
            as.integer(index_na.rm),
            # control #
            as.integer(2),
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

   } # for (iter in 1:maxiter)

   if (verbose) cat("\n", file=stderr())

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
      as.integer(1),
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

   return(list(vPath=vitC[[8]], loglik=C_call_2[[9]]/n, Q=Q,
      alpha=alpha, theta=theta, p=p, q=q,
      iterations=iter, blocks=blocks, index=index))

}
