simAnneal <- function (data, Q, alpha=1, beta, gammas, tol=1e-5) {

###############################################
#              OPTION PROCESSING              #
###############################################

   # Number of states is hard coded.
   m = 3

   n <- nrow(data)
   r <- ncol(data)-2
   y <- data[,2]
   data[which(y > quantile(y, .9999, na.rm=TRUE)),2] <- NA

   if (!missing(Q)) {
      if (!all.equal(dim(Q), c(m,m)))
         stop("Q must be an m by m matrix")
   }
   else {
      Q <- matrix(c(.99, .005, 0, .01, .99, 0.05, 0, .005, .95), nrow=3)
   }

   if (missing(beta)) beta <- mean(data[,2], na.rm=TRUE) / alpha
   if (missing(gammas)) {
      gammas <- matrix(NA, ncol=3, nrow=r)
      for (i in 1:r) gammas[i,] <- quantile(data[,i+2],
         prob=c(.3, .6, .9), na.rm=TRUE) / (alpha * beta)
   }



###############################################
#            VARIABLE DEFINITIONS             #
###############################################

   # Coerce 'y' and 'z' for C call.
   yz_as_matrix <- t(as.matrix(data[,-1]))

   # Tabulation.
   tab.z.y <- c(tabulate(rowSums(data[,-1]) + 1L), -1)

   blockSizes <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   index <- match(unique(data[,1]), names(blockSizes))
   blockSizes <- blockSizes[index]


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


   index <- as.integer(rep(-1L, n))

   oldll <- -Inf
   oldalpha <- alpha
   oldbeta <- beta
   oldQ <- Q
   oldgammas <- gammas

   schedule <- c(seq(from=1, to=1e-3, by=-1e-3), rep(1e-3, 100))

   for (T in schedule) {

      cat(paste("temperature:", T, "\r"))

      alpha <- oldalpha * exp(rnorm(1, sd=.1))
      beta <- oldbeta * exp(rnorm(1, sd=.1))
      Q <- oldQ * exp(rnorm(9, sd=.1)); Q <- Q / rowSums(Q)
      gammas <- oldgammas * exp(rnorm(3*r, sd=.1))

      initialProb <- steady_state_probs(Q)

      C_call_1 <- .C("compute_loglik",
         # input #
         as.integer(n),
         as.integer(r),
         as.integer(yz_as_matrix),
         # params #
         as.double(initialProb),
         as.double(Q),
         as.double(alpha),
         as.double(beta),
         as.double(gammas),
         # index #
         index,
         as.integer(tab.z.y),
         # output #
         double(1),                    # Probability ratio.
         # extra '.C()' arguments #
         NAOK = TRUE,
         DUP = FALSE,
         PACKAGE = "HummingBee"
      )

      ll <- C_call_1[[11]]/n

      if (is.na(ll)) next
      if (ll > oldll || runif(1) < exp((ll - oldll)/T)) {
         oldll <- ll
         oldalpha <- alpha
         oldbeta <- beta
         oldQ <- Q
         oldgammas <- gammas
      }

   } # for (iter in 1:maxiter)

   cat("\n")

   C_call_1 <- .C("compute_pratio",
      # input #
      as.integer(n),
      as.integer(r),
      as.integer(yz_as_matrix),
      # params #
      as.double(oldalpha),
      as.double(oldbeta),
      as.double(oldgammas),
      # index #
      index,
      # output #
      double(2*n),                    # Probability ratio.
      # extra '.C()' arguments #
      NAOK = TRUE,
      DUP = FALSE,
      PACKAGE = "HummingBee"
   )

   initialProb <- steady_state_probs(Q)

   pem <- cbind(t(matrix(C_call_1[[8]], nrow=2)), rep(1,n))
   pem[is.infinite(pem)] <- .Machine$double.xmax
   pem <- pem / rowSums(pem)
   log_p <- log(t(pem))
   log_p[log_p < -320] <- -320
   log_Q <- log(Q)
   log_Q[log_Q < -320] <- -320

   vitC <- .C("block_viterbi",
      # input #
      as.integer(m),
      as.integer(length(blockSizes)),
      as.integer(blockSizes),
      log(initialProb),
      log_p,
      log_Q,
      # output #
      integer(n)
   )

   return(list(loglik=oldll, Q=oldQ, alpha=oldalpha,
      beta=oldbeta, gamma=oldgammas, vPath=vitC[[7]],
      emissionProb=pem))

}
