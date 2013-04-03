HummingBee <- function (data, ...) {

   n <- nrow(data)
   r <- ncol(data)-2
   yz_as_matrix <- t(as.matrix(data[,-1]))

   simAnneal_call <- .C("simAnneal",
      # input #
      as.integer(n),
      as.integer(r),
      as.integer(yz_as_matrix),
      # output #
      double(12+3*r),
      NAOK = TRUE,
      DUP = FALSE,
      PACKAGE = "HummingBee"
   )

   # First fit with simulated annealing.
   SAfit <- simAnneal_call[[4]]
   loglik <- SAfit[1]
   Q <- matrix(SAfit[2:10], nrow=3, ncol=3)
   alpha <- SAfit[11]
   beta <- SAfit[12]
   gammas <- matrix(SAfit[13:(12+3*r)], ncol=3)

   # Compute emission probabilities (needed for Viterbi algorithm).
   compute_pratio_call <- .C("compute_pratio",
      # input #
      as.integer(n),
      as.integer(r),
      as.integer(yz_as_matrix),
      # params #
      as.double(alpha),
      as.double(beta),
      as.double(gammas),
      # index #
      as.integer(rep(-1L, n)),
      # output #
      double(2*n),                    # Probability ratio.
      # extra '.C()' arguments #
      NAOK = TRUE,
      DUP = FALSE,
      PACKAGE = "HummingBee"
   )

   # Precompute the log of emission and transition probabilities
   # for input of the Viterbi algorithm.
   pem <- cbind(t(matrix(compute_pratio_call[[8]], nrow=2)), rep(1,n))
   pem[is.infinite(pem)] <- .Machine$double.xmax
   pem <- pem / rowSums(pem)
   log_init <- log(rep(1/3,3))
   log_Q <- log(Q)
   log_p <- log(t(pem))

   # Set -320 as the lowest value (corresponds to the smallest positive
   # real on many machines.
   log_p[log_p < -320] <- -320
   log_Q[log_Q < -320] <- -320

   blocks <- tapply(X=rep(1,n), INDEX=as.character(data[,1]), FUN=sum)
   id <- match(unique(data[,1]), names(blocks))
   blocks <- blocks[id]

   # Viterbi algorithm.
   block_viterbi_call <- .C("block_viterbi",
      # input #
      as.integer(3),
      as.integer(length(blocks)),
      as.integer(blocks),
      log_init,
      log_p,
      log_Q,
      # output #
      integer(n)
   )

   return(list(loglik=loglik, Q=Q, alpha=alpha, beta=beta, gamma=gammas,
      vPath=block_viterbi_call[[7]], emissionProb=pem, blocks=blocks))

}
