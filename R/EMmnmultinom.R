EM.mnmultinom <- function(x, a, p, q, theta, tol=1e-6, maxiter=100) {
   # Does maximum likelihood estimation of negative binomial
   # parameters for sample 'x'.

   # Ignore NAs.
   n <- sum(!is.na(x))
   m <- mean(x, na.rm=TRUE)
   tab <- tabulate(x+1L)
   u <- 0:(length(tab)-1L)

   if (missing(a)) {
      a <- 1
      new.a <- a + 2*tol

      while (abs(new.a - a) > tol) {
         a <- new.a
         num <- sum(tab*digamma(a+u)) - n*digamma(a) + n*log(a/(m+a))
         denom <- sum(tab*trigamma(a+u)) - n*trigamma(a) + n*m/(a*(m+a))
         new.a <- a - num/denom
      }
   }

   # We have the initial value for a. Now give initial values to
   # the other parameters. The values of 'p' and 'q' are arbitrary.
   # Hands on experience shows that the algorithm may converge to
   # a local optimum if initial values are taken too close to each
   # other.
   if (missing(theta)) theta <- 0.5
   if (missing(p)) p <- .1
   if (missing(q)) q <- .5

   old_params <- c(theta, a, p, q)

   # Start the EM cycles.
   for (iter in 1:maxiter) {
      # Compute the weights in log space for numeric stability.
      # This is equivalent to w = L1 / (L1 + L2), where
      #    L1 = theta * p^a * (1-p)^u and
      #    L2 = (1-theta) * q^a * (1-q)^u,
      # avoiding the undefined case 0/0.
      logterm <- log(1-theta)-log(theta) + a*(log(q)-log(p)) +
         u*(log(1-q)-log(1-p))
      w <- 1 / (1 + exp(logterm))

      x1 <- sum(w*u*tab) / sum(w*tab)
      x2 <- sum((1-w)*u*tab) / sum((1-w)*tab)

      new.a <- a + 2*tol

      # Solve equation in alpha (a) numerically with
      # Newton-Raphson method.
      while (abs(new.a - a) > tol) {
         a <- new.a
         num <- -n*digamma(a) + 
            sum(tab*(digamma(a+u) - w*log(1+x1/a)-(1-w)*log(1+x2/a)))
         denom <- -n*trigamma(a) +
            sum(tab*(trigamma(a+u) + w*x1/(a*(x1+a))+(1-w)*x2/(a*(x2+a))))
         new.a <- a - num/denom
      }

      a <- new.a
      theta <- sum(w*tab) / n
      p <- a / (x1+a)
      q <- a / (x2+a)

      if (all(abs(c(theta, a, p, q) - old_params) < tol)) break
      old_params <- c(theta, a, p, q)

   }

   # Swap parameters so that 'p' is always smaller than 'q'.
   if (p > q) {
      tmp <- p; p <- q; q <- tmp;
      theta <- 1-theta
   }

   # The general log-likelihood equation is given by
   # 'sum(lgamma(x+a) - lgamma(a) - lgamma(x+1) + log(L1+L2))', where
   #    L1 = theta * p^a * (1-p)^x, and
   #    L2 = (1-theta) * q^a * (1-q)^x
   # We rewrite the last term of the equation as
   #    log(L1+L2) = log(1+L2/L1) + log(L1)
   # Because p < q, the term 'L2/L1' is bounded, which prevents
   # overflow. Because 'L1' is a product, the term 'log(L1)' will
   # not underflow.
   logL1L2 <- log(1 + (1-theta)/theta * (q/p)^a * ((1-q)/(1-p))^u) + 
      log(theta) + a*log(p) + u*log(1-p)
   loglik <- sum(tab*lgamma(a+u))/n - lgamma(a) -
      sum(tab*lgamma(u+1))/n + sum(tab*(logL1L2))/n

   return(c(loglik=loglik, theta=theta, a=a, p=p, q=q))

}
