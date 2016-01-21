# Functions for doing connectivity estimation with multiple sites using multinomial distributions

# Function for dirichlet distribution. Adapted from the ddirichlet function in gtools
ddirichlet <- function (x, alpha, log=FALSE) 
{
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- (alpha - 1) * log(x)
    s <- ifelse(alpha == 1 & x == 0, -Inf, s)
    return(sum(s) - logD)
  }
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    else x <- t(x)
  if (!is.matrix(alpha)) 
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                    byrow = TRUE)
  if (any(dim(x) != dim(alpha))) 
    stop("Mismatch between dimensions of 'x' and 'alpha'.")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
                                                         ])
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- -Inf
  pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- -Inf
  
  if (!log)
    pd=exp(pd)
  
  pd
}

d.rel.conn.multinomial.unnorm <- function(phis,ps,ks,n,log=FALSE,
                                          dirichlet.prior.alphas=1/(length(phis)+1)) {
  alphas = c(ks,n-sum(ks))+1
  x = c(phis*ps,1-sum(phis*ps))
  phis = c(phis,1-sum(phis))
  
  if (length(dirichlet.prior.alphas)==1)
    dirichlet.prior.alphas=rep(dirichlet.prior.alphas,length(alphas))
  
  if (log)
    d = ddirichlet(x,alphas,log=log)+ddirichlet(phis,dirichlet.prior.alphas,log=log)
  else
    d = ddirichlet(x,alphas,log=log)*ddirichlet(phis,dirichlet.prior.alphas,log=log)
  
  return(d)
}
