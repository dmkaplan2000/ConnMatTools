library(ConnMatTools)

k <- 10 # Number of marked settlers among sample
n.obs <- 87 # Number of settlers in sample
n.settlers <- 100 # Total size of settler pool

p <- 0.4 # Fraction of eggs that was marked
phi <- seq(0,1,length.out=101) # Values for relative connectivity

# Probability distribution assuming infinite settler pool and uniform prior
drc <- d.rel.conn.unif.prior(phi,p,k,n.obs)
prc <- p.rel.conn.unif.prior(phi,p,k,n.obs)
qrc <- q.rel.conn.unif.prior(c(0.025,0.975),p,k,n.obs) # 95% confidence interval

n.bootstraps <- 1000 # Increase to improve accuracy of brute force estimations

# Test with brute force finite settlement function and large (approx. infinite) settler pool
dbf <- d.rel.conn.finite.settlement(p,n.obs,7*n.obs,
                                    n.bootstraps=n.bootstraps)
dbf$d <- dbf$res[as.character(k),]

# Approximate quantiles - could be slow for large cohort sizes or numbers of bootstraps
dbf$q = quantile(do.call(c,sapply(1:length(dbf$d),function(i) rep(dbf$phi[i],dbf$d[i]))),
                 c(0.025,0.975))

# Normalize probability density
dbf$d <- dbf$d / sum(dbf$d)


# Finite settler pool
dfr <- d.rel.conn.finite.settlement(p,n.obs,n.settlers,
                                    n.bootstraps=n.bootstraps)
dfr$d <- dfr$res[as.character(k),]

# Approximate quantiles - could be slow for large cohort sizes or numbers of bootstraps
dfr$q = quantile(do.call(c,sapply(1:length(dfr$d),function(i) rep(dfr$phi[i],dfr$d[i]))),
                 c(0.025,0.975))

# Normalize probability density
dfr$d <- dfr$d / sum(dfr$d)

# Make a plot of different distributions
plot(phi,drc,type="l",main="Probability of relative connectivity values",
     xlab=expression(phi),ylab="Probability density")
lines(phi,prc,col="blue")
lines(dbf$phi,dbf$d/diff(dbf$phi[1:2]),col="black",lty="dashed")
lines(dfr$phi,dfr$d/diff(dfr$phi[1:2]),col="red",lty="dashed")
abline(v=qrc,col="black")
abline(v=dbf$q,col="black",lty="dashed")
abline(v=dfr$q,col="red",lty="dashed")
