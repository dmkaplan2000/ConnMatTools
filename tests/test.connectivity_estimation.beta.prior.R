library(ConnMatTools)

k <- 10 # Number of marked settlers among sample
n.obs <- 87 # Number of settlers in sample
n.settlers <- 100 # Total size of settler pool

p <- 0.4 # Fraction of eggs that was marked
phi <- seq(0.001,1-0.001,length.out=101) # Values for relative connectivity

# Probability distribution assuming infinite settler pool and uniform prior
drc <- d.rel.conn.unif.prior(phi,p,k,n.obs)
qrc <- q.rel.conn.unif.prior(c(0.025,0.975),p,k,n.obs) # 95% confidence interval

# Probability distribution assuming infinite settler pool and using reference/Jeffreys prior
drp <- d.rel.conn.beta.prior(phi,p,k,n.obs)
prp <- p.rel.conn.beta.prior(phi,p,k,n.obs)
qrp <- q.rel.conn.beta.prior(c(0.025,0.975),p,k,n.obs) # 95% confidence interval

n.bootstraps <- 1000 # Increase to improve accuracy of brute force estimations

# Finite settler pool using Jeffreys/Reference prior
dfr <- d.rel.conn.finite.settlement(phi,p,n.obs,n.settlers,
                                    n.bootstraps=n.bootstraps)
dfr <- dfr[as.character(k),]

# Weight by prior
dfr <- dfr * dbeta(phi,0.5,0.5)

# Normalize and calculate cumulative probability function
pfr <- cumsum(c(0,(dfr[1:length(dfr)-1]+dfr[2:length(dfr)])/2*diff(phi)))
dfr <- dfr / max(pfr)
pfr <- pfr / max(pfr)

# Approximate quantiles
qfr <- approx(pfr,phi,c(0.025,0.975))$y

# Make a plot of different distributions
plot(phi,drp,type="l",main="Probability of relative connectivity values",
     xlab=expression(phi),ylab="Probability density")
lines(phi,drc,col="blue")
lines(phi,dfr,col="red",lty="dashed")
abline(v=qrp,col="black")
abline(v=qrc,col="blue")
abline(v=qfr,col="red",lty="dashed")
