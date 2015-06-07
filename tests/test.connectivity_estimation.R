k <- 10 # Number of marked settlers among sample
n.obs <- 87 # Number of settlers in sample
n.settlers <- 100 # Total size of settler pool

p <- 0.4 # Fraction of eggs that was marked
phi <- seq(0,1,length.out=101) # Values for relative connectivity

drc <- d.relative.connectivity(phi,p,k,n.obs)
prc <- p.relative.connectivity(phi,p,k,n.obs)
qrc <- q.relative.connectivity(c(0.025,0.975),p,k,n.obs) # 95% confidence interval

n.bootstraps <- 1000 # Increase to improve accuracy of brute force estimations

# Test with brute force finite settlement function
dbf <- d.rel.conn.finite.settlement(phi,p,n.obs,100*n.obs,
                                    n.bootstraps=n.bootstraps)
dbf <- dbf[as.character(k),]

# Normalize and calculate cumulative probability function
pbf <- cumsum(c(0,(dbf[1:length(dbf)-1]+dbf[2:length(dbf)])/2*diff(phi)))
dbf <- dbf / max(pbf)
pbf <- pbf / max(pbf)

# Approximate quantiles
qbf <- approx(pbf,phi,c(0.025,0.975))$y

# Finite recruitment
dfr <- d.rel.conn.finite.settlement(phi,p,n.obs,n.settlers,
                                    n.bootstraps=n.bootstraps)
dfr <- dfr[as.character(k),]

# Normalize and calculate cumulative probability function
pfr <- cumsum(c(0,(dfr[1:length(dfr)-1]+dfr[2:length(dfr)])/2*diff(phi)))
dfr <- dfr / max(pfr)
pfr <- pfr / max(pfr)

# Approximate quantiles
qfr <- approx(pfr,phi,c(0.025,0.975))$y

# Make a plot of different distributions
plot(phi,drc,type="l",main="Probability of relative connectivity values",
     xlab=expression(phi),ylab="Probability density")
lines(phi,prc,col="blue")
lines(phi,dbf,col="black",lty="dashed")
lines(phi,dfr,col="red",lty="dashed")
abline(v=qrc,col="black")
abline(v=qbf,col="black",lty="dashed")
abline(v=qfr,col="red",lty="dashed")
