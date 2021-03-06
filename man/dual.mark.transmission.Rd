% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimating_marked_egg_fraction.R
\encoding{UTF-8}
\name{dual.mark.transmission}
\alias{dual.mark.transmission}
\title{Fraction of eggs marked for male and female mark transmission}
\usage{
dual.mark.transmission(p.female, p.male = p.female)
}
\arguments{
\item{p.female}{Fraction of all adult females genotyped in the source 
population}

\item{p.male}{Fraction of all adult males genotyped in the source 
population. Defaults to be equal to \code{p.female}}
}
\value{
A list with the following elements: \describe{\item{prob.matrix}{2x2 
  matrix with probabilities for producing offspring with male or female known
  or unknown parents}\item{p}{fraction of all eggs produced at source site 
  that will come from at least one genotyped 
  parent}\item{p.female.known}{Fraction of eggs with a single known female 
  parent among all eggs that have one or more known 
  parents}\item{p.male.known}{Fraction of eggs with a single known male 
  parent among all eggs that have one or more known 
  parents}\item{p.two.known.parents}{Fraction of eggs with two known parents
  among all eggs that have one or more known parents}}
}
\description{
Estimates the fraction of eggs produced at the source site that are the 
result of crossing parents, one or both of which have been genotyped. Based 
on the assumption that probability of breeding between pairs of individuals 
is completely independent of whether or not one or more of those individuals 
was genotyped.
}
\examples{
library(ConnMatTools)
data(damselfish.lods)

# Histograms of simulated LODs
l <- seq(-1,30,0.5)
h.in <- hist(damselfish.lods$in.group,breaks=l)
h.out <- hist(damselfish.lods$out.group,breaks=l)

# PDFs for marked and unmarked individuals based on simulations
d.marked <- stepfun.hist(h.in)
d.unmarked <- stepfun.hist(h.out)

# Fraction of adults genotyped at source site
p.adults <- 0.25

# prior.shape1=1 # Uniform prior
prior.shape1=0.5 # Jeffreys prior

# Fraction of eggs from one or more genotyped parents
p <- dual.mark.transmission(p.adults)$p

# PDF for relative connectivity
D <- d.rel.conn.dists.func(damselfish.lods$real.children,
                           d.unmarked,d.marked,p,
                           prior.shape1=prior.shape1)

# Estimate most probable value for relative connectivity
phi.mx <- optim.rel.conn.dists(damselfish.lods$real.children,
                                    d.unmarked,d.marked,p)$phi

# Estimate 95\% confidence interval for relative connectivity
Q <- q.rel.conn.dists.func(damselfish.lods$real.children,
                           d.unmarked,d.marked,p,
                           prior.shape1=prior.shape1)

# Plot it up
phi <- seq(0,1,0.001)
plot(phi,D(phi),type="l",
     xlim=c(0,0.1),
     main="PDF for relative connectivity",
     xlab=expression(phi),
     ylab="Probability density")

abline(v=phi.mx,col="green",lty="dashed")
abline(v=Q(c(0.025,0.975)),col="red",lty="dashed")
}
\references{
Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
  (in press) Uncertainty in empirical estimates of marine larval connectivity. 
  ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
}
\seealso{
Other connectivity estimation: 
\code{\link{d.rel.conn.beta.prior}()},
\code{\link{d.rel.conn.dists.func}()},
\code{\link{d.rel.conn.finite.settlement}()},
\code{\link{d.rel.conn.multinomial.unnorm}()},
\code{\link{d.rel.conn.multiple}()},
\code{\link{d.rel.conn.unif.prior}()},
\code{\link{optim.rel.conn.dists}()},
\code{\link{r.marked.egg.fraction}()}
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\concept{connectivity estimation}
