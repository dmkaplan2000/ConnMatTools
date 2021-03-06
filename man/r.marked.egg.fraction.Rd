% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimating_marked_egg_fraction.R
\encoding{UTF-8}
\name{r.marked.egg.fraction}
\alias{r.marked.egg.fraction}
\title{Estimates of fraction of eggs marked accounting for variability in 
reproductive output}
\usage{
r.marked.egg.fraction(
  n,
  n.females,
  n.marked.females = round(n.females * p.marked.females),
  mean.female = 1,
  cv.female,
  dual = FALSE,
  male.uncert = FALSE,
  n.males = n.females,
  n.marked.males = tryCatch(round(n.males * p.marked.males), error = function(e)
    n.marked.females),
  mean.male = mean.female,
  cv.male = cv.female,
  p.marked.females,
  p.marked.males = p.marked.females
)
}
\arguments{
\item{n}{Number of random values to estimates}

\item{n.females}{Total number of mature females in the population}

\item{n.marked.females}{Number of marked females in population}

\item{mean.female}{Mean egg production of each mature female. Defaults to 1.}

\item{cv.female}{Coefficient of variation of reproductive output of an 
individual mature female}

\item{dual}{Logical variable. If \code{TRUE}, then the fraction of marked 
eggs is calculated assuming dual (male and female) mark transmission. 
Defaults to \code{FALSE}.}

\item{male.uncert}{Logical variable. If \code{TRUE}, then variability in male
sperm output is also taken into account when estimating the number of 
marked eggs. Defaults to \code{FALSE}.}

\item{n.males}{Total number of mature males in the population. Only used if 
\code{dual=TRUE}. Defaults to being equal to \code{n.females}.}

\item{n.marked.males}{Number of marked males in population. Only used if 
\code{dual=TRUE}. Defaults to being equal to \code{n.marked.females}.}

\item{mean.male}{Mean sperm production of each mature male. Only used if 
\code{dual=TRUE} and \code{male.uncert=TRUE}. Defaults to being equal to 
\code{mean.female}.}

\item{cv.male}{Coefficient of variation of reproductive output of an 
individual mature male. Only used if \code{dual=TRUE} and 
\code{male.uncert=TRUE}. Defaults to being equal to \code{cv.female}.}

\item{p.marked.females}{Fraction of marked females in population. Can be 
supplied instead of \code{n.marked.females}. Ignored if 
\code{n.marked.females} is given.}

\item{p.marked.males}{Fraction of marked males in population. Can be supplied
instead of \code{n.marked.males}. Only used if \code{dual=TRUE}. Ignored if
\code{n.marked.males} is given.}
}
\value{
A list with the following elements: \describe{\item{p}{Vector of 
  length \code{n} with estimates for fraction of marked 
  eggs}\item{eggs}{Vector of length \code{n} with estimates for total egg 
  production}\item{marked.eggs}{Vector of length \code{n} with estimates for 
  total number of marked eggs produced}\item{sperm}{Only returned if 
  \code{dual=TRUE}. If \code{male.uncert=FALSE}, then a scalar equal to 
  \code{n.males}. Otherwise, a vector of length \code{n} with estimates for 
  total sperm production}\item{marked.sperm}{Only returned if 
  \code{dual=TRUE}. If \code{male.uncert=FALSE}, then a scalar equal to 
  \code{n.marked.males}. Otherwise, a vector of length \code{n} with
  estimates for total marked sperm production}}
}
\description{
This function estimates the fraction of eggs "marked" at a site (where the 
"mark" could be micro-chemical or genetic) taking into account uncertainty in
female (and potentially male in the case of dual genetic mark transmission) 
reproductive output. It generates a set of potential values for the fraction 
of eggs marked assuming that reproductive output of each marked or unmarked 
mature individual is given by a random variable drawn from a single 
probability distribution with known mean and standard deviation (or 
equivalently coefficient of variation) \strong{and} that the numbers 
of marked and unmarked individuals are large enough that the central limit 
theorem applies and, therefore, their collective reproductive outputs are 
reasonably well described by a gamma distribution whose mean and standard 
deviation are appropriately scaled based on the number of individual 
reproducers.  The function also returns the total egg production 
corresponding to each fraction of marked eggs, needed for estimating absolute
connectivity values (i.e., elements of the connectivity matrix needed for 
assessing population persistence).
}
\examples{
library(ConnMatTools)

n.females <- 500
n.marked.females <- 100
p.marked.females <- n.marked.females/n.females
mn <- 1
cv <- 1
# Numbers of males and marked males and variance in male sperm production
# assumed the same as values for females

# Random values from distribution of pure female mark transmission
F=r.marked.egg.fraction(1000,n.females=n.females,n.marked.females=n.marked.females,
                        mean.female=mn,cv.female=cv)

# Random values from distribution of dual female-male mark transmission, but
# fraction of marked eggs only depends on fraction of marked males
Fm=r.marked.egg.fraction(1000,n.females=n.females,n.marked.females=n.marked.females,
                        mean.female=mn,cv.female=cv,dual=TRUE,male.uncert=FALSE)

# Random values from distribution of dual female-male mark transmission, with
# fraction of marked eggs depending on absolute marked and unmarked sperm output
FM=r.marked.egg.fraction(1000,n.females=n.females,n.marked.females=n.marked.females,
                         mean.female=mn,cv.female=cv,dual=TRUE,male.uncert=TRUE)

# Plot of pure female mark transmission
hist(F$p,50,main="Female mark transmission",
     xlab="Fraction of marked eggs",
     ylab="Frequency")

# Female+male mark transmission, but no variability in male mark transmission
h <- hist(Fm$p,50,main="Female+male mark transmission, no male uncert.",
          xlab="Fraction of marked eggs",
          ylab="Frequency")
hh <- hist((1-p.marked.females)*F$p + p.marked.females,
           breaks=c(-Inf,h$breaks,Inf),plot=FALSE)
lines(hh$mids,hh$counts,col="red")

# Plot of pure female mark transmission
h <- hist(FM$p,50,plot=FALSE)
hh <- hist(Fm$p,
           breaks=c(-Inf,h$breaks,Inf),plot=FALSE)
plot(h,ylim=c(0,1.1*max(hh$counts,h$counts)),
     main="Female+Male mark transmission, male uncert.",
     xlab="Fraction of marked eggs",
     ylab="Frequency")
lines(hh$mids,hh$counts,col="red")
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
\code{\link{dual.mark.transmission}()},
\code{\link{optim.rel.conn.dists}()}
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\concept{connectivity estimation}
