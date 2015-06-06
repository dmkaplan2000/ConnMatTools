# Code for detecting the fraction of recaptured individuals that are marked when we
# know the distributions of marked and unmarked individuals

#' Create a probability density step function from a histogram object
#' 
#' This function creates a step function from the bars in a \code{histogram} 
#' object. By default, the step function will be normalized so that it 
#' integrates to 1.
#' 
#' @param h an object of type \code{histogram}
#' @param \dots Additional arguments for the default \code{\link{stepfun}} 
#'   function.
#' @param normalize Boolean indicating whether or not to normalize the output 
#'   stepfun so that it integrates to 1. Defaults to \code{TRUE}.  If 
#'   \code{FALSE}, then the function will integrate to \code{sum(h$counts)}
#'   
#' @return A function of class \code{stepfun}.  The height of the steps will be
#'   divided by the distance between breaks and possibly the
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
stepfun.hist <- function(h,...,normalize=TRUE) {
  n = 1
  if (normalize)
    n = sum(h$counts)
  
  return(stepfun(x=h$breaks,y=c(0,h$counts / diff(h$breaks) / n,0),...))
}

#' Returns probability density function (PDF) for a mix of marked and unmarked 
#' individuals
#' 
#' This function returns a probability density function (PDF) for scores for a 
#' mix of marked and unmarked individuals with known fraction of marked 
#' individuals. The distributions for marked individuals and for unmarked 
#' individuals must be known.
#' 
#' @param d.unmarked A function representing the PDF of unmarked individuals. 
#'   Must be normalized so that it integrates to 1 for the function to work 
#'   properly.
#' @param d.marked A function representing the PDF of marked individuals.  Must 
#'   be normalized so that it integrates to 1 for the function to work properly.
#'   
#' @return A function representing the PDF of observations drawn from the mixed 
#'   distribution of marked and unmarked individuals.  The function takes two 
#'   arguments: \code{p.marked}, the fraction of marked individuals in the
#'   distribution; and \code{obs}, a vector of observed score values.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
d.mix.dists.func <- function(d.unmarked,d.marked) {
  return(function(p.marked,obs) (1-p.marked)*d.unmarked(obs)+
           p.marked*d.marked(obs))
}

#' Returns probability a set of observations correspond to marked individuals
#' 
#' This function returns the probability each of a set of observations 
#' corresponds to a marked individual given the distribution of scores for 
#' unmarked and marked individuals and the fraction of individuals that are 
#' marked.
#' 
#' @param p.marked The overall fraction of marked individuals in the entire 
#'   population
#' @param obs A vector of score values for a random sample of (marked and 
#'   unmarked) individuals from the population
#' @inheritParams d.mix.dists.func
#'   
#' @return A vector of the same size as \code{obs} containing the probability
#'   that each individual is marked
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
prob.marked <- function(p.marked,obs,d.unmarked,d.marked) {
  p.marked*d.marked(obs) / ((1-p.marked)*d.unmarked(obs)+p.marked*d.marked(obs))
}

#' Logarithm of probability for a set of observed score values given a PDF
#' 
#' This function returns log-probability for a set of observed score values 
#' given a PDF for the distribution of scores.
#' 
#' @param p Parameter in the PDF distribution.  This should be a single scalar
#'   value, presumably the fraction of marked individuals in the population, but
#'   the function is generic.
#' @param obs A vector of score values
#' @param dfunc A function representing the PDF of scores.  Should be normalized
#'   so that it integrates to 1 for the function to return a true 
#'   log-probability.
#'   
#' @return The log-probability.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
log.prob <- function(p,obs,dfunc) {
  v=dfunc(p,obs)
  
  if (any(v==0))
    return(-Inf)
  
  sum(log(v))  
}

#' Calculate the optimal fraction of marked individuals
#' 
#' This function calculates the fraction of marked individuals that best fits a 
#' set of observed score values and a pair of distributions for marked and 
#' unmarked individuals.
#' 
#' @param obs Vector of observed score values for potentially marked individuals
#' @inheritParams d.mix.dists.func
#' @param par Initial fraction of marked individuals for \code{\link{optim}} 
#'   function. Defaults to 0.5.
#' @param method Method variable for \code{\link{optim}} function. Defaults to
#'   \code{"Brent"}.
#' @param lower Lower limit for search for fraction of marked individuals. 
#'   Defaults to 0.
#' @param upper Upper limit for search for fraction of marked individuals. 
#'   Defaults to 1.
#' @param \dots Additional arguments for the \code{\link{optim}} function.
#'   
#' @return A list with results of optimization. Optimal fraction of marked
#'   individuals is in \code{par} field. See \code{\link{optim}} for details.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval
#'   connectivity estimation
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
optimal.fraction.from.distributions <- function(obs,d.unmarked,d.marked,
                                                par=0.5,method="Brent",lower=0,upper=1,
                                                ...) {
  f = function(p.marked)
    -log.prob(p.marked,obs,d.mix.dists.func(d.unmarked,d.marked))
  
  o=optim(par=par,fn=f,method=method,lower=lower,upper=upper,...)
  
  return(o)
}

.d.marked.fraction.internal <- function(obs,d.unmarked,d.marked) {
  f = d.mix.dists.func(d.unmarked,d.marked)
  ff = function(p.marked,obs)
    log.prob(p.marked,obs,f)
  
  f0 = optimal.fraction.from.distributions(obs,d.unmarked,d.marked)$value
  
  return(function(p) exp(ff(p,obs)+f0))
}

.d.marked.fraction.N <- function(obs)
  max(100,min(5000,2*length(obs)))

#' Functions for examining the probability distribution for the fraction of 
#' marked individuals
#' 
#' These functions return functions that calculates the probability density 
#' function (\code{d.marked.fraction}), the probability distribution function (aka 
#' the cumulative distribution function; \code{p.marked.fraction}) and the quantile 
#' function (\code{q.marked.fraction}) given a set of observed score values and 
#' distributions for unmarked and marked individuals.
#' 
#' The normalization of the probability distribution is carried out using a 
#' simple, fixed-step trapezoidal integration scheme.  The number of steps 
#' between marked fractions 0 and 1 defaults to \code{2*length(obs)} so long as 
#' that number is comprised between \code{100} and \code{5000}.
#' 
#' @param obs Vector of observed score values for potentially marked individuals
#' @inheritParams d.mix.dists.func
#' @param N number of steps in fixed-step, trapezoidal integration scheme used 
#'   to normalize probability distributions.  Defaults to \code{2*length(obs)}
#'   so long as that number is comprised between \code{100} and \code{5000}.
#' @param \dots Additional arguments for the \code{\link{approxfun}} function.
#'   
#' @return A function that takes one argument (the fraction of marked 
#'   individuals for \code{d.marked.fraction} and \code{p.marked.fraction}; the quantile
#'   for \code{q.marked.fraction}) and returns the probability density, cumulative 
#'   probability distribution or score value, respectively.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval
#'   connectivity estimation
#'   
#' @describeIn d.marked.fraction Returns a function that is PDF for fraction of 
#'   marked individuals
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
d.marked.fraction <- function(obs,d.unmarked,d.marked,N=.d.marked.fraction.N(obs),...) {
  p.marked=seq(0,1,length.out=N+1)
  
  f = sapply(p.marked,.d.marked.fraction.internal(obs,d.unmarked,d.marked))
  fs = sum(f[1:N]+diff(f)/2)*(p.marked[2]-p.marked[1])
  f2 = f / fs
  
  return(approxfun(p.marked,f2,...))
}

#' @describeIn d.marked.fraction Returns a function that is cumulative probability
#'   distribution for fraction of marked individuals
#' @export
p.marked.fraction <- function(obs,d.unmarked,d.marked,N=.d.marked.fraction.N(obs),...) {
  p.marked=seq(0,1,length.out=N+1)
  
  f = sapply(p.marked,.d.marked.fraction.internal(obs,d.unmarked,d.marked))
  fs = c(0,cumsum(f[1:N]+diff(f)/2)*(p.marked[2]-p.marked[1]))
  f2 = fs / fs[length(fs)]
  browser()
  return(approxfun(p.marked,f2,...))
}

#' @describeIn d.marked.fraction Returns a function that is quantile function for
#'   fraction of marked individuals
#' @export
q.marked.fraction <- function(obs,d.unmarked,d.marked,N=.d.marked.fraction.N(obs),...) {
  p.marked=seq(0,1,length.out=N+1)
  
  f = sapply(p.marked,.d.marked.fraction.internal(obs,d.unmarked,d.marked))
  fs = c(0,cumsum(f[1:N]+diff(f)/2)*(p.marked[2]-p.marked[1]))
  f2 = fs / fs[length(fs)]
  
  return(approxfun(f2,p.marked,...))
}

