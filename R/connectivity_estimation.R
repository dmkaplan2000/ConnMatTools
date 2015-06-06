################################################################
# This code estimates analytic probability distribution for a 
# single type of marked individuals among a set of collected individuals (larvae) 
# from a source population with known fraction of marked individuals (eggs).
# 
# In the code below, I first define recursion functions to estimate
# the integrals in the probability distribution.  There is one such
# for recursion to k=n, and another to k=0.  I then define a function
# that chooses which one to use based on whether or not k>n/2.
#
# The rest of the functions use this integral to calculate different
# aspects of the probability distribution.
################################################################

# Utility function to determine if input is an integer
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

###########################################
# From k=n
###########################################

# Integral needed for calculating quantiles
.sr.int.n <- function(p,k,n) {
  #   # Check for whole number inputs.  Comment out to increase efficiency a bit.
  #   if (!is.wholenumber(k) | !is.wholenumber(n))
  #     stop("k and n must be integers")
  
  # Define base of series
  if (k>n) return(0)
  
  # For other values, calculate by induction
  return(p^(k+1) * (1-p)^(n-k) / (k+1) + 
           (n-k)/(k+1) * .sr.int.n(p,k+1,n))
}

###########################################
# From k=0
###########################################

# Integral needed for calculating quantiles
# Recursion from k=0
.sr.int.0 <- function(p,k,n) {
  #   # Check for whole number inputs.  Comment out to increase efficiency a bit.
  #   if (!is.wholenumber(k) | !is.wholenumber(n))
  #     stop("k and n must be integers")
  
  # Define base of series
  if (k<1) return((1-(1-p)^(n+1))/(n+1))
  
  # For other values, calculate by induction
  return(k/(n-k+1) * .sr.int.0(p,k-1,n) -
           p^(k) * (1-p)^(n-k+1) / (n-k+1))
}

#' Evaluates integral used to for PDF calculation for estimating relative 
#' transport between a pair of sites
#' 
#' This function analytically evaluates the integral from \code{0} to \code{p} 
#' of \code{x^k*(1-x)^(n-k)}.  This integral is used to normalize the PDF 
#' functions for estimating relative (to the all recruits to the recruitment 
#' site) larval transport between a pair of sites.
#' 
#' This is a universal version of the integral function that does error checking
#' and picks one of two specific functions that actually do the calculations by 
#' iteration from \code{k} to \code{n} or to \code{0}.  The choice between the 
#' two versions is based on \code{k>n/2} or not.  Also looks for simple cases so
#' they can be evaluated efficiently.
#' 
#' @param p Vector of upper limits for integral
#' @param k Number of marked settlers found in sample
#' @param n Total number of settlers collected
#'   
#' @return A vector of size \code{length(p)} with the values of the integrals.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval
#'   connectivity estimation
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
.sr.int <- function(p,k,n) {
  # Check for whole number inputs.  Comment out to increase efficiency a bit.
  if (!is.wholenumber(k) || !is.wholenumber(n))
    stop("k and n must be integers")
  
  if (k<0 || k>n)
    stop("k must satisfy 0<=k<=n")
  
  # Do simple cases
  I = p == 0
  J = p == 1
  K = !I & !J
  x = NA*p
  x[I] = 0
  x[J] = 1/choose(n,k)/(n+1)
  
  if (any(K)) {
    if (k>n/2)
      x[K] = .sr.int.n(p[K],k,n)
    else
      x[K] = .sr.int.0(p[K],k,n)    
  }
  
  return(x)
}


################################################
# Functions for actually calculating PDF
################################################

#' Functions for estimating the probability distribution of relative 
#' connectivity values
#' 
#' These functions calculate the probability density function 
#' (\code{d.relative.connectivity}), the probability distribution function (aka 
#' the cumulative distribution function; \code{p.relative.connectivity}) and the
#' quantile function (\code{q.relative.connectivity}) for the relative (to all 
#' settlers at the destination site) connectivity value for larval transport 
#' between a source and destination site given a known fraction of marked 
#' individuals (i.e., eggs) in the source population.
#' 
#' Estimations of the probability distribution are analytic, except that
#' quantile estimation is performed using \code{\link{approxfun}} to perform
#' reverse estimation based on the cumulative probability distribution function
#' estimated at a finite number of points.
#' 
#' @param phi Vector of fractions of individuals (i.e., eggs) from the source 
#'   population settling at the destination population
#' @param q Vector of quantiles
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param k Number of marked settlers found in sample
#' @param n Total number of settlers collected
#' @param N Number of points at which to estimate cumulative probability 
#'   function for reverse approximation of quantile distribution. Defaults to 
#'   \code{1000}.
#'   
#' @return Vector of probabilities or quantiles.
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @describeIn d.relative.connectivity Returns the probability density for 
#'   relative connectivity between a paire of sites
#' @export
d.relative.connectivity <- function(phi,p,k,n)
  p*(p*phi)^k*(1-p*phi)^(n-k)/.sr.int(p,k,n)

#' @describeIn d.relative.connectivity Returns the cumulative probability
#'   distribution for relative connectivity between a paire of sites
#' @export
p.relative.connectivity <- function(phi,p,k,n)
  .sr.int(p*phi,k,n)/.sr.int(p,k,n)

#' @describeIn d.relative.connectivity Returns a function to estimate quantiles
#'   for the probability distribution function for relative connectivity between
#'   a pair of sites
#' @export
q.relative.connectivity.func <- function(p,k,n,N=1000){
  phi = seq(0,1,length.out=N)
  q = psr(phi,p,k,n)
  return(approxfun(q,phi))
}

#' @describeIn d.relative.connectivity Estimates quantiles for the probability
#'   distribution function for relative connectivity between a pair of sites
#' @export
q.relative.connectivity <-function(q,p,k,n,N=1000)
  (q.relative.connectivity.func(p,k,n,N=N))(q)
