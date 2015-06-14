############################################################################
# Functions for calculating a weighted average of relative connectivity
# probability functions to take into account multiple possible p,k and n values
############################################################################

.rel.conn.multiple.pars.df <- function(ps,ks,ns,ws) {
  pars = data.frame(ps=ps,ks=ks,ns=ns,ws=ws)
  
  # Normalize weights
  pars$ws = pars$ws / sum(pars$ws)  
  
  return(pars)
}

#' Functions for estimating the probability distribution of relative 
#' connectivity values as a weighted sum over possible input parameters
#' 
#' These functions calculate the probability density function 
#' (\code{d.rel.conn.multiple}), the probability distribution function (aka the 
#' cumulative distribution function; \code{p.rel.conn.multiple}) and the 
#' quantile function (\code{q.rel.conn.multiple}) for the relative (to all 
#' settlers at the destination site) connectivity value for larval transport 
#' between a source and destination site. This version allows one to input 
#' multiple possible fractions of individuals (i.e., eggs) marked at the source 
#' site and multiple possible numbers of settlers collected and marked 
#' individuals observed in the sample.  This gives one the possibility to 
#' produce ensemble averages over different input parameter values with 
#' different probabilities of being correct.
#' 
#' If \code{ps}, \code{ks}, \code{ns} and \code{weights} can be scalars or
#' vectors of the same length (or lengths divisible into that of the largest
#' input parameter).  \code{weights} are normalized to sum to 1 before being
#' used to sum probabilities from each individual set of input parameters.
#' 
#' Estimations of the probability distributions are analytic, with the exception
#' of quantile estimation.  See \code{\link{q.relative.connectivity}} for more 
#' details.
#' 
#' @param phi Vector of fractions of individuals (i.e., eggs) from the source 
#'   population settling at the destination population
#' @param q Vector of quantiles
#' @param ps Vector of fractions of individuals (i.e., eggs) marked in the source 
#'   population
#' @param ks Vector of numbers of marked settlers found in sample
#' @param ns Vector of total numbers of settlers collected
#' @param weights Vector of weights for each set of p, k and n values
#' @param N Number of points at which to estimate cumulative probability 
#'   function for reverse approximation of quantile distribution. Defaults to 
#'   \code{1000}.
#'   
#' @return Vector of probabilities or quantiles, or a function in the case of 
#'   \code{\link{q.rel.conn.multiple.func}}
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @seealso See also \code{\link{d.relative.connectivity}}, 
#'   \code{\link{d.rel.conn.finite.settlement}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @describeIn d.rel.conn.multiple Estimates quantiles for the probability 
#'   distribution function for relative connectivity between a pair of sites for
#'   multiple possible \code{p}, \code{k} and \code{n} values.
#' @include connectivity_estimation.R
#' @example tests/test.connectivity_estimation.multiple.R
#' @export
d.rel.conn.multiple <- function (phi,ps,ks,ns,weights=1) {
  pars = .rel.conn.multiple.pars.df(ps,ks,ns,weights)

  d=rep(0,length(phi))
  for (M in 1:dim(pars)[1]) {
    d = d + 
      d.relative.connectivity(phi,pars$ps[M],pars$ks[M],pars$ns[M]) *
      pars$ws[M]
  }
  
  return(d)
}

#' @describeIn d.rel.conn.multiple  Estimates the cumulative probability
#'   distribution for relative connectivity between a paire of sites for
#'   multiple possible \code{p}, \code{k} and \code{n} values.
#' @include connectivity_estimation.R
#' @export
p.rel.conn.multiple <- function (phi,ps,ks,ns,weights=1) {
  pars = .rel.conn.multiple.pars.df(ps,ks,ns,weights)
  
  d=rep(0,length(phi))
  for (M in 1:dim(pars)[1]) {
    d = d + 
      p.relative.connectivity(phi,pars$ps[M],pars$ks[M],pars$ns[M]) *
      pars$ws[M]
  }
  
  return(d)
}

#' @describeIn d.rel.conn.multiple Returns a function to estimate quantiles for
#'   the probability distribution function for relative connectivity between a
#'   pair of sites for multiple possible \code{p}, \code{k} and \code{n} values.
#' @include connectivity_estimation.R
#' @export
q.rel.conn.multiple.func <- function(ps,ks,ns,weights=1,N=1000){
  phi = seq(0,1,length.out=N)
  q = p.rel.conn.multiple(phi,ps,ks,ns,weights)
  return(approxfun(q,phi))
}

#' @describeIn d.rel.conn.multiple Estimates quantiles for the probability 
#'   distribution function for relative connectivity between a pair of sites for
#'   multiple possible \code{p}, \code{k} and \code{n} values.
#' @include connectivity_estimation.R
#' @export
q.rel.conn.multiple <-function(q,ps,ks,ns,weights=1,N=1000)
  (q.rel.conn.multiple.func(ps,ks,ns,weights,N=N))(q)