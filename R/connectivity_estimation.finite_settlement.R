#' Frequency of the number of marked settlers for finite total settlement
#' 
#' This function calculates the relative frequency with which a number of marked
#' settlers would be observed as a function of relative connectivity, sample 
#' size and the total number of settlers to a site. Probable frequency of 
#' observation is calculated by randomly drawing \code{n.bootstraps} artificial 
#' samples and calculating for each the number of would-be marked settlers in 
#' the sample for all values of \code{n.origin}, the true number of settlers from
#' the marked site in the entire cohort. \code{n.origin} is related to relative 
#' connectivity, \code{phi}, by \code{phi=n.origin/n.settlers}.
#' 
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param n.obs Total number of settlers collected
#' @param n.settlers Total number of settlers at the destination site from which
#'   the \code{n.obs} (\code{<=n.settlers}) settlers are collected
#' @param n.bootstraps Number of simulated samples to generate to calculate 
#'   relative frequencies. Defaults to \code{1000}.
#' @param n.origin Vector of integers of possible numbers of settlers in the 
#'   cohort that originated at the site of marking. All values should be 
#'   integers \code{<=n.settlers}.  Defaults to \code{0:n.settlers}
#'   
#' @return A list with the following elements: \describe{\item{res}{matrix with 
#'   \code{n.obs+1} rows and \code{length(n.origin)} columns, the element of 
#'   which are counts indicating the number of times among the 
#'   \code{n.bootstraps} samples that \code{k} marked settlers are found among 
#'   the \code{n.obs} settlers collected from the total settlement pool of 
#'   \code{n.settlers} individuals. Results for \code{k} marked settlers are in 
#'   row \code{k+1} of the matrix, the first row being used for \code{k=0}.  The
#'   row names of the matrix reflect the value of \code{k}, whereas columns 
#'   names reflect the elements of \code{n.origin}.}\item{phi}{relative
#'   connectivity values corresponding to each column of \code{res}}}
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.connectivity_estimation.R
#' @export
#' @importFrom stats runif
d.rel.conn.finite.settlement <- function(p,n.obs,n.settlers,n.bootstraps=1000,
                                         n.origin=0:n.settlers) {
  res = matrix(0,nrow=n.obs+1,ncol=length(n.origin))
  rownames(res) = 0:n.obs
  colnames(res) = n.origin
  
  for (k in 1:n.bootstraps) {
    x = sort(sample(n.settlers,n.obs))
    y = as.numeric(runif(n.obs)<=p)
    i=j=0
    for (l in 1:length(n.origin)) {
      if (i < n.obs) {
        for (i in (i+1):n.obs) {
          if (x[i]>n.origin[l]) {
            i = i-1
            break            
          }
          j = j + y[i]
        }
      }
      res[j+1,l] = res[j+1,l]+1
    }
    
  }
  
  return(list(res=res,n.origin=n.origin,phi=n.origin/n.settlers))
}



################################################
# Alternative using Hypergeometric distribution to calculate probability
################################################
# Works, but possibly can be sped up significantly by pre-calculating all values
# of dhyper(k,K,n.settlers-K,n.obs)
# and then precalculating all values of f(k:n.settlers)
#
# Also, use isTRUE(all.equal(n.conns,round(n.conns))) to test for all integral values
d.rel.conn.finite.settlement2 <- function(n.conns,p,k,n.obs,n.settlers,
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  f = function(N) {
    r = sapply(N,function(NN) sum(dbinom(k:NN,NN,p)*dhyper(k,k:NN,n.settlers-(k:NN),n.obs) ))
    return(r)
  }
  
  f(n.conns) / sum(f(k:n.settlers))
}
