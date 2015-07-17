#' Frequency of the number of marked settlers for finite total settlement
#' 
#' This function calculates the relative frequency with which a number of marked
#' settlers would be observed as a function of relative connectivity, sample 
#' size and the total number of settlers to a site. Probable frequency of 
#' observation is calculated by randomly drawing \code{n.bootstraps} artificial 
#' samples and calculating for each the number of would-be marked settlers for 
#' all values of \code{phi}, the relative connectivity value.
#' 
#' @param phi Vector of fractions of settlers at the destination population that
#'   originated at the source population
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param n.obs Total number of settlers collected
#' @param n.settlers Total number of settlers at the destination site from which
#'   the \code{n.obs} (\code{<=n.settlers}) settlers are collected
#' @param n.bootstraps Number of simulated samples to generate to calculate 
#'   relative frequencies. Defaults to \code{1000}.
#'   
#' @return A matrix with \code{n.obs+1} rows and \code{length(phi)} columns, the
#'   element of which are counts indicating the number of times among the 
#'   \code{n.bootstraps} samples that \code{k} marked settlers are found among 
#'   the \code{n.obs} settlers collected from the total settlement pool of 
#'   \code{n.settlers} individuals. Results for \code{k} marked settlers are in 
#'   row \code{k+1} of the matrix, the first row being used for \code{k=0}.  The
#'   row names of the matrix reflect the value of \code{k}.
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
d.rel.conn.finite.settlement <- function(phi,p,n.obs,n.settlers,n.bootstraps=1000) {
  ss = phi * n.settlers
  
  res = matrix(0,nrow=n.obs+1,ncol=length(ss))
  row.names(res) = 0:n.obs
  
  for (k in 1:n.bootstraps) {
    x = sample(n.settlers,n.obs)
    y = runif(n.obs)<=p
    for (l in 1:length(ss)) {
      nn = sum(x<=ss[l] & y)
      res[nn+1,l] = res[nn+1,l]+1
    }
    
  }
  
  return(res)
}
