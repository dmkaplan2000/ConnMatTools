#' Frequency of the number of marked settlers for finite total settlement
#' 
#' This function calculates the relative frequency with which a number of marked
#' settlers would be observed as a function of relative connectivity, sample 
#' size and the total number of settlers to a site. Probable frequency of 
#' observation is calculated by randomly drawing \code{n.bootstraps} artificial 
#' samples and calculating for each the number of would-be marked settlers in 
#' the sample for all values of \code{n.conns}, the true number of settlers from
#' the marked site in the entire cohort. \code{n.conns} is related to relative 
#' connectivity, \code{phi}, by \code{phi=n.conns/n.settlers}.
#' 
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param n.obs Total number of settlers collected
#' @param n.settlers Total number of settlers at the destination site from which
#'   the \code{n.obs} (\code{<=n.settlers}) settlers are collected
#' @param n.bootstraps Number of simulated samples to generate to calculate 
#'   relative frequencies. Defaults to \code{1000}.
#' @param n.conns Vector of integers of possible numbers of settlers in the 
#'   cohort that originated at the site of marking. All values should be 
#'   integers \code{<=n.settlers}.  Defaults to \code{0:n.settlers}
#'   
#' @return A list with the following elements: \describe{\item{res}{matrix with 
#'   \code{n.obs+1} rows and \code{length(n.conns)} columns, the element of 
#'   which are counts indicating the number of times among the 
#'   \code{n.bootstraps} samples that \code{k} marked settlers are found among 
#'   the \code{n.obs} settlers collected from the total settlement pool of 
#'   \code{n.settlers} individuals. Results for \code{k} marked settlers are in 
#'   row \code{k+1} of the matrix, the first row being used for \code{k=0}.  The
#'   row names of the matrix reflect the value of \code{k}, whereas columns 
#'   names reflect the elements of \code{n.conns}.}\item{phi}{relative
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
d.rel.conn.finite.settlement <- function(p,n.obs,n.settlers,n.bootstraps=1000,
                                         n.conns=0:n.settlers) {
  res = matrix(0,nrow=n.obs+1,ncol=length(n.conns))
  rownames(res) = 0:n.obs
  colnames(res) = n.conns
  
  for (k in 1:n.bootstraps) {
    x = sort(sample(n.settlers,n.obs))
    y = as.numeric(runif(n.obs)<=p)
    i=j=0
    for (l in 1:length(n.conns)) {
      if (i < n.obs) {
        for (i in (i+1):n.obs) {
          if (x[i]>n.conns[l]) {
            i = i-1
            break            
          }
          j = j + y[i]
        }
      }
      res[j+1,l] = res[j+1,l]+1
    }
    
  }
  
  return(list(res=res,n.conns=n.conns,phi=n.conns/n.settlers))
}
