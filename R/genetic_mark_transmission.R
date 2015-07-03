#' Fraction of eggs marked for male and female mark transmission
#' 
#' Estimates the fraction of eggs produced at the source site that are the 
#' result of crossing parents, one or both of which have been genotyped. Based 
#' on the assumption that probability of breeding between pairs of individuals 
#' is completely independent of whether or not one or more of those individuals 
#' was genotyped.
#' 
#' @param p.female Fraction of all adult females genotyped in the source 
#'   population
#' @param p.male  Fraction of all adult males genotyped in the source 
#'   population. Defaults to be equal to \code{p.female}
#'   
#' @return A list with the following elements: \describe{\item{prob.matrix}{2x2 
#'   matrix with probabilities for producing offspring with male or female known
#'   or unknown parents}\item{p}{fraction of all eggs produced at source site 
#'   that will come from at least one genotyped 
#'   parent}\item{p.female.known}{Fraction of eggs with a single known female 
#'   parent among all eggs that have one or more known 
#'   parents}\item{p.male.known}{Fraction of eggs with a single known male 
#'   parent among all eggs that have one or more known 
#'   parents}\item{p.two.known.parents}{Fraction of eggs with two known parents
#'   among all eggs that have one or more known parents}}
#'   
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#'   
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.connectivity_estimation.distributions.R
#' @encoding UTF-8
#' @export
dual.mark.transmission <- function(p.female,p.male=p.female) {
  if (p.female>1 || p.female<0 || p.male<0 || p.male>1 )
    stop('Probabilities must lie between 0 and 1.')
  
  pf = c(1-p.female,p.female)
  pm = c(1-p.male,p.male)
  
  pp = pf %*% t(pm)
  
  rownames(pp) = paste("female",c("unknown","known"),sep=".")
  colnames(pp) = paste("male",c("unknown","known"),sep=".")
  
  p = 1 - pp["female.unknown","male.unknown"]
  
  p.male.known = pp["female.unknown","male.known"] / p
  p.female.known = pp["female.known","male.unknown"] / p
  p.two.known.parents = pp["female.known","male.known"] / p
  
  return(list(prob.matrix=pp,p=p,p.male.known=p.male.known,
              p.female.known=p.female.known,
              p.two.known.parents=p.two.known.parents))
}
