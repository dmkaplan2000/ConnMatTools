# Implements code in the following references:
# 
# @references Grüss, A., Kaplan, D. M., and Lett, C. 2012. Estimating local 
#   settler–recruit relationship parameters for complex spatially explicit 
#   models. Fisheries Research, 127–128: 34–39.
# @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal 
#   per recruit: An efficient method for assessing sustainability in marine 
#   reserve networks. Ecological Applications, 16: 2248–2263.
# @references White, J. W. 2010. Adapting the steepness parameter from 
#   stock-recruit curves for use in spatially explicit models. Fisheries 
#   Research, 102: 330–334.
# @references Grüss A, Kaplan DM, Hart DR (2011) Relative Impacts of Adult
#   Movement, Larval Dispersal and Harvester Movement on the Effectiveness of
#   Reserve Networks. PLoS ONE 6:e19960
# @references Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish
#   populations. H.M.S.O., London. 533 pp.
NULL

#' Correction slope of settler-recruit relationship
#' 
#' This function corrects the slope of the settler-recruit curve so that the 
#' collapse point of the spatially-explicit population model corresponding to 
#' the connectivity matrix agrees with that of the global non-spatially-explicit
#' model.  Uses the method in White (2010).
#' 
#' @param slope slope at the origin of the settler-recruit relationship.  Can be
#'   a vector of length = \code{dim(conn.mat)[2]}.
#' @param conn.mat a square connectivity matrix.
#' @param natural.LEP value of lifetime-egg-production (LEP), also known as 
#'   eggs-per-recruit, in the absence of fishing.  Can be a vector of length = 
#'   \code{dim(conn.mat)[2]}.  Defaults to 1.
#' @param critical.FLEP Fraction of natural.LEP at which collapse occurs. 
#'   Defaults to 0.35.
#' @param use.arpack Boolean determining if calculation is to be done with 
#'   arpack function from the igraph package. This is much quicker for large 
#'   matrices, but requires the igraph package. Defaults to TRUE, but will use 
#'   eigen instead if igraph is not found.
#'   
#' @return The slope argument corrected so that collapse happens with LEP is
#'   critical.FLEP * natural.LEP.
#'   
#' @references White, J. W. 2010. Adapting the steepness parameter from 
#'   stock-recruit curves for use in spatially explicit models. Fisheries 
#'   Research, 102: 330–334.
#'   
#' @author David Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
settler.recruit.slope.correction <- function(slope,conn.mat,natural.LEP=1,
                                             critical.FLEP=0.35,use.arpack=TRUE) {
  if (class(conn.mat) != "matrix")
    stop("Input conn.mat must be a matrix.")
  
  C = diag(slope) %*% conn.mat %*% (natural.LEP * critical.FLEP)
  
  if (use.arpack && require(igraph)) {
    ii = igraph.arpack.default
    ii$n = dim(conn.mat)[2]
    ff = function(x,extra=NULL) { conn.mat%*%x }
    v = arpack(ff, extra=NULL, sym=FALSE, options=ii)$values[1]
  } else {
    v = eigen(C)$values[1]
  }
  return(slope/Mod(v))
}

#' Beverton-Holt settler-recruit relationship
#' 
#' Calculates recruitment based on the settler-recruit relationship from 
#' Beverton & Holt (1957): \code{ slope * settlers / (1+slope*settlers/Rmax) }
#' 
#' \code{slope} and \code{Rmax} can both either be scalars or vectors of the
#' same length as \code{S}.
#' 
#' @param S a vector of settlement values, 1 for each site.
#' @param slope slope at the origin of the settler-recruit relationship.  Can be
#'   a vector of same length as \code{S}.
#' @param Rmax maximum recruitment value.
#'   
#' @return A vector of recruitment values.
#'   
#' @references Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish 
#'   populations. H.M.S.O., London. 533 pp.
#'   
#' @author David Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
beverton.holt <- function(S,slope,Rmax) {
  return( (slope * S) / (1 + slope * S / Rmax ) )
}

#' Hockey-stick settler-recruit relationship
#' 
#' Calculates recruitment based on a settler-recruit relationship that increases
#' linearly until it reaches a maximum values.
#' 
#' \code{slope} and \code{Rmax} can both either be scalars or vectors of the 
#' same length as \code{S}.
#' 
#' @param S a vector of settlement values, 1 for each site.
#' @param slope slope at the origin of the settler-recruit relationship.  Can be
#'   a vector of same length as \code{S}.
#' @param Rmax maximum recruitment value.
#'   
#' @return A vector of recruitment values.
#'   
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248–2263.
#'   
#' @author David Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
hockey.stick <- function(S,slope,Rmax) {
  R = slope * S
  R[ R>Rmax] = Rmax
  return( R )
}

#' Population dynamics model based on lifetime-egg-production
#' 
#' This function implements the marine population dynamics model described in
#' Kaplan et al. (2006).  This model is most appropriate for examining
#' equilibrium dynamics of age-structured populations or temporal dynamics of
#' semelparous populations.
#' 
#' @param LEP a vector of lifetime-egg-production (LEP; also known as
#'   eggs-per-recruit (EPR)) for each site.
#' @param conn.mat a square connectivity matrix.  \code{dim(conn.mat) =
#'   rep(length(LEP),2)}
#' @param recruits0 a vector of initial recruitment values for each site.
#' @param timesteps a vector of timesteps at which to record egg production,
#'   settlement and recruitment.
#' @param settler.recruit.func a function to calculate recruitment from the
#'   number of settlers at each site.  Defaults to \code{\link{beverton.holt}}.
#' @param \dots additional arguments to settler.recruit.func.
#'   
#' @return A list with the following elements:
#' \item{eggs}{egg production for the timesteps in \code{timesteps}}
#'
#' \item{settlers}{Similar for settlement}
#'   
#' \item{recruits}{Similar for recruitment}
#'   
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248–2263.
#'   
#' @author David Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
dispersal.per.recruit.model <- 
  function(LEP,conn.mat,
           recruits0=matrix(tryCatch(Rmax,finally=1),nrow=dim(conn.mat)[2],ncol=1),
           timesteps=10, settler.recruit.func=beverton.holt,...) {
    
    if (class(conn.mat) != "matrix")
      stop("Input conn.mat must be a matrix.")
    
    r = recruits0
    
    eggs = matrix(NA,nrow=dim(conn.mat)[2],ncol=length(timesteps))
    settlers = eggs
    recruits = settlers
    
    for (t in 1:max(timesteps)) {
      e = r * LEP
      s = conn.mat %*% e
      r = settler.recruit.func(s,...)
      
      I = which( t == timesteps )
      if (length(I)>0) {
        eggs[,I] = e
        settlers[,I] = s
        recruits[,I] = r        
      }
    }
    
    return(list(eggs=eggs,settlers=settlers,recruits=recruits))
  }
