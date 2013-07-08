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

beverton.holt <- function(S,slope,Rmax) {
  return( (slope * S) / (1 + slope * S / Rmax ) )
}

hockey.stick <- function(S,slope,Rmax) {
  R = slope * S
  R[ R>Rmax] = Rmax
  return( R )
}

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
