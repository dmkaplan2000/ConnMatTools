#' Local retention of a connectivity matrix
#' 
#' Local retention is defined as the diagonal elements of the connectivity matrix.
#' 
#' @param conn.mat A square connectivity matrix.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
localRetention <- function(conn.mat)
  diag(conn.mat)

#' Relative local retention of a connectivity matrix
#' 
#' Relative local retention is defined as the diagonal elements of the
#' connectivity matrix divided by the sum of the corresponding column of the
#' connectivity matrix.
#' 
#' @param conn.mat A square connectivity matrix.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
relativeLocalRetention <- function(conn.mat) {
  ss = apply(conn.mat,2,sum)
  ss[ ss == 0 ] = 1
  localRetention(conn.mat) / ss  
}

#' Self recruitment of a connectivity matrix
#' 
#' Self recruitment is defined as the diagonal elements of the
#' connectivity matrix divided by the sum of the corresponding row of the
#' connectivity matrix.
#' 
#' @param conn.mat A square connectivity matrix.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
selfRecruitment <- function(conn.mat) {
  ss = apply(conn.mat,1,sum)
  ss[ ss == 0 ] = 1
  localRetention(conn.mat) / ss  
}

