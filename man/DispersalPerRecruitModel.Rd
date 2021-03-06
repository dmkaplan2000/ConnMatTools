% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpr_model.R
\name{DispersalPerRecruitModel}
\alias{DispersalPerRecruitModel}
\title{Population dynamics model based on lifetime-egg-production}
\usage{
DispersalPerRecruitModel(
  LEP,
  conn.mat,
  recruits0,
  timesteps = 10,
  settler.recruit.func = hockeyStick,
  ...
)
}
\arguments{
\item{LEP}{a vector of lifetime-egg-production (LEP; also known as 
eggs-per-recruit (EPR)) for each site.}

\item{conn.mat}{a square connectivity matrix.  \code{dim(conn.mat) = 
rep(length(LEP),2)}}

\item{recruits0}{a vector of initial recruitment values for each site.}

\item{timesteps}{a vector of timesteps at which to record egg production, 
settlement and recruitment.}

\item{settler.recruit.func}{a function to calculate recruitment from the 
number of settlers at each site.  Defaults to \code{\link{hockeyStick}}.}

\item{\dots}{additional arguments to settler.recruit.func.  Typically
\code{Rmax} and \code{slope}.}
}
\value{
A list with the following elements: \item{eggs}{egg production for
  the timesteps in \code{timesteps}}
  
  \item{settlers}{Similar for settlement}
  
  \item{recruits}{Similar for recruitment}
}
\description{
This function implements the marine population dynamics model described in 
Kaplan et al. (2006).  This model is most appropriate for examining 
equilibrium dynamics of age-structured populations or temporal dynamics of 
semelparous populations.
}
\examples{
library(ConnMatTools)
data(chile.loco)

# Get appropriate collapse slope
# critical.FLEP=0.2 is just an example
slope <- settlerRecruitSlopeCorrection(chile.loco,critical.FLEP=0.2)

# Make the middle 20 sites a reserve
# All other sites: scorched earth
n <- dim(chile.loco)[2]
LEP <- rep(0,n)
nn <- round(n/2)-9
LEP[nn:(nn+19)] <- 1

Rmax <- 1

recruits0 <- rep(Rmax,n)

# Use DPR model
ret <- DispersalPerRecruitModel(LEP,chile.loco,recruits0,1:20,slope=slope,Rmax=Rmax,
                                settler.recruit.func=BevertonHolt)
image(1:n,1:20,ret$settlers,xlab="sites",ylab="timesteps",
      main=c("Settlement","click to proceed"))
locator(1)
plot(ret$settlers[,20],xlab="sites",ylab="equilibrium settlement",
     main="click to proceed")
locator(1)

# Same, but with a uniform Laplacian dispersal matrix and hockeyStick
cm <- laplacianConnMat(n,10,0,"circular")
ret <- DispersalPerRecruitModel(LEP,cm,recruits0,1:20,slope=1/0.35,Rmax=Rmax)
image(1:n,1:20,ret$settlers,xlab="sites",ylab="timesteps",
      main=c("Settlement","click to proceed"))
locator(1)
plot(ret$settlers[,20],xlab="sites",ylab="equilibrium settlement")
}
\references{
Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal
  per recruit: An efficient method for assessing sustainability in marine 
  reserve networks. Ecological Applications, 16: 2248-2263.
}
\seealso{
See also \code{\link{BevertonHolt}}, \code{\link{hockeyStick}}
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
