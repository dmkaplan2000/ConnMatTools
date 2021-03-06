% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\encoding{UTF-8}
\name{damselfish.lods}
\alias{damselfish.lods}
\title{Sample LOD score data for simulated and real parent-child pairs}
\format{A list with 3 elements: \describe{\item{in.group}{5000 maximum LOD 
  scores for simulated children from random crossing of real potential 
  parents}\item{out.group}{5000 maximum LOD scores for simulated children 
  from random crossing of artificial potential parents based on observed 
  allelic frequencies}\item{real.children}{Maximum LOD scores for 200 real
  juvenile fish}}}
\description{
This dataset contains both simulated and real 'log of the odds ratio' (LOD) 
scores for potential parent-child pairs of humbug damselfish (\emph{Dascyllus
aruanus}) from New Caledonia. Data was generated using 
\href{http://www.pierroton.inra.fr/genetics/labo/Software/Famoz/index.html}{FaMoz}.
In all cases, results are for the potential parent with the highest LOD score
for a given larval fish (child). Simulated data is based on artificial 
children generated from either real potential parent-pairs (the 'in' group) 
or artificial parents generated from observed allelic frequencies (the 'out' 
group).
}
\references{
Gerber S, Chabrier P, Kremer A (2003) FAMOZ: a software for 
  parentage analysis using dominant, codominant and uniparentally inherited 
  markers. Molecular Ecology Notes 3:479-481. 
  \href{http://dx.doi.org/10.1046/j.1471-8286.2003.00439.x}{doi:10.1046/j.1471-8286.2003.00439.x}

Kaplan et al. (submitted) Uncertainty in marine larval 
  connectivity estimation
}
\seealso{
See also \code{\link{d.rel.conn.dists.func}}
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\keyword{data}
