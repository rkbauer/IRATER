\name{read.par}
\alias{read.par}
\title{Read in parameters from an IRATE run}
\description{
  Reads in par-file from an IRATE run
}
\usage{
  read.par(fn)
}

\arguments{
  \item{fn}{
Character string indicating the name of the \code{.par}-file to be read.
}
}

\value{
a list containing the parameters from an IRATE run:

\describe{
   \item{\code{napram}}{number of parameters estimated}
   \item{\code{ofv}}{objective function value}
   \item{\code{max_grad_comp}}{maximum gradient component}
   \item{\code{e_F}}{fishing mortality estimates}
   \item{\code{e_M}}{natural mortality estimates}
   \item{\code{e_FA}}{tag mortality}
   \item{\code{LRR1}, \code{LRR2}}{reporting rate}
   \item{\code{PHI1R}, \code{PHI2R}}{phi}
   \item{\code{e_nmixh}, \code{e_nmixr}}{non-mixing}
   \item{\code{sel1}, \code{sel2}}{catch-at-age selectivity}
}



}

\seealso{\link{read.rep}} 

\author{
Robert K. Bauer
}
