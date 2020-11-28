\name{read.dat}
\alias{read.dat}
\title{
loads existing IRATE-model setup (incl. tag release & recapture data)
}
\description{
loads an existing IRATE-model setup (incl. tag release and recapture data) as a list}
\usage{
read.dat(setup.file.name,skip=0)
}
\arguments{
\item{setup.file.name}{
  character string defining the name of the existing IRATE-model input file (\emph{.dat}) to be loaded.
}

  \item{skip}{
Number of lines to skip in setup file (default = 0).
}

}

\seealso{\link{IRATE.examples}, \link{make.dat}, \link{remake.dat}} 


\author{
Robert K. Bauer
}

\examples{
### read a preinstalled IRATE setup example:
examples.folder <- system.file("IRATE.examples", package = "IRATER")
setwd(examples.folder)                       
fnames <- Sys.glob('*.dat')
read.dat(fnames[1])

}

