\name{run.IRATE.example}
\alias{run.IRATE.example}
\title{compile and run preinstalled example IRATE-model setups}
\description{
  compile and run preinstalled example IRATE-model setups
}
\usage{
  run.IRATE.example(setup.name, run.folder, ...)
}

\arguments{
  \item{setup.name}{
Character string defining preinstalled IRATE example setups/datasets to be run. If missing or not valid, the user will be asked to choose among valid examples.
Check IRATE.examples() for valid, preinstalled IRATE examples.
}

\item{run.folder}{
Character string defining the folder in which the selected example shall be run. If missing a folder, named after the example setup file, will be created in the current working directory.
}

\item{...}{
Further arguments to be passed to \link{run.IRATE}
}

}

%\value{
%}

\seealso{To check model setup files see: \link{read.dat} and \link{make.dat}.\cr For run results see: \link{read.rep} and \link{read.par}. \cr To delete run files see: \link{clean.IRATE}. \cr For more information on admb model compiliation and run prodcudres see: \link{compile_admb} and \link{run_admb} from the \link{R2admb}-package} 
\author{
Robert K. Bauer
}


\examples{
### setup example not specified, user will be prompted to select among exisisting example setups:
# run.IRATE.example() 

### example name specified
# run.IRATE.example("JiangADCR")
}


