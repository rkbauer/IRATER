\name{clean.admb}
\alias{clean.admb}
\alias{clean.IRATE}
\title{delete IRATE (ADMB) run files}
\description{
  delete IRATE (ADMB) run files (not model setup files).
}
\usage{
  clean.admb(setup.name)
  clean.IRATE(setup.name)
}

\arguments{
  \item{setup.name}{
Character string defining preinstalled IRATE example setups/datasets to be run. If missing or not valid, the user will be asked to choose among valid examples.
Check IRATE.examples() for valid, preinstalled IRATE examples.
}

}

%\value{
%}

\seealso{\link{compile_admb}, \link{run_admb}, \link{IRATE.examples}} 

\author{
Robert K. Bauer
}


\examples{
example.setup <- IRATE.examples()
new.setup <- new.setup.path <- example.setup[1] # select old IRATE setup to reparameterize
print(new.setup) # print setup name to be run
old.setup.path <- system.file(paste0("IRATE.examples/",new.setup,".dat"), package = "IRATER")

system(paste("mkdir -p ",new.setup.path)) # create run folder for new setup
system(paste("cp",old.setup.path, new.setup.path)) # copy old setup in new run folder
setwd(new.setup.path)

run.IRATE(new.setup) # compile and run setup
clean.IRATE(new.setup) # deletes run files (leaving only the setup file for reruns)
}


