\name{run.IRATE}
\alias{run.IRATE}
\title{compile and run NOAA-IRATE models from within R}
\description{
  compile and run \href{http://nft.nefsc.noaa.gov/IRATE.html}{NOAA}-IRATE models from within R
}
\usage{
  run.IRATE(setup.file, safe=F, re=F, verbose=T,
            admb_errors=c("stop","warn","ignore"),
            mcmc=F, mcmc.opts=mcmc.control(),
            profile=F, extra.args="")
}

\arguments{
  \item{setup.file}{
Character string defining IRATE setup to be run. Check \link{IRATE.examples} for preinstalled IRATE setup examples and \link{make.dat} to create new or \link{remake.dat} to alter existing setup files.
}

\item{safe}{
 (logical) Compile in safe mode? (default = \emph{FALSE})
}

\item{re}{
 (logical) Compile in random effects (ADMB-RE) mode? (default = \emph{FALSE})
}

\item{verbose}{
(logical) Verbose output? (default = \emph{TRUE})
}

\item{admb_errors}{
(character) how to handle compilation/linking errors?
}

\item{mcmc}{
(logical) run post-hoc MCMC? (default = \emph{FALSE})
}

\item{mcmc.opts}{
options for MCMC run (see \link{mcmc.control})
}

\item{profile}{
(logical) Run likelihood profiles? (default = \emph{FALSE})
}

\item{extra.args}{
(character) extra (ADMB-) arguments for IRATE run
}

}


\seealso{To check model setup files see: \link{read.dat} and \link{make.dat}.\cr To read run results see: \link{read.rep} and \link{read.par}. \cr For preinstalled example runs see: \link{run.IRATE.example}.  \cr To delete run files see: \link{clean.IRATE}. \cr For more information on admb model compiliation and run prodcudres see: \link{compile_admb} and \link{run_admb} from the \link{R2admb}-package} 

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
}

\references{
\url{http://nft.nefsc.noaa.gov/IRATE.html}
}
