\name{remake.dat}
\alias{remake.dat}
\title{
load and reparameterize an existing IRATE-model setup file
}
\description{
loads an existing IRATE-model setup file (incl. tag release and recapture data) and allows to reparameterize and save it as a new setup file for IRATE model runs. Note that reparameterization will depend on the tagging dataset inlcuded in the old setup. (e.g. a "harvest and catch-&-release"-serup  is not feasable with a former harvest-only data set.)}
\usage{
remake.dat(old.setup.name, new.setup.name, age.dependent=T,
           only.harvested=T,...)}
%         ,HM
%         ,eM=T,pM,sM,Mb # natural mortality settings
%         ,eF=T,pF,sF,Fb # fishing mortality settings
%         ,eT=T,pT,sT,Tb # tagging mortality settings
%         ,combine.Hrr_and_Rrr=F # combine Harvest & Release reporting rates
%         ,eHrr=F,pHrr,sHrr,Hrrb # Harvest reporting rates
%         ,eRrr=F,pRrr,sRrr,Rrrb # Release reporting rates
%         ,combine.HRS_and_RRS=F # combine Harvest & Release Retention Survival
%         ,eHRS=F,pHRS,sHRS,HRSb #  Harvest Retention Survival
%         ,eRRS=F,pRRS,sRRS,RRSb #  Release Retention Survival
%         ,incomplete.mix=F,pNonMix,sNonMix,NonMixb #  Non-Mixing
%         ,combine.HS_and_RS=T # combine Harvest & Release Selectivity
%         ,eHS=T,sHS,HSb # age-depdendent Harvest Selectivity
%         ,eRS=F,sRS,RSb # age-depdendent Release Selectivity
%         ,eMS=T,aMS,sMS # age-depdendent M Selectivity
%         ,use.L.const=T # Likelihood Constant
%
%)
%}
\arguments{
\item{old.setup.name}{
  character string defining the name of the existing IRATE-model setup file (\emph{.dat}) to be loaded.
}

\item{new.setup.name}{
  character string defining the name of the IRATE-model setup file (\emph{.dat}) to be created.
}
\item{age.dependent}{
  whether setup file should be created for an age dependent model (default = \emph{TRUE}).
}
\item{only.harvested}{
  optional, whether actual releases of recaptures should be treated as harvested or a release-recapture model should be defined. (default = \emph{TRUE}).
}

\item{...}{Additional arguments, see Details
}
}
\details{


\code{HM}{\cr
  hooking-mortality value for each recapture year (default = 0.09), only required for release-recapture model, otherwise being ignored.
}

\code{eM,pM,sM,Mb}{
  \cr\emph{eM}: whether to estimate natural mortality (default = \emph{TRUE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pM}: beginning years of different natural mortality periods, 
  \cr\emph{sM}: starting values of different natural mortality periods,
  \cr\emph{Mb}: natural mortality lower and upper boundary values for model runs.
}

\code{eF,pF,sF,Fb}{
  \cr\emph{eF}: whether to estimate fishing mortality (default = \emph{TRUE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pF}: beginning years of different fishing mortality periods, 
  \cr\emph{sF}: starting values of different fishing mortality periods,
  \cr\emph{Fb}: fishing mortality lower and upper boundary values for model runs.
}

\code{eT=T,pT,sT,Tb}{
  \cr\emph{eT}: whether to estimate tagging mortality (default = \emph{TRUE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pT}: beginning years of different tagging mortality periods, 
  \cr\emph{sT}: starting values of different tagging mortality periods,
  \cr\emph{Tb}: tagging mortality lower and upper boundary values for model runs.
}

\code{combine.Hrr_and_Rrr}{
  whether to combine Harvest & Release reporting rates (default = \emph{FALSE})
}

\code{eHrr=F,pHrr,sHrr,Hrrb}{
  \cr\emph{eHrr}: whether to estimate harvest reporting rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pHrr}: beginning years of different harvest reporting periods, 
  \cr\emph{sHrr}: starting values of different harvest reporting periods,
  \cr\emph{Hrrb}: harvest reporting lower and upper boundary values for model runs.
}

\code{eRrr=F,pRrr,sRrr,Rrrb}{
  \cr\emph{eRrr}: whether to estimate release reporting rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pRrr}: beginning years of different release reporting periods, 
  \cr\emph{sRrr}: starting values of different release reporting periods,
  \cr\emph{Rrrb}: release reporting lower and upper boundary values for model runs.
}

\code{combine.HRS_and_RRS}{
  whether to combine Harvest & Release Retention Survival (default = \emph{FALSE})
}

\code{eHRS=F,pHRS,sHRS,HRSb}{
  \cr\emph{eHRS}: whether to estimate harvest-retention survival rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pHRS}: beginning years of harvest-retention survival periods, 
  \cr\emph{sHRS}: starting values of different-harvest retention survival periods,
  \cr\emph{HRSb}: lower and upper boundary values of harvest-retention survival rates for model runs.
}

\code{eRRS=F,pRRS,sRRS,RRSb}{
  \cr\emph{eRRS}: whether to estimate release-retention survival rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pRRS}: beginning years of release-retention survival periods, 
  \cr\emph{sRRS}: starting values of different release-retention survival periods,
  \cr\emph{RRSb}: lower and upper boundary values of release-retention survival rates for model runs.
}

\code{incomplete.mix=F,pNonMix,sNonMix,NonMixb}{
  \cr\emph{incomplete.mix}: whether to assume incomplete mixing of tagged and wild population (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pNonMix}: beginning years of incomplete mixing periods, 
  \cr\emph{sNonMix}: starting values of incomplete mixing periods,
  \cr\emph{NonMixb}: of release-retention survival rates lower and upper boundary values for model runs.
}

\code{combine.HS_and_RS}{
  whether to combine Harvest & Release Selectivity (default = \emph{TRUE})
}

\code{eHS=F,pHS,sHS,HSb}{
  \cr\emph{eHS}: whether to estimate age-depdendent harvest selectivity rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pHS}: beginning years of age-depdendent harvest selectivity periods, 
  \cr\emph{sHS}: starting values of age-depdendent harvest selectivity periods,
  \cr\emph{HSb}: lower and upper boundary values of age-depdendent harvest selectivity rates for model runs.
}

\code{eRS=F,pRS,sRS,RSb}{
  \cr\emph{eRS}: whether to estimate age-depdendent release selectivity rates (default = \emph{FALSE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{pRS}: beginning years of age-depdendent release selectivity periods, 
  \cr\emph{sRS}: starting values of age-depdendent release selectivity periods,
  \cr\emph{RSb}: lower and upper boundary values of age-depdendent release selectivity rates for model runs.
}

\code{eMS=F,aMS,sMS}{
  \cr\emph{eMS}: whether to estimate age-depdendent natural mortality rates (default = \emph{TRUE}; only considered if \emph{age.dependent} is set \emph{TRUE}),
  \cr\bold{if set \emph{TRUE}:}
    \cr\emph{aMS}: beginning age classes of age-dependent natural mortality rates (age classes are counted from 1 onwards irrespective of their age-values), 
  \cr\emph{sMS}: starting values of (grouped) age-depdendent natural mortality rates,
}

\code{use.L.const}{
  whether to set likelihood function constant (default = \emph{TRUE})
}

}

\seealso{\link{IRATE.examples}, \link{make.dat}, \link{read.dat}} 


\author{
Robert K. Bauer
}

\examples{

# old.setup <- "JiangADCR" # select old IRATE setup to reparameterize 
# old.setup.path <- system.file(paste0("IRATE.examples/",old.setup,".dat"), package = "IRATER")

# new.setup <- new.setup.path <- "JiangAICR_manual"
# system(paste("mkdir -p ",new.setup.path)) # create run folder for new setup
# system(paste("cp",old.setup.path, new.setup.path)) # copy old setup in new run folder

# setwd(new.setup.path)
# remake.dat(old.setup,new.setup,age.dependent = F) # reparameterize old setup
# run.IRATE(setup.new) # run new setup

}

