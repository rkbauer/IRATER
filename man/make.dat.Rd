\name{make.dat}
\alias{make.dat}
\title{
writes and sets up an input data set (setup file) for IRATE model runs}
\description{
writes and sets up an input data set (setup file) for IRATE model runs. Obligatory arguments are:
\code{setup.file.name}, \code{releases} and \code{recapture}
}
\usage{
make.dat(setup.file.name, releases, recaptures, age.dependent=F,only.harvested=F,HM
         ,eM=T,pM,sM,Mb # natural mortality settings
         ,eF=T,pF,sF,Fb # fishing mortality settings
         ,eT=T,pT,sT,Tb # tagging mortality settings
         ,combine.Hrr_and_Rrr=F # combine Harvest & Release reporting rates
         ,eHrr=F,pHrr,sHrr,Hrrb # Harvest reporting rates
         ,eRrr=F,pRrr,sRrr,Rrrb # Release reporting rates
         ,combine.HRS_and_RRS=F # combine Harvest & Release Retention Survival
         ,eHRS=F,pHRS,sHRS,HRSb #  Harvest Retention Survival
         ,eRRS=F,pRRS,sRRS,RRSb #  Release Retention Survival
         ,incomplete.mix=F,pNonMix,sNonMix,NonMixb #  Non-Mixing
         ,combine.HS_and_RS=T # combine Harvest & Release Selectivity
         ,eHS=T,sHS,HSb # age-depdendent Harvest Selectivity
         ,eRS=F,sRS,RSb # age-depdendent Release Selectivity
         ,eMS=T,aMS,sMS # age-depdendent M Selectivity
         ,use.L.const=T # Likelihood Constant
         %,open.dat_file=T
)
}
\arguments{
\item{setup.file.name}{
  character string defining the name of the IRATE-model setup file (\code{.dat}) to be created.
}
\item{releases}{
  raw \code{data.frame} with release dates (\code{rel.date} in "\code{Date}"-format) as well as the \code{age} (numeric) of released individuals (if age dependent).
}
\item{recaptures}{
  raw \code{data.frame} with release & recpature dates (\code{rel.date}, \code{rec.date} in "\code{Date}"-format) as well as the \code{age} (numeric) of recaptured individuals and a vector \code{harvested} to indicate whether recaptured individuals were harvested or released (if age dependent).
}
\item{age.dependent}{
  whether setup file should be created for an age dependent model (default = \code{TRUE}).
}
\item{only.harvested}{
  optional, whether actual releases of recaptures should be treated as harvested or a release-recapture model should be defined. (default = \code{TRUE}).
}
\item{HM}{
  hooking-mortality value for each recapture year (default = 0.09), only required for release-recapture model, otherwise being ignored.
}

\item{eM,pM,sM,Mb}{
  \code{eM}: whether to estimate natural mortality (default = \code{TRUE}),
    if set \code{TRUE}:
  \cr\code{pM}: beginning years of different natural mortality periods, 
  \cr\code{sM}: starting values of different natural mortality periods,
  \cr\code{Mb}: natural mortality lower and upper boundary values for model runs.
}

\item{eF,pF,sF,Fb}{
  \code{eF}: whether to estimate fishing mortality (default = \code{TRUE}),
    if set \code{TRUE}:
  \cr\code{pF}: beginning years of different fishing mortality periods, 
  \cr\code{sF}: starting values of different fishing mortality periods,
  \cr\code{Fb}: fishing mortality lower and upper boundary values for model runs.
}

\item{eT,pT,sT,Tb}{
  \code{eT}: whether to estimate tagging mortality (default = \code{TRUE}),
    if set \code{TRUE}:
  \cr\code{pT}: beginning years of different tagging mortality periods, 
  \cr\code{sT}: starting values of different tagging mortality periods,
  \cr\code{Tb}: tagging mortality lower and upper boundary values for model runs.
}

\item{combine.Hrr_and_Rrr}{
   whether to combine Harvest & Release reporting rates (default = \code{FALSE})
}

\item{eHrr,pHrr,sHrr,Hrrb}{
  \code{eHrr}: whether to estimate harvest reporting rates (default = \code{FALSE}),
    if set \code{TRUE}:
  \cr\code{pHrr}: beginning years of different harvest reporting periods, 
  \cr\code{sHrr}: starting values of different harvest reporting periods,
  \cr\code{Hrrb}: harvest reporting lower and upper boundary values for model runs.
}

\item{eRrr,pRrr,sRrr,Rrrb}{
  \code{eRrr}: whether to estimate release reporting rates (default = \code{FALSE}),
    if set \code{TRUE}: 
  \cr\code{pRrr}: beginning years of different release reporting periods, 
  \cr\code{sRrr}: starting values of different release reporting periods,
  \cr\code{Rrrb}: release reporting lower and upper boundary values for model runs.
}

\item{combine.HRS_and_RRS}{
   whether to combine Harvest & Release Retention Survival (default = \code{FALSE})
}

\item{eHRS,pHRS,sHRS,HRSb}{
  \code{eHRS}: whether to estimate harvest-retention survival rates (default = \code{FALSE}),
  if set \code{TRUE}:
  \cr\code{pHRS}: beginning years of harvest-retention survival periods, 
  \cr\code{sHRS}: starting values of different-harvest retention survival periods,
  \cr\code{HRSb}: lower and upper boundary values of harvest-retention survival rates for model runs.
}

\item{eRRS,pRRS,sRRS,RRSb}{
  \code{eRRS}: whether to estimate release-retention survival rates (default = \code{FALSE}),
    if set \code{TRUE}:
  \cr\code{pRRS}: beginning years of release-retention survival periods, 
  \cr\code{sRRS}: starting values of different release-retention survival periods,
  \cr\code{RRSb}: lower and upper boundary values of release-retention survival rates for model runs.
}

\item{incomplete.mix,pNonMix,sNonMix,NonMixb}{
  \code{incomplete.mix}: whether to assume incomplete mixing of tagged and wild population (default = \code{FALSE}),
    if set \code{TRUE}:
  \cr\code{pNonMix}: beginning years of incomplete mixing periods, 
  \cr\code{sNonMix}: starting values of incomplete mixing periods,
  \cr\code{NonMixb}: of release-retention survival rates lower and upper boundary values for model runs.
}

\item{combine.HS_and_RS}{
   whether to combine Harvest & Release Selectivity (default = \code{TRUE})
}

\item{eHS,sHS,HSb}{
  \code{eHS}: whether to estimate age-depdendent harvest selectivity rates (default = \code{FALSE}),
    if set \code{TRUE}:%\cr\code{pHS}: beginning years of age-depdendent harvest selectivity periods, 
  \cr\code{sHS}: starting values of age-depdendent harvest selectivity periods,
  \cr\code{HSb}: lower and upper boundary values of age-depdendent harvest selectivity rates for model runs.
}

\item{eRS,sRS,RSb}{
  \code{eRS}: whether to estimate age-depdendent release selectivity rates (default = \code{FALSE}),
    if set \code{TRUE}:%  \cr\code{pRS}: beginning years of age-depdendent release selectivity periods, 
  \cr\code{sRS}: starting values of age-depdendent release selectivity periods,
  \cr\code{RSb}: lower and upper boundary values of age-depdendent release selectivity rates for model runs.
}

\item{eMS,aMS,sMS}{
  \code{eMS}: whether to estimate age-depdendent natural mortality rates (default = \code{TRUE}; only considered if \code{age.dependent} is set \code{TRUE}),
    if set \code{TRUE}:
  \cr\code{aMS}: beginning age classes of age-dependent natural mortality rates (age classes are counted from 1 onwards irrespective of their age-values), 
  \cr\code{sMS}: starting values of (grouped) age-depdendent natural mortality rates,
}

\item{use.L.const}{
    whether to set likelihood function constant (default = \code{TRUE})
}

%\item{open.dat_file}{
%    whether to open generated IRATE-model setup file (default = \code{TRUE})
%}


}

\seealso{To check model setup files in detail see: \link{read.dat}. \cr To alter existent model setup files see: \link{remake.dat}. \cr To read model run results see: \link{read.rep} and \link{read.par}. \cr For preinstalled example runs and setup data see: \link{run.IRATE.example}.  \cr To delete run files see: \link{clean.IRATE}. \cr For more information on admb model compiliation and run prodcudres see: \link{compile_admb} and \link{run_admb} from the \link{R2admb}-package} 


\author{
Robert K. Bauer
}

\references{
\url{http://nft.nefsc.noaa.gov/IRATE.html}
}

%\examples{
%}

