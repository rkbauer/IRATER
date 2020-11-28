\name{IRATER}
\alias{IRATER}
\docType{package}
\title{
IRATER - A R Interface for the Instantaneous RATEs (IRATE) Model
}
\description{
\href{http://nft.nefsc.noaa.gov/IRATE.html}{IRATE} is a software tool, developed by Dr. Gary Nelson, to analyze band recovery (conventional tagging) data. It can estimate fishing and natural mortality rates based on likewise estimated exploitation and survival rates of (continuous) fisheries. IRATE can handle both age-independent and age-dependent instantaneous rates models (Hoenig et al., 1998; Jiang et al., 2007) to be fitted to multi-year fish tag return data. IRATE allows model development with either age-dependent harvest-only or harvest and catch-release tag returns or similar age independent models. It also allows estimation of harvest reporting rates, catch and release reporting rates, and tag retention of harvested and/or released fish.  However, not all parameters in the model can be estimated simultaneously with tag data alone. Some parameters must be fixed and assumed known (usually reporting rate and tag loss) to obtain good estimates of remaining parameters.  Additionally, the model can account for non-mixing of the tagged fish in the first release year and adjust for harvest and M selectivity in the age-based models.   The negative log likelihood is used as the objective function to obtain maximum likelihood estimates of parameters.  Several model fit statistics are provided that can be used to select the best model formulation; these include the Akaike Information Criterion (AIC), c-hat (a measure of overdispersion) and standard residuals.
\cr \code{Source}: \url{http://nft.nefsc.noaa.gov/IRATE.html}

\code{IRATER} is an interface to setup IRATE model runs and execute them from within R, using the ADMB engine of IRATE.
It further provides functions to read and access the model results back into R.

\bold{Major functions:}\cr
\itemize{
\item{\link{make.dat}: writes and sets up an input data set (setup file) for IRATE model runs}
\item{\link{run.IRATE}: compile and run IRATE model from within R}
\item{\link{read.par} & \link{read.rep}: read run results}

}

\bold{Getting Started}\cr
Check out the help files of the principle functions, listed above.
}

\references{
\url{http://nft.nefsc.noaa.gov/IRATE.html}

Jiang, H.H., K.H. Pollock, C. Brownie, J.M. Hoenig, R.J. Latour, B.K. Wells, and J.E. Hightower. 2007. Tag return models allowing for harvest and catch and release: evidence of environmental and management impacts on striped bass fishing and natural mortality rates. N. Amer. J. Fish. Man. 27: 387-396. \href{http://dx.doi.org/10.1577/M06-089.1}{doi: 10.1577/M06-089.1}

Hoenig, J.M., N.J. Barrowman, W.S. Hearn, and K.H. Pollock. 1998. Multiyear tagging studies incorporating fishing effort data. Can. J. Fish. Aquat. Sci. 55: 1466-1476. \href{http://dx.doi.org/10.1139/f97-256}{doi: 10.1139/f97-256}

Brownie, C., D.R. Anderson, K.P. Burnham, and D.S. Robson. 1985. Statistical inference from band recovery data: a handbook. Resource publication 156 (2nd edition), US Fish and Wildlife Service. 305 p.

}


\author{
Robert K. Bauer
}
