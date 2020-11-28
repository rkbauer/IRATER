\name{IRATE.examples}
\alias{IRATE.examples}
\title{list preinstalled IRATE examples}
\description{
  list preinstalled IRATE examples
}
\usage{
  IRATE.examples()
}

\value{
\describe{

\item{\code{HoenigNonMixing}}{Age-dependent, harvest only (no catch-&-release) conventional tagging dataset from 1960-1964 of lake trout (Salvelinus namaycush), assuming an incomplete mixing of tagged fish withe the wild population after being released. Based on: Hoenig et al. 1998 (incomplete mixing Table 4 for more details).}

\item{\code{JiangADCR}}{Age-dependent, harvest and catch-&-release conventional tagging dataset of Striped bass (Morone saxatilis) from 1991-2003, assuming an incomplete mixing of tagged fish withe the wild population after being released (Similar to the \code{JiangADH} and \code{JiangAICR} data sets). Based on: Jiang et al., 2007 (Table 2 Model b).}

\item{\code{JiangADH}}{Age-dependent, harvest only (no catch-&-release) conventional tagging dataset of Striped bass (Morone saxatilis) from 1991-2003 (Similar to the \code{JiangADH} and \code{JiangAICR} data sets), assuming an incomplete mixing of tagged fish withe the wild population after being released. Based on: Jiang et al., 2007b (Table 4 Model a).}

\item{\code{JiangAICR}}{Age-independent, harvest and catch-&-release conventional tagging dataset of Striped bass (Morone saxatilis) from 1991-2003 (Similar to the \code{JiangADCR} and \code{JiangADH} data sets). Based on Jiang Dissertation (Page 81 Model Fy,F'y,M91-99, M00-03).}

}
}

\references{
Hoenig JM, Barrowman NJ, Pollock KH, Brooks EN, Hearn WS, Polacheck T. (1998) Models for tagging data that allow for incomplete mixing of newly tagged animals. Canadian Journal of Fisheries and Aquatic Sciences 55:1477-1483, doi: 10.1139/f97-258

Jiang H. (2005) Age-dependent tag return models for estimating fishing mortality, natural mortality, and selectivity. Ph.D. thesis, North Carolina State University, Raleigh

Jiang H, Pollock KH, Brownie C, Hightower JE, Hoenig JM, Hearn WS. (2007b) Age-dependent tag return models for estimating fishing mortality, natural mortality, and selectivity. Journal of Agricultural, Biological, and Environmental Statistics 12:177-194, doi: 10.1198/108571107X197382

Jiang H, Pollock KH, Brownie C, Hoenig JM, Latour RJ, Wells BK, Hightower JE. (2007) Tag return models allowing for harvest and catch and release: Evidence of environmental and management impacts on striped bass fishing and natural mortality rates. North American Journal of Fisheries Management 27:387-396, doi: 10.1577/M06-089.1

}

\examples{
IRATE.examples()
}


\seealso{\link{run.IRATE.example} and \link{remake.dat} to change the model configuration of an existing IRATE-setup} 

\author{
Robert K. Bauer
}
