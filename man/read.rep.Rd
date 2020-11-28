\name{read.rep}
\alias{read.rep}
\alias{read.short.rep}
\title{Read report-file from an IRATE run}
\description{
  Reads in the report-file from an IRATE run
}
\usage{
  read.rep(fn, skip=0, show.sm=F, short=F)
  read.short.rep(fn, skip=0, show.sm=F, short=T)
}

\arguments{
  \item{fn}{
Character string indicating the name of the \code{.rep}-file to be read.
}
  \item{skip}{
Number of lines to skip in report file (default = 0).
}
  \item{show.sm}{
Show model report summary while reading the file.
}
  \item{short}{
read short report or long report.
}

}

\value{
a list containing the model statistics and parameter estimates from an IRATE run:

for \code{read.short.rep} or \code{read.rep(short=T)}:

\code{Log_L}, \code{AIC}, \code{AICc}, \code{Effective_Sample_Size}, \code{Unpooled_Chi_square}, \code{Upooled_df}, \code{Unpooled_c_hat}, \code{Pooled_Chi_square}, \code{Pooled_df}, \code{Pooled_c_hat}

othervise the full report file will be returned, containing the following 


\describe{
\item{\code{rep_info}}{report file information}
\item{\code{Log_L}}{Log-Likelihood}
\item{\code{K}}{Number of parameters}
\item{\code{AIC}}{Akaike information criterion}
\item{\code{AICc}}{AIC corrected for finite sample sizes}
\item{\code{Effective_Sample_Size}}{}
\item{\code{Unpooled_Chi_square}}{}
\item{\code{Upooled_df}}{}
\item{\code{Unpooled_c_hat}}{}
\item{\code{Pooled_Chi_square}}{}
\item{\code{Pooled_df}}{}
\item{\code{Pooled_c_hat}}{}
\item{\code{Harvest_Residuals, Harvest_Residuals_of_Age_Class_*}}{}
\item{\code{Release_Residuals, Release_Residuals_of_Age_Class_*}}{}
\item{\code{Not_seen_Residuals, Not_seen_Residuals_of_Age_Class_*}}{}
\item{\code{Obs_Recoveries_of_harvest_fish}}{}
\item{\code{Obs_Recoveries_of_release_fish}}{}
\item{\code{Total_Released}}{}
\item{\code{Total_Recovered_Tags}}{}
\item{\code{s_matrix}}{matrix of survival rates of the tags (Sj,tags)}
\item{\code{S_prob_matrix}}{matrix of survival rates used in the calculation of expected cell probabilities}
\item{\code{Exploitation_Rate_of_harvested_fish}}{matrix of estimated exploitation rates for harvested fish (Uj,tags(kept))}
\item{\code{Exploitation_Rate_of_released_fish}}{matrix of estimated exploitation rates for released fish (Uj,tags(rels))}
\item{\code{Expected_Probability_of_harvested_fish}}{matrix of cell expected probabilities for harvested tag returns (Pij)}
\item{\code{Expected_Probability_of_released_fish}}{matrix of cell expected probabilities for released tag returns (P'ij)}
\item{\code{Not_Seen_Probability}}{matrix of cell expected probabilities for not seen tag returns (P-sum(Pij+P'ij))}
\item{\code{Expected_Number_of_harvested_fish}}{matrix of expected number of harvested tag returns (E(Rij)}
\item{\code{Expected_Number_of_released_fish}}{matrix of expected number of released tag returns (E(R'ij))}
\item{\code{Expected_Number_of_not_seen}}{matrix of expected number of not seen tag returns (N-sum(Rij+R'ij))}
\item{\code{Cell_Likelihoods_of_harvested_fish}}{individual cell likelihood values for harvested tag returns}
\item{\code{Cell_Likelihoods_of_released_fish}}{individual cell likelihood values for released tag returns}
\item{\code{Cell_Likelihoods_of_not_seen}}{individual cell likelihood values for not seen tag returns}
\item{\code{Unpooled_Chi_squares_of_Harvested_Fish}}{individual unpooled cell chi-square values for harvested tag returns}
\item{\code{Unpooled_Chi_squares_of_Released_Fish}}{individual unpooled cell chi-square values for released tag returns}
\item{\code{Chi_squares_of_Not_Seen}}{individual unpooled cell chi-square values for unseen tag returns}
\item{\code{Pooled_Cells_of_Harvested_Fish}}{matrix of pooled, observed harvested tag returns}
\item{\code{Pooled_Cells_of_Released_Fish}}{matrix of pooled, observed harvested tag returns}
\item{\code{Pooled_Expected_Cells_of_Harvested_Fish}}{matrix of pooled expected number of harvested tag returns}
\item{\code{Pooled_Expected_Cells_of_Released_Fish}}{matrix of pooled expected number of released tag returns}
\item{\code{Pooled_Chi_squares_of_Harvested_Fish}}{matrix of cell chi-square values for harvested tag returns }
\item{\code{Pooled_Chi_squares_of_Released_Fish}}{matrix of pooled cell chi-square values for released tag returns}
\item{\code{Pearson_Residuals_for_harvested_fish}}{matrix of Pearson residuals for harvested tag returns}
\item{\code{Pearson_Residuals_for_released_fish}}{matrix of Pearson residuals for released tag returns}
\item{\code{Pearson_Residuals_for_not_seen}}{matrix of Pearson residuals for not seen tag returns}

}


}

\seealso{\link{read.par}}

\author{
Robert K. Bauer
}
