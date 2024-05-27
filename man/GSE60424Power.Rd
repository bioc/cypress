
\name{GSE60424Power}

\alias{GSE60424Power}
\alias{quickPowerGSE60424}
\alias{GSE60424Power}
\alias{iad_param}

\docType{data}
\title{Power calculation results from immune-related disease (IAD) study (GSE60424)}
\format{
  A S4 object.
  \describe{
  \item{\code{ct_TDR_bio_smry}}{Cell-type-specific target TDR.}
  \item{\code{TDR_bio_smry}}{Average target TDR.}
  \item{\code{ct_PWR_bio_smry}}{Cell-type-specific target power.}
  \item{\code{PWR_bio_smry}}{Average target power.}
  \item{\code{PWR_strata_ct_bio_smry}}{Cell type specific target power by gene expression stratification.}
  \item{\code{PWR_strata_bio_smry}}{Average target power by gene expression stratification.}
  \item{\code{ct_FDC_bio_smry}}{Cell type specific target FDC.}
  \item{\code{FDC_bio_smry}}{Average target FDC.}

  }
}

\usage{
data('quickPowerGSE60424')
}
\description{
Pre-calculated power evaluation results from immune-related disease (IAD) study (GSE60424). The results
can be used to create plots using plot functions (\code{plotFDC}, \code{plotPower}, \code{plotTDR}).
}
\value{
One S4 object.
}
\examples{
data('quickPowerGSE60424')

}
\keyword{datasets}
