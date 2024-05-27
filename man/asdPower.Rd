\name{asd_nopropPower}

\alias{asd_nopropPower}
\alias{quickPowerASD}
\alias{asdPower}

\docType{data}
\title{Power calculation results From ASD data}
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
data('quickPowerASD')
}

\description{
Pre-calculated power evaluation results from  Autism Spectrum Disorder (ASD) study. The results
can be used to create plots using plot functions (\code{plotFDC}, \code{plotPower}, \code{plotTDR}).
}
\value{
One S4 object.
}
\examples{
data('quickPowerASD')
}
\keyword{datasets}
