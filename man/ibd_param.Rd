\docType{data}
\name{ibd_prop_param}

\alias{ibd_prop_param}
\alias{quickParaIBD}
\alias{ibd_param}

\title{Simulation parameters estimated from pediatric inflammatory bowel disease (IBD) study (GSE57945)}
\format{
Simulation parameters for \code{simFromParam} function.
  \describe{
  \item{\code{health_lm_mean}}{A numeric vector containing the log-mean parameter estimates for each cell type from healthy samples.}
  \item{\code{health_lm_mean_d}}{A matrix containing the variance-covariance estimates of log-mean parameters across cell types from healthy samples.}
  \item{\code{lod_m}}{A numeric vector containing the log-dispersion parameter estimates for each cell type from healthy samples.}
  \item{\code{lod_d}}{A matrix containing the variance-covariance estimates of log-dispersion parameters across cell types from healthy samples.}
  \item{\code{health_alpha}}{A numeric vector of the estimated alpha parameter used to simulate cell type proportions for healthy samples.}
  \item{\code{case_alpha}}{A numeric vector of the estimated alpha parameter used to simulate cell type proportions for case samples.
}
  }
}
\usage{
data('quickParaIBD')
}
\description{
  An S4 object that stores simulation parameter estimated from the inflammatory bowel disease (IBD) study (GSE57945). This object contains a variety of numerical vectors and matrices representing different statistical parameters used in the simulation. The patients were drawn from healthy subjects from the immune-related diseases study.
}
\value{
One S4 object.
}
\examples{
data('quickParaIBD')
}
\keyword{datasets}
