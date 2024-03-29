% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refitME_package.r
\docType{data}
\name{Corymbiaeximiadata}
\alias{Corymbiaeximiadata}
\title{The Corymbia eximia presence-only data set}
\format{
A data set that contains: 8 columns with 86,316 observations (or sites). The columns are defined as follows:
\describe{
\item{\code{X}}{Longitude coordinate.}
\item{\code{Y}}{Latitude coordinate.}
\item{\code{FC}}{Recorded number of fire counts for each site.}
\item{\code{MNT}}{Recorded minimum temperatures for each site.}
\item{\code{MXT}}{Recorded maximum temperature for each site.}
\item{\code{Rain}}{Recorded rainfall for each site.}
\item{\code{D.Main}}{Recorded distance from nearest major road.}
\item{\code{Y.obs}}{Presences for the plant species \emph{Corymbia eximia} for each site.}
}
}
\source{
See Renner and Warton (2013) for full details on the data and study.
}
\usage{
Corymbiaeximiadata
}
\description{
Data set consisting of presence-only records for the plant species \emph{Corymbia eximia}, site coordinates 5 covariates for each site.
}
\examples{
# Load the data.

data(Corymbiaeximiadata)
}
\references{
Renner, I. W. and Warton, D. I. (2013). Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. \emph{Biometrics}, \strong{69}, 274–281.
}
\keyword{datasets}
