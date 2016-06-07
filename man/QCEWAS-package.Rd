\name{QCEWAS-package}
\alias{QCEWAS-package}
\alias{QCEWAS}
\docType{package}
\title{Quality Control of Epigenome-Wide Association Study Results}
\description{
  Functions for automated and manual quality control of
  Epigenome-Wide Association Study results.
}
\details{
  \tabular{ll}{
    Package: \tab QCEWAS\cr
    Type: \tab Package\cr
    Version: \tab 1.0-7\cr
    Date: \tab 2016-05-26\cr
    License: \tab GPL (>= 3)\cr
  }
  When combining the results of multiple Epigenome-Wide
  Association Studies (EWAS) into a meta-analysis, it is
  recommended to run a quality check (QC) over the individual
  files first. This is not only to check if the individual
  results are valid, reliable, and of high quality, but also to
  ensure that they are comparable between cohorts. The
  \code{QCEWAS} package was developed to automate this process.
  
  The core function of \code{QCEWAS} is \code{\link{EWAS_QC}}.
  This function accepts the filename of a single EWAS results
  file and performs a thorough quality check. For more
  information, see the \code{\link{EWAS_QC}} help page. The
  function \code{\link{EWAS_series}} is a wrapper function that
  accepts multiple file names, and then runs \code{\link{EWAS_QC}}
  for all of them. It also generates a few additional plots to
  compare the results of the files, allowing the user to spot
  differences in effect-size distribution or standard errors.
  
  The functions \code{\link{EWAS_plots}} and
  \code{\link{P_correlation}} are subroutines of
  \code{\link{EWAS_QC}} that can also be called by the user to
  perform specific QC tasks.
}
\author{
P.J. van der Most, Leanne K. Kupers and Ilja M. Nolte

Maintainer: P.J. van der Most <p.j.van.der.most@umcg.nl>
}
\keyword{ package }
