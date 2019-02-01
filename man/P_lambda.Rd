\name{P_lambda}
\alias{P_lambda}
\title{Calculation of the Lambda value}
\description{
  The Lambda value represents the inflation of p-values
  compared to a normal distribution of p.
}
\usage{
P_lambda(p)
}
\arguments{
  \item{p}{a numeric vector of p-values}
}
\details{
  The function removes any missing values from \code{p}, and then
  returns:
  
  \code{median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)}
  
  The lambda value represents the inflation of the p-values
  compared to a normal distribution. In a
  genome-wide study, one would expect the results for the
  vast majority of CpG
  sites to accord with the null hypothesis, i.e. the p-values
  are random, and have a normal
  distribution. Only sites that are significantly associated
  with the phenotype of interest should lie outside of the
  normal distribution.
  
  Ideally the lambda value should be 1. Lambda represents the
  \emph{overall} difference with the expected distribution
  - so the presence of a few significant results (i.e. p-values
  that do not follow the normal distribution) does not bias it.
  
  However, if lambda is 2 or higher, it means that a
  substantial portion of your dataset is more significant than
  expected for a genome-wide study (i.e.
  oversignificance). This could mean your dataset has been
  filtered for low-significance markers. If this is not
  the case, you should consider doing a genomic control
  correction on the p-values, to correct the oversignificance.
  
  Similary,
  values of 0.8 or lower indicate that your results are less
  significant than would be expected from a random distribution
  of p-values.
}
\value{
  A single numeric value, the lambda value.
}
\examples{
  pvector <- ppoints(10000)
  P_lambda(pvector)
  # The lambda of a random distribution of p-values equals 1
  
  pvector[pvector > 0.9 & pvector < 0.91] <- NA
  P_lambda(pvector)
  # If low-significance results are removed (i.e. there are more
  # significant results than expected) lambda increases
}
