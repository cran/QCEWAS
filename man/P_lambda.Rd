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
  
  The lambda value represents the inflation of the p-values in a
  results file compared to a normal distribution. In a
  genome-wide study, one would expect the results for the
  majority of CpG
  sites to accord with the null hypothesis, i.e. the p-values
  are random, and combining them results in a normal
  distribution. Only sites that are significantly associated
  with the phenotype of interest should lie outside of the
  normal distribution.
  
  The lambda value represents the overall difference with the
  expected distribution - so the presence of a few significant
  results (i.e. p-values that do not follow the normal
  distribution) does not bias it. Ideally the lambda value
  should be 1. Values of 1.1 or higher are cause for concern,
  as it means that the p values are lower (more significant)
  than expected. Similary, but far less frequently,
  values of 0.9 or lower indicate that your results
  have a higher p-value than expected.
}
\value{
  A single numeric value, the lambda value.
}
\examples{
  # Create a distribution of probabilities:
  pvector <- ppoints(10000)
  # and introduce bias by removing a selection of high values:
  pvector[pvector > 0.9 & pvector < 0.91] <- NA
  P_lambda(pvector)
}
