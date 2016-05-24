\name{P_correlation}
\alias{P_correlation}
\title{Testing P-value distribution}
\description{
  A sub-function of \code{\link{EWAS_QC}}, \code{P_correlation}
  tests if the reported p-values match the p-value that can be
  derived from the effect-size and standard error values.
  Aberrations between these indicate that the p-values have been
  adjusted, or that there is some other problem with the data.
  It also creates a plot of reported vs. expected p-values that
  shows the correlation.
}
\usage{
  P_correlation(dataset,
                plot_correlation = TRUE,
                plot_if_threshold = FALSE,
                threshold_r = 0.99,
                save_name = "dataset",
                header_translations, ...)
}
\arguments{
  \item{dataset}{
    a data frame with the columns \code{BETA} (effect size),
    \code{SE} (standard error), and \code{P_VAL} (p value). If
    the column names differ from the above, the argument
    \code{header_translations} can be used to translate them.
}
  \item{plot_correlation}{
    logical, determines whether a graph is made of reported vs.
    expected p values.
}
  \item{plot_if_threshold}{
    logical. If \code{TRUE}, the plot is only generated if the
    p-value correlation is below the specified threshold.
}
  \item{threshold_r}{
    numeric. If the p-value correlation is below this, a warning
    is generated.
}
  \item{save_name}{
    character string used for the output file. Do not add an
    extension; \code{P_correlation} will do so automatically.
}
  \item{header_translations}{
    a translation table for the header of \code{dataset}. See
    \code{\link{translate_header}} for details.
}
  \item{\dots}{
    arguments passed to the generic R \code{\link{plot}} function.
}
}
\details{
  \code{P_correlation} is primarly a subfunction of
  \code{\link{EWAS_QC}}, but it can be used separately.
}
\value{
  \code{P_correlation} returns a single numeric value,
  representing the correlation between reported and expected
  p-values.
}
