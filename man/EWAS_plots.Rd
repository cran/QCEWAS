\name{EWAS_plots}
\alias{EWAS_plots}
\title{Manhattan and Quantile-Quantile plots for EWAS results files}
\description{
  This function is used by \code{\link{EWAS_QC}} to generate
  quantile-quantile (QQ) and Manhattan plots. It can also be
  called by users. Note that it does not generate the histogram
  or volcano plot - this is done by \code{\link{EWAS_QC}}
  itself.
}
\usage{
  EWAS_plots(dataset,
             plot_QQ = TRUE,
             plot_Man = TRUE,
             plot_cutoff_p = 0.05,
             plot_QQ_bands = FALSE,
             high_quality_plots = FALSE,
             save_name = "dataset",
             header_translations)
}
\arguments{
  \item{dataset}{
    either a vector of p-values, or a data frame containing the
    columns \code{CHR} (chromosome number), \code{MAPINFO}
    (base-pair position), and \code{P_VAL} (p-value). If
    different columnnames are used, the
    \code{header_translations} argument can be used to translate
    these. \code{CHR} and \code{MAPINFO} are only required for
    generating a Manhattan plot. Note that, unlike
    \code{EWAS_QC}, this function does not accept filenames,
    only data frames or vectors.
}
  \item{plot_QQ,
        plot_Man}{
    logicals determining whether a QQ and Manhattan plot are
    made.
}
  \item{plot_cutoff_p}{
    numeric: the threshold of p-values to be shown in the QQ and
    Manhattan plots. Higher (less significant) p-values are
    excluded from the plot. The default setting is 0.05, which
    excludes 95\% of data-points. It's \emph{NOT} recommended to
    increase the value above 0.05, as this may dramatically
    increase running time and memory usage.
}
  \item{plot_QQ_bands}{
    logical, if \code{TRUE}, probability bands are added to the
    QQ plot.
}
  \item{high_quality_plots}{
    logical. Setting this to TRUE will save the graphs as 
    high-resolution tiff images.
}
  \item{save_name}{
    character string, the name used for the plot files (do not
    add an extension: \code{EWAS_plots} will do this
    automatically).
}
  \item{header_translations}{
    a table that translates the column names of \code{dataset}
    to the standard names. See \code{\link{translate_header}}
    for details.
}
}
\details{
  \code{EWAS_plots} is a fairly straightforward function. It
  accepts a data table or a vector of p-values, and generates
  QQ and (when chromosome and position data are included)
  Manhattan plots from these.
}
\value{
  \code{EWAS_plots}' most important output are the two graphs.
  However, it also returns a single, invisible, numeric value,
  representing the lambda calculated over the p-values.
}