\name{EWAS_series}
\alias{EWAS_series}
\title{
  Quality Control and Comparison of multiple EWAS results files
}
\description{
  This function runs a QC (via the function \code{\link{EWAS_QC}})
  over multiple files and generates additional graphs to
  comparing the results of these files.
}
\usage{
EWAS_series(EWAS_files,
            output_files,
            map,
            N,
            header_translations,
            save_final_dataset = TRUE,
            gzip_final_dataset = TRUE,
            high_quality_plots = FALSE,
            N_plot_beta = 500000L,
            ...)
}
\arguments{
  \item{EWAS_files}{
    a character vector containing the filenames of the EWAS
    results to be QC'ed.
}
  \item{output_files}{
    a character vector containing the filenames of the output
    files. Do not add an extension; \code{\link{EWAS_QC}} does
    so automatically.
}
  \item{map}{
    a data frame with chromosome and position values of the CpGs
    in \code{data}, or the name of a file containing the same.
    See \code{\link{EWAS_QC}} for details. This argument is
    optional: if not specified, \code{\link{EWAS_QC}} will not
    generate a Manhattan plot and no filter for X and Y markers
    can be performed.
}
  \item{N}{
    a data frame containing the filenames (as listed in the
    \code{EWAS_files} argument) and sample sizes of the
    datasets, or the name of a file containing the same. The
    data frame must contain the columns \code{file} and
    \code{N}, with those exact names. All files listed in
    the \code{EWAS_files} argument must be included in the 
    \code{file} column. This argument is optional: if not
    specified, \code{EWAS_series} will not generate a precision
    plot.
}
  \item{header_translations}{
    a translation table for the column names of the EWAS files,
    or the name of a file containing the same. See
    \code{\link{translate_header}} for details.
}
  \item{save_final_dataset,
  gzip_final_dataset,
  high_quality_plots}{
    logical values. See \code{\link{EWAS_QC}} for details.
}
  \item{N_plot_beta}{
    integer specifying how many beta values per file should
    be used in the effect-size comparison plot. Set this to
    a value larger than the number of markers in the datasets
    to include all markers.
}
  \item{\dots}{
    arguments passed to \code{\link{EWAS_QC}}.
}
}
\details{
 \code{QCEWAS} includes a Quick-Start guide in the doc
  folder of the library. This guide will explain how to
  run a QC and how to interpret the results.
  The start-up message when loading
  \code{QCEWAS} will indicate where it can be found on your
  computer. In brief,
  \code{EWAS_series} works by calling \code{\link{EWAS_QC}} for
  every filename given in \code{EWAS_files}. After all files
  have been processed, it will generate two additional graphs: a
  precision plot (provided \code{N} was specified) and a
  beta-distribution plot. The former shows the distribution of
  precision (1 / median standard error) against the
  square root of the sample size of the results file. Normally,
  one expects to see a roughly positive correlation (i.e. the
  cohorts ought to cluster around the linear diagonal from the
  lower left to the upper right). The presence of outliers means
  that the outlying cohort(s) have a far higher/lower
  uncertainty in their estimates that can be expected from their
  sample size. This could indicate a different method, a
  different measure (check the effect-size distribution plot) or
  possibly over- or undersignificance of their estimates (check
  the QQ plot and lambda value).

  The effect-size distribution plot allows comparison of the
  effect-size scale of different files. One  expects the
  distribution to become somewhat narrower as sample size
  increases. However, large differences in scale suggest that
  the files used different units for their measurements.
  
  As of version 1.2-0, the effect-size distribution plot shows
  a random (rather than proportional) selection of effect-sizes
  from the cohort. As a consequence, rerunning QC over a
  dataset may result in a slightly different distribution plot
  in each run. This is only a cosmetic issue (as the default
  sample size is sufficiently large to include the majority of a
  normally-sized EWAS dataset) and can be averted entirely by
  changing the \code{N_plot_beta} argument to a value exceeding
  the number of markers in the dataset(s).

  Both plots use numbers rather than names to identify files.
  The full filenames and corresponding numbers are listed in the
  EWAS_QC_legend.txt file that is generated after \code{EWAS_series}
  completes.
}
\value{
  The main output of \code{EWAS_series} are the cleaned results
  files, logs and graphs. The function also returns an invisible
  data frame (also saved as \code{EWAS_QC_legend.txt}), listing
  the input file names, file numbers,
  whether they passed a complete QC (note that this merely
  indicates that the QC was completed, not that there were no
  problems), the standard error and, if specified, the sample
  size.
}
\seealso{
  \code{\link{EWAS_QC}}
}
\examples{
# For use in this example, the 4 sample files in the
# extdata folder of the QCEWAS library will be copied
# to your current R working directory. Running the QC
# generates several files in your working directory:
# consult the Quick-Start Guide for more information
# on how to interpret these.
\dontrun{
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "sample_map.txt.gz"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "sample1.txt.gz"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "sample2.txt.gz"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "translation_table.txt"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)

sample_list <- c("sample1.txt.gz", "sample2.txt.gz")
sample_N <- data.frame(file = sample_list,
                       N = c(77, 79),
                       stringsAsFactors = FALSE)
                       


QC_results <- EWAS_series(EWAS_files = sample_list,
                          output_files = c("sample_output1", "sample_output2"),
                          map = "sample_map.txt.gz",
                          N = sample_N,
                          header_translations = "translation_table.txt",
                          save_final_dataset = FALSE,
                          threshold_outliers = c(-20, 20),
                          exclude_outliers = FALSE,
                          exclude_X = TRUE, exclude_Y = FALSE)
}
}