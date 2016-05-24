\name{EWAS_QC}
\alias{EWAS_QC}
\title{Automated Quality Control of EWAS results files}
\description{
    The main function of the \code{\link{QCEWAS}} package.
    \code{EWAS_QC} accepts a single EWAS results file and runs a
    thorough quality check (including data integrity),
    optionally applies various filters and generates QQ, Volcano
    and Manhattan plots. The function \code{\link{EWAS_series}}
    can be used to process multiple results files sequentially.
}
\usage{
  EWAS_QC(data,
          map,
          outputname,
          header_translations,
          threshold_outliers = c(NA, NA),
          exclude_outliers = FALSE,
          exclude_X = FALSE, exclude_Y = FALSE,
          save_final_dataset = TRUE, gzip_final_dataset = TRUE,
          header_final_dataset = "standard",
          return_beta = FALSE, N_return_beta = 500000L,
          ...)
}
\arguments{
  \item{data}{
    a data frame with EWAS results, or the name of a file
    containing the same. The table must include the columns
    \code{PROBEID}, \code{BETA}, \code{SE}, and \code{P_VAL}.
    Other columns may be included but will be ignored. If the
    column names differ from the above, the argument
    \code{header_translations} can be used to translate them.
    If a filename is entered in this argument, it
    will be imported via the \code{\link{read.table}} function.
    \code{\link{read.table}} can handle a variety of formats,
    including files compressed in the .gz format. \code{EWAS_QC}
    will pass any named, unknown arguments to
    \code{\link{read.table}}, so you can specify the column
    separator and NA string in \code{EWAS_QC} with the usual
     \code{\link{read.table}} arguments. (Note that this only
     applied to importing the EWAS results, and not the map or
     translation files.)
}
  \item{map}{
    a data frame with chromosome and position values of the
    probes, or the name of a file containing the same. This
    argument is optional: if no map is specified,
    \code{EWAS_QC} will skip the Manhattan plot and chromosome
    filters. \code{map}	must include the columns \code{PROBEID},
    \code{CHR} (chromosome), and \code{POS} (position), using
    those exact names. Other columns may be included but will be
    ignored. f a filename is entered in this argument, it
    will be imported via the \code{\link{read.table}} function.
    \code{\link{read.table}} can handle a variety of formats,
    including files compressed in the .gz format.
}
  \item{outputname}{
    a character string specifying the intended filename for the
    output. This includes not only the cleaned results file and
    the log, but also any graphs created. Do not include an
    extension; \code{EWAS_QC} adds these automatically.
}
  \item{header_translations}{
    a translation table for the column names of the input file,
    or the name of a file containing the same. This argument is
    optional: if not specified, \code{EWAS_QC} assumes the
    default column names are used. See
    \code{\link{translate_header}} for information on the
    format.
}
  \item{threshold_outliers}{
    a numeric string of length two. This defines which effect
    sizes will be treated as outliers. The first value specifies
    the lower limit (i.e. markers with effect sizes below this
    value are considered outliers), the second the upper limit.
    The check for low or high outliers is skipped if the
    respective value is set to \code{NA}. To skip the check
    entirely, set this argument to \code{c(NA, NA)}.
}
  \item{exclude_outliers}{
    a logical value determining how outliers are treated. If
    \code{TRUE}, they are excluded from the final dataset. If
    \code{FALSE}, they are merely counted.
}
  \item{exclude_X, exclude_Y}{
    logical values determining whether markers at the X and Y
    chromosome respectively are excluded from the final dataset.
    This requires a \code{map} to be specified.
}
  \item{save_final_dataset}{
    logical determining whether the cleaned dataset will be
    saved.
}
  \item{gzip_final_dataset}{
    logical determining whether the saved dataset will be
    compressed in the .gz format.
}
  \item{header_final_dataset}{
    either a character vector or a table determining the header
    names used in the final dataset, or the name of a file
    containing the same. If \code{"original"}, the
    final dataset will use the same column names as the original
    input file. If \code{"standard"}, it will use the default
    \code{EWAS_QC} column names. If a table, it will be passed
    to \code{\link{translate_header}} to convert the column
    names. If a table, the default column names (\code{PROBEID},
    \code{BETA}, \code{SE}, and \code{P_VAL}) must be in the
    second column, and the desired column names in the first.
}
  \item{return_beta, N_return_beta}{
    arguments used by \code{\link{EWAS_series}}. These are not
    important for users and can be ignored. For the sake of
    completeness: \code{return_beta} is a logical value; if
    \code{TRUE}, the function return value includes a vector of
    effect sizes. \code{N_return_beta} defines the length of the
    vector.
}
  \item{\dots}{
    arguments passed to \code{\link{read.table}} for importing
    the EWAS results file.
}
}
\details{
  The Quality Control has the following 5 stages:
  
    \itemize{
    \item Checking data integrity:
    
    The values inside the EWAS results are tested for validity.
    If impossible p-values, effect-sizes, etc. are encountered,
    \code{EWAS_QC} generates a warning in the R console and sets
    them to \code{NA}.
    
    \item Filter for outliers and sex-chromosomes (optional)
    
    Counts the number of outlying markers, as well as chromosome
    X and Y markers, and deletes them if specified.
    
    \item Generating QC plots
    
    A histogram of beta and standard error distribution is
    plotted.
    
    The p-values are checked by correlating and plotting them
    against p-values calculated from the effect size and
    standard error.
    
    A QQ plot is generated to test for over/undersignificance.
    
    A Manhattan plot is generated to see where the signals (if
    any) are located.
    
    A Volcano plot is generated to check the distribution of
    effect sizes vs. p values.
    
    \item Creating a QC log
    
    The log contains notes about any problems encountered during
    the QC, as well as several tables describing the data.
    
    \item Saving the cleaned dataset (optional)
  }
}
\value{
  The main output of \code{EWAS_QC} are the cleaned results
  file, log file and QC graphs. However, the function also
  returns a list with 9 elements:
  
  \item{data_input}{the file name of the input file, if loaded
  from a file. If not, this will be an empty character string.}
  \item{file}{the filename of the cleaned results file.}
  \item{QC_success}{logical, indicates whether \code{EWAS_QC} 
  was able to run a full QC on the file. Note that a \code{TRUE}
  value does not mean that no problems where encountered,
  merely that the full QC was executed.}
  \item{lambda}{the lambda value of reported p-values in the
  cleaned dataset.}
  \item{p_cor}{the correlation between reported and expected
  (based on effect size and standard error) p values.}
  \item{N}{a named integer vector reporting how many markers
  were in the original dataset, how many had missing values,
  how many were on chromosomes X and Y, how many were outliers,
  how many were removed and how many are in the final, cleaned
  dataset. Has no relation to the \code{N} argument of
  \code{\link{EWAS_series}}.}
  \item{SE_median}{a numeric value: the median of the standard
  errors in the cleaned dataset.}
  \item{mean_methylation}{a \code{NULL}: this functionality
  has not been implemented yet.}
  \item{effect_size}{if \code{return_beta} is \code{TRUE}, this
  is a numeric vector of length \code{N_return_beta},
  containing a representative selection of effect sizes from the
  filtered dataset. If \code{FALSE}, this will be \code{NULL}.}
}
\note{
  The function will return a warning if it encounters p-values
  < 1e-300, as this is close to the smallest number that R can
  process correctly. Various functions in the \code{QCEWAS}
  package will set these values to 1e-300 to ensure proper
  handling.
}
\seealso{
  See \code{\link{EWAS_series}} for running a QC over multiple
  files.
  
  See \code{\link{EWAS_plots}} and \code{\link{P_correlation}}
  for carrying out specific steps of the QC.
}
\examples{
# For use in this example, the 2 sample files in the
# extdata folder of the QCEWAS library will be copied
# to your current R working directory. Running the QC
# generates 7 new files in your working directory:
# a cleaned, post-QC dataset, a log file, and 5 graphs.
# Consult the manual for more information on how to
# interpret these.
  
\dontrun{

file.copy(from = file.path(system.file("extdata", package = "QCEWAS"), "sample_map.txt.gz"),
               to = getwd(), overwrite = FALSE, recursive = FALSE)
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"), "sample1.txt.gz"),
               to = getwd(), overwrite = FALSE, recursive = FALSE)

QC_results <- EWAS_QC(data = "sample1.txt.gz",
                      map = "sample_map.txt.gz",
                      outputname = "sample_output",
                      threshold_outliers = c(-20, 20),
                      exclude_outliers = FALSE,
                      exclude_X = TRUE, exclude_Y = FALSE,
                      save_final_dataset = TRUE, gzip_final_dataset = FALSE)

}
}
