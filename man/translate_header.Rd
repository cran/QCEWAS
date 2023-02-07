\name{translate_header}
\alias{translate_header}
\title{Translate column names into standard names}
\description{
  This function is used to translate non-standard column names
  into the standard ones used by \code{\link{EWAS_QC}} and other
  functions.
}
\usage{
  translate_header(header,
                   standard = c("PROBEID","BETA","SE","P_VAL"),
                   alternative)
}
\arguments{
  \item{header}{character vector; the header to be translated.}
  \item{standard}{character vector; the names \code{header}
    should be translated into.}
  \item{alternative}{translation table; see below for more information.}
}
\details{
  The function takes the entries in \code{standard} one by one, and
  checks them against the translation table for alternatives.
  It will report any missing standard headers, as well as
  duplicate ones.
}
\value{
  \code{translate_header} returns an object of class 'list' with 6 components:
  \item{header_h}{character vector; the translated header.
    Unknown columns are included under their old names.}
  \item{missing_h}{character vector; the standard column names
    that were not found. If none, this returns \code{NULL}.}
  \item{unknown_h}{character vector; column names that could not
    be converted to a standard name. Note that these columns
    are also included in \code{header_h}. If none, this returns
    \code{NULL}.}
  \item{header_N, missing_N, unknown_N}{integer; the lengths of
    the above three vectors}
}
\section{Translation Table}{
  The translation table must meet the following requirements:
  \itemize{
    \item 2 columns, with the default column names (i.e. the ones
    in the standard argument) in the first column, and the
    alternatives in the second.
    \item Multiple alternatives are allowed for a single standard
    name, but every alternative name must be in a separate row.
    \item The alternatives must be capitalized.
    \item No duplicate alternatives are allowed.
    \item A header line is not required, and will be ignored if present.
  }
}
\note{
  The function will automatically capitalize the elements of
  the \code{header} argument (so the alternatives in the translation
  table must also be capitalized). Also, elements that are not
  in \code{standard} will not be translated, even if they are present
  in the translation table.
}
\examples{
# For use in this example, the 2 sample files in the
# extdata folder of the QCEWAS library will be copied
# to your current R working directory
\dontrun{
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "sample2.txt.gz"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)
file.copy(from = file.path(system.file("extdata", package = "QCEWAS"),
                           "translation_table.txt"),
          to = getwd(), overwrite = FALSE, recursive = FALSE)

sample_ewas <- read.table("sample2.txt.gz", header = TRUE,
                          stringsAsFactors = FALSE, nrow = 10)
colnames(sample_ewas)


translation_table <- read.table("translation_table.txt", header = TRUE,
                                stringsAsFactors = FALSE)
sample_translation <- translate_header(header = colnames(sample_ewas),
                                       alternative = translation_table)
sample_translation                                         

colnames(sample_ewas) <- sample_translation$header_h

colnames(sample_ewas)
}
}