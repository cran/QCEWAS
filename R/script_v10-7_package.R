# QCEWAS v1.0-7
# Created by Peter van der Most, April 2015 - May 2016 with samples
# Based on code for the QCGWAS package


# Remember to adjust version number here and in output log, and activate onAttach functie
# don't export the zf_testLogical function


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("QCEWAS library, version 1.0-7")
}


# zf_testLogical is an internal function to test if arguments are logical
zf_testLogical <- function(param, notNA = TRUE){
  OK <- is.logical(param) & is.vector(param) & length(param) == 1L
  if(OK & notNA) { OK <- !is.na(param) }
  return(OK)    
}


# Function P_correlation is no longer the same check_p as in QCGWAS. The save name has been changed, the log entries & HQ subset removed and the colnames have been updated
P_correlation <- function(dataset,
                    plot_correlation = TRUE, plot_if_threshold = FALSE, threshold_r = 0.99,
                    save_name = "dataset", header_translations, ...) {
  
  if(!missing(header_translations)) {
    header_test <- translate_header(header = colnames(dataset), standard = c("P_VAL", "BETA", "SE"), alternative = header_translations)
    if(any(duplicated(header_test$header_h))) stop("cannot translate header - duplicate column names")
    if(header_test$missing_N > 0L) stop(paste("Cannot identify data column(s):", paste(header_test$missing_h, collapse = ", ")))
    colnames(dataset) <- header_test$header_h
  }
  
  goodOnes <- !(is.na(dataset$P_VAL) | is.na(dataset$BETA) | is.na(dataset$SE))
  if(sum(goodOnes) > 10L) {
    p_obs <- -log10(dataset$P_VAL[goodOnes])
    p_exp <- -log10(pchisq((dataset$BETA[goodOnes]/dataset$SE[goodOnes])^2, 1, lower.tail=FALSE))
    if(any(p_obs > 300)) {
      p_obs[p_obs > 300] <- 300
      warning("Extreme p-values in dataset - temporarly set to 1e-300 for correlation calculation and plot")
    }
    if(any(p_exp > 300)) {
      p_exp[p_exp > 300] <- 300
      warning("Expected p-values capped at 1e-300")
    }
    p_cor <- cor(p_obs, p_exp, use = "everything")
    if(is.na(p_cor)){
      warning("Unable to calculate p-value correlation")
      return(NA) # uses return here so that the if loops below are skipped      
    }
    if(p_cor < threshold_r) {
      print(paste0("Warning: reported p-values correlate poorly to expected values(r = ", round(p_cor, digits = 3), ")"), quote = FALSE)
    }
    if(plot_correlation & (p_cor < threshold_r | !plot_if_threshold )) {
      p_max <- max(p_obs, p_exp)
      png(paste0(save_name, ".graph_p.png"),
          width = 720, height = 720, res = 144)
      plot(x = 0, y = 0, col = "white",
           xlim = c(0, p_max), ylim = c(0, p_max),
           main="P-value correlation", xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
           sub = save_name, cex.sub = 1.3, ...)
      lines(x = c(0, 1.1 * p_max), y = c(0, 1.1 * p_max))
      points(p_exp, p_obs, col = "black", pch = 20, cex = 0.8)
      text(0.0 * p_max, 0.95 * p_max, paste("r =", round(p_cor, digits = 3)), pos = 4, cex=1.0, col = ifelse(p_cor < threshold_r, "red", "black") )
      dev.off()
    }
  } else {
    print("p-test aborted: less than 10 entries with the data required to calculate corelation", quote = FALSE)
    p_cor <- NA
  }
  return(p_cor)
}


# Function EWAS_plots is an edited version of QC_plots from QCGWAS.
# EWAS_plots was edited to prevent a bug from occuring, and to suppress messages
# as well as suppresing an underline in the QQ plot
# also removed a HQ filter call
# lateron, we removed all code relating to the filters

P_lambda <- function(p){
  p <- p[!is.na(p)]
  if(length(p) < 2L) stop("'p' does not contain sufficient non-missing values to calculate lambda") 
  return(median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1))
}

EWAS_plots <- function(dataset, plot_QQ = TRUE, plot_Man = TRUE,
                     plot_cutoff_p = 0.05,
                     plot_QQ_bands = FALSE,
                     save_name = "dataset",
                     header_translations) {
  
  if(is.vector(dataset)) { dataset <- data.frame(P_VAL = dataset) }
  colnames(dataset) <- toupper(colnames(dataset))
  header_std <- c("P_VAL", "CHR", "POS")[c(TRUE, plot_Man, plot_Man)]
  
  if(missing(header_translations)) {
    header_missing <- header_std[!header_std %in% colnames(dataset)]
    if(length(header_missing) > 0L) { stop(paste0("Missing columns: ", paste0(header_missing, collapse = ", "))) }
  } else {
    header_test <- translate_header(header = colnames(dataset), standard = header_std, alternative = header_translations)
    if(any(duplicated(header_test$header_h))) stop(paste0("Cannot translate header - duplicate columns: ", paste0(header_test$header_h[duplicated(header_test$header_h)], collapse = ", ") ))
    if(header_test$missing_N > 0L) stop(paste("Cannot identify data column(s):", paste(header_test$missing_h, collapse = ", ")))
    colnames(dataset) <- header_test$header_h
  }
    
   
  # Stage 1: creating filters & calculating lambda
 
  # Calculating lambda
  # First the script creates a short-list by throwing out NA's (!) to calculate the lambda's,
  #	and a p-cutoff filter is created (!). The script checks whether there are still more than
  #	10 datapoints available at both exclamation marks. If yes, the script uses the shortlist
  #	to generate the unfiltered expected and observed QQ plots; and the long list plus the 
  #	addQQplot function to generate the filtered plots.
  #	WARNING: note that the content of QQ_obs_short changes: it starts out as the p-values
  #	minus NA (*), then -log10 values (#), then -log10 values minus p > 0.05 (***). Lambda needs
  #	to be calculated at the (*). QQ_exp is calculated at (#), but (*) works too; the contents
  #	of the vector used to calculate QQ_exp do not matter, only its length. QQ_exp changes
  #	as well: a copy without > 0.05 is made at (***: QQ_exp_short). The copy is used for
  #	the graphs, but QQ_exp is necessary for the probability clouds. If so, it's shortened
  #	to 1000 points to prevent the polygon functions from being overloaded.
  
  QQ_obs_p	 <- dataset$P_VAL
  QQ_obs_short <- subset(QQ_obs_p, !is.na(QQ_obs_p))
  QQ_obs_N	 <- length(QQ_obs_short)	# N of non-missing p's: the N of all p's is dataN
  
  if(QQ_obs_N > 10L) {
    lambda <- P_lambda(QQ_obs_short)
    
    if(any(QQ_obs_short < 1e-300)) {
      warning("Extreme p-values in dataset - temporarly set to 1e-300 for QQ & Manhattan plots")
      QQ_obs_short <- ifelse(QQ_obs_short < 1e-300, 1e-300, QQ_obs_short)
      dataset$P_VAL <- ifelse(dataset$P_VAL < 1e-300 & !is.na(dataset$P_VAL), 1e-300, dataset$P_VAL)
    }
    
    QQ_obs_short<- sort(-log10(QQ_obs_short))
    QQ_incl	<- QQ_obs_short >= -log10(plot_cutoff_p)
  } else {
    lambda	<- NA
    QQ_incl	<- 0
  }
  
  # Stage 2: creating QQ plots
  if(plot_QQ & sum(QQ_incl) > 10L ) {
    QQ_exp	<- sort(-log10(ppoints(QQ_obs_N)))
    QQ_exp_short<- QQ_exp[QQ_incl]
    QQ_obs_short<- QQ_obs_short[QQ_incl]
    QQ_exp_min	<- QQ_exp_short[1]
    QQ_exp_max	<- QQ_exp_short[length(QQ_exp_short)]
    QQ_obs_min	<- QQ_obs_short[1]
    QQ_obs_max	<- QQ_obs_short[length(QQ_obs_short)]
    
    if(plot_QQ_bands) {
      temp <- (1:QQ_obs_N)
      i1000 <- c(1, (1:1000) * floor(QQ_obs_N / 1000), QQ_obs_N)
      QQ_band_upper <- sort(-log10(qbeta( 1 - 0.05 / 2, temp, QQ_obs_N - temp + 1 ) ) )[i1000]
      QQ_band_lower <- sort(-log10(qbeta(		 0.05 / 2, temp, QQ_obs_N - temp + 1 ) ) )[i1000]
      QQ_exp <- QQ_exp[i1000]
    }
    
    png(paste0(save_name, ".graph_QQ.png"),
        width = 720, height = 720, res = 144)
    plot(c(QQ_exp_min, QQ_exp_max), c(QQ_obs_min, QQ_obs_max), xlim = c(0, QQ_exp_max), ylim = c(0, QQ_obs_max),
         main = "QQ plot", xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
         pch = 1, col = "black", sub = save_name, cex.sub = 1.3)
    if(plot_QQ_bands) {
      polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_upper, rev(QQ_exp)), col="grey", border = NA )
      polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_lower, rev(QQ_exp)), col="grey", border = NA )
    }
    points(QQ_exp_short, QQ_obs_short, pch = 1, col = "black")
    abline(0,1)
    text(0, 0.98 * QQ_obs_max, substitute(paste(lambda, " = ", x), list(x=round(lambda, digits = 3))),
         pos = 4, col = ifelse(lambda > 1.5, "red", "black"))
    dev.off()
  } else { if(plot_QQ) warning("Insufficient significant p-values to create QQ plot") }
  
  # Stage 3: creating Manhattan plot
  if(plot_Man & sum(QQ_incl) > 10L ) {
    
    man_set <- dataset[!is.na(dataset$POS) & !is.na(dataset$CHR) & !is.na(dataset$P_VAL) & dataset$P_VAL <= plot_cutoff_p, c("CHR", "POS", "P_VAL")]
    
    if(!is.numeric(man_set$CHR)) {
      if(is.character(man_set$CHR)){
        man_set$CHR <- toupper(man_set$CHR)
        man_set$CHR[man_set$CHR == "X"] <- 23L
        man_set$CHR[man_set$CHR == "Y"] <- 24L
        man_set$CHR[man_set$CHR == "XY"] <- 25L
        man_set$CHR[man_set$CHR %in% c("M", "MT")] <- 26L
      }
      man_set$CHR <- as.integer(man_set$CHR)
      if(any(is.na(man_set$CHR))) {
        warning(paste0(sum(is.na(man_set$CHR)), " entries with untranslatable chromosome values"))
        man_set <- man_set[!is.na(man_set$CHR), ]
      }
    }
    
    if(any(!man_set$CHR %in% 1:26)) {
      warning(paste0(sum(!man_set$CHR %in% 1:26), " entries with chromosome values outside of the normal range"))
      man_set <- man_set[man_set$CHR %in% 1:26, ]				
    }
    manhattanN <- nrow(man_set)
    
    
    if(manhattanN > 9L) { 		# testing if there are sufficient p-values < 0.05
      chr_size <- data.frame(chromosome = 1:27, #	chr2				chr3			chr4				chr5				chr6				chr7				chr8				chr9			chr10				chr11			 chr12			chr13				chr14				chr15				chr16			chr17			chr18				chr19			chr20				chr21			chr22			X23					Y24		XY25,M26, 27 = end of M
                             size = c(247249719, 242951149, 199501827, 191273063, 180857866,	170899992,	158821424,	146274826,	140273252,	135374737,	134452384,	132349534,	114142980,	106368585,	100338915,	 88827254,	 78774742,	 76117153,	 63811651,	 62435964,	 46944323,	 49691432,	154913754,	 57772954,	0,	0,	0 ),
                             start= c(	 500000, 248249719, 491700868, 691702695, 883475758, 1064833624, 1236233616, 1395555040, 1542329866, 1683103118, 1818977855, 1953930239, 2086779773, 2201422753, 2308291338, 2409130253, 2498457507, 2577732249, 2654349402, 2718661053, 2781597017, 2829041340, 2879232772, 3034646526, NA, NA, NA))
      
      new_pos <- integer(length = manhattanN)
      for (mi in 1:24 ) new_pos <- ifelse(man_set$CHR==chr_size$chromosome[mi], chr_size[mi, 3] + man_set$POS, new_pos)
      
      use_Y	<- any(man_set$CHR == 24L)
      use_XY <- any(man_set$CHR == 25L)
      use_M	<- any(man_set$CHR == 26L)
      man_label <- c(1:22, "X")
      at_label	<- chr_size$start[2:24] - 250000 - ( chr_size$start[2:24] - chr_size$start[1:23] ) / 2
      
      if(use_Y) {
        man_label <- c(man_label, "Y")
        at_label <- c(at_label, chr_size$start[25] - 250000 - ( chr_size$start[25] - chr_size$start[24] ) / 2)
        chr_size$start[25] <- chr_size$start[24] + chr_size$size[24] + 500000
      } else { chr_size$start[25] <- chr_size$start[24] }
      
      if(use_XY) {
        man_label <- c(man_label, "XY")
        at_label <- c(at_label, chr_size$start[26] - 250000 - ( chr_size$start[26] - chr_size$start[25] ) / 2)
        chr_size$size[25] <- max(man_set$POS[man_set$CHR == 25])
        chr_size$start[26] <- chr_size$start[25] + chr_size$size[25] + 500000
      } else { chr_size$start[26] <- chr_size$start[25] }
      
      if(use_M ) {
        man_label <- c(man_label, "M")
        at_label <- c(at_label, chr_size$start[27] - 250000 - ( chr_size$start[27] - chr_size$start[26] ) / 2)
        chr_size$size[26] <- max(man_set$POS[man_set$CHR == 26])
        chr_size$start[27] <- chr_size$start[26] + chr_size$size[26] + 500000
      } else { chr_size$start[27] <- chr_size$start[26] }
      
      manMax <- ceiling(-log10(min(man_set$P_VAL)))
      if(manMax < 10L) manMax <- 10L
      
      png(paste0(save_name, ".graph_M.png"),
          width = 960, height = 480)
      par(mfrow = c(1,1), mgp = c(2.5, 0.9, 0))
      palette(c("red", "green3", "blue", "cyan"))
      plot(new_pos, -log10(man_set$P_VAL), pch = 20, xaxs = "i", xaxt = "n", ylim = c(1, manMax),
           xlim = c(0, chr_size$start[27]), col = man_set$CHR, cex.lab = 1.8, cex.axis = 1.4,
           xlab = "Chromosome", ylab = "Observed -log10(p-value)", main = "Manhattan plot", sub = save_name, cex.sub = 1.3, col.sub = "red")
      abline(v = chr_size[2:26, 3] - 250000, lty = 2)
      abline(h = -log10(0.05/sum(QQ_incl)),  #  -log10(5e-8),
             lty = 3, col="red")
      axis(1, at = at_label, labels = man_label, cex.axis = 1.5)
      dev.off()
    } else warning("Insufficient SNPs remain to create Manhattan plot")
  } else {
    if(plot_Man) warning("Insufficient SNPs remain to create Manhattan plot")
    #manhattanN	<- NA
  }
  return(invisible(lambda))
}

# Function translate_header has the standard argument altered & the identify_column function made internal
translate_header <- function(header, standard = c("PROBEID","BETA","SE","P_VAL"), alternative){
  if(any(duplicated(alternative[ ,2]))) stop("duplicated elements in alternative, column 2")
  capitalized <- toupper(header)
  unknowns  <- !logical(length = length(header))
  missings	<- logical(length = length(standard))
  identify_column <- function(std_name, alt_names, header) { return(which(header %in% c(alt_names[ which(alt_names[ ,1]==std_name), 2], std_name))) }
  for(forI in 1:length(standard)) {
    column_no <- identify_column(std_name = standard[forI], alt_names = alternative, header = capitalized)
    if(length(column_no) == 0L) {
      missings[forI] <- TRUE
    } else {
      header[column_no] <- standard[forI]
      unknowns[column_no] <- FALSE
    }
  }	
  return(list(header_N = length(header), header_h = header,
              missing_N = sum(missings), missing_h = if(sum(missings) == 0L) NULL else standard[missings],
              unknown_N = sum(unknowns), unknown_h = if(sum(unknowns) == 0L) NULL else header[unknowns]	) )
}

EWAS_QC <- function(data, # datatable with ewas results, or filename of the same
                    map, # datatable with chr & pos of CpGs in idata or filename of the same
                    outputname, # characterstring to identify dataset in output
                    header_translations, # translation table for file headers
                    threshold_outliers = c(NA, NA),
                    exclude_outliers = FALSE,
                    exclude_X = FALSE,
                    exclude_Y = FALSE,
                    save_final_dataset = TRUE,
                    gzip_final_dataset = TRUE,
                    header_final_dataset = "standard",
                    return_beta = FALSE,
                    N_return_beta = 500000L,
                    ...){
  ### Stage 0: checking input
  # 0: Prepation: creating functions, checking params & output template
  zv_startime <- date()
  
  use_map <- !missing(map)
  if(use_map) { use_map <- !is.null(map) }
  zf_emergencyExit <- function(value, message){ # removed the warning to avoid duplicates (and it's less than informative) since the name of zf_EmergencyExit is returned
    # The function uses "file" from value to determine filename
    # it returns 'value' so that it can be passed to the return that aborts EWAS_QC
    #warning(paste0("Error in file ", value$file, " : ", message))
    print(paste0("Error in file ", value$file, " : ", message), quote = FALSE)
    return(value)
  }

  outcome_QC <- integer(8L)
  names(outcome_QC) <- c("input", "unusable", "chrX", "chrY", "below_lower_limit", "above_upper_limit", "removed", "final")
  outcome_QC <- list(data_input = "",
                     file = "",
                     QC_success = FALSE,
                     lambda = -1,
                     p_cor = -1,
                     N = outcome_QC,
                     SE_median = -1,
                     mean_methylation = NULL,
                     effect_size = NULL)
  
  
  # 0: Checking exclude & threshold outliers
  if(!zf_testLogical(exclude_outliers)) stop("'exclude_outliers' is not a single logical value")
  if(!(is.vector(threshold_outliers) & length(threshold_outliers) == 2L)) {
    stop("'threshold_outliers' is not a numeric vector of length 2") }
  
  use_outliers <- !all(is.na(threshold_outliers))
  if(use_outliers){
    if(threshold_outliers[1] > threshold_outliers[2] & !is.na(threshold_outliers[1]) & !is.na(threshold_outliers[2])) {
      threshold_outliers <- threshold_outliers[2:1]
    }
  } else {
    if(exclude_outliers) stop("'exclude_outliers' is set to TRUE, but no thresholds have been specified")
  }
  
  # 0: Checking exclude_X & exclude_Y
  if(!zf_testLogical(exclude_X)) stop("'exclude_X' is not a single logical value")
  if(!zf_testLogical(exclude_Y)) stop("'exclude_Y' is not a single logical value")
  if((exclude_X | exclude_Y) & !use_map) stop("No map has been specified to exclude X or Y chromosome markers")

  # 0: Checking save_final_dataset & return_beta
  if(!zf_testLogical(save_final_dataset)) stop("'save_final_dataset' is not a single logical value")
  if(!zf_testLogical(return_beta)) stop("'return_beta' is not a single logical value")
  if(return_beta) stopifnot(is.integer(N_return_beta), is.vector(N_return_beta), length(N_return_beta) == 1L, N_return_beta > 50L)
  
  # 0: checking outputname & gzip_final_dataset
  if(!(is.character(outputname) & length(outputname) == 1L)) stop("'outputname' is not a valid filename")
  if(nchar(outputname) == 0L) stop("'outputname' is not a valid filename")
  outcome_QC$file <- outputname # will be changed if save & gzip are TRUE
  
  if(save_final_dataset){
    if(!zf_testLogical(gzip_final_dataset)) stop("'gzip_final_dataset' is not a single logical value")
    outcome_QC$file <- if(gzip_final_dataset) paste0(outputname, ".txt.gz") else paste0(outputname, ".txt")
    if(file.exists(outcome_QC$file)) stop(paste0("Cannot save output file: ", outcome_QC$file, " already exists"))
    
    if(length(header_final_dataset) == 1L) {
      stopifnot(is.character(header_final_dataset), nchar(header_final_dataset) > 2L)
      if(header_final_dataset %in% c("original", "standard")) {
        use_outheader <- header_final_dataset
      } else {
        if(!file.exists(header_final_dataset)) stop("'header_final_dataset' is not a standard name nor a filename in the working directory")
        header_final_dataset <- read.table(header_final_dataset, header = F, stringsAsFactors = FALSE)
        use_outheader <- "table"
      }
    } else {
      if(!is.matrix(header_final_dataset) & !is.data.frame(header_final_dataset)) stop("'header_final_dataset' is not a table, filename or standard name")
      use_outheader <- "table"
    }
    
    if(use_outheader == "table"){
      if(is.matrix(header_final_dataset) | is.data.frame(header_final_dataset)) {
        if(ncol(header_final_dataset) != 2L) { stop("'header_final_dataset' does not have 2 columns") }
        if(any(is.na(header_final_dataset))) { stop("'header_final_dataset' contains missing values") }
        if(!is.character(header_final_dataset[,1])) { stop("'header_final_dataset' contains non-character values") }
        if(!is.character(header_final_dataset[,2])) { stop("'header_final_dataset' contains non-character values") }
        if(any(duplicated(header_final_dataset[ ,1]),duplicated(header_final_dataset[ ,2]))) { stop("'header_final_dataset' contains duplicated names") }
      } else {
        stop("'header_final_dataset' is not a table, filename or standard name")
      }
    }
  } else { gzip_final_dataset <- FALSE } # this is necessary for log
  
  # 0: checking header translation table
  use_translation <- !missing(header_translations)
  if(use_translation) { use_translation <- !is.null(header_translations) }
  if(use_translation) {
    if(length(header_translations) == 1L) {
      stopifnot(is.character(header_translations), nchar(header_translations) > 2L, file.exists(header_translations))
      header_translations <- read.table(header_translations, stringsAsFactors = FALSE)
    }
    if(!is.data.frame(header_translations)) stop("'header_translations' is not dataframe or a filename leading to such")
    if(ncol(header_translations) != 2L) stop("'header_translations' does not have two columns")
    if(any(duplicated(header_translations[ ,2]))) stop("'header_translations' contains duplicated elements in column 2")
  }
  
  
  #### Stage 1: loading data & map & checking column names
  
  print("", quote = FALSE)
  print(paste0(" *** QC'ing ", outputname), quote = FALSE)
  print("", quote = FALSE)
  flush.console()
  
  # 1: map
  if(use_map){
    if(is.character(map) & length(map) == 1L) {
      if(!file.exists(map)) stop(paste0("Cannot find file: ", map))
      map <- read.table(map, header = T, stringsAsFactors = FALSE, comment.char = "")
    }
    if(!is.data.frame(map)) stop("'map' is not a dataframe or a filename leading to such")
    
    colnames(map) <- toupper(colnames(map))
    
    if(!all(c("PROBEID", "CHR", "POS") %in% colnames(map))) {
      stop(paste0("map does not have columns:  ",
                  paste0(c("PROBEID", "CHR", "POS")[!c("PROBEID", "CHR", "POS") %in% colnames(map)],
                         collapse = ", ")))
    }
    # We cannot allow duplicate or missing probeIDs to be present in map,
    # because otherwise temp_map won't be of equal length to data.
    if(any(is.na(map$PROBEID))) stop("'map' contains missing probeIDs")
    if(any(duplicated(map$PROBEID))) stop("'map' contains duplicate probeIDs")
  }
  
  # 1: data
  if(is.character(data) & length(data) == 1L) {
    if(!file.exists(data)) stop(paste0("Cannot find file: ", data))
    outcome_QC$data_input <- data
    data <- try(read.table(data, header = T, stringsAsFactors = FALSE, ...))
    if(class(data) == "try-error") return(zf_emergencyExit(value = outcome_QC,
                                                           message = paste0("cannot load file ", outcome_QC$data_input, "; check warnings() for details")))
  } else { outcome_QC$data_input <- "user-supplied table" }
  if(!is.data.frame(data)) stop("'data' is not a dataframe or a filename leading to such")
  # The above error won't occur for EWAS_series, as that loads all data via read.table
  # hence we use stop isntead of the return function
  
  
  header_orig <- colnames(data)
  header_std <- c("PROBEID", "BETA", "SE", "P_VAL")
  if(use_translation){
    header_info <- translate_header(header = toupper(header_orig),
                                    standard = header_std,
                                    alternative = header_translations)
    if(any(duplicated(header_info$header_h))) {
      return(zf_emergencyExit(value = outcome_QC,
                              message = paste0("Duplicate columns: ", paste0(header_info$header_h[duplicated(header_info$header_h)], collapse = ", ") )))
    }
    colnames(data) <- header_info$header_h
  } else { colnames(data) <- toupper(header_orig) }
  
  if(!all(header_std %in% colnames(data))) {
    return(zf_emergencyExit(value = outcome_QC,
                            message = paste0("Columns: ",
                                             paste0(header_std[!header_std %in% colnames(data)], collapse = ", "),
                                             " are missing from data. Data does contain (but could not identify) column(s) ", paste0(colnames(data)[!colnames(data) %in% header_std], collapse = ", "), "."
                                             )))
  }
  
  
  #### Stage 2: data integrity
  outcome_QC$N["input"] <- nrow(data)
  if(outcome_QC$N["input"] == 0L) return(zf_emergencyExit(value = outcome_QC,
                                                          message = "'data' contains no entries"))
  vec_unusable <- logical(outcome_QC$N["input"])
  
  # 2: probeID - checking for duplicates
  vec_NA_ID <- is.na(data$PROBEID)
  N_NA_ID <- sum(vec_NA_ID)
  if(N_NA_ID > 0L) { vec_unusable <- vec_unusable | vec_NA_ID }
  
  vec_dupli <- duplicated(data$PROBEID)
  N_dupli <- sum(vec_dupli)
  if(N_dupli > 0L) {
    temp_duplis <- data$PROBEID[vec_dupli]
    N_dupli <- sum(data$PROBEID %in% temp_duplis)
    rm(temp_duplis)
    warning(paste0(N_dupli, " entries with duplicated IDs in file: ", outputname, "."))
  }

  # 2: Beta
  vec_NA_beta <- is.na(data$BETA)
  N_NA_beta <- sum(vec_NA_beta)
  if(N_NA_beta == outcome_QC$N["input"]) return(zf_emergencyExit(value = outcome_QC,
                                                                 message = "'data' contains no BETA values"))
  if(!is.numeric(data$BETA)) return(zf_emergencyExit(value = outcome_QC,
                                                     message = "BETA column inside 'data' does not contain numeric values"))
  if(N_NA_beta > 0L) { vec_unusable <- vec_unusable | vec_NA_beta }
  
  # 2: SE - checking for negative values
  vec_temp <- is.na(data$SE)
  N_NA_SE <- sum(vec_temp)
  if(N_NA_SE == outcome_QC$N["input"]) return(zf_emergencyExit(value = outcome_QC,
                                                               message = "'data' contains no SE values"))
  if(!is.numeric(data$SE)) return(zf_emergencyExit(value = outcome_QC,
                                                   message = "SE column inside 'data' does not contain numeric values"))
  if(N_NA_SE > 0L) { vec_unusable <- vec_unusable | vec_temp }
  
  vec_temp <- data$SE <= 0 & !vec_temp  # vec_temp represent missing values up until this point
  N_bad_SE <- sum(vec_temp)
  if(N_bad_SE > 0L) {
    warning(paste0(N_bad_SE, " entries with a negative or zero standard-error in file: ", outputname, ". These values have been set to NA."))
    data$SE[vec_temp] <- NA
    vec_unusable <- vec_unusable | vec_temp
  }
  
  # 2: Pval - checking for missing, bad en p = 0
  vec_temp <- is.na(data$P_VAL)
  N_NA_p <- sum(vec_temp)
  if(N_NA_p == outcome_QC$N["input"]) return(zf_emergencyExit(value = outcome_QC,
                                                              message = "'data' contains no p values"))
  if(!is.numeric(data$P_VAL)) return(zf_emergencyExit(value = outcome_QC,
                                                      message = "P_VAL column inside 'data' does not contain numeric values"))
  if(N_NA_p > 0L) { vec_unusable <- vec_unusable | vec_temp }
  
  vec_temp <- (data$P_VAL < 0 | data$P_VAL > 1) & !vec_temp  # vec_temp represent missing values up until this point
  N_bad_p <- sum(vec_temp)
  if(N_bad_p > 0L) {
    warning(paste0(N_bad_p, " entries with p-values outside or the range 0 - 1 in file: ", outputname, ". These values have been set to NA."))
    data$P_VAL[vec_temp] <- NA
    vec_unusable <- vec_unusable | vec_temp
  }

  N_low_p <- sum(data$P_VAL < 1e-300 & !is.na(data$P_VAL))
  if(N_low_p > 0L){
    N_zero_p <- sum(data$P_VAL == 0 & !is.na(data$P_VAL))
    if(N_zero_p > 0L) warning(paste0(N_zero_p, " entries with p = 0 in file: ", outputname, ". This indicates a rounding error."))
    warning(paste0(N_low_p, " entries with p < 1e-300 in file: ", outputname, ". These values will be temporarily set to 1e-300 for graph creation."))
  } else { N_zero_p <- 0L }
  
  # 2: final check
  rm(vec_temp)
  outcome_QC$N["unusable"] <- sum(vec_unusable)
  if(outcome_QC$N["unusable"] == outcome_QC$N["input"]) return(zf_emergencyExit(value = outcome_QC,
                                                                  message = "All entries in 'data' contain missing or invalid values"))
  #print(paste0(if(outcome_QC$N["unusable"] == 0L) "No" else outcome_QC$N["unusable"], " entries with missing values"), quote = FALSE)
  print(paste0(outcome_QC$N["unusable"], " entries with missing values"), quote = FALSE)

  
  #### Stage 3: data cleaning
  # generate matching map for use in Manhattan plot

  if(use_map) {
    N_map <- sum(data$PROBEID %in% map$PROBEID)
    if(N_map == 0L){
      use_map <- FALSE
      exclude_X <- FALSE
      exclude_Y <- FALSE
      warning(paste0("Cannot matchs probeIDs in file ", outputname, " with map. Check map or file probeIDs for problems."))
    } else {
      if(N_map < (outcome_QC$N["input"]-N_NA_ID) ) { warning(paste0((outcome_QC$N["input"] - N_NA_ID) - N_map, " entries in file ", outputname, " do not appear in the map.")) }
      
      temp_map <- data.frame(data[, c("PROBEID", "P_VAL")],
                             order = 1:outcome_QC$N["input"],
                             stringsAsFactors = FALSE)
      temp_map <- merge(temp_map, map[, c("PROBEID", "CHR", "POS")], by = "PROBEID",
                        all.x = TRUE, all.y = FALSE, sort = FALSE)
      temp_map <- temp_map[order(temp_map$order), -3]
      if(!all(temp_map$PROBEID == data$PROBEID)) stop("DEBUG ERROR - temp_map does not match data-probeIDs")
      temp_map$CHR <- ifelse(is.na(temp_map$CHR), -9, temp_map$CHR) # so we don't have to check for missing values in the following logic test
    }
    rm(map)
  } else { N_map <- NA }
  
  # 3: counting X & Y chr
  vec_remove <- logical(outcome_QC$N["input"])
  
  if(use_map){ # the use_map test is so that use_map can be disabled above
    vec_temp <- temp_map$CHR == "X" |  temp_map$CHR == 23
    outcome_QC$N["chrX"] <- sum(vec_temp)
    if(outcome_QC$N["chrX"] > 0L & exclude_X) vec_remove <- vec_remove | vec_temp
    
    vec_temp <- temp_map$CHR == "Y" |  temp_map$CHR == 24
    outcome_QC$N["chrY"] <- sum(vec_temp)
    if(outcome_QC$N["chrY"] > 0L & exclude_Y) vec_remove <- vec_remove| vec_temp
  } else {
    outcome_QC$N["chrX"] <- NA
    outcome_QC$N["chrY"] <- NA
  }
  
  # 3: counting extremes
  if(is.na(threshold_outliers[1])){
    outcome_QC$N["below_lower_limit"] <- 0L
  } else {
    vec_temp <- data$BETA < threshold_outliers[1] & !vec_NA_beta & !vec_remove # this excludes chr X and Y, if they have been excluded
    outcome_QC$N["below_lower_limit"] <- sum(vec_temp)
    if(outcome_QC$N["below_lower_limit"] > 0L & exclude_outliers) vec_remove <- vec_remove| vec_temp
  }
  
  if(is.na(threshold_outliers[2])){
    outcome_QC$N["above_upper_limit"] <- 0L
  } else {
    vec_temp <- data$BETA > threshold_outliers[2] & !vec_NA_beta & !vec_remove # this excludes chr X and Y and low_outliers, which don't matter
    outcome_QC$N["above_upper_limit"] <- sum(vec_temp)
    if(outcome_QC$N["above_upper_limit"] > 0L & exclude_outliers) vec_remove <- vec_remove| vec_temp
  }
  
  
  # 3: removing extremes & X & Y
  if(use_map | use_outliers) rm(vec_temp)
  outcome_QC$N["removed"] <- sum(vec_remove)
  
  
  if(use_outliers){
    print(paste0(outcome_QC$N["below_lower_limit"] + outcome_QC$N["above_upper_limit"], " outliers ", if(exclude_outliers) "excluded" else "found"), quote = FALSE)
  }
  if(exclude_X) print(paste0(outcome_QC$N["chrX"], " X-chromosome probes excluded"), quote = FALSE)
  if(exclude_Y) print(paste0(outcome_QC$N["chrY"], " Y-chromosome probes excluded"), quote = FALSE)
  if(sum(exclude_outliers, exclude_X, exclude_Y) > 1L) print(paste0(outcome_QC$N["removed"], " total probes excluded"), quote = FALSE)
  flush.console()
  
  if(outcome_QC$N["removed"] > 0L){
    if(outcome_QC$N["removed"] == outcome_QC$N["input"]) return(zf_emergencyExit(value = outcome_QC,
                                                                                 message = "All entries have been removed from 'data'"))
    
    data <- data[!vec_remove, ]
    if(use_map) { temp_map <- temp_map[!vec_remove, ] }
    vec_unusable <- vec_unusable[!vec_remove]
    vec_NA_ID <- vec_NA_ID[!vec_remove]
    vec_dupli <- vec_dupli[!vec_remove]
    vec_NA_beta <- vec_NA_beta[!vec_remove]
    
    if(all(vec_unusable)) return(zf_emergencyExit(value = outcome_QC,
                                                  message = "All entries in 'data' contain missing or invalid values"))
  }
  rm(vec_remove)
  outcome_QC$N["final"] <- nrow(data)
  
  
  #### Phase 4: generating plots
  
  # 4: histograms
  png(paste0(outputname, ".histo.png"),  width = 960, height = 480, res = 108)
  par(mfrow = c(1, 2))
  
  hist(data$BETA,
       breaks = 30, freq = FALSE, col = "blue", plot = TRUE,
       main = "Beta", xlab = "Beta", cex.main = 1.5)
  
  hist(data$SE,
       breaks = 30, freq = FALSE, col = "red", plot = TRUE,
       main = "Standard error", xlab = "Standard error", cex.main = 1.5,
       sub = outputname, cex.sub = 1.3)
  
  dev.off()
  
  # 4: P correlation plot
  outcome_QC$p_cor <- P_correlation(dataset = data,
                              plot_correlation = TRUE, plot_if_threshold = FALSE, threshold_r = 0.99,
                              save_name = outputname)
  
  # 4: Manhattan & QQ plots plot - note that we run this even if there is no map in order to calculate lamda
  outcome_QC$lambda <- EWAS_plots(dataset = if(use_map) temp_map else data$P_VAL,
                                  plot_QQ = TRUE, plot_Man = use_map, plot_cutoff_p = 0.05,
                                  plot_QQ_bands = FALSE,
                                  save_name = outputname)
  if(use_map) rm(temp_map)
  
  
  # 4: Volcano plot - Effect/Beta op X-as, y-as = -log10 p-waarde
  temp_p <- if(N_low_p > 0L) ifelse(data$P_VAL < 1e-300, 1e-300, data$P_VAL) else data$P_VAL
  png(paste0(outputname, ".graph_volcano.png"),
      width = 720, height = 720, res = 144)
  plot(data$BETA, -log10(temp_p),
       main = "Volcano plot", sub = outputname, cex.sub = 1.3,
       xlab = "Effect Size", ylab = "-log10(p value)")
  dev.off()
  rm(temp_p)
  
  # 4: Effect size against methylation beta (beta vs. beta-value) - ook even in de ijskast
  

  
  #### 5 : generate output
  # 5: update outcome_QC - this needs to occur before creation of the log file
  
  outcome_QC$SE_median <- median(data$SE, na.rm = TRUE)
  #outcome_QC$mean_methylation <- NULL # not implemented yet
  if(return_beta) {
    #N_return_beta <- 50000L
    N_final_beta <- sum(!vec_NA_beta)
    if(N_final_beta < N_return_beta) {
      outcome_QC$effect_size <- c(data$BETA[!vec_NA_beta], rep(NA, N_return_beta - N_final_beta))
    } else {
      outcome_QC$effect_size <- sort(data$BETA[!vec_NA_beta])[c(1:N_return_beta)*(N_final_beta %/% N_return_beta)]
    }
  }
  outcome_QC$QC_success <- TRUE
  
  
  # 5: create log file - this needs to occur before the file headers are changed
  
  print("", quote = FALSE)
  print("", quote = FALSE)
  write.table(c(
    "****************************************",
    "\t\tEWAS_QC LOG FILE",
    "****************************************",
    ""), paste0(outputname, ".log"), append = FALSE,
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  zc_log <- file(description = paste0(outputname, ".log"), open = "a")
  
  write.table(data.frame(
    v1 = c("File", "QC Start Time", "QC End Time", "EWAS_QC Version"),
    v2 = "",
    v3 = c(outcome_QC$file, zv_startime, date(), "1.0-7"),
    stringsAsFactors = FALSE),
    zc_log, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(c("",
                "",
                "****************************************",
                "\t1. Number of entries in file",
                "****************************************",
                ""), zc_log, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  write.table(t(c("", "", "N", "", "", "", "N")),
              zc_log, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(
    v1 = c("Start", "ID Missing", "ID Duplicated", "Beta Missing", "SE Missing", "SE Invalid", "P Missing", "P Invalid", "Total missing & invalid"),
    v2 = "",
    v3 = c(outcome_QC$N["input"], N_NA_ID, N_dupli, N_NA_beta, N_NA_SE, N_bad_SE, N_NA_p, N_bad_p, outcome_QC$N["unusable"]),
    v4 = "",
    v5 = c("Start", "Not found in map", "Found in map", "> X", "> Y", "Outliers High", "Outliers Low", "Total Removed", "Final"),
    v6 = "",
    v7 = c(outcome_QC$N["input"], outcome_QC$N["input"] - N_map, N_map, outcome_QC$N["chrX"], outcome_QC$N["chrY"], outcome_QC$N["below_lower_limit"], outcome_QC$N["above_upper_limit"], outcome_QC$N["removed"], outcome_QC$N["final"]),
    stringsAsFactors = FALSE),
    zc_log, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(c("",
                "",
                "****************************************",
                "\t2. Summary statistics",
                "****************************************",
                ""), zc_log, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  write.table(data.frame(
    v1 = c("p-value correlation", "Lambda", "SE median"),
    v2 = "",
    v3 = c(outcome_QC$p_cor, outcome_QC$lambda, outcome_QC$SE_median),
    stringsAsFactors = FALSE),
    zc_log, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)
  
  zf_saveStats <- function(stat, statname){ # we have to add name, otherwise c() will convert all integers to numeric
    Nall <- length(stat)
    NNA <- sum(is.na(stat))
    PNA <- 100*NNA/Nall
    if(PNA == 100) return(c(statname, Nall - NNA, NNA, round(PNA, digits = 1), rep(NA, 6)))
    #return(c(statname, Nall - NNA, NNA, round(PNA, digits = 1),
    #         round(c(mean(stat, na.rm = TRUE), quantile(stat, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE))[c(2,3,4,1,5,6)], digits = 2)))
    return(c(statname, Nall - NNA, NNA, round(PNA, digits = 1),
             c(mean(stat, na.rm = TRUE), quantile(stat, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE))[c(2,3,4,1,5,6)]))
  }
  outcome_stats <- data.frame(stat = c("Beta", "StdErr", "P"), Nfinal = integer(3), NNA = integer(3), PNA = integer(3),
                              min = numeric(3), Q25 = numeric(3), median = numeric(3), mean = numeric(3), Q75 = numeric(3), max = numeric(3),
                              stringsAsFactors = FALSE)
  outcome_stats[1, ] <- zf_saveStats(stat = data$BETA , statname = "Beta")
  outcome_stats[2, ] <- zf_saveStats(stat = data$SE   , statname = "StdErr")
  outcome_stats[3, ] <- zf_saveStats(stat = data$P_VAL, statname = "P")
  
  write.table("", zc_log, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(t(c("" ,"N final", "N missing", "% missing", "min", "25%", "median", "mean", "75%", "max")),
              zc_log, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(outcome_stats,
              zc_log, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  close(zc_log)

  # Saving final dataset

  if(save_final_dataset) {
    if(use_outheader == "original") colnames(data) <- header_orig
    if(use_outheader == "table") {
      temp_header <- translate_header(header = colnames(data), standard = header_final_dataset[ ,1], alternative = header_final_dataset)
      colnames(data) <- temp_header$header_h
    }
    write.table(data,
                if(gzip_final_dataset) gzfile(outcome_QC$file) else outcome_QC$file,
                row.names = FALSE, quote = FALSE)
  }

  rm(data, vec_unusable, vec_NA_ID, vec_dupli, vec_NA_beta)
  gc(FALSE, FALSE)
  return(outcome_QC)
}

EWAS_series <- function(EWAS_files, # datatable with ewas results, or filename of the same
                        output_files,
                        map, # datatable with chr & pos of CpGs in idata or filename of the same
                        N, # sample sizes for precision plot
                        header_translations, # translation table for file headers
                        save_final_dataset = TRUE,
                        gzip_final_dataset = TRUE,
                        ...){
  # all other arguments will be checked by the first EWAS_QC, so no need to check 'em here
  zf_loadAndReturn <- function(filename, filedesc, ...){
    
    if(is.data.frame(filename) | is.matrix(filename)) return(filename)
    
    if(!is.character(filename) | length(filename) != 1L) stop(paste0("argument '", filedesc, "' is neither a dataset nor a filename"))
    
    if(!file.exists(filename)) stop(paste0("Cannot find file '", filename, "' in the working directory"))
    
    return(read.table(file = filename, stringsAsFactors = FALSE, ...))
  }
  
  use_map <- !missing(map)
  if(use_map) { use_map <- !is.null(map) }
  use_translation <- !missing(header_translations)
  if(use_translation) { use_translation <- !is.null(header_translations) }
  use_N <- !missing(N)
  if(use_N) { use_N <- !is.null(N) }
  
  # checking input & output filenames, as well as save commands
  stopifnot(is.vector(EWAS_files), is.character(EWAS_files))
  if(any(nchar(EWAS_files) == 0L)) stop("'EWAS_files' contains invalid filenames")
  if(any(duplicated(EWAS_files))) stop("duplicate filenames in 'EWAS_files'")
  if(any(!file.exists(EWAS_files))) stop(paste0("files: ", paste0(EWAS_files[!file.exists(EWAS_files)], collapse = ", "), " do not exist"))
  N_file <- length(EWAS_files)
  
  
  if(missing(output_files)) { output_files <- paste0("QC_", EWAS_files)
  } else {
    stopifnot(is.vector(output_files), is.character(output_files))
    if(N_file != length(output_files)) stop("vectors 'EWAS_files' (input) and 'output_files' are of unequal length")
  }
  # removing .gz
  output_files <- ifelse(substr(output_files, nchar(output_files) - 2, nchar(output_files)) == ".gz",
                         substr(output_files, 1, nchar(output_files) - 3),
                         output_files)
  if(any(nchar(output_files) == 0L)) stop("'output_files' contains valid filename")
  if(any(duplicated(output_files))) stop("duplicate filenames in 'output_files'")
  
  if(save_final_dataset){
    if(!zf_testLogical(gzip_final_dataset)) stop("'gzip_final_dataset' is not a single logical value")
    #if(!zf_testLogical(save_standard_header)) stop("'save_standard_header' is not a single logical value")
    if(gzip_final_dataset){
      if(any(file.exists(paste0(output_files, ".gz")))) stop("cannot save output: one of the specified 'output_files' already exists")
    } else {
      if(any(file.exists(output_files))) stop("cannot save output: one of the specified 'output_files' already exists")
    }
  }
  
  # checking N & preparing precision plot
  dat_prec <- data.frame(file = EWAS_files, no = 1:N_file,
                         use = logical(N_file), SE = numeric(N_file),
                         stringsAsFactors = FALSE)
  if(use_N){
    if(is.character(N) & length(N) == 1L) { N <- zf_loadAndReturn(filename = N, filedesc = "N", header = TRUE) }
    if(!is.data.frame(N)) stop("'N' is not dataframe or a filename leading to such")
    if(any(!colnames(N) %in% c("file", "N"))) stop("'N' does not contain the specified columns 'file' and/or 'N'")
    if(any(!EWAS_files %in% N$file)) stop("Cannot find matching entries in 'N' for all files in 'EWAS_files'")
    
    dat_prec <- merge(dat_prec, N[ , c("file", "N")], by = "file", all.x = TRUE, all.y = FALSE, sort = FALSE)
    rm(N)
    dat_prec <- dat_prec[order(dat_prec$no), ]
    
    if(!is.numeric(dat_prec$N)) stop("column N in table 'N' does not contain numeric values")
    if(any(dat_prec$N < 1, na.rm = TRUE)) stop("column N in table 'N' does not contain valid numeric values")
    if(sum(!is.na(dat_prec$N)) < 2) {
      warning("insufficient non-missing N values to generate precision plot")
      use_N <- FALSE
    }
  }
  
  # Checking header translation table
  if(use_translation) {
    header_translations <- zf_loadAndReturn(filename = header_translations, filedesc = "header_translations")
    #if(!is.data.frame(header_translations)) stop("'header_translations' is not dataframe or a filename leading to such")
    if(ncol(header_translations) != 2L) stop("'header_translations' does not have two columns")
    if(any(duplicated(header_translations[ ,2]))) stop("'header_translations' contains duplicated elements in column 2")
  } else { header_translations <- NULL }
  
  # Checking map
  if(use_map){
    map <- zf_loadAndReturn(filename = map, filedesc = "map", header = TRUE)
    #if(!is.data.frame(map)) stop("'map' is not a dataframe or a filename leading to such")
    
    colnames(map) <- toupper(colnames(map))
    if(!all(c("PROBEID", "CHR", "POS") %in% colnames(map))) {
      stop(paste0("map does not have columns:  ",
                  paste0(c("PROBEID", "CHR", "POS")[!c("PROBEID", "CHR", "POS") %in% colnames(map)],
                         collapse = ", ")))
    }
    # We cannot allow duplicate or missing probeIDs to be present in map,
    # because otherwise temp_map won't be of equal length to data.
    if(any(is.na(map$PROBEID))) stop("'map' contains missing probeIDs")
    if(any(duplicated(map$PROBEID))) stop("'map' contains duplicate probeIDs")
  } else { map <- NULL }
  
  
  #### Starting the actual QC
  
  N_beta <- 500000L
  dat_beta <- matrix(data = numeric(N_beta * N_file), nrow = N_beta, ncol = N_file)
  for(ci in 1:N_file){
    #print(paste0(" *** QC'ing file ", ci, " of ", N_file), quote = FALSE)
    #flush.console() # unnecessary as EWAS_QC will already flush
    temp_out <- EWAS_QC(data = EWAS_files[ci],
                        map = map,
                        outputname = output_files[ci],
                        header_translations = header_translations,
                        save_final_dataset = save_final_dataset, gzip_final_dataset = gzip_final_dataset,
                        return_beta = TRUE, N_return_beta = N_beta, ...)
    if(temp_out$QC_success){
      dat_prec$use[ci] <- TRUE
      dat_prec$SE[ci] <- temp_out$SE_median
      # mean methylation is not used, for the moment
      dat_beta[, ci] <- temp_out$effect_size
    }
    rm(temp_out)
  }
  print("", quote = FALSE)
  
  # output:
  N_pass <- sum(dat_prec$use)
  if(N_pass > 1L){
    #Precision plot
    dat_temp <- dat_prec[dat_prec$use, ]
    
    if(use_N){
      temp_prec <- 1/dat_temp$SE
      temp_N <- sqrt(dat_temp$N)
      temp_ticks <- pretty(temp_prec) # using pretty() will hopefully ensure neat values on both y-axes
      
      png("EWAS_QC_graph_precision.png")
      par(mar = c(5, 4, 4, 4.5) + 0.1) # add +2 to right margin, so to leave space for label
      plot(x = temp_N, y = temp_prec, yaxt = "n", # suppresses the y-axis
           main = "Precision by Sample Size", xlab = "sqrt(sample size)", ylab = "Precision ( = 1 / median SE)")
      axis(side = 2, at = temp_ticks, labels = TRUE)
      axis(side = 4, at = temp_ticks, labels = 1/temp_ticks)
      mtext(text = "median SE", side = 4, line = 3) # adds label to right margin
      text(temp_N, temp_prec, labels = dat_temp$no, pos = 4)
      dev.off()
      rm(temp_N, temp_prec, temp_ticks)
    }
    
    #mean_methylation density plot skipped
    
    #effect size boxplot
    
    dat_beta <- dat_beta[ , dat_prec$use]
    if(use_N){
      dat_beta <- dat_beta[ , order(dat_temp$N, na.last = FALSE)]
      dat_temp <- dat_temp[order(dat_temp$N, na.last = FALSE), ]
    }
    
    dat_quant <- matrix(data = 0.0, ncol = 3, nrow = N_pass,
                        dimnames = list(NULL, c("quantile25", "median", "quantile75")))
    for(ti in 1:ncol(dat_beta)) dat_quant[ti, ] <- quantile(dat_beta[, ti], names = FALSE, na.rm = TRUE)[2:4]
    
    png("EWAS_QC_graph_effectsize.png", width = 400 + 80 * N_pass, height = 480)
    boxplot(dat_beta, names = if(use_N) paste0(dat_temp$no, ", N=", dat_temp$N) else dat_temp$no,
            las = 0, main = "Effect-size distribution")
    abline(h = c(median(dat_quant[ ,1]), median(dat_quant[ ,2]), median(dat_quant[ ,3])), lty = 3)
    dev.off()
    rm(dat_temp)
  } else { print("Insufficient files passed the QC to generate graphs", quote = FALSE)}
  
  rm(dat_beta, map)
  write.table(dat_prec, "EWAS_QC_legend.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  gc(FALSE, FALSE)
  return(invisible(dat_prec))
}
