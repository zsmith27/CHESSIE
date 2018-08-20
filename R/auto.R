#==============================================================================
#==============================================================================
# Author: Zachary M. Smith
# Maintainer: Zachary M. Smith
# Organization: ICPRB
# Created: December 2016
# Updated: 3-14-2017
# Purpose: Functions to automate the development of the Chessie BIBI.
#==============================================================================
#==============================================================================
#'Automate Spatial Classification
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

auto_bio <- function(my.data, index_res, bioregion = NULL){
  my.data$BIOREGION <- ifelse(my.data$ICPRB_BIOREGION_ID %in% "NGV", "UNP",
                              as.character(my.data$ICPRB_BIOREGION_ID))
  if(index_res %in% "BASIN") my.data$BIOREGION <- "BASIN"
  if(index_res %in% "COAST"){
    #my.data <- my.data[my.data$BIOREGION %in% c("MAC", "SEP"), ]
    my.data <- my.data[my.data$ICPRB_BIOREGION_ID %in% c("MAC", "SEP"), ]
    my.data$BIOREGION <- "COAST"
  }
  if(index_res %in% "INLAND"){
    #my.data <- my.data[!my.data$BIOREGION %in% c("MAC", "SEP"), ]
    my.data <- my.data[my.data$ICPRB_BIOREGION_ID %in% c("BLUE", "CA", "LNP",
                                                         "NAPU", "NCA", "NRV",
                                                         "PIED", "SGV", "SRV",
                                                         "UNP"), ]
    my.data$BIOREGION <- "INLAND"
  }
  if (!index_res %in% c("BASIN", "INLAND", "COAST")) {
    my.data <- my.data[my.data$BIOREGION %in% bioregion, ]
  }
  return(my.data)
}

#==============================================================================
#'Automate Metric Calculation
#'
#'@param all.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = 
#'@param taxon.rank =
#'@param todays_date =
#'@return 
#'@export
#'
auto_metrics <- function(all.data, index.res, bioregion, taxon.rank,
                         todays_date = format(Sys.time(), "%m_%d_%y")){
  #==============================================================================
  # Prepare the environmental data.
  env.df <- prep_env(Prep.data = all.data)
  # Assign the disturbance gradient classes (Reference, Degraded, etc.)
  test.class <- prep_class_multi(env.df)
  # Assign the appropriate spatial resolution name.
  all.data <- auto_bio(all.data, index.res, bioregion)
  #names(all.data)[names(all.data) %in% "ICPRB_BIOREGION_ID"] <- "BIOREGION"
  # Keep only the rows where the BIOREGION column is NOT NA.
  bio.data <- all.data[!is.na(all.data$BIOREGION), ]
  #table(test.class$BIOREGION, test.class$CATEGORY) # Test proper subset.
  #==============================================================================
  # Create a standard data frame to collect output
  env.class <- unique(test.class[, c("EVENT_ID", "STATION_ID", "DATE",
                                     "AGENCY_CODE", "SAMPLE_NUMBER", "CATEGORY")])
  #==============================================================================
  # Calculate Standard Metrics.
  bio.data.metrics <- all_metrics(master, bio.data, taxa.rank =  taxon.rank,
                                  rare = "DIV_RARE", pct_un = 10, bibi.standard = TRUE,
                                  seed = TRUE)
  # Apply the 2016-1017 taxa standardization to the master taxa list. 
  master2 <- clean_taxa(master)
  
  #==============================================================================
  # Calculate the percent of each taxon in the data base by Sequencing through
  # each taxon rank
  #==============================================================================
  # Identify the appropriate columns to replace with NAs.
  # Taxonomic resolution greater than the specified taxon.rank should not be
  # calculated; therefore, replace all values in these columns with NA.
  if(taxon.rank %in% "ORDER") ranks.list <- c("SUBORDER", "FAMILY", "SUBFAMILY",
                                              "TRIBE", "GENUS", "SPECIES")
  if(taxon.rank %in% "FAMILY") ranks.list <- c("SUBFAMILY", "TRIBE", "GENUS",
                                               "SPECIES")
  if(taxon.rank %in% "GENUS") ranks.list <- c("SPECIES")
  master2[, ranks.list] <- NA
  # Sequence through each taxonomic rank and calculate the percentage of each
  # unique taxon.
  bio.data.seq <- seq_pct_taxa(bio.data, master2)
  # Check for duplicate EVENT_IDs.
  test <- bio.data.seq[duplicated(bio.data.metrics$EVENT_ID), ]
  #==============================================================================
  # Merge the Standard Metrics with the Sequenced Metrics.
  bio.data.merge <- merge(bio.data.metrics, bio.data.seq,
                          by = c("EVENT_ID", "STATION_ID", "DATE",
                                 "SAMPLE_NUMBER", "AGENCY_CODE"))
  # Remove any column names containing "PCT_UNIDENTIFIED" or ending with ".y".
  bio.data.merge <- bio.data.merge[, !grepl("PCT_UNIDENTIFIED", names(bio.data.merge))]
  bio.data.merge <- bio.data.merge[, !grepl(".y", names(bio.data.merge))]
  # Remove ".x" from the trailing end of column names.
  names(bio.data.merge) <- gsub(".x", "", names(bio.data.merge))
  #==============================================================================
  # Check Metric Min and Max
  #==============================================================================
  # Create a new data frame with all of the metric names listed in the first column.
  #new <- data.frame(names(bio.data.merge[, 7:ncol(bio.data.merge)]))
  # Identify the min and max of each metric.
  #new$MIN <- round(apply(bio.data.merge[, 7:ncol(bio.data.merge)], 2, min))
  #new$MAX <- round(apply(bio.data.merge[, 7:ncol(bio.data.merge)], 2, max), 0)
  #==============================================================================
  # Merge the Metrics with the Site Classes
  m.bio.data <- merge(env.class, bio.data.merge,
                      by = c("EVENT_ID", "STATION_ID", "DATE",
                             "AGENCY_CODE", "SAMPLE_NUMBER"), all.y = TRUE)
  #==============================================================================
  # REMOVE??? Seems redundant.
  #==============================================================================
  # Double check that all column names containing "PCT_UNIDENTIFIED" or 
  # ending with ".y" have been removed.
  #if("PCT_UNIDENTIFIED"  %in% names(m.bio.data)) m.bio.data <- m.bio.data[, !grepl("PCT_UNIDENTIFIED", names(m.bio.data))]
  #if(any(grepl(".y", names(m.bio.data)))) m.bio.data <- m.bio.data[, !grepl(".y", names(m.bio.data))]
  #if(any(grepl(".x", names(m.bio.data)))) names(m.bio.data) <- gsub(".x", "", names(m.bio.data))
  #m.bio.data$CATEGORY <- factor(m.bio.data$CATEGORY, levels =c("REF", "MIN", "MOD", "SEV", "MIX"))
  #table(m.bio.data$CATEGORY)
  #test <- data.frame(names(m.bio.data[, 7:ncol(m.bio.data)]))
  #==============================================================================
  # Re-assign the the bioregion to the newest data frame.
  merged <- merge(unique(bio.data[, c("EVENT_ID", "BIOREGION")]), m.bio.data,
                  by = "EVENT_ID", all.y = TRUE)
  #write.csv(merged, "Test2_Metrics_11_7_16.csv", row.names = FALSE)
  #==============================================================================
  # Change REF+ to REF and DEG+ to SEV.
  merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF+", "REF", 
                            ifelse(merged$CATEGORY %in% c("DEG+", "DEG"), "MOD", 
                                   ifelse(merged$CATEGORY %in% c("DEG2"), "SEV",
                                          merged$CATEGORY)))
  #==============================================================================
  # Coastal Reference sampling events were specifically selected by C. Buchanan
  # and Z. Smith and must be re-assigned.  All of the EVENT_IDs that were 
  # classified as Reference by the R-scripts but were not in the imported list,
  # were changed to the "MIN" classification (Minimally Degraded).
  #***************Changed 5/8/2017***************************************************************************
  if(index.res %in% "COAST2"){
    #Set Working directory (The folder files are imported and exported to)
    #setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
    coastal_ref <- read.csv("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016/new_coastal_ref_11_3_16.csv")
    merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF" & merged$EVENT_ID %in% coastal_ref$EVENT_ID, "REF",
                              ifelse(merged$CATEGORY %in% "REF" & !merged$EVENT_ID %in% coastal_ref$EVENT_ID, "MIN", as.character(merged$CATEGORY)))
    
  }
  
  #
  
  #merged$BIOREGION <- ifelse(merged$BIOREGION %in% c("MAC", "SEP"), "COAST", as.character(merged$BIOREGION))
  #==============================================================================
  # Test Strict Classification
  #merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF", "REF", 
  #                          ifelse(merged$CATEGORY %in% c("DEG+", "DEG"), "MOD",
  #                                 ifelse(merged$CATEGORY %in% c("DEG2"), "SEV",
  #                                        ifelse(merged$CATEGORY %in% "REF+", "REF", merged$CATEGORY))))
  
  
  #==============================================================================
  # Change NGV to UNP.
  # NGV was an originally treated as an independent bioregion but was later
  # subsumed into the UNP bioregion.
  merged$BIOREGION <- ifelse(merged$BIOREGION %in% c("NGV"), "UNP",
                             as.character(merged$BIOREGION))
  #==============================================================================
  # Export the metrics table as a csv.
  write.csv(merged, paste(bioregion, "_", taxon.rank, "_raw_metrics_",
                          todays_date, ".csv", sep = ""), row.names = F)
  #==============================================================================
  return(merged)
}


#==============================================================================
#'Automate Scoring
#'
#'@param raw.data =
#'@param bioregion = 
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param taxon.rank =
#'@param todays_date =
#'@return 
#'@export
#'

auto_score <- function(raw.data, bioregion, taxon_rank, prod_date, metric.freq,
                       master.df = master, metric.class = m.c, 
                       todays_date = format(Sys.time(), "%m_%d_%y"),
                       num_itr= "10i",
                       previous.sub.dir){
  #============================================================================
  
  redundant <- T
  range.var <- T
  zero.inflate <- F
  fam.group <- F
  metric.groups <- F
  ibi.method <- "H"
  #jack.iterate = 2
  #metric.class <- m.c
  #============================================================================
  # Identify the most frequently occuring metrics from the output of the 
  # iterative metric selection process.
  #============================================================================
  # I do NOT think this script is necessary. Update below.
  #if(!bioregion %in% "BASIN"){
  #  fm <- read.csv(paste0(previous.sub.dir, "/",
  #                        paste(bioregion, taxon_rank, num_itr, "Metric_Selection",
  #                       paste(prod_date, ".csv", sep = ""), sep = "_")))
  #}else{
  #  fm <- read.csv(paste0(previous.sub.dir, "/",
  #                        paste("BASIN", taxon_rank, num_itr, "Metric_Selection",
  #                       paste(prod_date, ".csv", sep = ""), sep = "_")))
  #}
  #============================================================================
  # UPDATE.
  fm <- read.csv(paste0(previous.sub.dir, "/",
                        paste(bioregion, taxon_rank, num_itr, "Metric_Selection",
                              paste(prod_date, ".csv", sep = ""), sep = "_")))
  
  #============================================================================
  # Select frequently occuring metrics.
  #============================================================================
  # We are aiming to create indices with a minimum of five metrics.  However,
  # there are instansous when it is not possile to select a minimum of five 
  # metrics becuase the too many of the metrics were not sensitive enough or
  # were redundant with one another.
  # The first if statement below checks for at least five metrics that met our
  # frequency expectations (metric.freq). If there are fewer than five metrics
  # that meet our frequency expectations but there are more than 5 metric
  # choices, then sort the data frame by FREQUENCY in descending order
  # and select the top five metrics.  If there are fewer than five metrics
  # reported, then select all available metrics.
  
  if (nrow(fm[fm$PERCENT >= metric.freq, ]) >= 5){
    final.metrics <- fm[fm$PERCENT >= metric.freq, "METRICS"]
  } else {
    if (nrow(fm) >= 5) {
      fm <- fm[order(-fm$PERCENT), ] 
      metric.freq2 <- fm[5, "PERCENT"]
      final.metrics <- fm[fm$PERCENT >= metric.freq2, "METRICS"]
    } else {
      fm <- fm[order(-fm$PERCENT), ] 
      final.metrics <- fm$METRICS
    }
    
  }
  final.metrics <- as.character(final.metrics)
  #============================================================================
  # Columns to keep: Sample information and the metrics selected to represent
  # the index.
  keep.cols <- c("EVENT_ID", "BIOREGION", "CATEGORY","STATION_ID", "SAMPLE_NUMBER",
                 "AGENCY_CODE", "DATE", as.character(final.metrics))
  #============================================================================
  # Select only data from the bioregion of interest
  new.df <- raw.data[raw.data$BIOREGION %in% bioregion, ]
  new.df <- new.df[new.df$ABUNDANCE > 70, ]
  new.df <- new.df[, names(new.df) %in% keep.cols]
  # Samples must contain more than 70 observed individuals.
  # Smaller sample sizes may result in odd percentage metric values.
  #============================================================================
  if(!bioregion %in% "BASIN"){
    summary.metrics <- read.csv(paste0(previous.sub.dir, "/", 
                                       paste(bioregion, taxon_rank, num_itr, "Metric_Summary",
                                      paste(prod_date, ".csv", sep = ""), sep = "_")))
  }else{
    summary.metrics <- read.csv(paste0(previous.sub.dir, "/", 
                                       paste("BASIN", taxon_rank, num_itr, "Metric_Summary",
                                      paste(prod_date, ".csv", sep = ""), sep = "_")))
  }
  #============================================================================
  redund.test <- redundancy(new.df[, !names(new.df) %in% "BIOREGION"], "wilcox",
                            summary.metrics, upper.coef = 0.8, lower.coef = -0.8,
                            upper.class = "REF", lower.class = "SEV")
  #new.test <- new.df[, !names(new.df) %in% "BIOREGION"]
  #new.test2 <- new.test[complete.cases(new.test[, 7:ncol(new.test)]), ]
  #corr.df <- cor(new.test2[, 7:ncol(new.test2)], method = "spearman")
  #corr.df[upper.tri(corr.df, diag = "TRUE")] <- NA
  #corr.df <- data.frame(corr.df)
  #corr.df$Metrics <- rownames(corr.df)
  #final.df <- tidyr::gather(corr.df, METRICS, COEF, -Metrics)
  #final.df <- final.df[!is.na(final.df$COEF), ]
  keep.metrics <- unique(as.character(redund.test$METRICS))
  metrics.removed <- final.metrics[!final.metrics %in% keep.metrics]
  if (length(metrics.removed) > 0) {
    warning(paste("The following metric(s) were removed due to redundancy:",
                  paste(metrics.removed, collapse = ", ")))
    new.df <- new.df[, !names(new.df) %in% metrics.removed]
  }
  #============================================================================
  disturb.react <- unique(summary.metrics[, c("METRICS", "DISTURBANCE")])
  if (any(duplicated(disturb.react$METRICS))) {
    dups <- disturb.react[duplicated(disturb.react$METRICS), ]
    fine.df <- disturb.react[!disturb.react$METRICS %in% dups$METRICS, ]
    freqy <- data.frame(table(summary.metrics[summary.metrics$METRICS %in% dups$METRICS, "DISTURBANCE"],
                              summary.metrics[summary.metrics$METRICS %in% dups$METRICS, "METRICS"]))
    freqy <- freqy[freqy$Freq > 0, ]
    freqy$METRICS <- factor(freqy$Var2)
    alt.df <- do.call(rbind,lapply(split(freqy, freqy$METRICS), function(chunk) {
      chunk[which.max(chunk$Freq),]
    } ))
    alt.df <- alt.df[, c(2, 1)]
    names(alt.df) <- c("METRICS", "DISTURBANCE")
    disturb.react <- rbind(fine.df, alt.df)
  }
  
  
  if(!bioregion %in% "BASIN"){
    ms <- read.csv(paste0(previous.sub.dir, "/", 
                          paste(bioregion, taxon_rank, num_itr, "Mean_Thresholds",
                         paste(prod_date, ".csv", sep = ""), sep = "_")))
  }else{
    ms <- read.csv(paste0(previous.sub.dir, "/", 
                          paste("BASIN", taxon_rank, num_itr, "Mean_Thresholds",
                         paste(prod_date, ".csv", sep = ""), sep = "_")))
  }
  
  names(ms)[names(ms) %in% c("MEAN_CEILING", "MEAN_FLOOR")] <- c("REF_MEDIAN", "BOUND_BI_CMA")
  ms <- merge(disturb.react, ms, by = "METRICS")
  ms <- ms[ms$METRICS %in% names(new.df), ]
  #============================================================================
  sub.new <- new.df[new.df$CATEGORY %in% c("REF", "SEV"), ]
  ms$METRICS <- as.character(ms$METRICS)
  ms$SENSITIVITY <- lapply(unique(ms$METRICS), function(metric) {
    print(metric)
    sub.df <- sub.new[, c("CATEGORY", as.character(metric))]
    ref.df <- sub.df[sub.df$CATEGORY %in% "REF", ]
    deg.df <- sub.df[sub.df$CATEGORY %in% "SEV", ]
    sub.ms <- ms[ms$METRICS %in% metric, ]
    
    if (sub.ms$DISTURBANCE %in% "DECREASE") {
      sub.ms$BSP <- sub.ms$REF_MEDIAN - (abs(sub.ms$REF_MEDIAN - sub.ms$BOUND_BI_CMA) / 2)
      deg.pct <- (nrow(deg.df[deg.df[, metric] < sub.ms$BSP, ]) / nrow(deg.df)) * 100
      ref.pct <- (nrow(ref.df[ref.df[, metric] >= sub.ms$BSP, ]) / nrow(ref.df)) * 100
    }
    
    if (sub.ms$DISTURBANCE %in% "INCREASE") {
      sub.ms$BSP <- sub.ms$REF_MEDIAN  + (abs(sub.ms$REF_MEDIAN - sub.ms$BOUND_BI_CMA) / 2)
      ref.pct <- (nrow(ref.df[ref.df[, metric] <= sub.ms$BSP, ]) / nrow(ref.df)) * 100
      deg.pct <- (nrow(deg.df[deg.df[, metric] > sub.ms$BSP, ]) / nrow(deg.df)) * 100
    }
    final.value <- mean(ref.pct, deg.pct)
  })
  ms$SENSITIVITY <- as.numeric(as.character(ms$SENSITIVITY))
  #============================================================================
  
  #ms2 <- metrics_summary(new.df, bioregion, zero = zero.inflate)
  #ms <- merge(ms, ms2[, c("METRICS", "SENSITIVITY")], by = "METRICS")
  
  write.csv(ms, paste(bioregion, "_", taxon_rank, "_ms_", todays_date, ".csv", sep = ""), row.names = F)
  #============================================================================
  
  scored <- old_scoring3(new.df, bioregion, zero_null = zero.inflate, bound.lim = TRUE, metric.summary = ms)
  if (ncol(scored) > 7) {
    scored[, 7:ncol(scored)] <- apply(scored[, 7:ncol(scored)], 2, function(x) as.numeric(as.character(x)))
    scored$FINAL_SCORE <- apply(scored[, 7:ncol(scored)], 1, mean, na.rm = TRUE)
  } else {
    scored[, 7] <- as.numeric(as.character(scored[, 7]))
    scored$FINAL_SCORE <- scored[, 7]
  }
  
  #============================================================================
  auto_plot(scored, new.df, taxon_rank, bioregion, todays_date)
  #============================================================================
  scored$BIOREGION <- bioregion
  scored <- scored[, c(ncol(scored), 1:(ncol(scored) - 1))]
  write.csv(scored, paste(bioregion, "_", taxon_rank, "_scored_metrics_", todays_date, ".csv", sep = ""), row.names = F)
  
  return(scored)
  
}

#==============================================================================
#'Automate Plots
#'
#'@param raw.data =
#'@param bioregion = 
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param taxon.rank =
#'@param todays_date =
#'@return 
#'@export
#'

auto_plot <- function(plot.score, plot.value, taxon_rank, bioregion,
                      todays_date = format(Sys.time(), "%m_%d_%y")){
  
  #============================================================================
  # Create a pdf of the final index distributions
  pdf(paste(bioregion, "_", taxon_rank, "_Index_", todays_date, ".pdf", sep = ""))
  #png(paste(bioregion, "_", taxon_rank,   "_Index_", todays_date, ".png", sep = ""),
  #   units = "in", res = 720,  width = 6.5, height = 7)
  plot.score <- plot.score[!plot.score$CATEGORY %in% "MIX", ]
  plot.score$CATEGORY <- factor(plot.score$CATEGORY,
                                levels = c("REF", "MIN", "MOD", "SEV"))
  boxplot(plot.score[, "FINAL_SCORE"] * 100 ~ plot.score$CATEGORY,
          a <- paste(taxon_rank, bioregion),
          data = plot.score,
          col="white", #las = 2,
          main = substitute(paste(a)))
  
  dev.off()
  
  #============================================================================
  # Create a pdf of all of the individual metric used in the final index
  pdf(paste(bioregion, "_", taxon_rank, "_raw_metrics_", todays_date, ".pdf", sep = ""))
  par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(2,2,2,2))
  
  plot.me <- plot.value[,names(plot.value) %in% names(plot.score)]
  plot.me <- plot.me[!plot.me$CATEGORY %in% "MIX", ]
  col.names <- names(plot.me)
  col.names <- col.names[!col.names %in% c("div", "FFG", "HABIT",
                                           "tol", "comp", "fam", "FINAL_SCORE")]
  for (j in col.names[!col.names %in% c("EVENT_ID", "CATEGORY", "STATION_ID", 
                                        "SAMPLE_NUMBER", "AGENCY_CODE", "DATE")]){
    a <- j
    #c_list <- list('REF', 'NEAR', 'MIN', 'MOD', 'SEV', 'MIX')
    plot.me$CATEGORY <- factor(plot.me$CATEGORY,
                               levels = c("REF", "MIN", "MOD", "SEV"))
    
    boxplot(plot.me[, j] ~ plot.me$CATEGORY,
            #data = raw.data,
            #names = c_list,
            col="white", las = 0,
            main = substitute(paste(a)))
  }
  dev.off()
  #============================================================================
  # Create a pdf of all of the individual metric used in the final index
  pdf(paste(bioregion, "_", taxon_rank, "_scored_metrics_", todays_date, ".pdf", sep = ""))
  par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(2,2,2,2))
  plot.me <- plot.score
  #***
  plot.me <- plot.me[!plot.me$CATEGORY %in% "MIX", ]
  col.names <- names(plot.me)
  col.names <- col.names[!col.names %in% c("div", "FFG", "HABIT",
                                           "tol", "comp", "fam", "FINAL_SCORE")]
  for (j in col.names[!col.names %in% c("EVENT_ID", "CATEGORY", "STATION_ID", 
                                        "SAMPLE_NUMBER", "AGENCY_CODE", "DATE")]){
    a <- j
    #c_list <- list('REF', 'NEAR', 'MIN', 'MOD', 'SEV', 'MIX')
    plot.me$CATEGORY <- factor(plot.me$CATEGORY,
                               levels = c("REF", "MIN", "MOD", "SEV"))
    
    boxplot(plot.me[, j] ~ plot.me$CATEGORY,
            #data = raw.data,
            #names = c_list,
            col="white", las = 0,
            main = substitute(paste(a)))
  }
  dev.off()
  #============================================================================
  
  
  
}



#==============================================================================
#'Combine Multiple Automation Functions to Create an Index
#'
#'@param raw.data =
#'@param bioregion = 
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param taxon.rank =
#'@param todays_date =
#'@return 
#'@export
#'


auto_index <- function(my.data, bioregion, taxon_rank,
                       calc.date = "11_09_16", jack.runs = 100, 
                       todays_date = format(Sys.time(), "%m_%d_%y"), seed = TRUE,
                       num.itr = "10i",
                       main.dir = "D:/ZSmith/Projects/Chessie_BIBI/Output/March_2017",
                       previous.subdir = NULL){
  #==============================================================================
  # Prep for subdirectory
  if (bioregion %in% "BASIN") spat.res <- "BASIN"
  if (bioregion %in% c("COAST", "INLAND")) spat.res <- "REGION"
  if (bioregion %in% c("BLUE", "CA", "LNP", "MAC", "NAPU", "NCA",
                       "NRV", "PIED", "SEP", "SGV", "SRV", "UNP")) spat.res <- "BIOREGION"
  # Create/specify the appropriate subdirectory
  new.dir <- create_subdir(main.dir, spat.res, taxon_rank)
  setwd(new.dir)
  if (is.null(previous.subdir)) previous.subdir <- new.dir
  #==============================================================================
  # Rename todays_date to avoid errors when specifying todays date in the 
  # functions below. An error will occur if todays_date = todays_date in any of
  # the functions below.
  td <- todays_date
  #==============================================================================
  # If seed is TRUE then used the randomly selected seed of 62.
  # This should eliminate the variability caused by rarefaction.
  if (seed == TRUE) set.seed(62)
  #==============================================================================
  # Calculate all of the metrics.
  metrics.df <- auto_metrics(my.data, bioregion, bioregion, taxon_rank,
                             todays_date = td)
  metrics.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  #==============================================================================
  # Changed frequency to 20 on 4/17/2017
  scores.df <- auto_score(metrics.df, bioregion, taxon_rank, calc.date, metric.freq = 20, 
                          todays_date = td, num_itr = num.itr,
                          previous.sub.dir = previous.subdir)
  #==============================================================================
  if (seed == TRUE) set.seed(62)
  jack.df <- bibi_jackknife(metrics.df[, names(metrics.df) %in% names(scores.df)],
                            jack.runs, 0.75, bioregion, m.c, Fam = F, Master = master,
                            zero.null = F, metric_types = F, redund.df = F,
                            method = "H", seed = TRUE)
  #==============================================================================
  t.bde <- bde(scores.df, bioregion)
  write.csv(t.bde, paste(bioregion, "_", taxon_rank, "_bde_", todays_date, ".csv", sep = ""), row.names = F)
  #==============================================================================
  sub.bde <- t.bde[1, ]
  jack.df$ORIGINAL_CE <-  sub.bde$CE
  jack.df$ORGINAL_BSP <- sub.bde$THRESHOLD
  write.csv(jack.df, paste(bioregion, "_", taxon_rank, "_jackknife_", todays_date, ".csv", sep = ""), row.names = F)
  #==============================================================================
  if(seed == TRUE) set.seed(sample(1:10000, 1))
  #==============================================================================
}























