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
    my.data <- my.data[my.data$BIOREGION %in% c("MAC", "SEP"), ]
    my.data$BIOREGION <- "COAST"
  }
  if(index_res %in% "INLAND"){
    my.data <- my.data[!my.data$BIOREGION %in% c("MAC", "SEP"), ]
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
  
  env.df <- prep_env(Prep.data = all.data)
  test.class <- prep_class_multi(env.df)
  all.data <- auto_bio(all.data, index.res, bioregion)
  #names(all.data)[names(all.data) %in% "ICPRB_BIOREGION_ID"] <- "BIOREGION"
  bio.data <- all.data[!is.na(all.data$BIOREGION), ]
  table(test.class$BIOREGION, test.class$CATEGORY)
  env.class <- unique(test.class[, c("EVENT_ID", "STATION_ID", "DATE",
                                     "AGENCY_CODE", "SAMPLE_NUMBER", "CATEGORY")])
  #==============================================================================
  # Calculate Standard Metrics
  
  bio.data.metrics <- all_metrics(master, bio.data, taxa.rank =  taxon.rank,
                                  rare = "DIV_RARE", pct_un = 10, bibi.standard = TRUE,
                                  seed = TRUE)
  master2 <- clean_taxa(master)
  #==============================================================================
  # Calculate the percent of each taxon in the data base by Sequencing through each taxon rank
  if(taxon.rank %in% "ORDER") ranks.list <- c("SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  if(taxon.rank %in% "FAMILY") ranks.list <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  if(taxon.rank %in% "GENUS") ranks.list <- c("SPECIES")
  master2[, ranks.list] <- NA
  bio.data.seq <- seq_pct_taxa(bio.data, master2)
  test <- bio.data.seq[duplicated(bio.data.metrics$EVENT_ID), ]
  #==============================================================================
  # Merge the Standard Metrics with the Sequenced Metrics
  bio.data.merge <- merge(bio.data.metrics, bio.data.seq, by = c("EVENT_ID", "STATION_ID",
                                                                 "DATE", "SAMPLE_NUMBER", "AGENCY_CODE"))
  bio.data.merge <- bio.data.merge[, !grepl("PCT_UNIDENTIFIED", names(bio.data.merge))]
  bio.data.merge <- bio.data.merge[, !grepl(".y", names(bio.data.merge))]
  
  new <- data.frame(names(bio.data.merge[, 7:ncol(bio.data.merge)]))
  new$MIN <- round(apply(bio.data.merge[, 7:ncol(bio.data.merge)], 2, min))
  new$MAX <- round(apply(bio.data.merge[, 7:ncol(bio.data.merge)], 2, max), 0)
  #==============================================================================
  # Merge the Metrics with the Site Classes
  m.bio.data <- merge(env.class, bio.data.merge,
                      by = c("EVENT_ID", "STATION_ID", "DATE",
                             "AGENCY_CODE", "SAMPLE_NUMBER"), all.y = TRUE)
  if("PCT_UNIDENTIFIED"  %in% names(m.bio.data)){
    m.bio.data <- m.bio.data[, !grepl("PCT_UNIDENTIFIED", names(m.bio.data))]
  }
  
  if(any(grepl(".y", names(m.bio.data)))){
    m.bio.data <- m.bio.data[, !grepl(".y", names(m.bio.data))]
  }
  
  if(any(grepl(".x", names(m.bio.data)))){
    names(m.bio.data) <- gsub(".x", "", names(m.bio.data))
  }
  #m.bio.data$CATEGORY <- factor(m.bio.data$CATEGORY, levels =c("REF", "MIN", "MOD", "SEV", "MIX"))
  table(m.bio.data$CATEGORY)
  test <- data.frame(names(m.bio.data[, 7:ncol(m.bio.data)]))
  #==============================================================================
  # Re-assign the the bioregion to the newest data frame
  merged <- merge(unique(bio.data[, c("EVENT_ID", "BIOREGION")]), m.bio.data, by = "EVENT_ID", all.y = TRUE)
  #write.csv(merged, "Test2_Metrics_11_7_16.csv", row.names = FALSE)
  
  #==============================================================================
  # Change REF+ to REF and DEG+ to SEV
  merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF+", "REF", 
                            ifelse(merged$CATEGORY %in% c("DEG+", "DEG"), "MOD", 
                                   ifelse(merged$CATEGORY %in% c("DEG2"), "SEV", merged$CATEGORY)))
  if(index.res %in% "COAST"){
    #Set Working directory (The folder files are imported and exported to)
    setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
    coastal_ref <- read.csv("new_coastal_ref_11_3_16.csv")
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
  # Change NGV to UNP
  merged$BIOREGION <- ifelse(merged$BIOREGION %in% c("NGV"), "UNP", as.character(merged$BIOREGION))
  #==============================================================================
  write.csv(merged, paste(bioregion, "_", taxon.rank, "_raw_metrics_", todays_date, ".csv", sep = ""), row.names = F)
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
                       todays_date = format(Sys.time(), "%m_%d_%y")){
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
  # Identify the most frequently occuring metrics...
  if(!bioregion %in% "BASIN"){
    fm <- read.csv(paste(bioregion, taxon_rank, "10i", "Metric_Selection",
                         paste(prod_date, ".csv", sep = ""), sep = "_"))
  }else{
    fm <- read.csv(paste("BASIN", taxon_rank, "10i", "Metric_Selection",
                         paste(prod_date, ".csv", sep = ""), sep = "_"))
  }
  
  if (nrow(fm[fm$FREQUENCY >= metric.freq, ]) >= 5){
    final.metrics <- fm[fm$FREQUENCY >= metric.freq, "METRICS"]
  } else {
    if (nrow(fm) >= 5) {
      fm <- fm[order(-fm$FREQUENCY), ] 
      metric.freq2 <- fm[5, 2]
      final.metrics <- fm[fm$FREQUENCY >= metric.freq2, "METRICS"]
    } else {
      fm <- fm[order(-fm$FREQUENCY), ] 
      final.metrics <- fm$METRICS
    }
    
  }
  
  #============================================================================
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
    summary.metrics <- read.csv(paste(bioregion, taxon_rank, "10i", "Metric_Summary",
                                      paste(prod_date, ".csv", sep = ""), sep = "_"))
  }else{
    summary.metrics <- read.csv(paste("BASIN", taxon_rank, "10i", "Metric_Summary",
                                      paste(prod_date, ".csv", sep = ""), sep = "_"))
  }
  
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
    ms <- read.csv(paste(bioregion, taxon_rank, "10i", "Mean_Thresholds",
                         paste(prod_date, ".csv", sep = ""), sep = "_"))
  }else{
    ms <- read.csv(paste("BASIN", taxon_rank, "10i", "Mean_Thresholds",
                         paste(prod_date, ".csv", sep = ""), sep = "_"))
  }
  
  names(ms)[names(ms) %in% c("MEAN_CEILING", "MEAN_FLOOR")] <- c("REF_MEDIAN", "BOUND_BI_CMA")
  ms <- merge(disturb.react, ms, by = "METRICS")
  ms <- ms[ms$METRICS %in% names(new.df), ]
  
  ms2 <- metrics_summary(new.df, bioregion, zero = zero.inflate)
  ms <- merge(ms, ms2[, c("METRICS", "SENSITIVITY")], by = "METRICS")
  
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
                       todays_date = format(Sys.time(), "%m_%d_%y"), seed = TRUE){
  
  #==============================================================================
  td <- todays_date
  #==============================================================================
  if (seed == TRUE) set.seed(62)
  #==============================================================================
  metrics.df <- auto_metrics(my.data, bioregion, bioregion, taxon_rank,
                             todays_date = td)
  metrics.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  #==============================================================================
  scores.df <- auto_score(metrics.df, bioregion, taxon_rank, calc.date, metric.freq = 8, 
                          todays_date = td)
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























