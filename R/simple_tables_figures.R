#==============================================================================
#==============================================================================
# Author: Zachary M. Smith
# Maintainer: Zachary M. Smith
# Organization: ICPRB
# Created: December 2016
# Updated: 3-15-2017
# Purpose: Simplify the development of tables and figures for BIBI report.
#==============================================================================
#==============================================================================
#'Group Me 2
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

group.me2 <- function(spatial, taxon_rank, data_set, calc.date, final_score = FALSE){
  keep.cols <- c("BIOREGION", "EVENT_ID", "STATION_ID", "SAMPLE_NUMBER", "DATE",
                 "AGENCY_CODE", "CATEGORY", "FINAL_SCORE")
  
  if(spatial %in% "BIOREGION"){
    bioregions <- c("BLUE", "CA", "LNP", "MAC", "NAPU", "NCA",
                    "NRV", "PIED", "SEP", "SGV", "SRV", "UNP")
    data.list <- list()
    for(i in bioregions){
      csv.file <- paste(paste(i, taxon_rank, data_set, calc.date, sep = "_"), ".csv", sep = "")
      if(file_test("-f", csv.file)){
        bio.df <- read.csv(csv.file)
        bio.df$BIOREGION <- i
        if (final_score == FALSE) {
          data.list[[i]] <- bio.df
        } else {
          data.list[[i]] <- bio.df[, keep.cols]
        }
        
      }
    }
    bound <- do.call(rbind, data.list)
  }
  
  if(spatial %in% "REGION"){
    regions <- c("COAST", "INLAND")
    data.list <- list()
    for(i in regions){
      csv.file <- paste(paste(i, taxon_rank,  data_set, calc.date, sep = "_"), ".csv", sep = "")
      if(file_test("-f", csv.file)){
        bio.df <- read.csv(csv.file)
        bio.df$BIOREGION <- i
        if (final_score == FALSE) {
          data.list[[i]] <- bio.df
        } else {
          data.list[[i]] <- bio.df[, keep.cols]
        }
      }
    }
    bound <- do.call(rbind, data.list)
  }
  
  if(spatial %in% "BASIN"){
    if (file.exists(paste(paste("ALL", taxon_rank,  data_set, calc.date, sep = "_"), ".csv", sep = ""))){
      bio.df <- read.csv(paste(paste("ALL", taxon_rank,  data_set, calc.date, sep = "_"), ".csv", sep = ""))
    } else {
      bio.df <- read.csv(paste(paste("BASIN", taxon_rank,  data_set, calc.date, sep = "_"), ".csv", sep = ""))
    }
    
    bio.df$BIOREGION <- "BASIN"
    if (final_score == FALSE) {
      bound <- bio.df
    } else {
      bound <- bio.df[, keep.cols]
    }
  }
  
  return(bound)
}

#==============================================================================
#'Metric Table
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

metric_table <- function(spatial, taxon_rank, data_set, calc.date,  m.c, 
                         todays_date = format(Sys.time(), "%m_%d_%y")){
  #grouped <- group.me2("H_","FAM_10_13_16", "_ms_")
  grouped <- group.me2(spatial, taxon_rank, data_set, calc.date)
  
  
  metric.sum <- grouped[, c("BIOREGION", "METRICS", "DISTURBANCE",
                            "REF_MEDIAN", "BOUND_BI_CMA", "SENSITIVITY")]
  
  metric.sum <- merge(m.c, metric.sum, by = "METRICS", all.y = T)
  metric.sum$METRIC_CLASS <- ifelse(is.na(metric.sum$METRIC_CLASS), "Composition",
                                    as.character(metric.sum$METRIC_CLASS))
  
  #test <- metric.sum[!metric.sum$METRICS %in% m.c$METRICS,]
  
  metric.sum <- metric.sum[order(metric.sum$BIOREGION, metric.sum$METRIC_CLASS, metric.sum$METRICS), ]
  metric.sum2 <- metric.sum[, c( "BIOREGION", "METRICS", "METRIC_CLASS", "DISTURBANCE",
                                 "REF_MEDIAN", "BOUND_BI_CMA", "SENSITIVITY")]
  metric.sum2$METRIC_CLASS <- ifelse(metric.sum2$METRIC_CLASS %in% "DIVERSITY", "Richness/Diversity",
                                     ifelse(metric.sum2$METRIC_CLASS %in% "COMPOSITION", "Composition",
                                            ifelse(metric.sum2$METRIC_CLASS %in% "HABIT", "Habit",
                                                   ifelse(metric.sum2$METRIC_CLASS %in% "TOLERANCE", "Tolerance",
                                                          as.character(metric.sum2$METRIC_CLASS)))))
  metric.sum2$DISTURBANCE <- as.character(metric.sum2$DISTURBANCE)
  metric.sum2[, c("DISTURBANCE")] <- paste(toupper(substr(metric.sum2[, c("DISTURBANCE")], 1, 1)),
                                           tolower(substr(metric.sum2[, c("DISTURBANCE")], 2,
                                                          nchar(metric.sum2[, c("DISTURBANCE")]))), sep="")
  names(metric.sum2) <- c("Bioregion", "Metric", "Metric Class", "Influence of Disturbance", "Reference Median",
                          "Bound", "Metric BDE")
  metric.sum2[, 5:ncol(metric.sum2)] <- round(metric.sum2[, 5:ncol(metric.sum2)], 2)
  metric.sum2$Bound <- ifelse(metric.sum2$Bound > 100, 100,
                              ifelse(metric.sum2$Bound < 0, 0, as.numeric(metric.sum2$Bound)))
  #============================================================================
  #metric.sum2$NUM <- ave(metric.sum2$'Metric BDE', metric.sum2$Bioregion, FUN = seq_along)
  metric.sum2 <- metric.sum2[, c("Bioregion", "Metric", "Metric Class",
                                 "Influence of Disturbance", "Reference Median",
                                 "Bound", "Metric BDE")]
  #============================================================================
  if(spatial %in% "BASIN"){
    write.csv(metric.sum2, paste("BASIN", "_", taxon_rank, "_metrics_table_",
                                 todays_date, ".csv", sep = ""), row.names = F)
  }
  if(spatial %in% "REGION"){
    write.csv(metric.sum2, paste("REGION", "_", taxon_rank, "_metrics_table_",
                                 todays_date, ".csv", sep = ""), row.names = F)
  }
  if(spatial %in% "BIOREGION"){
    write.csv(metric.sum2, paste("BIOREGION", "_", taxon_rank, "_metrics_table_",
                                 todays_date, ".csv", sep = ""), row.names = F)
  }
  
}

#==============================================================================
#'Index Table
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

index_table <- function(spatial, taxon_rank, data_set, calc.date,
                        todays_date = format(Sys.time(), "%m_%d_%y")){
  grouped <- group.me2(spatial, taxon_rank, data_set, calc.date)
  #============================================================================
  original.df <- unique(grouped[, c("BIOREGION", "ORGINAL_BSP", "ORIGINAL_CE")])
  names(original.df) <- c("Bioregion", "Index BSP", "Index CE")
  #original.df$NUM <- seq_along(original.df[, 1])
  #original.df <- original.df[, c("NUM", "Bioregion", "Index BSP", "Index CE")]
  original.df <- original.df[, c("Bioregion", "Index BSP", "Index CE")]
  write.csv(original.df, paste(spatial, "_", taxon_rank, "_Index_Parameters_",
                               todays_date, ".csv", sep = ""), row.names = F)
}

#==============================================================================
#'Jackknife Table
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

jack_table <- function(spatial, data_set, calc.date,
                       todays_date = format(Sys.time(), "%m_%d_%y"),
                       maindir = main.dir,
                       prevdir = previous.subdir){
  data.list <- list()
  for(i in c("ORDER", "FAMILY", "GENUS")){
    #==============================================================================
    if (is.null(prevdir)) {
    month.year <- format(Sys.Date(), format="%B_%Y")
    todays.date <- format(Sys.Date(), format="%m_%d_%Y")
    check.dir <- paste0(maindir, "/", month.year)
    check.dir <- paste0(check.dir, "/", todays.date)
    check.dir <- paste0(check.dir, "/", spatial)
    check.dir <- paste0(check.dir, "/", i)
    if (!dir.exists(check.dir)) next
    # Create/specify the appropriate subdirectory
    new.dir <- create_subdir(maindir, spatial, i)
    setwd(new.dir)
    } else {
      check.dir <- paste0(prevdir, "/", i)
      if (!dir.exists(check.dir)) {
        print("skip")
        next
      } 
      setwd(check.dir)
    }
    #==============================================================================
    grouped <- group.me2(spatial, i, data_set, calc.date)
    #==========================================================================
    jv <- grouped[, c("BIOREGION", "MEASURE", "ORIGINAL_CE", "ORGINAL_BSP", "MEAN_SIM_VALUE", "RMSE")]
    jv$ORIGINAL_VALUE <- ifelse(jv$MEASURE %in% c("VAL_CE", "TRAIN_CE"), jv$ORIGINAL_CE,
                                ifelse(jv$MEASURE %in% "BSP", jv$ORGINAL_BSP, "ERROR"))
    jv <- jv[, c("BIOREGION", "MEASURE", "ORIGINAL_VALUE", "MEAN_SIM_VALUE", "RMSE")]
    #jv <- jv[jv$MEASURE %in% c("VAL_CE", "TRAIN_CE"), ]
    names(jv) <- c("Bioregion", "Measure", "Original Value", "Mean Simulated Value", "RMSE")
    jv$'Taxonomic Tier' <- paste(substring(i, 1, 1), tolower(substring(i,2)), sep = "")
    data.list[[i]] <- jv
  }
  ord <- data.list$ORDER
  fam <- data.list$FAMILY
  gen <- data.list$GENUS
  jack.df <- rbind(ord, fam, gen)
  
  jack.df$Measure <- ifelse(jack.df$Measure %in% "TRAIN_CE", "Training CE",
                            ifelse(jack.df$Measure %in% "VAL_CE", "Validation CE",
                                   ifelse(jack.df$Measure %in% "BSP", "BSP", "ERROR")))
  jack.df$Measure <- factor(jack.df$Measure, levels = c("Training CE", "Validation CE", "BSP"))
  jack.df$`Taxonomic Tier` <- factor(jack.df$`Taxonomic Tier`, levels = c("Order", "Family", "Genus"))
  jack.df <- jack.df[order(jack.df$Bioregion, jack.df$'Taxonomic Tier', jack.df$Measure), ]
  jack.df <- jack.df[, c(1, 6, 2:5)]
  write.csv(jack.df, paste(spatial, "_Jackknife_Table_",
                           todays_date, ".csv", sep = ""), row.names = F)
  #============================================================================
  
  #data.list <- list()
  #for(i in c("ORDER", "FAMILY", "GENUS")){
  # grouped <- group.me2(spatial, i, data_set, calc.date)
  #==========================================================================
  #  jv <- grouped[, c("BIOREGION", "MEASURE", "ORGINAL_BSP", "MEAN_SIM_VALUE", "RMSE")]
  #  jv <- jv[jv$MEASURE %in% c("BSP"), ]
  #  names(jv) <- c("Bioregion", "Measure", paste(i, "Original BSP"),
  #                 paste(i, "Mean Simulated BSP"), 
  #                 paste(i, "RMSE BSP"))
  #  data.list[[i]] <- jv
  #}
  #ord <- data.list$ORDER
  #fam <- data.list$FAMILY
  #gen <- data.list$GENUS
  #jack.bsp <- plyr::join_all(list(ord, fam, gen), by = c("Bioregion", "Measure"), match = "all")
  #write.csv(jack.bsp, paste(spatial, "_Jackknife_BSP_",
  #                         todays_date, ".csv", sep = ""), row.names = F)
  
}
#==============================================================================
#'Plot Me2
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

plot.me2 <- function(scores.df, spatial, taxon_rank, todays_date){
  library(ggplot2)
  library(stringr)
  #scores.df <- scores.df[scores.df$BIOREGION %in% bioregion, ]
  scores.df <- scores.df[!scores.df$CATEGORY %in% "MIX", ]
  
  scores.df$CATEGORY <- factor(scores.df$CATEGORY,
                               #levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
                               levels = c("REF", "MIN", "MOD", "SEV"))
  #----------------------------------------------------------------------------
  agg.df <- aggregate(scores.df[, "SCORE"] ~ scores.df[, "CATEGORY"],
                      data = scores.df, FUN = max)
  names(agg.df) <- c("CATEGORY", "MAX")
  len.df <- aggregate(scores.df[, "SCORE"] ~ scores.df[, "CATEGORY"],
                      data = scores.df, FUN = length)
  names(len.df) <- c("CATEGORY", "N")
  n.size <- merge(agg.df, len.df, by = "CATEGORY")
  #----------------------------------------------------------------------------
  
  p <- ggplot(data = scores.df, aes(CATEGORY, SCORE)) +
    stat_boxplot(geom ='errorbar', width = 0.5) + 
    geom_boxplot(notch = FALSE) + #, aes(fill = factor(CATEGORY))) +
    geom_hline(aes(yintercept = ORIGINAL_BSP), size = 1, color = "red2", linetype = "dashed") +
    #scale_fill_manual( values = c("forestgreen", "gold", "goldenrod2", "red", "red4")) +
    scale_x_discrete(labels = c("Reference", "Minimally\nDegraded",
                                #"Mixed", 
                                "Moderately\nDegraded", "Degraded")) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          # text = element_text(size = 20),
          #axis.text.x = element_text(size = 15),
          #axis.title.y = element_text(margin= margin(0, 10, 0, 0)),
          #axis.title.x = element_text(margin= margin(10, 0, 0, 0)),
          #plot.margin  = unit(c(0.5, 1, 0.5, 0.5), "cm"),
          text = element_text(size = 10)) +
    #ylim(0, 100) +
    ylim(0, max(n.size$MAX) + 4) +
    ylab("Final Score") +
    xlab("Class") +
    geom_text(data = n.size, aes(y = MAX + 1, label = N), vjust = 0)
  return(p)
}

#==============================================================================
#'Plot Me
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

plot.me <- function(scores.df, spatial, taxon_rank, todays_date){
  library(ggplot2)
  library(stringr)
  #scores.df <- scores.df[scores.df$BIOREGION %in% bioregion, ]
  scores.df <- scores.df[!scores.df$CATEGORY %in% "MIX", ]
  
  scores.df$CATEGORY <- factor(scores.df$CATEGORY,
                               #levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
                               levels = c("REF", "MIN", "MOD", "SEV"))
  
  
  plot.list <-lapply(unique(scores.df$BIOREGION), function(x){
    sub.scores <- scores.df[scores.df$BIOREGION %in% x, ]
    #----------------------------------------------------------------------------
    agg.df <- aggregate(sub.scores[, "SCORE"] ~ sub.scores[, "CATEGORY"],
                        data = sub.scores, FUN = max)
    names(agg.df) <- c("CATEGORY", "MAX")
    len.df <- aggregate(sub.scores[, "SCORE"] ~ sub.scores[, "CATEGORY"],
                        data = sub.scores, FUN = length)
    names(len.df) <- c("CATEGORY", "N")
    n.size <- merge(agg.df, len.df, by = "CATEGORY")
    #----------------------------------------------------------------------------
    plot.grob <- ggplot(data = sub.scores, aes(CATEGORY, SCORE)) +
      stat_boxplot(geom ='errorbar', width = 0.5) + 
      geom_boxplot(notch = FALSE) + #, aes(fill = factor(CATEGORY))) +
      geom_hline(aes(yintercept = ORIGINAL_BSP), size = 1, color = "red2", linetype = "dashed") +
      #scale_fill_manual( values = c("forestgreen", "gold", "goldenrod2", "red", "red4")) +
      scale_x_discrete(labels = c("Reference", "Minimally\nDegraded",
                                  #"Mixed", 
                                  "Moderately\nDegraded", "Degraded")) +
      theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
      #labs(title = x) +
      ggtitle(x) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.line.x = element_line(color="black", size = 0.5),
            axis.line.y = element_line(color="black", size = 0.5),
            # text = element_text(size = 20),
            #axis.text.x = element_text(size = 15),
            #axis.title.y = element_text(margin= margin(0, 10, 0, 0)),
            #axis.title.x = element_text(margin= margin(10, 0, 0, 0)),
            #plot.margin  = unit(c(0.5, 1, 0.5, 0.5), "cm"),
            text = element_text(size = 10)) +
      #ylim(0, 100) +
      #ylim(0, max(n.size$MAX) + 8) +
      scale_y_continuous(limits = c(-2, 110), breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
      ylab("Final Score") +
      xlab("Class") +
      #geom_text(data = n.size, aes(y = MAX + 1, label = N), vjust = 0)
      geom_text(data = n.size, aes(y = 101, label = N), vjust = 0)
    
  })
  #----------------------------------------------------------------------------
  if (spatial %in% "BASIN") {
    png(paste(paste(spatial, taxon_rank, todays_date, sep = "_"),
              ".png", sep = ""),
        width     = 3.75,
        height    = 3,
        units     = "in",
        res       = 900)
    gridExtra::grid.arrange(plot.list[[1]])
    dev.off()
  }
  #----------------------------------------------------------------------------
  if (spatial %in% "REGION") {
    png(paste(paste(spatial, taxon_rank, todays_date, sep = "_"),
              ".png", sep = ""),
        width     = 6.2,
        height    = 3,
        units     = "in",
        res       = 900)
    gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]],
                            ncol = 2)
    dev.off()
  }
  #----------------------------------------------------------------------------
  if (spatial %in% "BIOREGION") {
    name.num <- if ("CA" %in% scores.df$BIOREGION) 1 else 2
    png(paste(paste(spatial, taxon_rank, name.num, todays_date, sep = "_"),
              ".png", sep = ""),
        width     = 6.2,
        height    = 6.5,
        units     = "in",
        res       = 900)
    gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]],
                            plot.list[[4]], plot.list[[5]], plot.list[[6]],
                            ncol = 2)
    dev.off()
    
  }
  
  
}
#==============================================================================
#'Plot Me Old
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

plot.me.old <- function(scores.df, spatial, taxon_rank, todays_date){
  library(ggplot2)
  library(stringr)
  #scores.df <- scores.df[scores.df$BIOREGION %in% bioregion, ]
  scores.df <- scores.df[!scores.df$CATEGORY %in% "MIX", ]
  
  scores.df$CATEGORY <- factor(scores.df$CATEGORY,
                               #levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
                               levels = c("REF", "MIN", "MOD", "SEV"))
  #----------------------------------------------------------------------------
  agg.df <- aggregate(scores.df[, "SCORE"] ~ scores.df[, "CATEGORY"],
                      data = scores.df, FUN = max)
  names(agg.df) <- c("CATEGORY", "MAX")
  len.df <- aggregate(scores.df[, "SCORE"] ~ scores.df[, "CATEGORY"],
                      data = scores.df, FUN = length)
  names(len.df) <- c("CATEGORY", "N")
  n.size <- merge(agg.df, len.df, by = "CATEGORY")
  #----------------------------------------------------------------------------
  
  p <- ggplot(data = scores.df, aes(CATEGORY, SCORE)) +
    stat_boxplot(geom ='errorbar', width = 0.5) + 
    geom_boxplot(notch = FALSE) + #, aes(fill = factor(CATEGORY))) +
    geom_hline(aes(yintercept = ORIGINAL_BSP), size = 1, color = "red2", linetype = "dashed") +
    #scale_fill_manual( values = c("forestgreen", "gold", "goldenrod2", "red", "red4")) +
    scale_x_discrete(labels = c("Reference", "Minimally\nDegraded",
                                #"Mixed", 
                                "Moderately\nDegraded", "Degraded")) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          # text = element_text(size = 20),
          #axis.text.x = element_text(size = 15),
          #axis.title.y = element_text(margin= margin(0, 10, 0, 0)),
          #axis.title.x = element_text(margin= margin(10, 0, 0, 0)),
          #plot.margin  = unit(c(0.5, 1, 0.5, 0.5), "cm"),
          text = element_text(size = 10)) +
    #ylim(0, 100) +
    #ylim(0, max(n.size$MAX) + 4) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, max(n.size$MAX) + 4)) +
    ylab("Final Score") +
    xlab("Class") +
    geom_text(data = n.size, aes(y = MAX + 1, label = N), vjust = 0)
  

  #labs(title = scores.df$BIOREGION) +
  if (spatial %in% "BASIN") {
    p #+ facet_wrap(~ BIOREGION) + theme(strip.text.x = element_text(size = 10,
      #                                                              face = "bold"),
      #                                  strip.background = element_blank())
    ggsave(file = paste(paste(spatial, taxon_rank, todays_date, sep = "_"), ".png", sep = ""),
           units = "in", 
           width = 3.75, 
           height = 3)
  }
  if (spatial %in% "REGION") {
    p #+ facet_wrap(~ BIOREGION, scales = "free") + theme(panel.spacing.x = unit(2, "lines"),
      #                                                   strip.text.x = element_text(size = 10,
      #                                                                               face = "bold"),
      #                                                   strip.background = element_blank())
    ggsave(file = paste(paste(spatial, taxon_rank, todays_date, sep = "_"), ".png", sep = ""),
           units = "in", 
           width = 6.2, 
           height = 3)
  }
  if (spatial %in% "BIOREGION") {
    p #+ facet_wrap(~ BIOREGION, nrow = 3, scales = "free") + theme(panel.spacing.x = unit(2, "lines"),
      #                                                             panel.spacing.y = unit(1, "lines"),
      #                                                             strip.text.x = element_text(size = 10,
      #                                                                                         face = "bold"),
      #                                                             strip.background = element_blank())
    name.num <- if ("CA" %in% scores.df$BIOREGION) 1 else 2
    ggsave(file=paste(paste(spatial, taxon_rank, name.num, todays_date, sep = "_"), ".png", sep = ""),
           units = "in", 
           width = 6.2, 
           height = 7)
  }
  
  
}

#==============================================================================
#'Plot Jackknife Results
#'
#'@param my.data = data.frame.
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

plot_jack <- function(spatial, data_set, calc.date, taxon_rank, todays_date,
                      rmse = TRUE, maindir = main.dir, prevdir = previous.subdir){
  library(ggplot2)
  library(stringr)
  data.list <- list()
  for(i in c("ORDER", "FAMILY", "GENUS")){
    #==============================================================================
    #==============================================================================
    if (is.null(prevdir)) {
      month.year <- format(Sys.Date(), format="%B_%Y")
      todays.date <- format(Sys.Date(), format="%m_%d_%Y")
      check.dir <- paste0(maindir, "/", month.year)
      check.dir <- paste0(check.dir, "/", todays.date)
      check.dir <- paste0(check.dir, "/", spatial)
      check.dir <- paste0(check.dir, "/", i)
      if (!dir.exists(check.dir)) next
      # Create/specify the appropriate subdirectory
      new.dir <- create_subdir(maindir, spatial, i)
      setwd(new.dir)
    } else {
      check.dir <- paste0(prevdir, "/", i)
      if (!dir.exists(check.dir)) {
        print("skip")
        next
      } 
      setwd(check.dir)
    }
    #==============================================================================
    #==============================================================================
    grouped <- group.me2(spatial, i, data_set, calc.date)
    #==========================================================================
    jv <- grouped[, c("BIOREGION", "MEASURE", "ORIGINAL_CE", "MEAN_SIM_VALUE", "RMSE")]
    jv <- jv[jv$MEASURE %in% c("VAL_CE", "TRAIN_CE"), ]
    #names(jv) <- c("Bioregion", "Measure", "Original CE", "Mean Simulated", "RMSE")
    jv$RANK <- i
    data.list[[i]] <- jv
  }
  ord <- data.list$ORDER
  fam <- data.list$FAMILY
  gen <- data.list$GENUS
  jack.df <- rbind(ord, fam, gen)
  jack.df$RANK <- factor(jack.df$RANK, levels = c("ORDER", "FAMILY", "GENUS"))
  #============================================================================
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  for (j in unique(jack.df$MEASURE)) {
    sub.jack <- jack.df[jack.df$MEASURE %in% j, ]
    
    
    
    if (spatial %in% "BASIN") {
      p <- ggplot(sub.jack, aes(x = RANK, y = MEAN_SIM_VALUE)) + 
        geom_point(size = 2, colour = "blue") + 
        #ggtitle(bioregion) + 
        xlab("Taxonomic Tier") + ylab("Classification Efficiency") + 
        scale_colour_manual(values = cbPalette) +
        
        geom_errorbar(aes(ymax = if(rmse == TRUE){
          ifelse(MEAN_SIM_VALUE + RMSE > 100, 100, MEAN_SIM_VALUE + RMSE) 
        }else{
          UP.CI
        }  , ymin = if(rmse == TRUE){
          MEAN_SIM_VALUE - RMSE
        }else{
          LOW.CI
        } )) +
        
        theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5),
              text = element_text(size = 10)) +
        ylim(40, 100) +
        geom_point(data = sub.jack, aes(x = RANK, y = ORIGINAL_CE, colour = Z), size = 2, shape = 17, colour = "red")
      
      p + facet_wrap(~ BIOREGION) + theme(strip.text.x = element_text(size = 10,
                                                                      face = "bold"),
                                          strip.background = element_blank())
      ggsave(file = paste(paste(spatial,   j, todays_date, sep = "_"), ".png", sep = ""),
             units = "in", 
             width = 3.75, 
             height = 3)
    }
    if (spatial %in% "REGION") {
      p <- ggplot(sub.jack, aes(x = RANK, y = MEAN_SIM_VALUE)) + 
        geom_point(size = 2, colour = "blue") + 
        #ggtitle(bioregion) + 
        xlab("Taxonomic Tier") + ylab("Classification Efficiency") + 
        scale_colour_manual(values = cbPalette) +
        
        geom_errorbar(aes(ymax = if(rmse == TRUE){
          ifelse(MEAN_SIM_VALUE + RMSE > 100, 100, MEAN_SIM_VALUE + RMSE) 
        }else{
          UP.CI
        }  , ymin = if(rmse == TRUE){
          MEAN_SIM_VALUE - RMSE
        }else{
          LOW.CI
        } )) +
        
        theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5),
              text = element_text(size = 10)) +
        ylim(40, 100) +
        geom_point(data = sub.jack, aes(x = RANK, y = ORIGINAL_CE, colour = Z), size = 2, shape = 17, colour = "red")
      
      p + facet_wrap(~ BIOREGION, scales = "free") + theme(panel.spacing.x = unit(2, "lines"),
                                                           strip.text.x = element_text(size = 10,
                                                                                       face = "bold"),
                                                           strip.background = element_blank())
      ggsave(file = paste(paste(spatial,  j,  todays_date, sep = "_"), ".png", sep = ""),
             units = "in", 
             width = 6.2, 
             height = 3)
    }
    if (spatial %in% "BIOREGION") {
      sub.jack1 <- sub.jack[sub.jack$BIOREGION %in% c("NAPU", "NCA", "CA", "NRV", "SGV", "SRV"), ]
      p1 <- ggplot(sub.jack1, aes(x = RANK, y = MEAN_SIM_VALUE)) + 
        geom_point(size = 2, colour = "blue") + 
        #ggtitle(bioregion) + 
        xlab("Taxonomic Tier") + ylab("Classification Efficiency") + 
        scale_colour_manual(values = cbPalette) +
        
        geom_errorbar(aes(ymax = if(rmse == TRUE){
          ifelse(MEAN_SIM_VALUE + RMSE > 100, 100, MEAN_SIM_VALUE + RMSE) 
        }else{
          UP.CI
        }  , ymin = if(rmse == TRUE){
          MEAN_SIM_VALUE - RMSE
        }else{
          LOW.CI
        } )) +
        
        theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5),
              text = element_text(size = 10)) +
        ylim(40, 100) +
        geom_point(data = sub.jack1, aes(x = RANK, y = ORIGINAL_CE, colour = Z), size = 2, shape = 17, colour = "red")
      
      p1 + facet_wrap(~ BIOREGION, nrow = 3, scales = "free") + theme(panel.spacing.x = unit(2, "lines"),
                                                                      panel.spacing.y = unit(1, "lines"),
                                                                      strip.text.x = element_text(size = 10,
                                                                                                  face = "bold"),
                                                                      strip.background = element_blank())
      name.num <- if ("CA" %in% sub.jack1$BIOREGION) 1 else 2
      ggsave(file=paste(paste(spatial,  j, name.num, todays_date, sep = "_"), ".png", sep = ""),
             units = "in", 
             width = 6.2, 
             height = 7)
      
      sub.jack2 <- sub.jack[sub.jack$BIOREGION %in% c("BLUE", "UNP", "LNP", "PIED", "SEP", "MAC"), ]
      p2 <- ggplot(sub.jack2, aes(x = RANK, y = MEAN_SIM_VALUE)) + 
        geom_point(size = 2, colour = "blue") + 
        #ggtitle(bioregion) + 
        xlab("Taxonomic Tier") + ylab("Classification Efficiency") + 
        scale_colour_manual(values = cbPalette) +
        
        geom_errorbar(aes(ymax = if(rmse == TRUE){
          ifelse(MEAN_SIM_VALUE + RMSE > 100, 100, MEAN_SIM_VALUE + RMSE) 
        }else{
          UP.CI
        }  , ymin = if(rmse == TRUE){
          MEAN_SIM_VALUE - RMSE
        }else{
          LOW.CI
        } )) +
        
        theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5),
              text = element_text(size = 10)) +
        ylim(40, 100) +
        geom_point(data = sub.jack2, aes(x = RANK, y = ORIGINAL_CE, colour = Z), size = 2, shape = 17, colour = "red")
      
      p2 + facet_wrap(~ BIOREGION, nrow = 3, scales = "free") + theme(panel.spacing.x = unit(2, "lines"),
                                                                      panel.spacing.y = unit(1, "lines"),
                                                                      strip.text.x = element_text(size = 10,
                                                                                                  face = "bold"),
                                                                      strip.background = element_blank())
      name.num <- if ("CA" %in% sub.jack2$BIOREGION) 1 else 2
      ggsave(file=paste(paste(spatial,  j, name.num, todays_date, sep = "_"), ".png", sep = ""),
             units = "in", 
             width = 6.2, 
             height = 7)
    }
    
    
  }
  
  
  
}



#==============================================================================
#'BIBI Figures
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

bibi_figures <- function(spatial, taxon_rank, calc.date,
                         todays_date = format(Sys.time(), "%m_%d_%y")){
  grouped.score <- group.me2(spatial, taxon_rank, "scored_metrics", calc.date, final_score = T)
  grouped.score$SCORE <- grouped.score$FINAL_SCORE * 100
  #============================================================================
  grouped.jack <- group.me2(spatial, taxon_rank, "jackknife", calc.date)
  if(any(names(grouped.jack) %in% "ORGINAL_BSP")){
    names(grouped.jack)[names(grouped.jack) %in% "ORGINAL_BSP"] <- "ORIGINAL_BSP"
  } 
  grouped.jack.sub <- unique(grouped.jack[, c("BIOREGION", "ORIGINAL_CE", "ORIGINAL_BSP")])
  #============================================================================
  plot.this <- merge(grouped.jack.sub, grouped.score, by = "BIOREGION")
  
  #============================================================================
  if (spatial %in% "BASIN") plot.me(plot.this, spatial, taxon_rank, todays_date)

  if (spatial %in% "REGION") plot.me(plot.this, spatial, taxon_rank, todays_date)

  
  if (spatial %in% "BIOREGION") {
    plot.this.1 <- plot.this[plot.this$BIOREGION %in% c("NAPU", "NCA", "CA", "NRV", "SGV", "SRV"), ]

    plot.me(plot.this.1, spatial, taxon_rank, todays_date)

    #==========================================================================
    plot.this.2 <- plot.this[plot.this$BIOREGION %in% c("BLUE", "UNP", "LNP", "PIED", "SEP", "MAC"), ]

    plot.me(plot.this.2, spatial, taxon_rank, todays_date)

  }
  
}


#==============================================================================
#'Group Me 3
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

group.me3 <- function(front, proc_date, data_set){
  tidy.me <- function(bioregion) tidyr::gather(bioregion, METRIC, SCORE, 8:ncol(bioregion))
  blue <- read.csv(paste(front, "BLUE", data_set, proc_date, ".csv", sep = ""))
  blue$X <- "BLUE"
  blue <- tidy.me(blue)
  ca <- read.csv(paste(front, "CA", data_set, proc_date, ".csv", sep = ""))
  ca$X <- "CA"
  ca <- tidy.me(ca)
  lnp <- read.csv(paste(front, "LNP", data_set, proc_date, ".csv", sep = ""))
  lnp$X <- "LNP"
  lnp <- tidy.me(lnp)
  mac <- read.csv(paste(front, "MAC", data_set, proc_date, ".csv", sep = ""))
  mac$X <- "MAC"
  mac <- tidy.me(mac)
  napu <- read.csv(paste(front, "NAPU", data_set, proc_date, ".csv", sep = ""))
  napu$X <- "NAPU"
  napu <- tidy.me(napu)
  nca <- read.csv(paste(front, "NCA", data_set, proc_date, ".csv", sep = ""))
  nca$X <- "NCA"
  nca <- tidy.me(nca)
  nrv <- read.csv(paste(front, "NRV", data_set, proc_date, ".csv", sep = ""))
  nrv$X <- "NRV"
  nrv <- tidy.me(nrv)
  pied <- read.csv(paste(front, "PIED", data_set, proc_date, ".csv", sep = ""))
  pied$X <- "PIED"
  pied <- tidy.me(pied)
  sep <- read.csv(paste(front, "SEP", data_set, proc_date, ".csv", sep = ""))
  sep$X <- "SEP"
  sep <- tidy.me(sep)
  sgv <- read.csv(paste(front, "SGV", data_set, proc_date, ".csv", sep = ""))
  sgv$X <- "SGV"
  sgv <- tidy.me(sgv)
  srv <- read.csv(paste(front, "SRV", data_set, proc_date, ".csv", sep = ""))
  srv$X <- "SRV"
  srv <- tidy.me(srv)
  unp <- read.csv(paste(front, "UNP", data_set, proc_date, ".csv", sep = ""))
  unp$X <- "UNP"
  unp <- tidy.me(unp)
  
  bound <- rbind(blue, ca, lnp, mac, napu, nca, nrv, pied, sep, sgv, srv, unp)
  
  if("ORGINAL_VALUE" %in% names(bound) | "ORGINAL_THRESHOLD" %in% names(bound)){
    names(bound)[names(bound) %in% "ORGINAL_VALUE"|
                   names(bound) %in% "ORGINAL_THRESHOLD"] <- "ORIGINAL_VALUE"
  } 
  
  if("MEAN_SIM_THRESHOLD" %in% names(bound)){
    names(bound)[names(bound) %in% "MEAN_SIM_THRESHOLD"] <- "MEAN_SIM_VALUE"
  } 
  return(bound)
}

#==============================================================================
#'Precentile Summary
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

percentile_summary <- function(fam.final, vp_thresh, outlier = TRUE){
  
  
  fam.ref <- fam.final[fam.final$CATEGORY %in% "REF", ]
  
  if (outlier == TRUE) {
    outlier <- function(scores.df, mult.val, job = "REF"){
      datalist = list()
      for (i in unique(scores.df$BIOREGION)) {
        bioregion.df <- scores.df[scores.df$BIOREGION %in% i, ]
        
        quantile.vec <- quantile(bioregion.df$FINAL_SCORE, c(0.25, 0.75))
        if (job %in% "REF"){
          outlier.thresh <- quantile.vec[1] - (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
        }
        if (job %in% "SEV"){
          outlier.thresh <- quantile.vec[2] + (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
        }
        
        outlier.df <- data.frame(BIOREGION = i)
        outlier.df$THRESH <- outlier.thresh
        datalist[[i]] <- outlier.df # add it to your list
      }
      final.df <- do.call(rbind, datalist)
      return(final.df)
    }
    
    out.df <- outlier(fam.ref, 1.5, job = "REF")
    
    fam.ref<- merge(fam.ref, out.df, by = "BIOREGION")
    fam.ref$OUTLIER <- ifelse(fam.ref$FINAL_SCORE >= fam.ref$THRESH, "NO", "YES")
    fam.ref <- fam.ref[fam.ref$OUTLIER %in% "NO", ]
    
    quant.df <- data.frame(do.call("rbind", tapply(fam.ref$FINAL_SCORE, fam.ref$BIOREGION,
                                                   quantile, c(0.1, 0.25, 0.5))))
  } else {
    quant.df <- tapply(fam.ref$FINAL_SCORE, fam.ref$BIOREGION,
           quantile, c(0.1, 0.25, 0.5))
  }
  
  
  library(plyr)
  #do.call("rbind", tapply(fam.ref$FINAL_SCORE, fam.ref$BIOREGION, quantile))
  
  
  names(quant.df) <- c("%10", "%25", "%50")
  quant.df$BIOREGION <- row.names(quant.df)
  quant.df <- quant.df[, c(ncol(quant.df), 1:ncol(quant.df[, -1]))]
  
  
  
  median.calc <- function(scores.df){
    datalist = list()
    for(i in unique(scores.df$BIOREGION)){
      bioregion.df <- scores.df[scores.df$BIOREGION %in% i, ]
      deg.df <- data.frame(BIOREGION = i)
      out.deg <- outlier(bioregion.df, 1.5, job = "SEV")
      
      fam.deg <- merge(bioregion.df, out.deg, by = "BIOREGION")
      fam.deg$OUTLIER <- ifelse(fam.deg$FINAL_SCORE <= fam.deg$THRESH, "NO", "YES")
      fam.deg <- fam.deg[fam.deg$OUTLIER %in% "NO", ]
      
      deg.df$MEDIAN <- quantile(fam.deg$FINAL_SCORE, c(0.5))
      deg.df$COUNT <- nrow(fam.deg)
      datalist[[i]] <- deg.df # add it to your list
    }
    final.df <- do.call(rbind, datalist)
    return(final.df)
  }
  if(vp_thresh %in% "DEG_50"){
    deg <- fam.final[fam.final$CATEGORY %in% "SEV", ]
    med.deg <- median.calc(deg)
    final.df <- merge(quant.df, med.deg[, -3], by = "BIOREGION")
  }
  
  if(vp_thresh %in% "half_10"){
    quant.df$DEG_50 <- quant.df$`%10` / 2
    final.df <- quant.df
  }
  names(final.df) <- c("BIOREGION", "REF_10", "REF_25", "REF_50", "DEG_50")
  
  
  return(final.df)
}

#==============================================================================
#'Rate Me 2
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

rate.me2 <- function(scores.df, percentiles, taxa.rank = "FAMILY", napu = TRUE){
  final.df <- merge(scores.df, percentiles, by = "BIOREGION", all.x = T)
  if(taxa.rank %in% c("ORDER", "FAMILY", "GENUS")){
    final.df$RATING <- ifelse (final.df$FINAL_SCORE < final.df$DEG_50, "VeryPoor",
                               ifelse (final.df$FINAL_SCORE < final.df$REF_10, "Poor",
                                       ifelse (final.df$FINAL_SCORE < final.df$REF_25, "Fair",
                                               ifelse (final.df$FINAL_SCORE < final.df$REF_50, "Good",
                                                       ifelse (final.df$FINAL_SCORE >= final.df$REF_50 &
                                                                 final.df$FINAL_SCORE <= 100,
                                                               "Excellent", "ERROR")))))
  }
  
  #if(taxa.rank %in% c("ORDER")){
  #  final.df$RATING <- ifelse (final.df$FINAL_SCORE < final.df$DEG_50, "Poor",
  #                             ifelse (final.df$FINAL_SCORE < final.df$REF_50, "Fair",
  #                                     ifelse (final.df$FINAL_SCORE >= final.df$REF_50 &
  #                                               final.df$FINAL_SCORE <= 100,
  #                                             "Good", "ERROR")))
  #}
  
  
  if (napu == FALSE) final.df$RATING <- ifelse(final.df$BIOREGION %in% "NAPU", "TBD", as.character(final.df$RATING))
  return(final.df)
}

#==============================================================================
#'Random Events
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

random_events <- function(tab.events){
  # Returns string w/o trailing whitespace.
  tab_event$SITE_TYPE_CODE <- trimws(tab_event$SITE_TYPE_CODE)
  # Select only the sampling events designated random sampling locations.
  rand_event <- tab_event[tab_event$SITE_TYPE_CODE %in% c("R", "RR", "TS"), ]
  return(rand_event)
}

#==============================================================================
#'Group Me Low Resolution
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

group.me.low.res <- function(front, proc_date, data_set, job){
  #==============================================================================
  tidy.me <- function(bioregion) tidyr::gather(bioregion, METRIC, SCORE, 8:ncol(bioregion))
  #==============================================================================
  
  bioregions <- c("COAST", "NON_COAST") #, "BASIN")
  
  datalist <- list()
  for (i in bioregions){
    bio <- read.csv(paste(front, data_set, i, proc_date, ".csv", sep = ""))
    bio <- bio[, !names(bio) %in% "X"]
    bio$BIOREGION <- i
    bio <- bio[, c(ncol(bio), (1:(ncol(bio) - 1)))]
    if(job %in% "METRICS") bio <- tidy.me(bio)
    datalist[[i]] <- bio # add it to your list
  }
  low.res <- unique(do.call(rbind, datalist))
  #==============================================================================
  bound <- low.res
  bound$BIOREGION <- ifelse(bound$BIOREGION %in% "NON_COATAL", "INLAND", as.character(bound$BIOREGION))
  #==============================================================================
  fam.final <- bound[bound$METRIC %in% "FINAL_SCORE", ]
  names(fam.final)[1] <- "BIOREGION"
  fam.final$SCORE <- fam.final$SCORE * 100
  return(fam.final)
}

#==============================================================================
#'Category Counts
#'
#'@param spatial = "BASIN", "REGION", or "BIOREGION".
#'@param taxon.rank = "ORDER", "FAMILY", or "GENUS".
#'@param todays.date = todays date.
#'@return Provides a table of Category counts (Ref, Min, Mix, Mod, Deg) by
#' bioregion.
#'@export

category_counts <- function(spatial, taxon.rank, todays.date) {
  #setwd(working.dr)
  file.name <- paste(spatial, taxon.rank, todays.date,
                     "All_Event_Ratings.csv", sep = "_")
  my.df <- read.csv(file.name, stringsAsFactors = FALSE)
  #----------------------------------------------------------------------------
  new.df <- as.data.frame(table(my.df$BIOREGION, my.df$CATEGORY))
  new.df <- tidyr::spread(new.df, Var2, Freq)
  names(new.df)[1] <- "Bioregion"
  new.df <- new.df[, c("Bioregion", "REF", "MIN", "MIX", "MOD", "SEV")]
  names(new.df)[names(new.df) %in% "SEV"] <- "DEG"
  new.df$'Total Count' <- rowSums(new.df[, 2:ncol(new.df)])
  #----------------------------------------------------------------------------
  pct.df <- new.df[, c("Bioregion", "Total Count")]
  cat.cols <- c("REF", "MIN", "MIX", "MOD", "DEG")
  pct.df[, cat.cols] <- lapply(cat.cols, function(x){
    new.df[, x] / new.df$`Total Count` * 100
  })
  
  pct.df[, cat.cols] <- apply(pct.df[, cat.cols], 2, function(x) {
    paste0("(", round(x, 1), "%)")
  }) 
  
  pct.df <- pct.df[, c("Bioregion", "REF", "MIN", "MIX", "MOD", "DEG", "Total Count")]
  
  #----------------------------------------------------------------------------
  final.df <- rbind(new.df, pct.df)
  final.df <- final.df[order(final.df$Bioregion), ]
  #----------------------------------------------------------------------------
  file.name <- paste(spatial, taxon.rank, todays.date,
                     "Category_Count.csv", sep = "_")
  write.csv(final.df, file.name, row.names = FALSE)
}
#==============================================================================
#'Year Counts
#'
#'@param spatial = "BASIN", "REGION", or "BIOREGION".
#'@param taxon.rank = "ORDER", "FAMILY", or "GENUS".
#'@param todays.date = todays date.
#'@return Provides a a bar plot figure of sample counts by year.
#'@export

year_counts_barplot <- function(spatial, taxon.rank, todays.date) {
  file.name <- paste(spatial, taxon.rank, todays.date,
                     "All_Event_Ratings.csv", sep = "_")
  my.df <- read.csv(file.name, stringsAsFactors = FALSE)
  #----------------------------------------------------------------------------
  my.df$DATE <- as.Date(my.df$DATE, "%m/%d/%Y")
  my.df$YEAR <- format(my.df$DATE, "%Y")
  #----------------------------------------------------------------------------
  final.df <- aggregate(EVENT_ID ~ YEAR, data = my.df, length)
  names(final.df) <- c("YEAR", "COUNT")
  #----------------------------------------------------------------------------
  ggplot() +
    geom_bar(data = final.df, aes(x = YEAR, y = COUNT), stat="identity", fill = "#0072B2") +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          # text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),,
          #plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
          #axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
          plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
    xlab("Year") +
    ylab("Sample Count") +
    scale_y_continuous(limits = c(0, 2000), expand = c(0, 0))
  ggsave(file = paste0(paste(spatial, taxon.rank, todays.date, sep = "_"), "_Sample_Count.png"),
         units = "in", 
         width = 6.2, 
         height = 4)
}

#==============================================================================
#' 2011 Ratings vs. 2016 Ratings
#'
#'@param curr.dir = the directory to the 2016 data.
#'@param todays.date = todays date.
#'@return Provides two tables: 1) Coastal rating comparisons and 
#'2) Inland rating comparisons.
#'@export

bibi_2011_vs_2016 <- function(curr.dir, spatial, taxon.rank, todays.date) {
  setwd("D:/ZSmith/Projects/Chessie_BIBI/2011_BIBI_Data")
  bibi_2011 <- read.csv("ALL_AVAILABLE_SITE_BIBI_DATA.csv",
                        stringsAsFactors = FALSE)
  bibi_2011 <- unique(bibi_2011[, c("EVENT_ID", "cbp_RANK")])
  names(bibi_2011) <- c("EVENT_ID", "RANK_2011")
  table(bibi_2011$RANK_2011)
  bibi_2011$RANK_2011 <- ifelse(bibi_2011$RANK_2011 %in% "EXCELLENT", "Excellent",
                                ifelse(bibi_2011$RANK_2011 %in% "GOOD", "Good",
                                       ifelse(bibi_2011$RANK_2011 %in% "FAIR", "Fair",
                                              ifelse(bibi_2011$RANK_2011 %in% "POOR", "Poor",
                                                     ifelse(bibi_2011$RANK_2011 %in% "VERY_POOR",
                                                            "Very Poor", "ERROR")))))
  #------------------------------------------------------------------------------
  setwd(curr.dir)
  file.2016 <- paste(spatial, taxon.rank, todays.date, "All_Event_Ratings.csv",
                     sep = "_")
  bibi_2016 <- read.csv(file.2016, stringsAsFactors = FALSE)
  bibi_2016 <- unique(bibi_2016[, c("EVENT_ID", "SAMPLE_NUMBER", "RATING", "BIOREGION")])
  bibi_2016$RATING <- ifelse(bibi_2016$RATING %in% "VeryPoor", "Very Poor", bibi_2016$RATING)
  table(bibi_2016$RATING)
  #==============================================================================
  merged.df <- merge(bibi_2016, bibi_2011, by = "EVENT_ID")
  merged.df$RATING <- factor(merged.df$RATING, levels = c("Excellent", "Good",
                                                          "Fair", "Poor",
                                                          "Very Poor"))
  merged.df$RANK_2011 <- factor(merged.df$RANK_2011, levels = c("Excellent", "Good",
                                                                "Fair", "Poor",
                                                                "Very Poor"))
  if (spatial %in% "REGION") {
    coast.df <- merged.df[merged.df$BIOREGION %in% "COAST", ]
    inland.df <- merged.df[merged.df$BIOREGION %in% "INLAND", ]
  }  
  if (spatial %in% "BIOREGION") {
    coast.df <- merged.df[merged.df$BIOREGION %in% c("MAC", "SEP"), ]
    inland.df <- merged.df[!merged.df$BIOREGION %in% c("MAC", "SEP"), ]
  }
  #==============================================================================
  wide.table <- function(my.df) {
    agg.df <- aggregate(EVENT_ID ~ RATING + RANK_2011, data = my.df, length)
    
    wide.df <- tidyr::spread(agg.df, RATING, EVENT_ID)
    wide.df[is.na(wide.df)] <- 0
    wide.df[, 2:ncol(wide.df)] <- lapply(wide.df[, 2:ncol(wide.df)], function(x){
      paste0(round(x / sum(wide.df[, 2:ncol(wide.df)]) * 100, 1), "%")
    })
    
    names(wide.df)[1] <- "RATING"
    return(wide.df)
  }
  #==============================================================================
  coast.table <- wide.table(coast.df)
  inland.table <- wide.table(inland.df)
  #==============================================================================
  write.csv(coast.table, paste(spatial, taxon.rank, todays.date,
                               "coast_2011_vs_2016.csv", sep = "_"),
            row.names = FALSE)
  write.csv(inland.table, paste(spatial, taxon.rank, todays.date,
                                "inland_2011_vs_2016.csv", sep = "_"),
            row.names = FALSE)
}

#==============================================================================
#'All Tables
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

all_tables <- function(spatial, calc.date, m.c,
                       all.data, tab.event,
                       todays_date = format(Sys.time(), "%m_%d_%y"),
                       random_event, main.dir,
                       previous.subdir, map.agg.col){
  rate_stacks("REGION",  sub("/[^/]*$", "", previous.subdir),
              calc.date, todays_date, very.poor.thresh = "half_10",
              random_event)
  rate_stacks("BIOREGION",  sub("/[^/]*$", "", previous.subdir), 
              calc.date, todays_date, very.poor.thresh = "half_10",
              random_event)
  
  jack_table(spatial, "jackknife", calc.date,
             todays_date = format(Sys.time(), "%m_%d_%y"), maindir = main.dir,
             prevdir = previous.subdir)
  plot_jack(spatial, "jackknife",
            calc.date, todays_date = format(Sys.time(), "%m_%d_%y"), maindir = main.dir,
            prevdir = previous.subdir)

  for(i in c("ORDER", "FAMILY", "GENUS")){
    print(i)
    #==============================================================================
    if (is.null(previous.subdir)) {
      month.year <- format(Sys.Date(), format="%B_%Y")
      todays.date <- format(Sys.Date(), format="%m_%d_%Y")
      check.dir <- paste0(main.dir, "/", month.year)
      check.dir <- paste0(check.dir, "/", todays.date)
      check.dir <- paste0(check.dir, "/", spatial)
      check.dir <- paste0(check.dir, "/", i)
      if (!dir.exists(check.dir)) {
        print("skip")
        next
      } 
      # Create/specify the appropriate subdirectory
      new.dir <- create_subdir(main.dir, spatial, i)
      setwd(new.dir)
    } else {
      new.dir <- paste0(previous.subdir, "/", i)
      if (!dir.exists(new.dir)) {
        print("skip")
        next
      } 
      setwd(new.dir)
    }
    
    #==============================================================================
    metric_table(spatial, i, "ms", calc.date, m.c,
                 todays_date = format(Sys.time(), "%m_%d_%y"))
    index_table(spatial, i, "jackknife", calc.date,
                todays_date = format(Sys.time(), "%m_%d_%y"))
    bibi_figures(spatial, taxon_rank = i, calc.date, todays_date)
    map_this(spatial, calc.date, taxon.rank = i, all.data, tab.event, todays_date,
             title.me = paste(i, "-level ", spatial, " Ratings (Random)", sep = ""),
             very.poor.thresh = "half_10", random_event, sub.dir = new.dir,
             map.agg.col)
    category_counts(spatial, i, todays_date)
    year_counts_barplot(spatial, i, todays_date)
    if(i %in% "FAMILY" & spatial %in% c("REGION", "BIOREGION")) {
      bibi_2011_vs_2016(curr.dir = new.dir, spatial, i, todays_date)
    }
  }
  

  
}



