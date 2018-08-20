#==============================================================================
#'Prep Rate Stacks
#'
#'@param spatial = "BASIN", "REGION", or "BIOREGION".
#'@param previous.subdir = specify the location of the previously calculated data.
#'@return Prepares the data for the stacked barplots of the ratings.
#'@export


prep_rate_stacks <- function(spatial, previous.subdir, calc.date,
                             very.poor.thresh, random_event){
  taxon.list <- lapply(c("ORDER", "FAMILY", "GENUS"), function(x){
    #test <- lapply(c("BASIN", "REGION", "BIOREGION"), function(spat.x){
    spat.x <- spatial
    if (is.null(previous.subdir)) {
      month.year <- format(Sys.Date(), format="%B_%Y")
      todays.date <- format(Sys.Date(), format="%m_%d_%Y")
      check.dir <- paste0(main.dir, "/", month.year)
      check.dir <- paste0(check.dir, "/", todays.date)
      check.dir <- paste0(check.dir, "/", spat.x)
      check.dir <- paste0(check.dir, "/", i)
      if (!dir.exists(check.dir)) {
        print("skip")
        next
      } 
      # Create/specify the appropriate subdirectory
      new.dir <- create_subdir(main.dir, spat.x, x)
      setwd(new.dir)
    } else {
      new.dir <- paste0(previous.subdir, "/", spat.x, "/", x)
      if (!dir.exists(new.dir)) {
        print("skip")
        next
      } 
      setwd(new.dir)
    }
    # Import the appropriate data based on the specified spatial resolution (spatial),
    # taxonomic resolution (taxon.rank), and the date the scores were calculated (calc.date).
    final.scores <- group.me2(spat.x, x, "scored_metrics", calc.date, final_score = TRUE)
    #==============================================================================
    # Establish the rating thresholds based on the Reference distribution.
    # Outliers in this case reflect ONLY the Reference distribution.
    pct_sum <- percentile_summary(final.scores, very.poor.thresh, outlier = TRUE)
    #==============================================================================
    # random_events was built specifically to select only randomly sampled events
    # according to the SITE_TYPE_CODE in the Chessie BIBI database.
    # It appears that there are currently (3/20/17) inaccuracies associated with
    # the SITE_TYPE_CODE; therefore, random_events should alwasy be set to FALSE
    # until the randomly selected sampling events can be verified.
    if (random_event == TRUE) rand.events <- random_events(tab.events)
    #==============================================================================
    # Apply the rating scheme.
    #all.scores <- rate.me2(final.scores, pct_sum, x) # Do NOT use rate.me2
    all.scores <- percentile_summary(final.scores, "half_10", outlier = TRUE)
    
    #==============================================================================
    #ggplot(all.scores, aes(BIOREGION, c(0:100))) +
    sub.scores <- unique(all.scores[, c("BIOREGION", "REF_10", "REF_25",
                                        "REF_50", "DEG_50")])
    rate.thresh <- tidyr::gather(sub.scores, THRESH_NAME, THRESH, REF_10:DEG_50)
    #rate.thresh$THRESH_NAME <- factor(rate.thresh$THRESH_NAME)
    #levels(rate.thresh$THRESH_NAME) <- c("DEG_50", "REF_10", "REF_25", "REF_50")
    rate.thresh$THRESH <- rate.thresh$THRESH * 100
    #==============================================================================
    todays.date <- format(Sys.time(), "%m_%d_%y")
    write.csv(rate.thresh, paste0(paste(x, "Threshold", spatial, todays.date, sep = "_"),
                                  ".csv"), row.names = FALSE)
    
    #==============================================================================
    rate.list <- lapply(unique(rate.thresh$BIOREGION), function(x){
      sub.thresh <- rate.thresh[rate.thresh$BIOREGION %in% x, ]
      sub.thresh <- rbind(sub.thresh, c(x, "MAX", 100))
      sub.thresh$THRESH_NAME <- factor(sub.thresh$THRESH_NAME,
                                       levels = c("DEG_50", "REF_10", "REF_25", "REF_50", "MAX"))
      sub.thresh <- sub.thresh[order(sub.thresh$THRESH_NAME), ]
      sub.thresh$THRESH <- as.numeric(sub.thresh$THRESH)
      sub.thresh$DIFFS <- c(sub.thresh[sub.thresh$THRESH_NAME %in% "DEG_50", "THRESH"], diff(sub.thresh$THRESH))
      sub.thresh$RATING <- c("Very Poor", "Poor", "Fair", "Good", "Excellent")
      final.thresh <- sub.thresh[, c("BIOREGION", "RATING", "DIFFS")]
    })
    
    rate.df <- do.call(rbind, rate.list)
  })
  #=============================================================================
  ord.df <- taxon.list[[1]]
  ord.mean <- aggregate(DIFFS ~ RATING, data = ord.df, mean)
  ord.mean$BIOREGION <- "MEAN"
  ord.mean <- ord.mean[, c("BIOREGION", "RATING", "DIFFS")]
  ord.df <- rbind(ord.df, ord.mean)
  #----------------------------------------------------------------------------
  fam.df <- taxon.list[[2]]
  fam.mean <- aggregate(DIFFS ~ RATING, data = fam.df, mean)
  fam.mean$BIOREGION <- "MEAN"
  fam.mean <- fam.mean[, c("BIOREGION", "RATING", "DIFFS")]
  fam.df <- rbind(fam.df, fam.mean)
  #----------------------------------------------------------------------------
  gen.df <- taxon.list[[3]]
  gen.mean <- aggregate(DIFFS ~ RATING, data = gen.df, mean)
  gen.mean$BIOREGION <- "MEAN"
  gen.mean <- gen.mean[, c("BIOREGION", "RATING", "DIFFS")]
  gen.df <- rbind(gen.df, gen.mean)
  #============================================================================= 
  final.list <- list(ord.df, fam.df, gen.df)
}

#==============================================================================
#'Rate Stacks
#'
#'@param spatial = "BASIN", "REGION", or "BIOREGION".
#'@param previous.subdir = specify the location of the previously calculated data.
#'@return Creates stacked barplots of the ratings.
#'@import ggplot2
#'@export

rate_stacks <- function(spatial, previous.subdir, calc.date, todays_date,
                        very.poor.thresh, random_event){
  my.list <- prep_rate_stacks(spatial, previous.subdir, calc.date,
                              very.poor.thresh, random_event)
  #if(taxon.rank %in% "ORDER") rate.df <- my.list[[1]]
  #if(taxon.rank %in% "FAMILY") rate.df <- my.list[[2]]
  #if(taxon.rank %in% "GENUS") rate.df <- my.list[[3]]
  ord.df <- my.list[[1]]
  ord.df$RANK <- "Order"
  fam.df <- my.list[[2]]
  fam.df$RANK <- "Family"
  gen.df <- my.list[[3]]
  gen.df$RANK <- "Genus"
  rate.df <- rbind(ord.df, fam.df, gen.df)
  rate.df$RATING <- factor(rate.df$RATING, levels = c("Excellent", "Good", "Fair", "Poor", "Very Poor"))
  #============================================================================
  bioregions <- c("BASIN", "COAST", "INLAND", "CA", "NAPU", "NCA", "NRV",
                  "UNP", "BLUE", "LNP", "PIED", "SGV", "SRV", "SEP", "MAC",
                  "MEAN")
  rate.df$BIOREGION <- factor(ord.df$BIOREGION,levels = bioregions)
  #ord.df$BIOREGION <- factor(ord.df$BIOREGION,levels = bioregions)
  #fam.df$BIOREGION <- factor(fam.df$BIOREGION,levels = bioregions)
  #gen.df$BIOREGION <- factor(gen.df$BIOREGION,levels = bioregions)
  #============================================================================
  rate.df$RANK <- factor(rate.df$RANK,levels = c("Order", "Family", "Genus"))
  #============================================================================
  cols <- c("darkgreen", "green3", "yellow2", "orange2", "red3")
  #============================================================================
  p <- ggplot() + geom_bar(aes(y = DIFFS, x = BIOREGION, fill = RATING),
                           width = 0.9,
                           data = rate.df,
                           stat="identity",
                           color = "black") +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          text = element_text(size = 12),
          #axis.text.x = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.title.x = element_blank(),
          #axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
          #plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
          legend.position = "bottom") +
    guides(fill = guide_legend(title = element_blank())) +
    #xlab("Scoring Method") +
    ylab("Index Score") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = cols)
  #============================================================================  
  setwd(previous.subdir)
  #============================================================================  
  if (spatial %in% "REGION") {
    p  + facet_grid(. ~RANK, scale="free") + 
      theme(strip.text.y = element_text(size = 10, face = "bold"), 
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank())# +
    #annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    #annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
    
    ggsave(file = paste(paste(spatial, "Stacked_Ratings", todays_date, sep = "_"),
                        ".png", sep = ""),
           units = "in", 
           width = 6.3, 
           height = 3.5)
  }
  #============================================================================  
  if (spatial %in% "BIOREGION") {
    p  + facet_grid(RANK ~ ., scale="free") + 
      theme(strip.text.y = element_text(size = 10, face = "bold"), 
            panel.spacing = unit(1.3, "lines"),
            strip.background = element_blank()) +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
    
    ggsave(file = paste(paste(spatial, "Stacked_Ratings", todays_date, sep = "_"),
                        ".png", sep = ""),
           units = "in", 
           width = 6.3, 
           height = 7.2)
  }
}


