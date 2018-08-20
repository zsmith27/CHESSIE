#==============================================================================
#==============================================================================
# Author: Zachary M. Smith
# Maintainer: Zachary M. Smith
# Organization: ICPRB
# Created: December 2016
# Updated: 3-14-2017
# Purpose: Iterative development of an index.
#==============================================================================
#==============================================================================
#'Iterative Index Development
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return 
#'@export

iterate_dev <- function(all.data, index.res, taxon.rank, runs = 10,
                        todays_date = format(Sys.time(), "%m_%d_%y"),
                        unique.station = FALSE, standard_count = FALSE,
                        main.dir = "D:/ZSmith/Projects/Chessie_BIBI/Output"){
  #==============================================================================
  # Prep for subdirectory
  if (index.res %in% "BASIN") spat.res <- "BASIN"
  if (index.res %in% c("COAST", "INLAND")) spat.res <- "REGION"
  if (index.res %in% c("BLUE", "CA", "LNP", "MAC", "NAPU", "NCA",
                       "NRV", "PIED", "SEP", "SGV", "SRV", "UNP",
                       "BIOREGION")) spat.res <- "BIOREGION"
  # Create the appropriate subdirectory.
  new.dir <- create_subdir(main.dir, spat.res, taxon.rank)
  # Set the working directory to the appropriate subdirectory.
  setwd(new.dir)
  #==============================================================================
  prep.bio <- function(my.data, index_res){

    my.data$BIOREGION <- ifelse(my.data$ICPRB_BIOREGION_ID %in% "NGV", "UNP",
                                as.character(my.data$ICPRB_BIOREGION_ID))
    if(index_res %in% "BASIN") my.data$BIOREGION <- "ALL"
    if(index_res %in% "COAST"){
      #my.data <- my.data[my.data$BIOREGION %in% c("MAC", "SEP"), ]
      my.data <- my.data[my.data$ICPRB_BIOREGION_ID %in% c("MAC", "SEP"), ]
      my.data$BIOREGION <- "COAST"
    }
    if(index_res %in% "INLAND"){
      #my.data <- my.data[my.data$BIOREGION %in% c("CA", "NRV", "UNP", "LNP",
      #                                            "SRV", "NCA", "NAPU", "SGV",
      #                                            "BLUE", "PIED"), ]
      my.data <- my.data[my.data$ICPRB_BIOREGION_ID %in% c("CA", "NRV", "UNP", "LNP",
                                                  "SRV", "NCA", "NAPU", "SGV",
                                                  "BLUE", "PIED"), ]
      my.data$BIOREGION <- "INLAND"
    }
    return(my.data)
  }
  
  #============================================================================
  test.df <- data.frame(unique(all.data$EVENT_ID))
  env.df <- prep_env(Prep.data = all.data)
  test.class <- prep_class_multi(env.df)
  all.data <- prep.bio(all.data, index.res)
  #names(all.data)[names(all.data) %in% "ICPRB_BIOREGION_ID"] <- "BIOREGION"
  all.data <- all.data[!is.na(all.data$BIOREGION), ]
  #table(test.class$BIOREGION, test.class$CATEGORY)
  env.class <- unique(test.class[, c("EVENT_ID", "STATION_ID", "DATE",
                                     "AGENCY_CODE", "SAMPLE_NUMBER", "CATEGORY")])
  #==========================================================================
  # Keep only the REF and SEV rows to reduce computation time.
  # These rows are the only rows necessary for the selection of metrics.
  initial.merge <- merge(all.data, env.class, by = c("EVENT_ID", "STATION_ID",
                                                   "DATE", "AGENCY_CODE",
                                                   "SAMPLE_NUMBER"),
                         all.x = TRUE)
  # The following categories are not necessary for metric selection and
  # are removed to speed up the script.
  ref.sev.df <- initial.merge[!initial.merge$CATEGORY %in% c("MIN", "MIX", "MOD",
                                                          "DEG", "DEG+"), ]
  # Remove the CATEGORY column so that it does not become dublicated during
  # subsequent merges.
  ref.sev.df <- ref.sev.df[, !names(ref.sev.df) %in% "CATEGORY"]
  #==========================================================================
  

  iterate_this <- function(bio.df, index.res, taxon.rank, runs, todays_date, env.class,
                           unique.station = FALSE, standard_count = FALSE){

    data.list <- list()
    data.list2 <- list()
    for(j in seq(1, runs, 1)){
      print(paste("Start Iteration:", j))
      #gc()
      
      # Prepare data for site classification
      bio.data <- prep.bio(bio.df, "BIOREGION")
      if (unique.station == TRUE) {
        one_event_per_station <- function(my.df){
          
          my.df$EVENT_ID <- paste(my.df$EVENT_ID, my.df$SAMPLE_NUMBER, sep = "_")
          my.df <- unique(my.df[, c("EVENT_ID", "STATION_ID", "SAMPLE_NUMBER", "AGENCY_CODE", "DATE")])
          
          new.list <- split(my.df$EVENT_ID, my.df$STATION_ID)
          
          new.short <- new.list[lengths(new.list) == 1]
          #test <- data.frame(unlist(new.short))
          new.long <- new.list[lengths(new.list) > 1]
          rand.long <- lapply(new.long, function(x) sample(x, 1))
          
          ns <- data.frame(EVENT_ID = unlist(new.short))
          nr <- data.frame(EVENT_ID = unlist(rand.long))
          final.df <- rbind(ns, nr)
          
          return(final.df)
        }
        
        keep.events <- one_event_per_station(bio.data)
        
        bio.data$EVENT_ID <- paste(bio.data$EVENT_ID, bio.data$SAMPLE_NUMBER, sep = "_")
        bio.data <- bio.data[bio.data$EVENT_ID %in% keep.events$EVENT_ID, ]
        #bio.data <- bio.data[!bio.data$BIOREGION %in% c("MAC", "SEP"), ]
        
        #names(bio.data)[names(bio.data) %in% "ICPRB_BIOREGION_ID"] <- "BIOREGION"
        
        #table(test.class$BIOREGION, test.class$CATEGORY)
        env.class <- unique(test.class[, c("EVENT_ID", "STATION_ID", "DATE",
                                           "AGENCY_CODE", "SAMPLE_NUMBER", "CATEGORY")])
        env.class$EVENT_ID <- paste(env.class$EVENT_ID, env.class$SAMPLE_NUMBER, sep = "_")
        env.class <- env.class[env.class$EVENT_ID %in% bio.data$EVENT_ID, ]
        table(env.class$CATEGORY)
        
        env.merge <- merge(bio.data, env.class, by = c("EVENT_ID", "STATION_ID", "DATE",
                                                       "AGENCY_CODE", "SAMPLE_NUMBER"))
        env.merge <- unique(env.merge[, c("EVENT_ID", "STATION_ID", "DATE",
                                          "AGENCY_CODE", "SAMPLE_NUMBER",
                                          "CATEGORY", "BIOREGION")])
        
        env.merge$CATEGORY <- ifelse(env.merge$CATEGORY %in% "REF+", "REF", 
                                     ifelse(env.merge$CATEGORY %in% c("DEG+", "DEG"), "MOD", 
                                            ifelse(env.merge$CATEGORY %in% c("DEG2"), "SEV",
                                                   as.character(env.merge$CATEGORY))))
        class.counts <- data.frame(table(env.merge$BIOREGION, env.merge$CATEGORY))
        names(class.counts) <- c("BIOREGION", "CATEGORY", "COUNT")
        class.counts <- class.counts[class.counts$CATEGORY %in% c("REF", "SEV"), ]
        
        if (standard_count == TRUE) {
          random_sub <- function(my.df){
            
            my.df$EVENT_ID <- paste(my.df$EVENT_ID, my.df$SAMPLE_NUMBER, sep = "_")
            my.df <- unique(my.df[, c("EVENT_ID", "STATION_ID", "SAMPLE_NUMBER",
                                      "AGENCY_CODE", "DATE", "CATEGORY", "BIOREGION")])
            data.list <- list()
            for (i in unique(class.counts$BIOREGION)) {
              sub.class <- class.counts[class.counts$BIOREGION %in% i, ]
              bio.sub <- env.merge[env.merge$BIOREGION %in% i, ]
              bio.ref <- bio.sub[bio.sub$CATEGORY %in% "REF", ]
              bio.sev <- bio.sub[bio.sub$CATEGORY %in% "SEV", ]
              if (sub.class[sub.class$CATEGORY %in% "REF", "COUNT"] > 50 & index.res %in% "INLAND") {
                rand.ref.list <- sample(unique(bio.ref$EVENT_ID), 50)
                rand.ref <- bio.ref[bio.ref$EVENT_ID %in% rand.ref.list, ]
              } else {
                rand.ref <- bio.ref
              }
              
              if (sub.class[sub.class$CATEGORY %in% "SEV", "COUNT"] > 50) {
                rand.sev.list <- sample(unique(bio.sev$EVENT_ID), 50)
                rand.sev <- bio.sev[bio.sev$EVENT_ID %in% rand.sev.list, ]
              } else {
                rand.sev <- bio.sev
              }
              
              rand.df <- rbind(rand.ref, rand.sev)
              data.list[[i]] <- rand.df
              
            }
            
            final.df <- do.call(rbind, data.list)
            table(final.df$BIOREGION, final.df$CATEGORY)
            
            return(final.df)
          }
          
          
          rand.ref.sev <- random_sub(env.merge)
          table(rand.ref.sev$BIOREGION, rand.ref.sev$CATEGORY)
          bio.ref.sev <- bio.data[bio.data$EVENT_ID %in% rand.ref.sev$EVENT_ID, ]
          category.other <- env.merge[!env.merge$CATEGORY %in% c("REF", "SEV"), ]
          bio.other <- bio.data[bio.data$EVENT_ID %in% category.other$EVENT_ID, ]
          
          bio.data <- rbind(bio.ref.sev, bio.other)
          env.class <- env.class[env.class$EVENT_ID %in% bio.data$EVENT_ID, ]
          #env.class$EVENT_ID <- gsub("_.", "", env.class$EVENT_ID)
        }
        
        bio.data <- bio.data[!is.na(bio.data$BIOREGION), ]
      } 
      
      #==============================================================================
      # Calculate Standard Metrics
      system.time(
      bio.data.metrics <- BIBI::specific_metrics(master, bio.data,
                                                 taxa.rank =  taxon.rank,
                                                 rare = "DIV_RARE", pct_un = 10,
                                                 bibi.standard = TRUE,
                                                 metrics.vec = "ALL")
      )
      master2 <- clean_taxa(master)
      #==============================================================================
      # Calculate the percent of each taxon in the data base by Sequencing through each taxon rank
      if(taxon.rank %in% "ORDER") ranks.list <- c("SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
      if(taxon.rank %in% "FAMILY") ranks.list <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
      if(taxon.rank %in% "GENUS") ranks.list <- c("SPECIES")
      master2[, ranks.list] <- NA
      bio.data.seq <- seq_pct_taxa(bio.data, master2)
      #test <- bio.data.seq[duplicated(bio.data.metrics$EVENT_ID), ]
      #test2 <- bio.data.seq[bio.data.seq$EVENT_ID %in% test$EVENT_ID, ]
      #==============================================================================
      # Merge the Standard Metrics with the Sequenced Metrics
      bio.data.merge <- merge(bio.data.metrics, bio.data.seq, by = c("EVENT_ID", "STATION_ID",
                                                                     "DATE", "SAMPLE_NUMBER",
                                                                     "AGENCY_CODE"))
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
      #table(m.bio.data$CATEGORY)
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
      #table(m.bio.data$CATEGORY)
      #test <- data.frame(names(m.bio.data[, 7:ncol(m.bio.data)]))
      #==============================================================================
      # Re-assign the the bioregion to the newest data frame
      merged <- merge(unique(bio.data[, c("EVENT_ID", "BIOREGION")]), m.bio.data,
                      by = "EVENT_ID", all.y = TRUE)
      #write.csv(merged, "Test2_Metrics_11_7_16.csv", row.names = FALSE)
      #table(merged$CATEGORY)
      #==============================================================================
      # Change REF+ to REF and DEG+ to SEV
      merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF+", "REF", 
                                ifelse(merged$CATEGORY %in% c("DEG+", "DEG"), "MOD", 
                                       ifelse(merged$CATEGORY %in% c("DEG2"), "SEV",
                                              as.character(merged$CATEGORY))))
      table(merged$CATEGORY, merged$BIOREGION)
      #**********************************************************************************************
      # Changed on 4/12/17
      #**********************************************************************************************
      if(index.res %in% "COAST2"){
        #Set Working directory (The folder files are imported and exported to)
        #setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
        coastal_ref <- read.csv("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016/new_coastal_ref_11_3_16.csv")
        merged$EVENT_ID2 <- gsub( "_.*$", "", merged$EVENT_ID)
        coastal_ref$EVENT_ID[!coastal_ref$EVENT_ID %in% merged$EVENT_ID2]
        merged$CATEGORY <- ifelse(merged$CATEGORY %in% "REF" & merged$EVENT_ID2 %in% coastal_ref$EVENT_ID, "REF",
                                  ifelse(merged$CATEGORY %in% "REF" &
                                           !merged$EVENT_ID2 %in% coastal_ref$EVENT_ID, "MIN", as.character(merged$CATEGORY)))
        coastal_ref[!coastal_ref$EVENT_ID %in% merged$EVENT_ID2, ]
        merged <- merged[, !names(merged) %in% "EVENT_ID2"]
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
      merged$BIOREGION <- ifelse(merged$BIOREGION %in% c("NGV"), "UNP",
                                 as.character(merged$BIOREGION))
      table(merged$CATEGORY, merged$BIOREGION)
      #==============================================================================
      # Change REF+ to REF and DEG+ to SEV
      m.bio.data$CATEGORY <- ifelse(m.bio.data$CATEGORY %in% "REF+", "REF", 
                                    ifelse(m.bio.data$CATEGORY %in% c("DEG+", "DEG", "DEG2"),
                                           "SEV", m.bio.data$CATEGORY))
      table(m.bio.data$CATEGORY)
      #==============================================================================
      if(taxon.rank %in% "ORDER"){
        remove.order.calc <- c("RICH_BURROW", "RICH_CLIMB", "RICH_CLING",
                               "RICH_COLLECT", "RICH_FILTER", "RICH_GATHER",
                               "RICH_INTOL", "RICH_MODTOL", "RICH_PREDATOR",
                               "RICH_SCRAPE", "RICH_sHRED", "RICH_SPRAWL",
                               "RICH_sWIM", "RICH_TOL", "PCT_COLLECT",
                               "PCT_FILTER", "PCT_GATHER", "PCT_PREDATOR",
                               "PCT_SCRAPE", "PCT_SHRED", "PCT_BURROW",
                               "PCT_CLIMB", "PCT_CLING", "PCT_SPRAWL",
                               "PCT_SWIM", "ASPT_MOD", "BECKS_V1",
                               "BECKS_V2", "HBI", "PCT_ATI",
                               "PCT_INTOL_0_3", "PCT_INTOL_0_4",
                               "PCT_MOD_TOL_4_6", "PCT_TOLERANT_5_10",
                               "PCT_TOLERANT_7_10", "PCT_URBAN_INTOL")
        merged <- merged[, !names(merged) %in% remove.order.calc]
      }
      #============================================================================
      # Inputs that change with each method type
      #todays_date <- "10_14_16"
      redundant <- TRUE
      range.var <- TRUE
      zero.inflate <- FALSE
      fam.group <- FALSE
      metric.groups <- FALSE
      ibi.method <- "H"
      #jack.iterate = 2
      master.df <- master
      metric.class <- m.c
      #============================================================================
      #if(index.res %in% "BASIN"){
      #  raw.data <- merged
      #  raw.data$BIOREGION <- "ALL"
      #} 
      #if(index.res %in% "COAST"){
      #  raw.data <- merged[merged$BIOREGION %in% c("MAC", "SEP"), ]
      #  raw.data$BIOREGION <- "COAST"
      #}
      #if(index.res %in% "INLAND"){
      #  raw.data <- merged[!merged$BIOREGION %in% c("MAC", "SEP"), ]
      #  raw.data$BIOREGION <- "INLAND"
      #}
      #if(index.res %in% "INLAND"){
      #  raw.data <- merged
      #}
      raw.data <- merged
      if (index.res %in% "COAST") raw.data$BIOREGION <- "COAST"
      if (index.res %in% "INLAND") raw.data$BIOREGION <- "INLAND"
      if (index.res %in% "BASIN") raw.data$BIOREGION <- "BASIN"
      table(raw.data$CATEGORY, raw.data$BIOREGION)
      
      for(i in unique(raw.data$BIOREGION)){
        #============================================================================
        # Select only data from the bioregion of interest
        new.df <- raw.data[raw.data$BIOREGION %in% i, ]
        # Samples must contain more than 70 observed individuals.
        # Smaller sample sizes may result in odd percentage metric values.
        new.df <- new.df[new.df$ABUNDANCE > 70, ]
        
        remove.cols.1 <- c("EFFECTIVE_RICH_SHANNON", "EFFECTIVE_RICH_SIMPSON",
                           "PCT_UNIDENTIFIED", "NO_MATCH")
        rc <- names(new.df)[colSums(new.df[, 8:ncol(new.df)]) == 0]
        remove.cols.2 <- c(remove.cols.1, rc)
        new.df <- new.df[, !names(new.df) %in% remove.cols.2]
        
        if(ibi.method %in% c("B", "C", "D", "E")){
          if(ibi.method %in% c("C", "D", "E")){
            reg.metrics <- c("PCT_EPT_RICH_NO_TOL", "RICH", "RICH_CLING",
                             "RICH_INTOL", "SIMPSONS", "BECKS_V3", "HBI",
                             "PCT_DOM3", "PCT_INTOL_0_3", "PCT_COLLECT",
                             "PCT_PREDATOR", "PCT_SCRAPE", "PCT_SHRED",
                             "PCT_BURROW", "PCT_CLIMB", "PCT_CLING",
                             "PCT_SPRAWL", "PCT_SWIM", "PCT_NON_INSECT",
                             "PCT_ANNELID_CHIRO", "PCT_COTE", "PCT_DIPTERA",
                             "PCT_EPHEMEROPTERA_NO_BAETID", "PCT_EPT")
            if(ibi.method %in% c("D", "E")){
              fam.metrics <- na.omit(unique(master.df$FAMILY))
              fam.metrics <- paste("PCT_", fam.metrics, sep = "")
              keep.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY",
                                "STATION_ID", "SAMPLE_NUMBER",
                                "AGENCY_CODE", "DATE",
                                reg.metrics, fam.metrics)
            }else{
              keep.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY",
                                "STATION_ID", "SAMPLE_NUMBER",
                                "AGENCY_CODE", "DATE",
                                reg.metrics)
            }
            
          }else{
            standard.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY","STATION_ID", "SAMPLE_NUMBER",
                                  "AGENCY_CODE", "DATE","PCT_INTOL_0_3", "PCT_MOD_TOL_4_6", "PCT_TOLERANT_7_10",
                                  "PCT_GATHER", "PCT_FILTER", "PCT_PREDATOR", "PCT_SCRAPE", "PCT_SHRED",
                                  "PCT_BURROW", "PCT_CLIMB", "PCT_CLING", "PCT_SPRAWL", "PCT_SWIM")
            div_comp.metrics <- as.character(m.c[m.c$METRIC_CLASS %in% c("DIVERSITY", "COMPOSITION"), "METRICS"])
            keep.metrics <- c(standard.metrics, div_comp.metrics)
            
          }
          
          new.df <- new.df[, names(new.df) %in% keep.metrics]
        }
        #============================================================================
        # Metric Assessment
        #keep.cols <- new.df[, 8:ncol(new.df)]
        ms <- metrics_summary(new.df, i, zero = zero.inflate)
        system.time(
        os <- old_scoring3(new.df, i, zero_null = zero.inflate,
                           bound.lim = TRUE) #, metric.summary = ms)
        )
        remove.cols <- c("ABUNDANCE", "EFFECTIVE_RICH_SHANNON", "EFFECTIVE_RICH_SIMPSON",
                         "PCT_UNIDENTIFIED", "NO_MATCH")
        os <- os[, !grepl(paste0(remove.cols, collapse = "|"), names(os))]
        if(ibi.method %in% c("A", "B", "H")){
          if(ibi.method %in% "H"){
            #Remove metrics that cannot identify disturbance better than a random flip of a coin
            new.df <- new.df[, !names(new.df) %in% ms[ms$SENSITIVITY < 50, "METRICS"]]
            ms <- ms[ms$SENSITIVITY >= 50, ]
            os <- os[, !names(os) %in% ms[ms$SENSITIVITY < 50, "METRICS"]]
            #===================================================================================
            b.metrics <- select_best(metrics.df = new.df, ms, m.c = metric.class,
                                     range_var = range.var,
                                     redund = redundant)
            bd <- best.distribution.H(new.df, ms, os, m.c = metric.class,
                                      range_var = range.var,
                                      redund = redundant,
                                      best.metrics = b.metrics)
            
            
            metrics.vec2 <- bd$METRICS
            
          }
          if(ibi.method %in% "A"){
            bd <- best.distribution2(metrics.df = new.df,
                                     metric.summary = ms,
                                     scores.df = os,
                                     m.c = metric.class,
                                     Family = fam.group,
                                     master = master.df,
                                     range_var = range.var)
            metrics.vec2 <- bd$METRICS
            
          }
          if(ibi.method %in% "B"){
            
            bd <- best.distribution3(new.df, ms, os, m.c= metric.class,
                                     Family = fam.group, 
                                     master <- master.df,
                                     range_var = range.var)
            metrics.vec2 <- bd$METRICS
          }
          
          
          #write.csv(bd, paste(ibi.method, "_", i, "_ms_", todays_date, ".csv", sep = ""))
          #============================================================================
          # Creat a vectorof  metrics with sensitivity values (DE) >= the best 
          # threshold provided by the best.distribution function.
          #metrics.vec <- ms[ms$SENSITIVITY >= bd[1, 7], "METRICS"]
          
          
          
          final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                           "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                           as.character(metrics.vec2)))
          final.metrics.vec2 <- final.metrics.vec[!final.metrics.vec %in% remove.cols]
          if(redundant == TRUE & length(final.metrics.vec) > 7){
            if(ibi.method %in% c("A", "B")){
              R_test <- metric_class_redund(new.df[, names(new.df) %in% final.metrics.vec],
                                            ms[ms$METRICS %in% final.metrics.vec,],
                                            m.c, method = ibi.method)
              final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                               "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                               as.character(R_test$METRICS)))
              #============================================================================
              # Select only the metrics selected for the final index from the os (contains scored values)
              scored <- os[, final.metrics.vec]
              scored[, 7:ncol(scored)] <- apply(scored[, 7:ncol(scored)], 2, function(x) as.numeric(as.character(x)))
              scored$FINAL_SCORE <- apply(scored[, 7:ncol(scored)], 1, mean, na.rm = TRUE)
              scored <- scoring2(metrics.df = new.df[, names(new.df) %in% final.metrics.vec],
                                 scores = os[, names(os) %in% final.metrics.vec],
                                 ms,
                                 m.c = metric.class,
                                 Family = fam.group,
                                 master = master.df,
                                 redund = redundant, 
                                 method = ibi.method)
              
            }else{
              R_test <- redundancy(new.df[, names(new.df) %in% final.metrics.vec],
                                   analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% final.metrics.vec, ],
                                   lower.class = "SEV", upper.class = "REF",
                                   lower.coef = -0.85, upper.coef = 0.85)
              final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                               "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                               as.character(R_test$METRICS)))
              #============================================================================
              # Select only the metrics selected for the final index from the os (contains scored values)
              scored <- os[, final.metrics.vec]
              scored[, 7] <- as.numeric(scored[, 7])
              scored$FINAL_SCORE <- scored[, 7:ncol(scored)]
              scored <- scoring2(metrics.df = new.df[, names(new.df) %in% final.metrics.vec],
                                 scores = os[, names(os) %in% final.metrics.vec],
                                 ms,
                                 m.c = metric.class,
                                 Family = fam.group,
                                 master = master.df,
                                 redund = redundant, 
                                 method = ibi.method)
            }
            bde.df <- bde(scored, i)
            data.list2[[j]] <- bde.df[1, ]
          }
        }
        final.df <- ms[ms$METRICS %in% final.metrics.vec[7:length((final.metrics.vec))], ]
        final.df$ITERATION <- j
        final.df <- final.df[, c(ncol(final.df), 1:(ncol(final.df) - 1))]
      }
      final.df$BIOREGION <- i
      final.df <- final.df[, c(1,ncol(final.df), 2:(ncol(final.df) - 1))]
      data.list[[j]] <- final.df
      print(paste("End Iteration:", j))
    }
    #========================================================================
    bound <- do.call(rbind, data.list)
    #========================================================================
    write.csv(bound, paste(as.character(unique(bound$BIOREGION)),
                           taxon.rank, paste(runs, "i", sep = ""), "Metric_Summary",
                           paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
    #========================================================================
    mean_ceiling <- aggregate(REF_MEDIAN ~ METRICS, data = bound, FUN = "mean")
    sd_ceiling <- aggregate(REF_MEDIAN ~ METRICS, data = bound, FUN = sd)
    min_ceiling <- aggregate(REF_MEDIAN ~ METRICS, data = bound, FUN = min)
    max_ceiling <- aggregate(REF_MEDIAN ~ METRICS, data = bound, FUN = max)
    mean_ceiling <- plyr::join_all(list(mean_ceiling, sd_ceiling, min_ceiling, max_ceiling), by = "METRICS")
    names(mean_ceiling) <- c("METRICS", "MEAN_CEILING", "SD_CEILING", "MIN_CEILING", "MAX_CEILING")
    
    mean_floor <- aggregate(BOUND_BI_CMA ~ METRICS, data = bound, FUN = mean)
    sd_floor <- aggregate(BOUND_BI_CMA ~ METRICS, data = bound, FUN = sd)
    min_floor <- aggregate(BOUND_BI_CMA ~ METRICS, data = bound, FUN = min)
    max_floor <- aggregate(BOUND_BI_CMA ~ METRICS, data = bound, FUN = max)
    mean_floor <- plyr::join_all(list(mean_floor, sd_floor, min_floor, max_floor), by = "METRICS")
    names(mean_floor) <- c("METRICS", "MEAN_FLOOR", "SD_FLOOR", "MIN_FLOOR", "MAX_FLOOR")
    
    mean_thresholds <- merge(mean_ceiling, mean_floor, by = "METRICS")
    
    write.csv(mean_thresholds, paste(as.character(unique(bound$BIOREGION)), taxon.rank, paste(runs, "i", sep = ""), "Mean_Thresholds",
                                     paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
    
    #========================================================================
    final.df <- data.frame(table(bound$METRICS))
    final.df <- final.df[order(-final.df$Freq), ]
    names(final.df) <- c("METRICS", "FREQUENCY")
    final.df <- final.df[final.df$FREQUENCY > 0, ]
    final.df$PERCENT <- (final.df$FREQUENCY / runs) * 100
    write.csv(final.df, paste(as.character(unique(bound$BIOREGION)), taxon.rank, paste(runs, "i", sep = ""), "Metric_Selection",
                              paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
    #========================================================================
    freq_thresh <- merge(final.df[, c("METRICS", "PERCENT")], mean_thresholds, by = "METRICS")
    write.csv(freq_thresh, paste(as.character(unique(bound$BIOREGION)), taxon.rank, paste(runs, "i", sep = ""),
                                 "Freq_Thresholds_Summary",
                                 paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
    #========================================================================
    
    ce.df <- do.call(rbind, data.list2)
    write.csv(ce.df, paste(as.character(unique(bound$BIOREGION)), taxon.rank, paste(runs, "i", sep = ""),
                           "Raw_CE",
                           paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
    ce.summary <- data.frame(Index = taxon.rank)
    ce.summary$Bioregion <- unique(ce.df$BIOREGION)
    ce.summary$Mean <- mean(ce.df$CE)
    ce.summary$'Standard Deviation' <- sd(ce.df$CE)
    ce.summary$Minimum <- min(ce.df$CE)
    ce.summary$Maximum <- max(ce.df$CE)
    write.csv(ce.summary, paste(as.character(unique(bound$BIOREGION)), taxon.rank, paste(runs, "i", sep = ""),
                                "CE_Summary",
                                paste(todays_date, ".csv", sep = ""), sep = "_"), row.names = F)
  }
  
  if(length(unique(all.data$BIOREGION)) > 1){
    for(b in unique(all.data$BIOREGION)){
      print("################################################################")
      print(paste("Start Bioregion:", b))
      print("################################################################")
      bio.data <- ref.sev.df[ref.sev.df$BIOREGION %in% b, ]
      #bio.data <- all.data[all.data$BIOREGION %in% b, ]
      #bio.data <- bio.data[, !names(bio.data) %in% "BIOREGION"]
      iterate_this(bio.data, index.res, taxon.rank, runs, todays_date, env.class,
                   unique.station, standard_count)
      print("################################################################")
      print(paste("End Bioregion:", b))
      print("################################################################")
    }
  }else{
    bio.data <- ref.sev.df
    #bio.data <- bio.data[, !names(bio.data) %in% "BIOREGION"]
    iterate_this(bio.data, index.res, taxon.rank, runs, todays_date, env.class,
                 unique.station, standard_count)
  }
}