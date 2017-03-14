#==============================================================================
#RIDGES
#==============================================================================
#'Ridge: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

ridge_clean_taxa <- function(Long){
  Long <- Long[!Long$TSN %in% Chessie::exclusions_2011$TSN, ]
  phy <- wide(Long, "PHYLUM")
  
  nematoda <- if("NEMATODA" %in% names(phy)){
    (phy$NEMATODA / rowSums(phy[,6:ncol(phy)])) * 100
  }else{
    0
  }
  cat("Nematoda  were excluded from the data.\
      The maximum percentage of Nematoda  found in a sample was",
      round(max(nematoda ), 2),"%.")
  
  nemertea <- if("NEMERTEA" %in% names(phy)){
    (phy$NEMERTEA / rowSums(phy[,6:ncol(phy)])) * 100
  }else{
    0
  }
  cat("Nemertea  were excluded from the data.\
      The maximum percentage of Nemertea  found in a sample was",
      round(max(nemertea ), 2),"%.")
  
  
  subphy <- wide(Long, "SUBPHYLUM")
  mites <- if("CHELICERATA" %in% names(subphy)){
    (subphy$CHELICERATA / rowSums(subphy[,6:ncol(subphy)])) * 100
  }else{
    0
  }
  
  cat("Mites were excluded from the data.\
      The maximum percentage of Mites found in a sample was", round(max(mites), 2),"%.")
  
  phy_spp <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
               "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$PHYLUM %in% c("NEMATODA", "NEMERTEA"), phy_spp] <- "UNIDENTIFIED"  
  Long[Long$SUBPHYLUM %in% c("CHELICERATA"), phy_spp] <- "UNIDENTIFIED"
  
  fam_spp <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  #Long[Long$CLASS %in% c("OLIGOCHAETA", "HIRUDINEA",
  #                       "GASTROPODA"), fam_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "OLIGOCHAETA", fam_spp] <- Long[Long$CLASS %in% "OLIGOCHAETA", "FAMILY"]
  Long[Long$CLASS %in% "HIRUDINEA", fam_spp] <- Long[Long$CLASS %in% "HIRUDINEA", "FAMILY"]
  Long[Long$CLASS %in% "GASTROPODA", fam_spp] <- Long[Long$CLASS %in% "GASTROPODA", "FAMILY"]
  Long[Long$FAMILY %in% c("CHIRONOMIDAE"), fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}

#==============================================================================
#'Ridges Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Ridges bioregion with a
#'Strahler Stream order <= 4
#'@export
prep_ridges <- function(Long) {
  sub.ridges <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                         c("66a", "66A", "66b", "66B", "67d", "67D", "67c", "67C",
                           "67i", "67I", "67h", "67H", "69a", "69A", "69b", "69B",
                           "70c", "70C"))
  sub.method <- prep_subset(sub.ridges)
  agg.ridges<- bioregion_agg(sub.method)
  final.df <- agg.ridges
  #final.df <- ridge_clean_taxa(agg.ridges)
  return(agg.ridges)
}

#==============================================================================
#'Ridges metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param habit_col = Specify Habit column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the ridges ecoregion
#'@export
metrics_ridges <- function(master, Long, ffg_col, habit_col) {
  Long <- Long[!Long$TSN %in% Chessie::exclusions_2011$TSN, ]
  #============================================================================
  # These should be used for taxa attribute related metrics
  
  
  taxa <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")
  
  #master2 <- fill_taxa(master)
  #master.fill <- unique(master2[, c("TSN_R", taxa)])
  
  #test <- (master.fill[duplicated(master.fill$TSN_R), ])
  long.fill <- Long
  #long.fill <- clean_taxa(long.fill)
  long.fill[long.fill == "UNIDENTIFIED"] <- NA
  long.fill <- long.fill[!is.na(long.fill$PHYLUM), ]
  
  long.sub <- fill_taxa(unique(long.fill[, c("TSN", taxa)]))
  long.fill <- long.fill[, !names(long.fill) %in% taxa]
  
  long.fill <- merge(long.fill, long.sub, by = "TSN", all.x = T)
  fam.fill <- wide(long.fill, "FAMILY")
  #============================================================================
  rare.long.fill <- BIBI::prep_rare(long.fill, master, "FAMILY", 100, NULL, FALSE)
  rare.long.fill[, taxa] <- fill_taxa(rare.long.fill[, taxa])
  rare.fam.fill <- BIBI::wide(rare.long.fill, "FAMILY")
  #============================================================================
  Family <- BIBI::wide(Long, "FAMILY")
  names(Family) <- toupper(colnames(Family))
  Order <- BIBI::wide(Long, "ORDER")
  #Rarefied <- BIBI::prarefy(Long)
  rare.long <- BIBI::prep_rare(Long, master, "FAMILY", 100, NULL, FALSE)
  rare.long[, taxa] <- fill_taxa(rare.long[, taxa])
  Rarefied <- BIBI::wide(rare.long, "FAMILY")

  metrics <- data.frame(Family[, 1:5])
  metrics$BECKS_100 <- BIBI::becks(rare.fam.fill, "FAMILY", master, beck.version = 1)
  metrics$PCT_EPHEMEROPTERA <- BIBI::pct_ephemeroptera(Order)
  metrics$PCT_SCRAPER <- BIBI::pct_attribute(fam.fill, master, ffg_col, "SC", "FAMILY")
  metrics$PCT_SWIMMER <- BIBI::pct_attribute(fam.fill, master, habit_col, "SW", "FAMILY")
  metrics$SW <- BIBI::shannon(Family)
  return(metrics)
}
#==============================================================================
#'Score Ridges Metrics
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@param habit_col = Specify Habit column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Metric scores for sites within the Ridges bioregion
#'@export
score_ridges <- function(Info, Long, ffg_col, habit_col, scoring) {

  metrics <- metrics_ridges(Info, Long, ffg_col, habit_col)

  ridges.thresh <- matrix(c(10, 13.51, 3.37, 4.59, 1.99, 13, 26.83, 11.06, 10.73,
                            2.26), nrow = 5, ncol = 2)
  rownames(ridges.thresh) <- c("Becks Index 100", "%EPHEMEROPTERA", "%SCRAPERS",
                               "%SWIMMERS", "SW")
  colnames(ridges.thresh) <- c("XT", "XM")
  thresh <- data.frame(ridges.thresh)

  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])
  if(scoring %in% "DISCRETE"){
    score$BECKS_100  <- score_1_3_5(metrics, thresh[1, ], "BECKS_100", "DECREASE")
    score$PCT_EPHEMEROPTERA <- score_1_3_5(metrics, thresh[2, ], "PCT_EPHEMEROPTERA", "DECREASE")
    score$PCT_SCRAPER <- score_1_3_5(metrics, thresh[3, ], "PCT_SCRAPER", "DECREASE")
    score$PCT_SWIMMER <- score_1_3_5(metrics, thresh[4, ], "PCT_SWIMMER", "DECREASE")
    score$SW <- score_1_3_5(metrics, thresh[5, ], "SW", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      score$BECKS_100 <- ifelse (metrics$BECKS_100 > thresh$XT[1] &
                                   metrics$BECKS_100 < thresh$XM[1],
                                 ((metrics$BECKS_100 - thresh$XT[1]) /
                                    (thresh$XM[1] - thresh$XT[1])) * 100,
                                 ifelse (metrics$BECKS_100 <= thresh$XT[1], 0,
                                         ifelse (metrics$BECKS_100 >= thresh$XM[1], 100, "ERROR")))
      
      score$PCT_EPHEMEROPTERA <- ifelse (metrics$PCT_EPHEMEROPTERA > thresh$XT[2] &
                                           metrics$PCT_EPHEMEROPTERA < thresh$XM[2],
                                         ((metrics$PCT_EPHEMEROPTERA - thresh$XT[2]) /
                                            (thresh$XM[2] - thresh$XT[2])) * 100,
                                         ifelse (metrics$PCT_EPHEMEROPTERA <= thresh$XT[2],  0,
                                                 ifelse (metrics$PCT_EPHEMEROPTERA >= thresh$XM[2],
                                                         100, "ERROR")))
      
      score$PCT_SCRAPER <- ifelse (metrics$PCT_SCRAPER > thresh$XT[3] &
                                     metrics$PCT_SCRAPER < thresh$XM[3],
                                   ((metrics$PCT_SCRAPER-thresh$XT[3]) /
                                      (thresh$XM[3] - thresh$XT[3])) * 100,
                                   ifelse (metrics$PCT_SCRAPER <= thresh$XT[3], 0,
                                           ifelse (metrics$PCT_SCRAPER >= thresh$XM[3], 100,
                                                   "ERROR")))
      
      score$PCT_SWIMMER <- ifelse (metrics$PCT_SWIMMER > thresh$XT[4] &
                                     metrics$PCT_SWIMMER < thresh$XM[4],
                                   ((metrics$PCT_SWIMMER  -thresh$XT[4]) /
                                      (thresh$XM[4] - thresh$XT[4])) * 100,
                                   ifelse (metrics$PCT_SWIMMER <= thresh$XT[4], 0,
                                           ifelse (metrics$PCT_SWIMMER >= thresh$XM[4], 100,
                                                   "ERROR")))
      
      
      score$SW <- ifelse(metrics$SW > thresh$XT[5] & metrics$SW < thresh$XM[5],
                         ((metrics$SW - thresh$XT[5]) /
                            (thresh$XM[5] - thresh$XT[5])) * 100,
                         ifelse(metrics$SW <= thresh$XT[5], 0,
                                ifelse(metrics$SW >= thresh$XM[5], 100, "ERROR")))
      
    }
  }

  
  return(prep_score(score, metrics))
}

#==============================================================================
#'Ridges
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param habit_col = Specify Habit column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the Ridges bioregion
#'@export
ridges <- function(Info, Long, ffg_col, habit_col, scoring) {
  prep <- prep_ridges(Long)
  r.score <- score_ridges(Info, prep, ffg_col, habit_col, scoring)
  return(r.score)
}

