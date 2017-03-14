#==============================================================================
# North Central Appalachians
#==============================================================================
#'NCA: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

nca_clean_taxa <- function(Long){
  Long <- Long[!Long$TSN %in% Chessie::exclusions_2011$TSN, ]
  phy <- wide(Long, "PHYLUM")
  
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
  Long[Long$PHYLUM %in% c("NEMERTEA"), phy_spp] <- "UNIDENTIFIED"  
  Long[Long$SUBPHYLUM %in% c("CHELICERATA"), phy_spp] <- "UNIDENTIFIED"
  
  
  subclass_spp <- c("SUBCLASS", "ORDER", "SUBORDER", "FAMILY",
                    "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  #Long[Long$CLASS %in% c("OLIGOCHAETA", "HIRUDINEA"), subclass_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "OLIGOCHAETA", subclass_spp] <- Long[Long$CLASS %in% "OLIGOCHAETA", "CLASS"]
  Long[Long$CLASS %in% "HIRUDINEA", subclass_spp] <- Long[Long$CLASS %in% "HIRUDINEA", "CLASS"]
  
  fam_spp <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$CLASS %in% c("GASTROPODA"), fam_spp] <- Long[Long$CLASS %in% c("GASTROPODA"), "FAMILY"]
  
  Long[Long$FAMILY %in% c("CHIRONOMIDAE"), fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}

#==============================================================================
#'NCA Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the NCA bioregion with a
#'Strahler Stream order <= 4.
#'@export
prep_nca <- function(Long) {
  sub.nca <- subset(Long, Long$ECOREGION_LEVEL_4 %in% c("62a", "62A", "62b",
                                                        "62B", "62c", "62C",
                                                        "62d", "62D"))
  sub.method <- prep_subset(sub.nca)
  agg.nca <- bioregion_agg(sub.method)
  final.df <- agg.nca
  #final.df <- nca_clean_taxa(agg.nca)
  return(final.df)
}

#==============================================================================
#'NCA metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the north central Appalachians (NCA) ecoregion
#'@export
metrics_nca <- function(master, Long, tol_col, ffg_col) {
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
  Family <- BIBI::wide(Long, "FAMILY")
  Order <- BIBI::wide(Long, "ORDER")
  #Rarefied <- BIBI::prarefy(Long)
  rare.long <- BIBI::prep_rare(Long, master, "FAMILY", 100, NULL, FALSE)
  rare.long[, taxa] <- fill_taxa(rare.long[, taxa])
  Rarefied <- BIBI::wide(rare.long, "FAMILY")

  metrics <- data.frame(Family[, 1:5])
  metrics$EPT_RICH_NO_TOL <- BIBI::ept_rich_no_tol(rare.long, "FAMILY", master, tolerance_value = tol_col)
  metrics$PCT_SCRAPER <- BIBI::pct_attribute(fam.fill, master, ffg_col, "SC", "FAMILY")
  metrics$TAXA_RICH_100 <- vegan::specnumber(Rarefied[, 6:ncol(Rarefied)])
  metrics$PCT_EPHEMEROPTERA <- BIBI::pct_ephemeroptera(Order)
  metrics$SW <- BIBI::shannon(Family)
  return(metrics)
}

#==============================================================================
#'Score NCA Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Metric scores for sites within the NCA bioregion
#'@export
score_nca <- function(Info, Long, tol_col, ffg_col, scoring) {

  metrics <- metrics_nca(Info, Long, tol_col, ffg_col)

  nca.thresh <- data.frame(XT = c(9, 8.45, 15, 24.26, 2.19),
                           XM = c(11, 13.3, 17, 40.22, 2.38))
  rownames(nca.thresh) <- c("EPT_TAXA_COUNT_NO_TOL_100", "%SCRAPER",
                            "%TAXA_RICH_100", "%EPHEM", "Shannon")
  thresh <- data.frame(nca.thresh)

  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])
  if(scoring %in% "DISCRETE"){
    score$EPT_RICH_NO_TOL<- score_1_3_5(metrics, thresh[1, ], "EPT_RICH_NO_TOL", "DECREASE")
    score$PCT_SCRAPER <- score_1_3_5(metrics, thresh[2, ], "PCT_SCRAPER", "DECREASE")
    score$TAXA_RICH_100 <- score_1_3_5(metrics, thresh[3, ], "TAXA_RICH_100", "DECREASE")
    score$PCT_EPHEMEROPTERA   <- score_1_3_5(metrics, thresh[4, ], "PCT_EPHEMEROPTERA", "DECREASE")
    score$SW <- score_1_3_5(metrics, thresh[5, ], "SW", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      score$EPT_RICH_NO_TOL <- ifelse (metrics$EPT_RICH_NO_TOL > thresh$XT[1] &
                                         metrics$EPT_RICH_NO_TOL < thresh$XM[1],
                                       ((metrics$EPT_RICH_NO_TOL - thresh$XT[1]) /
                                          (thresh$XM[1] - thresh$XT[1])) * 100,
                                       ifelse (metrics$EPT_RICH_NO_TOL <= thresh$XT[1], 0,
                                               ifelse (metrics$EPT_RICH_NO_TOL >= thresh$XM[1],
                                                       100, "ERROR")))
      
      score$PCT_SCRAPER <- ifelse (metrics$PCT_SCRAPER > thresh$XT[2] &
                                     metrics$PCT_SCRAPER < thresh$XM[2],
                                   ((metrics$PCT_SCRAPER - thresh$XT[2]) /
                                      (thresh$XM[2] - thresh$XT[2])) * 100,
                                   ifelse (metrics$PCT_SCRAPER <= thresh$XT[2], 0,
                                           ifelse (metrics$PCT_SCRAPER >= thresh$XM[2], 100,
                                                   "ERROR")))
      
      score$TAXA_RICH_100 <- ifelse (metrics$TAXA_RICH_100 > thresh$XT[3] &
                                       metrics$TAXA_RICH_100 < thresh$XM[3],
                                     ((metrics$TAXA_RICH_100 - thresh$XT[3]) /
                                        (thresh$XM[3] - thresh$XT[3])) * 100,
                                     ifelse (metrics$TAXA_RICH_100 <= thresh$XT[3], 0,
                                             ifelse (metrics$TAXA_RICH_100 >= thresh$XM[3], 100,
                                                     "ERROR")))
      
      score$PCT_EPHEMEROPTERA <- ifelse (metrics$PCT_EPHEMEROPTERA > thresh$XT[4] &
                                           metrics$PCT_EPHEMEROPTERA <
                                           thresh$XM[4],
                                         ((metrics$PCT_EPHEMEROPTERA -
                                             thresh$XT[4]) /
                                            (thresh$XM[4] - thresh$XT[4])) * 100,
                                         ifelse (metrics$PCT_EPHEMEROPTERA <=
                                                   thresh$XT[4], 0,
                                                 ifelse (metrics$PCT_EPHEMEROPTERA >=
                                                           thresh$XM[4],
                                                         100, "ERROR")))
      
      score$SW <- ifelse (metrics$SW > thresh$XT[5] & metrics$SW < thresh$XM[5],
                          ((metrics$SW - thresh$XT[5]) /
                             (thresh$XM[5] - thresh$XT[5])) * 100,
                          ifelse (metrics$SW <= thresh$XT[5], 0,
                                  ifelse (metrics$SW >= thresh$XM[5], 100, "ERROR")))
    }
  }

 

  return(prep_score(score, metrics))
}

#==============================================================================
#'NCA
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the NCA bioregion
#'@export
nca <- function(Info, Long, tol_col, ffg_col, scoring) {
  prep <- prep_nca(Long)
  n.score <- score_nca(Info, prep, tol_col, ffg_col, scoring)
  return(n.score)
}


