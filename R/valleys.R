#==============================================================================
#VALLEYS
#==============================================================================
#'Valleys: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

val_clean_taxa <- function(Long){

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
  
  
  phy_spp <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
               "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$PHYLUM %in% c("NEMATODA", "NEMERTEA"), phy_spp] <- "UNIDENTIFIED"  
  
  subclass_spp <- c("SUBCLASS", "ORDER", "SUBORDER", "FAMILY",
                    "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$SUBPHYLUM %in% c("CHELICERATA"), subclass_spp] <- Long[Long$SUBPHYLUM %in% c("CHELICERATA"), "CLASS"]
  #Long[Long$CLASS %in% c("OLIGOCHAETA","HIRUDINEA"), subclass_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "OLIGOCHAETA", subclass_spp] <- Long[Long$CLASS %in% "OLIGOCHAETA", "CLASS"]
  Long[Long$CLASS %in% "HIRUDINEA", subclass_spp] <- Long[Long$CLASS %in% "HIRUDINEA", "CLASS"]
  
  fam_spp <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$FAMILY %in% c("CHIRONOMIDAE"), fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}


#==============================================================================
#'Valley's Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Valley's bioregion with a
#'Strahler Stream order <= 4
#'@export
prep_valleys <- function(Long) {
  sub.valleys <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                          c("67e", "67E", "67a", "67A", "67b", "67B", "67f",
                            "67F","67g", "67G"))
  sub.methods <- prep_subset(sub.valleys)
  agg.valleys <- bioregion_agg(sub.methods)
  final.df <- agg.valleys
  #final.df <- val_clean_taxa(agg.valleys)
  return(final.df)
}

#==============================================================================
#'Valleys-All metrics
#'
#'@param Info = Taxonomic Information
#'@param Family = Familial taxonomic counts
#'@param Order = Ordinal taxonomic counts
#'@param Long_Data = Taxonomic data in long format
#'@param Level = Taxonomic level ("FAMILY" or "GENUS")
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the all of the valley ecoregion
#'@export
metrics_valleys <- function(master, Long, ffg_col) {
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
  
  rare.long.fill <- BIBI::prep_rare(long.fill, master, "FAMILY", 100, NULL, FALSE)
  rare.long.fill[, taxa] <- fill_taxa(rare.long.fill[, taxa])
  rare.fam.fill <- BIBI::wide(rare.long.fill, "FAMILY")
  #============================================================================
  Family <- BIBI::wide(Long, "FAMILY")
  names(Family) <- toupper(colnames(Family))
  Order <- BIBI::wide(Long, "ORDER")


  metrics <- data.frame(Family[, 1:5])
  metrics$BECKS_100 <- BIBI::becks(rare.fam.fill,  "FAMILY", master, beck.version = 1)
  metrics$PCT_EPHEMEROPTERA <- BIBI::pct_ephemeroptera(Order)
  metrics$PCT_EPT_TAXA_RICH <- BIBI::pct_ept_rich(Long, "FAMILY")
  metrics$PCT_SCRAPER <- BIBI::pct_attribute(fam.fill, master, ffg_col, "SC", "FAMILY")
  metrics$SW <- BIBI::shannon(Family)
  return(metrics)
}
#==============================================================================
#'Score valley's Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Metric scores for sites within the Valley's bioregion
#'@export

score_valleys <- function(Info, Long, ffg_col, scoring) {

  metrics <- metrics_valleys(Info, Long, ffg_col)

  valleys.thresh <- matrix(c(6.9, 18.97, 50, 7.02, 1.88, 9, 29.33, 55.56, 14.41, 2.09),
                           nrow = 5, ncol = 2)
  rownames(valleys.thresh) <- c("Becks Index 100", "%EPHEMEROPTERA", "%EPT TAXA RICH",
                                "%SCRAPERS", "SW")
  colnames(valleys.thresh) <- c("XT", "XM")
  thresh <- data.frame(valleys.thresh)

  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])
  
  if(scoring %in% "DISCRETE"){
    score$BECKS_100 <- score_1_3_5(metrics, thresh[1, ], "BECKS_100", "DECREASE")
    score$PCT_EPHEMEROPTERA <- score_1_3_5(metrics, thresh[2, ], "PCT_EPHEMEROPTERA", "DECREASE")
    score$PCT_EPT_TAXA_RICH <- score_1_3_5(metrics, thresh[3, ], "PCT_EPT_TAXA_RICH", "DECREASE")
    score$PCT_SCRAPER  <- score_1_3_5(metrics, thresh[4, ], "PCT_SCRAPER", "DECREASE")
    score$SW <- score_1_3_5(metrics, thresh[5, ], "SW", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      score$BECKS_100 <- ifelse (metrics$BECKS_100 > thresh$XT[1] &
                                   metrics$BECKS_100 < thresh$XM[1],
                                 ((metrics$BECKS_100 - thresh$XT[1]) /
                                    (thresh$XM[1] - thresh$XT[1])) * 100,
                                 ifelse (metrics$BECKS_100 <= thresh$XT[1], 0,
                                         ifelse (metrics$BECKS_100 >= thresh$XM[1], 100,
                                                 "ERROR")))
      
      
      score$PCT_EPHEMEROPTERA <- ifelse (metrics$PCT_EPHEMEROPTERA > thresh$XT[2] &
                                           metrics$PCT_EPHEMEROPTERA < thresh$XM[2],
                                         ((metrics$PCT_EPHEMEROPTERA - thresh$XT[2]) /
                                            (thresh$XM[2] - thresh$XT[2])) * 100,
                                         ifelse (metrics$PCT_EPHEMEROPTERA <= thresh$XT[2], 0,
                                                 ifelse (metrics$PCT_EPHEMEROPTERA >= thresh$XM[2],
                                                         100, "ERROR")))
      
      score$PCT_EPT_TAXA_RICH <- ifelse (metrics$PCT_EPT_TAXA_RICH > thresh$XT[3] &
                                           metrics$PCT_EPT_TAXA_RICH < thresh$XM[3],
                                         ((metrics$PCT_EPT_TAXA_RICH - thresh$XT[3]) /
                                            (thresh$XM[3] - thresh$XT[3])) * 100,
                                         ifelse (metrics$PCT_EPT_TAXA_RICH <= thresh$XT[3], 0,
                                                 ifelse (metrics$PCT_EPT_TAXA_RICH >= thresh$XM[3],
                                                         100, "ERROR")))
      
      
      score$PCT_SCRAPER <- ifelse (metrics$PCT_SCRAPER >thresh$XT[4] &
                                     metrics$PCT_SCRAPER < thresh$XM[4],
                                   ((metrics$PCT_SCRAPER - thresh$XT[4]) /
                                      (thresh$XM[4] - thresh$XT[4])) * 100,
                                   ifelse (metrics$PCT_SCRAPER <= thresh$XT[4], 0,
                                           ifelse (metrics$PCT_SCRAPER >= thresh$XM[4], 100,
                                                   "ERROR")))
      
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
#'Valleys
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the valleys bioregion
#'@export
valleys <- function(Info, Long, ffg_col, scoring) {
  prep <- prep_valleys(Long)
  v.score <- score_valleys(Info, prep, ffg_col, scoring)
  return(v.score)
}

