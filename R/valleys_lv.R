#==============================================================================
#VALLEYS LIMESTONE PRESENT
#==============================================================================
#'Valleys_LV PREP
#'
#'@param Long= Taxonomic data in long format
#'@return 2011 Chessie BIBI metrics and scores for the ridge ecoregion
#'@export
prep_valleys_lv <- function(Long) {
  sub.valleys <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                          c("67e", "67E", "67a", "67A", "67b", "67B", "67f",
                            "67F", "67g", "67G") &
                          !Long$ICPRB_KARST_ID %in% NA)
  sub.methods <- prep_subset(sub.valleys)
  agg.valleys <- bioregion_agg(sub.methods)
  return(agg.valleys)
}

#==============================================================================
#'Valleys_LV metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param aspt_col = Specify ASPT column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the all of the valley ecoregion
#'@export
metrics_valleys_lv <- function(master, Long, aspt_col, ffg_col) {
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
  Family <- wide(Long, "FAMILY")
  names(Family) <- toupper(colnames(Family))
  Order <- wide(Long, "ORDER")

  metrics <- data.frame(Family[, 1:5])
  metrics$ASPT_MOD <- tol_index(long.fill, master, aspt_col, "FAMILY")
  metrics$PCT_EPHEMEROPTERA <- pct_ephemeroptera(Order)
  metrics$PCT_EPT_TAXA_RICH <- pct_ept_rich(Long, "FAMILY")
  metrics$PCT_SCRAPER <- BIBI::pct_attribute(fam.fill, master, ffg_col, "SC", "FAMILY")
  metrics$SW <- shannon(Family)
  return(metrics)
}
#==============================================================================
#'Valleys_LV
#'
#'@param Info = Taxonomic Information
#'@param Long= Taxonomic data in long format
#'@param aspt_col = Specify ASPT column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return 2011 Chessie BIBI metrics and scores for the valleys ecoregion
#'@export
score_valleys_lv <- function(Info, Long, aspt_col, ffg_col, scoring){

  metrics <- metrics_valleys_lv(Info, Long, aspt_col, ffg_col)

  thresh <- matrix(c(4.45, 20.83, 50, 11.54, 1.81, 4.2, 38.51, 57.14, 17.54, 2.14), nrow=5, ncol=2)
  rownames(thresh) <- c("ASPT_MOD", "%EPHEMEROPTERA", "%EPT TAXA RICH",
                                "%SCRAPERS", "SW")
  colnames(thresh) <- c("XT", "XM")
  thresh <- data.frame(thresh)

  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])
  
  if(scoring %in% "DISCRETE"){
    score$ASPT_MOD <- score_1_3_5(metrics, thresh[1, ], "ASPT_MOD", "INCREASE")
    score$PCT_EPHEMEROPTERA <- score_1_3_5(metrics, thresh[2, ], "PCT_EPHEMEROPTERA", "DECREASE")
    score$PCT_EPT_TAXA_RICH <- score_1_3_5(metrics, thresh[3, ], "PCT_EPT_TAXA_RICH", "DECREASE")
    score$PCT_SCRAPER  <- score_1_3_5(metrics, thresh[4, ], "PCT_SCRAPER", "DECREASE")
    score$SW <- score_1_3_5(metrics, thresh[5, ], "SW", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      score$ASPT_MOD <- ifelse (metrics$ASPT_MOD < thresh$XT[1] &
                                  metrics$ASPT_MOD > thresh$XM[1],
                                ((thresh$XT[1] - metrics$ASPT_MOD) /
                                   (thresh$XT[1] - thresh$XM[1])) * 100,
                                ifelse (metrics$ASPT_MOD >= thresh$XT[1], 0,
                                        ifelse (metrics$ASPT_MOD <= thresh$XM[1], 100, "ERROR")))
      
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
#'Valleys_LV
#'
#'@param Info = Taxonomic Information
#'@param Long_Data = Taxonomic data in long format
#'@param aspt_col = Specify ASPT column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return 2011 Chessie BIBI metrics and scores for the valleys ecoregion with no limestone
#'@export
valleys_lv <- function(Info, Long, aspt_col, ffg_col, scoring) {
  prep <- prep_valleys_lv (Long)
  lv.score <- score_valleys_lv(Info, prep, aspt_col, ffg_col, scoring)
  return(lv.score)
}

