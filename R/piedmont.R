#==============================================================================
# Piedmont
#==============================================================================
#'Piedmont: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

pied_clean_taxa <- function(Long){
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
  Long[Long$PHYLUM %in% c("NEMATODA"), phy_spp] <- "UNIDENTIFIED"  
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
#'Piedmont Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Piedmont bioregion with a
#'Strahler Stream order <= 4
#'@export
prep_pied <- function(Long) {
  sub.pied <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                       c("45c", "45C", "45e", "45E", "45f", "45F", "45g",
                         "45G", "58h", "58H", "64d", "64D", "64c", "64C",
                         "64b", "64B", "64a", "64A"))
  sub.method <- prep_subset(sub.pied)
  agg.pied <- bioregion_agg(sub.method)
  final.df <- agg.pied
  #final.df <- pied_clean_taxa(agg.pied)
  return(final.df)
}

#==============================================================================
#'Piedmont metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the Piedmont ecoregion
#'@export
metrics_pied <- function(master, Long, taxa.rank, tol_col, ffg_col) {
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
  names(Family) <- toupper(colnames(Family))
  Order <- BIBI::wide(Long, "ORDER")
  pied_metrics <- data.frame(Family[, 1:5])
  #pied_metrics$FBI <- BIBI::tol_index(Long, Info, Index = "F_HILSENHOFF",
  #                              Level = "FAMILY")
  pied_metrics$FBI <- tol_index(long.fill, master, tol_col, taxa.rank)
  #pied_metrics$PCT_COLLECT <- BIBI::pct_group(Family, Info,
  #                                      "GUILD", c("CG", "CF"), "FAMILY")
  pied_metrics$PCT_COLLECT <- pct_attribute(fam.fill, master, ffg_col, c("CG", "CF"), taxa.rank)
  pied_metrics$PCT_DIPTERA <- BIBI::pct_diptera(Order)
  pied_metrics$PCT_EPT <- BIBI::pct_ept(Order)
  pied_metrics$SW <- BIBI::shannon(Family)
  return(pied_metrics)
}


#==============================================================================
#'Score Piedmont Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Metric scores for sites within the Piedmont bioregion
#'@export

score_pied <- function(Info, Long, taxa.rank, tol_col, ffg_col, scoring) {

  metrics <- metrics_pied(Info, Long, taxa.rank, tol_col, ffg_col)

  Piedmont_Thres <- matrix(c(4.54, 71.02, 11.72, 48.12, 1.92, 3.63, 52.71, 6.6,
                             72.24, 2.17), nrow = 5, ncol = 2)
  rownames(Piedmont_Thres) <- c("FBI", "%Col", "%Dip", "%EPT", "Shannon")
  colnames(Piedmont_Thres) <- c("XT", "XM")
  thresh <- data.frame(Piedmont_Thres)

  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])

  if(scoring %in% "DISCRETE"){
    score$FBI <- score_1_3_5(metrics, thresh[1, ], "FBI", "INCREASE")
    score$PCT_COLLECT <- score_1_3_5(metrics, thresh[2, ], "PCT_COLLECT", "INCREASE")
    score$PCT_DIPTERA <- score_1_3_5(metrics, thresh[3, ], "PCT_DIPTERA", "INCREASE")
    score$PCT_EPT <- score_1_3_5(metrics, thresh[4, ], "PCT_EPT", "DECREASE")
    score$SW <- score_1_3_5(metrics, thresh[5, ], "SW", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      score$FBI <- ifelse (metrics$FBI < thresh$XT[1] & metrics$FBI > thresh$XM[1],
                           ((thresh$XT[1] - metrics$FBI) / (thresh$XT[1] - thresh$XM[1])) * 100,
                           ifelse (metrics$FBI >= thresh$XT[1], 0,
                                   ifelse (metrics$FBI <= thresh$XM[1], 100, "ERROR")))
      
      score$PCT_COLLECT <- ifelse (metrics$PCT_COLLECT < thresh$XT[2] &
                                     metrics$PCT_COLLECT > thresh$XM[2],
                                   ((thresh$XT[2]-metrics$PCT_COLLECT) /
                                      (thresh$XT[2] - thresh$XM[2])) * 100,
                                   ifelse (metrics$PCT_COLLECT >= thresh$XT[2],  0,
                                           ifelse (metrics$PCT_COLLECT <= thresh$XM[2], 100, "ERROR")))
      
      score$PCT_DIPTERA <- ifelse (metrics$PCT_DIPTERA < thresh$XT[3] &
                                     metrics$PCT_DIPTERA > thresh$XM[3],
                                   ((thresh$XT[3] - metrics$PCT_DIPTERA) /
                                      (thresh$XT[3] - thresh$XM[3])) * 100,
                                   ifelse (metrics$PCT_DIPTERA >= thresh$XT[3], 0,
                                           ifelse (metrics$PCT_DIPTERA <= thresh$XM[3], 100,
                                                   "ERROR")))
      
      score$PCT_EPT <- ifelse (metrics$PCT_EPT > thresh$XT[4] & metrics$PCT_EPT < thresh$XM[4],
                               ((metrics$PCT_EPT - thresh$XT[4]) /
                                  (thresh$XM[4] - thresh$XT[4])) * 100,
                               ifelse (metrics$PCT_EPT <= thresh$XT[4],  0,
                                       ifelse (metrics$PCT_EPT >= thresh$XM[4], 100, "ERROR")))
      
      
      score$SW <- ifelse (metrics$SW > thresh$XT[5] & metrics$SW < thresh$XM[5],
                          ((metrics$SW - thresh$XT[5]) / (thresh$XM[5] - thresh$XT[5])) * 100,
                          ifelse(metrics$SW <= thresh$XT[5], 0,
                                 ifelse(metrics$SW >= thresh$XM[5], 100, "ERROR")))
      
    }
  }
  
  return(prep_score(score, metrics))
}

#==============================================================================
#'Piedmont
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the Piedmont bioregion
#'@export

pied <- function(Info, Long, taxa.rank, tol_col, ffg_col, scoring) {
  prep <- prep_pied(Long)
  p.score <- score_pied(Info, prep, taxa.rank, tol_col, ffg_col, scoring)
  return(p.score)
}

