#==============================================================================
#Middle Atlantic Coastal Plain
#==============================================================================
#'MAC: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

mac_clean_taxa <- function(Long){
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
 # Long[Long$CLASS %in% c("OLIGOCHAETA", "HIRUDINEA",
  #                       "GASTROPODA"), fam_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "OLIGOCHAETA", fam_spp] <- Long[Long$CLASS %in% "OLIGOCHAETA", "FAMILY"] 
  Long[Long$CLASS %in% "HIRUDINEA", fam_spp] <- Long[Long$CLASS %in% "HIRUDINEA", "FAMILY"]
  Long[Long$CLASS %in% "GASTROPODA", fam_spp] <- Long[Long$CLASS %in% "GASTROPODA", "FAMILY"]
  Long[Long$FAMILY %in% c("CHIRONOMIDAE"), fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}

#==============================================================================
#'Middle Atlantic Coastal Plain Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Middle Atlantic Coastal Plain bioregion
#'with a Strahler Stream order <= 4 and data collected with a kicknet.
#'@export
prep_mac <- function(Long) {
  sub.ma_coast <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                           c("63b", "63c", "63d", "63e", "63f"))
  sub.method <- prep_subset(sub.ma_coast)
  agg.ma_coast <- bioregion_agg(sub.method)
  final.df <- agg.ma_coast
  #final.df <- mac_clean_taxa(agg.ma_coast)
  return(final.df)
}

#==============================================================================
#'Middle Atlantic Coastal Plain Metrics
#'
#'@param Info = Taxonomic Information
#'@param Family = Familial taxonomic counts
#'@param Order = Ordinal taxonomic counts
#'@param Long_Data = Taxonomic data in long format
#'@param Level = Taxonomic level ("FAMILY" or "GENUS")
#'@param tol_col = Specify Tolerance value column.
#'@param habit_col = Specify Habit column.
#'@return 2011 Chessie BIBI metrics for the northern appalachian plateau (Middle Atlantic Coastal Plain) ecoregion
#'@export
metrics_mac <- function(master, Long, taxa.rank, tol_col, habit_col) {
  Long <- Long[!Long$TSN %in% Chessie::exclusions_2011$TSN, ]
  #============================================================================
  # These should be used for taxa attribute related metrics
  
  taxa <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")
  
  master.fill <- fill_taxa(master)
  master.fill <- unique(master.fill[, c("TSN_R", taxa)])
  #test <- (master.fill[duplicated(master.fill$TSN_R), ])
  long.fill <- Long
  long.fill <- long.fill[, !names(long.fill) %in% taxa]
  long.fill <- merge(long.fill, master.fill, by.x = "TSN", by.y = "TSN_R", all.x = T)
  fam.fill <- wide(long.fill, "FAMILY")
  rare.long.fill <- BIBI::prep_rare(long.fill, master, taxa.rank, 100, NULL, FALSE)
  rare.long.fill[, 8:18] <- fill_taxa(rare.long.fill[, 8:18])
  rare.fam.fill <- BIBI::wide(rare.long.fill, "FAMILY")
  #============================================================================

  rare.long <- BIBI::prep_rare(Long, master, taxa.rank, 100, NULL, FALSE)
  rare.long[, 8:18] <- fill_taxa(rare.long[, 8:18])
  rare.fam <- BIBI::wide(rare.long, "FAMILY")
  rare.ord <- BIBI::wide(rare.long, "ORDER")

  metrics <- data.frame(rare.fam[, 1:5])

  metrics$RICHNESS <- vegan::specnumber(rare.fam[, 6:ncol(rare.fam)])
  metrics$EPT_RICHNESS <- rich_ept(rare.long)
  metrics$PCT_EPHEMEROPTERA <- pct_ephemeroptera(rare.ord)
  metrics$FBI <- tol_index(long.fill, master, tol_col, taxa.rank)
  metrics$PCT_CLING <- pct_attribute(rare.fam.fill, master, habit_col, "CN", rare.fam.fill)

  return(metrics)
}


#==============================================================================
#'Score Middle Atlantic Coastal Plain Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param Level = Taxonomic level ("FAMILY" or "GENUS")
#'@param tol_col = Specify Tolerance value column.
#'@return Metric scores for sites within the Middle Atlantic Coastal Plain bioregion
#'@export

score_mac <-function(Info, Long, taxa.rank, tol_col, habit_col) {

  metrics <- metrics_mac(Info, Long, taxa.rank, tol_col, habit_col)

  ma_coast.thresh <- data.frame(XT = c(4.88, 56.72, 1.96, 13.4, 3.08),
                                XM = c(4.57, 49.33, 3.71, 15, 4.26))
  rownames(ma_coast.thresh) <- c("FBI","%GATHER", "%PLECOP", "%TAXA_RICH_100",
                                 "%TRICHOP_NO_HYDROPSYCHIDAE")
  thresh <- data.frame(ma_coast.thresh)
  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                  "AGENCY_CODE")])

  score$FBI <- ifelse (metrics$FBI < thresh$XT[1] & metrics$FBI > thresh$XM[1],
                       ((thresh$XT[1] - metrics$FBI) /
                          (thresh$XT[1] - thresh$XM[1])) * 100,
                       ifelse (metrics$FBI >= thresh$XT[1], 0,
                               ifelse (metrics$FBI <= thresh$XM[1], 100, "ERROR")))

  score$PCT_GATHER <- ifelse (metrics$PCT_GATHER < thresh$XT[2] &
                                metrics$PCT_GATHER > thresh$XM[2],
                              ((thresh$XT[2] - metrics$PCT_GATHER) /
                                 (thresh$XT[2] - thresh$XM[2])) * 100,
                              ifelse (metrics$PCT_GATHER >= thresh$XT[2], 0,
                                      ifelse (metrics$PCT_GATHER <= thresh$XM[2], 100,
                                              "ERROR")))

  score$PCT_PLECOPTERA <- ifelse (metrics$PCT_PLECOPTERA > thresh$XT[3] &
                                    metrics$PCT_PLECOPTERA < thresh$XM[3],
                                  ((metrics$PCT_PLECOPTERA - thresh$XT[3]) /
                                     (thresh$XM[3] - thresh$XT[3])) * 100,
                                  ifelse (metrics$PCT_PLECOPTERA <= thresh$XT[3], 0,
                                          ifelse (metrics$PCT_PLECOPTERA >= thresh$XM[3],
                                                  100, "ERROR")))

  score$TAXA_RICH_100 <- ifelse (metrics$TAXA_RICH_100 > thresh$XT[4] &
                                   metrics$TAXA_RICH_100 < thresh$XM[4],
                                 ((metrics$TAXA_RICH_100 - thresh$XT[4]) /
                                    (thresh$XM[4] - thresh$XT[4])) * 100,
                                 ifelse (metrics$TAXA_RICH_100 <= thresh$XT[4], 0,
                                         ifelse (metrics$TAXA_RICH_100 >= thresh$XM[4], 100,
                                                 "ERROR")))

  score$PCT_TRICHOPTERA_NO_TOL <- ifelse (metrics$PCT_TRICHOPTERA_NO_TOL >
                                            thresh$XT[5] &
                                            metrics$PCT_TRICHOPTERA_NO_TOL <
                                            thresh$XM[5],
                                          ((metrics$PCT_TRICHOPTERA_NO_TOL -
                                              thresh$XT[5]) /
                                             (thresh$XM[5] - thresh$XT[5])) *
                                            100,
                                          ifelse (metrics$PCT_TRICHOPTERA_NO_TOL <=
                                                    thresh$XT[5], 0,
                                                  ifelse (metrics$PCT_TRICHOPTERA_NO_TOL >=
                                                            thresh$XM[5], 100, "ERROR")))

  return(prep_score(score, metrics))
}
#==============================================================================
#'Middle Atlantic Coastal Plain
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@param tol_col = Specify Tolerance value column.
#'@param habit_col = Specify Habit column.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the Middle Atlantic Coastal Plain bioregion
#'@export
mac <- function(Info, Long, taxa.rank, tol_col, habit_col) {
  prep <- prep_mac(Long)
  n.score <- score_mac(Info, prep, taxa.rank, tol_col, habit_col)
  return(n.score)
}
