#==============================================================================
#Southeastern Plains
#==============================================================================
#'SEP: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

sep_clean_taxa <- function(Long){
  
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
  #Long[Long$CLASS %in% c("OLIGOCHAETA", "HIRUDINEA"), fam_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "OLIGOCHAETA", fam_spp] <- Long[Long$CLASS %in% "OLIGOCHAETA", "FAMILY"]
  Long[Long$CLASS %in% "HIRUDINEA", fam_spp] <- Long[Long$CLASS %in% "HIRUDINEA", "FAMILY"]
  Long[Long$FAMILY %in% c("CHIRONOMIDAE"), fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}

#==============================================================================
#'Southeastern Plains Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Southeastern Plains bioregion
#'with a Strahler Stream order <= 4 and data collected with a kicknet.
#'@export
prep_sep <-function(Long) {
  sub.sep <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                            c("65m", "65n"))
  sub.method <- prep_subset(sub.sep)
  agg.sep <- bioregion_agg(sub.method)
  final.df <- sep_clean_taxa(agg.sep)
  return(final.df)
}

#==============================================================================
#'Southeastern Plains Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@return 2011 Chessie BIBI metrics for the northern appalachian plateau
#' (Southeastern Plains) ecoregion
#'@export
metrics_sep <- function(master, Long) {
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
  #============================================================================
  Family <- wide(Long, "FAMILY")
  Order <- wide(Long, "ORDER")
  Rarefied <- prep_rare2(Long)
  metrics <- data.frame(Family[, 1:5])
  metrics$FBI <- tol_index(long.fill, master, Index = "F_HILSENHOFF",
                           Level = "FAMILY")
  metrics$PCT_GATHER <- pct_group(fam.fill, master, "GUILD", "CG", "FAMILY")
  metrics$PCT_PLECOPTERA <- pct_plecoptera(Order)
  metrics$TAXA_RICH_100 <- specnumber(Rarefied[, 6:ncol(Rarefied)])
  metrics$PCT_TRICHOPTERA_NO_TOL <- pct_trichoptera_no_tol(Family, Order)
  return(metrics)
}

#==============================================================================
#'Score Southeastern Plains Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@return Metric scores for sites within the Southeastern Plains bioregion
#'@export

score_sep<-function(Info, Long) {

  metrics <- metrics_sep(Info, Long)

  sep.thresh <- data.frame(XT = c(4.88, 56.72, 1.96, 13.4, 3.08),
                                 XM = c(4.57, 49.33, 3.71, 15, 4.26))
  rownames(sep.thresh) <- c("FBI","%GATHER", "%PLECOP", "%TAXA_RICH_100",
                                  "%TRICHOP_NO_HYDROPSYCHIDAE")
  thresh <- data.frame(sep.thresh)
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
#'Southeastern Plains
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the Southeastern Plains bioregion
#'@export
sep <- function(Info, Long) {
  prep <- prep_sep(Long)
  n.score <- score_sep(Info, prep)
  return(n.score)
}
