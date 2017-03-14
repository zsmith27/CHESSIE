#==============================================================================
#NAPU
#==============================================================================
#'NAPU: Standarize Taxonomic Resolution
#'
#'@param Long = Taxonomic data in long format
#'@return  Removes and adjusts taxonomic resolution to make data sets from
#'multiple sources comparable.
#'@export

napu_clean_taxa <- function(Long){
  
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
  
  class <- wide(Long, "CLASS")
  gerridae <- if("GERRIDAE" %in% names(class)){
    (class$GERRIDAE / rowSums(class[,6:ncol(class)])) * 100
  }else{
    0
  }
  cat("Gerridae were excluded from the data.\
The maximum percentage of Gerridae found in a sample was", round(max(gerridae), 2),"%.")
  
  mesoveliidae <- if("MESOVELIIDAE" %in% names(class)){
    (class$MESOVELIIDAE / rowSums(class[,6:ncol(class)])) * 100
  }else{
    0
  }
  cat("Mesoveliidae were excluded from the data.\
The maximum percentage of Mesoveliidae found in a sample was",
      round(max(mesoveliidae), 2),"%.")
  
  veliidae <- if("VELIIDAE" %in% names(class)){
    (class$VELIIDAE / rowSums(class[,6:ncol(class)])) * 100
  }else{
    0
  }
  cat("Veliidae were excluded from the data.\
The maximum percentage of Veliidae found in a sample was",
      round(max(veliidae), 2),"%.")
  
  
  phy_spp <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
               "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$PHYLUM %in% c("NEMATODA"), phy_spp] <- "UNIDENTIFIED"  
  Long[Long$SUBPHYLUM %in% c("CHELICERATA"), phy_spp] <- "UNIDENTIFIED"
  Long[Long$FAMILY %in% c("GERRIDAE", "MESOVELIIDAE",
                          "VELIIDAE"), phy_spp] <- "UNIDENTIFIED"
  
  subphy_spp <- c("SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
                  "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$PHYLUM %in% c("NEMERTEA"), subphy_spp] <- "NEMERTEA"
  
  subclass_spp <- c("SUBCLASS", "ORDER", "SUBORDER", "FAMILY",
                    "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  Long[Long$CLASS %in% c("OLIGOCHAETA"), subclass_spp] <- Long[Long$CLASS %in% c("OLIGOCHAETA"), "CLASS"]
  
  subord_spp <- c("SUBORDER", "FAMILY", "SUBFAMILY",
                  "TRIBE", "GENUS", "SPECIES")
  
  #Long[Long$ORDER %in% c("CUMACEA", "LEPIDOPTERA"), subord_spp] <- "UNIDENTIFIED"
  if("CUMACEA" %in% Long$ORDER){
    Long[Long$ORDER %in% "CUMACEA", subord_spp] <- Long[Long$ORDER %in% "CUMACEA", "ORDER"]
  }
  
  Long[Long$ORDER %in% "LEPIDOPTERA", subord_spp] <- Long[Long$ORDER %in% "LEPIDOPTERA", "ORDER"]
  
  fam_spp <- c("SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  #Long[Long$CLASS %in% c("GASTROPODA", "PELECYPODA"), fam_spp] <- "UNIDENTIFIED"
  Long[Long$CLASS %in% "GASTROPODA", fam_spp] <- Long[Long$CLASS %in% "GASTROPODA", "FAMILY"]
  
  if("PELECYPODA" %in% Long$CLASS){
    Long[Long$CLASS %in% "PELECYPODA", fam_spp] <- Long[Long$CLASS %in% "PELECYPODA", "FAMILY"]
  }
  
  Long[Long$ORDER %in% c("DECAPODA"), fam_spp] <- Long[Long$ORDER %in% c("DECAPODA"), "FAMILY"]
  
  #Long[Long$FAMILY %in% c("ANTHURIDAE", "IDOTEIDAE", "OEDICEROTIDAE",
  #                        "CORIXIDAE", "PELTOPERLIDAE", "SISYRIDAE",
  #                        "PSYCHODIDAE", "PTYCHOPTERIDAE", "DIXIDAE", "DRYOPIDAE",
  #                        "CERATOPOGONIDAE", "TABANIDAE", "DOLICHOPODIDAE",
  #                        "STRATIOMYIDAE", "EPHYDRIDAE", "MUSCIDAE",
  #                        "ANTHOMYIIDAE", "SCATHOPHAGIDAE", "CHIRONOMIDAE"),
  #     fam_spp] <- "UNIDENTIFIED"
  Long[Long$FAMILY %in% "ANTHURIDAE", fam_spp] <- "ANTHURIDAE"
  Long[Long$FAMILY %in% "IDOTEIDAE", fam_spp] <- "IDOTEIDAE"
  Long[Long$FAMILY %in% "OEDICEROTIDAE", fam_spp] <- "OEDICEROTIDAE"
  Long[Long$FAMILY %in% "CORIXIDAE", fam_spp] <- "CORIXIDAE"
  Long[Long$FAMILY %in% "PELTOPERLIDAE", fam_spp] <- "PELTOPERLIDAE"
  Long[Long$FAMILY %in% "SISYRIDAE", fam_spp] <- "SISYRIDAE"
  Long[Long$FAMILY %in% "PSYCHODIDAE", fam_spp] <- "PSYCHODIDAE"
  Long[Long$FAMILY %in% "PTYCHOPTERIDAE", fam_spp] <- "PTYCHOPTERIDAE"
  Long[Long$FAMILY %in% "DIXIDAE", fam_spp] <- "DIXIDAE"
  Long[Long$FAMILY %in% "DRYOPIDAE", fam_spp] <- "DRYOPIDAE"
  Long[Long$FAMILY %in% "CERATOPOGONIDAE", fam_spp] <- "CERATOPOGONIDAE"
  Long[Long$FAMILY %in% "TABANIDAE", fam_spp] <- "TABANIDAE"
  Long[Long$FAMILY %in% "DOLICHOPODIDAE", fam_spp] <- "DOLICHOPODIDAE"
  Long[Long$FAMILY %in% "STRATIOMYIDAE", fam_spp] <- "STRATIOMYIDAE"
  Long[Long$FAMILY %in% "EPHYDRIDAE", fam_spp] <- "EPHYDRIDAE"
  Long[Long$FAMILY %in% "MUSCIDAE", fam_spp] <- "MUSCIDAE"
  Long[Long$FAMILY %in% "ANTHOMYIIDAE", fam_spp] <- "ANTHOMYIIDAE"
  Long[Long$FAMILY %in% "SCATHOPHAGIDAE", fam_spp] <- "SCATHOPHAGIDAE"
  Long[Long$FAMILY %in% "CHIRONOMIDAE", fam_spp] <- "CHIRONOMIDAE"
  
  return(Long)
}

#==============================================================================
#'NAPU Data PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the NAPU bioregion with a
#'Strahler Stream order <= 4
#'@export
prep_napu<-function(Long) {
  sub.napu <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                       c("60a", "60b", "60d", "60e", "83f"))
  sub.method <- prep_subset(sub.napu)
  agg.napu <- bioregion_agg(sub.method)
  final.df <- agg.napu
  #final.df <- napu_clean_taxa(agg.napu)
  return(final.df)
}

#==============================================================================
#'NAPU Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@return 2011 Chessie BIBI metrics for the northern appalachian plateau
#' (NAPU) ecoregion.
#'@export
metrics_napu <- function(master, Long, taxa.rank, tol_col, ffg_col) {
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
  
  #test <- long.fill[!long.fill$TSN %in% master.fill$TSN_FINAL, ]
  #test <- Family[!Family$EVENT_ID %in% fam.fill$EVENT_ID, ]
  #============================================================================
  Family <- BIBI::wide(Long, "FAMILY")
  Order <- BIBI::wide(Long, "ORDER")
  #Rarefied <- prep_rare2(Long)
  rare.long <- BIBI::prep_rare(Long, master, taxa.rank, 100, NULL, FALSE)
  rare.long[, taxa] <- fill_taxa(rare.long[, taxa])
  Rarefied <- BIBI::wide(rare.long, "FAMILY")
  
  metrics <- data.frame(Family[, 1:5])
  
  #test <- dplyr::anti_join(Rarefied, metrics)
  #test <- metrics[duplicated(metrics),]
  #metrics$FBI <- tol_index(Long, Info, Index = "F_HILSENHOFF",
  #                         Level = "FAMILY")
  #metrics$PCT_GATHER <- pct_group(Family, Info, "GUILD", "CG", "FAMILY")
  
  metrics$FBI <- tol_index(long.fill, master, tol_col, taxa.rank)
  metrics$PCT_GATHER <- pct_attribute(fam.fill, master, ffg_col, "CG", taxa.rank)
  metrics$PCT_PLECOPTERA <- pct_plecoptera(Order)
  metrics$TAXA_RICH_100 <- vegan::specnumber(Rarefied[, 6:ncol(Rarefied)])
  metrics$PCT_NON_HYDROP_TRICHOPTERA <- pct_non_hydrop_trichoptera(Order, Family)

 
  return(metrics)
}

#==============================================================================
#'Score NAPU Metrics
#'
#'@param Info = Taxonomic Information
#'@param Long = Taxonomic data in long format
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Metric scores for sites within the NAPU bioregion
#'@export

score_napu<-function(Info, Long, taxa.rank, tol_col, ffg_col, scoring) {

  metrics <- metrics_napu(Info, Long, taxa.rank, tol_col, ffg_col)

  napu.thresh <- data.frame(XT = c(4.88, 56.72, 1.96, 13.4, 3.08),
                            XM = c(4.57, 49.33, 3.71, 15, 4.26))
  rownames(napu.thresh) <- c("FBI","%GATHER", "%PLECOP", "%TAXA_RICH_100",
                             "%TRICHOP_NO_HYDROPSYCHIDAE")
  thresh <- data.frame(napu.thresh)
  score <- data.frame(metrics[, c("EVENT_ID", "STATION_ID", "DATE",
                                  "SAMPLE_NUMBER", "AGENCY_CODE")])
  if(scoring %in% "DISCRETE"){
    score$FBI <- score_1_3_5(metrics, napu.thresh[1, ], "FBI", "INCREASE")
    score$PCT_GATHER <- score_1_3_5(metrics, napu.thresh[2, ], "PCT_GATHER", "INCREASE")
    score$PCT_PLECOPTERA <- score_1_3_5(metrics, napu.thresh[3, ], "PCT_PLECOPTERA", "DECREASE")
    score$TAXA_RICH_100 <- score_1_3_5(metrics, napu.thresh[4, ], "TAXA_RICH_100", "DECREASE")
    score$PCT_NON_HYDROP_TRICHOPTERA <- score_1_3_5(metrics, napu.thresh[5, ], "PCT_NON_HYDROP_TRICHOPTERA", "DECREASE")
  }else{
    if(scoring %in% "GRADIENT"){
      
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
      
      score$PCT_NON_HYDROP_TRICHOPTERA <- ifelse (metrics$PCT_NON_HYDROP_TRICHOPTERA>
                                                   thresh$XT[5] &
                                                   metrics$PCT_NON_HYDROP_TRICHOPTERA<
                                                   thresh$XM[5],
                                                 ((metrics$PCT_NON_HYDROP_TRICHOPTERA-
                                                     thresh$XT[5]) /
                                                    (thresh$XM[5] - thresh$XT[5])) *
                                                   100,
                                                 ifelse (metrics$PCT_NON_HYDROP_TRICHOPTERA<=
                                                           thresh$XT[5], 0,
                                                         ifelse (metrics$PCT_NON_HYDROP_TRICHOPTERA>=
                                                                   thresh$XM[5], 100, "ERROR")))
    }
  }



  return(prep_score(score, metrics))
}
#==============================================================================
#'NAPU
#'
#'@param Long = Taxonomic data in long format
#'@param Info = Taxonomic Information
#'@param tol_col = Specify Tolerance value column.
#'@param ffg_col = Specify Funtional Feeding Group column.
#'@param scoring = If scoring is set to "DISCRETE" the 1-3-5 method will be used.
#'If scoring is set to "GRADIENT" a continuous scoring procedure will be used.
#'@return Prepare the taxonomic data and score all metrics for each sites
#' within the NAPU bioregion
#'@export
napu <- function(Info, Long, taxa.rank, tol_col, ffg_col, scoring) {
  prep <- prep_napu(Long)
  n.score <- score_napu(Info, prep, taxa.rank, tol_col, ffg_col, scoring)
  return(n.score)
}
