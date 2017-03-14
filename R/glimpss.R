#==============================================================================
#'Calculate the metrics from GLIMPSS
#'
#'@param Long = a long data frame of taxonomic data
#'@param index = specify the regional and season index. 
#'Requires: Mountain Spring ("MT_SP"), Mountain Summer ("MT_SU"), 
#'Plateau Spring ("PL_SP"), or Plateau Summer ("PL_SU").
#'@return Calculates the metrics from the Genus Level Index of most 
#'Probable Stream Status (GLIMPSS) (POND et al. 2013).
#'@export

glimpss_metrics <- function(Long, index){
  ord <- wide(Long, "ORDER")
  fam <- wide(Long, "FAMILY")
  gen <- wide(Long, "GENUS")
  
  metrics <- fam[, 1:5]
  
  if(index  %in% "MT_SP"){
    metrics$RICH_INTOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 0, 4)
    metrics$RICH_EPHEM <- BIBI::rich_ephemeroptera(Long, "GENUS")
    metrics$RICH_PLECOP <- BIBI::rich_plecoptera(Long, "GENUS")
    metrics$RICH_TRICHOP <- BIBI::rich_trichoptera(Long, "GENUS")
    metrics$RICH_CLING <- BIBI::rich_attribute(gen, BIBI::master,
                                               "BIBI_HABIT", "CN", "GENUS")
    metrics$HBI <- BIBI::tol_index(Long, BIBI::master, Level = "GENUS")
    metrics$PCT_DOM5 <- BIBI::pct_dom5(gen)
    metrics$PCT_EPHEM <- BIBI::pct_ephemeroptera(ord)
    metrics$PCT_ORTHOCLADIINAE <- BIBI::pct_orthocladiinae(Long)
    metrics$RICH_SCRAPE <- BIBI::rich_attribute(gen, BIBI::master,
                                                "BIBI_FFG", "SC", "GENUS")
  }
  
  if(index %in% "MT_SU"){
    metrics$RICH <- vegan::specnumber(gen[, 6:ncol(gen)])
    metrics$RICH_INTOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 0, 4)
    metrics$RICH_EPHEM <- BIBI::rich_ephemeroptera(Long, "GENUS")
    metrics$RICH_PLECOP <- BIBI::rich_plecoptera(Long, "GENUS")
    metrics$RICH_CLING <- BIBI::rich_attribute(gen, BIBI::master,
                                               "BIBI_HABIT", "CN", "GENUS")
    metrics$HBI <- BIBI::tol_index(Long, BIBI::master, Level = "GENUS")
    metrics$PCT_DOM5 <- BIBI::pct_dom5(gen)
    metrics$PCT_EPT_CHEUMATO <- BIBI::pct_ept_cheumatopsyche(Long, ord, gen)
    metrics$PCT_ORTHOCLADIINAE <- BIBI::pct_orthocladiinae(Long)
    metrics$RICH_SHRED <- BIBI::rich_attribute(gen, BIBI::master,
                                                "BIBI_FFG", "SH", "GENUS")
    
  }
  
  if(index %in% "MT_SU_60"){
    metrics$RICH_INTOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 0, 4)
    metrics$RICH_EPT <- BIBI::rich_ept(Long, "GENUS")
    metrics$RICH_CLING <- BIBI::rich_attribute(gen, BIBI::master,
                                               "BIBI_HABIT", "CN", "GENUS")
    metrics$HBI <- BIBI::tol_index(Long, BIBI::master, Level = "GENUS")
    metrics$PCT_DOM5 <- BIBI::pct_dom5(gen)
    metrics$PCT_EPT_CHEUMATO <- BIBI::pct_ept_cheumatopsyche(Long, ord, gen)
    metrics$PCT_CHIRO <- BIBI::pct_chironomidae(fam)
  }
  
  if(index %in% "PL_SP"){
    phy <- wide(Long, "PHYLUM")
    metrics$RICH_INTOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 0, 4)
    metrics$RICH_EPHEM <- BIBI::rich_ephemeroptera(Long, "GENUS")
    metrics$RICH_PLECOP <- BIBI::rich_plecoptera(Long, "GENUS")
    metrics$RICH_CLING <- BIBI::rich_attribute(gen, BIBI::master,
                                               "BIBI_HABIT", "CN", "GENUS")
    metrics$HBI <- BIBI::tol_index(Long, BIBI::master, Level = "GENUS")
    metrics$RICH_TOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 6, 10)
    metrics$PCT_EPT_CHEUMATO <- BIBI::pct_ept_cheumatopsyche(Long, ord, gen)
    metrics$PCT_CHIRO_ANNELID <- BIBI::pct_chiro_annelid(fam, phy)
  }
  
  if(index %in% "PL_SU"){
    metrics$RICH <- vegan::specnumber(gen[, 6:ncol(gen)])
    metrics$RICH_INTOL <- BIBI::rich_tolerance(gen, master, "GENUS",
                                               "BIBI_TV", 0, 3)
    metrics$RICH_EPHEM <- BIBI::rich_ephemeroptera(Long, "GENUS")
    metrics$RICH_CLING <- BIBI::rich_attribute(gen, BIBI::master,
                                               "BIBI_HABIT", "CN", "GENUS")
    metrics$HBI <- BIBI::tol_index(Long, BIBI::master, Level = "GENUS")
    metrics$PCT_DOM5 <- BIBI::pct_dom5(gen)
    metrics$PCT_EPT_CHEUMATO <- BIBI::pct_ept_cheumatopsyche(Long, ord, gen)
    metrics$PCT_CHIRO <- BIBI::pct_chironomidae(fam)
    metrics$RICH_SCRAPE <- BIBI::rich_attribute(gen, BIBI::master,
                                                "BIBI_FFG", "SC", "GENUS")
  }
  
  return(metrics)
}

#==============================================================================
#'Calculate and score the metrics from GLIMPSS
#'
#'@param Long = a long data frame of taxonomic data
#'@param index = specify the regional and season index. 
#'Requires: Mountain Spring ("MT_SP"), Mountain Summer ("MT_SU"), 
#'Plateau Spring ("PL_SP"), or Plateau Summer ("PL_SU").
#'@return Calculates and scores the metrics from the Genus Level Index of most 
#'Probable Stream Status (GLIMPSS) (POND et al. 2013).
#'@export

glimpss <- function(Long, index){
  metrics <- glimpss_metrics(Long, index)
  
  if(index %in% "MT_SP"){
    metrics$SC_RICH_INTOL <- cf_score(metrics$RICH_INTOL, "DECREASE", 19, 1)
    metrics$SC_RICH_INTOL <- ifelse(metrics$SC_RICH_INTOL > 100, 100, metrics$SC_RICH_INTOL)
    metrics$SC_RICH_EPHEM <- cf_score(metrics$RICH_EPHEM, "DECREASE", 10, 1)
    metrics$SC_RICH_EPHEM <- ifelse(metrics$SC_RICH_EPHEM > 100, 100, metrics$SC_RICH_EPHEM)
    metrics$SC_RICH_PLECOP <- cf_score(metrics$RICH_PLECOP, "DECREASE", 8, 0)
    metrics$SC_RICH_PLECOP <- ifelse(metrics$SC_RICH_PLECOP > 100, 100, metrics$SC_RICH_PLECOP)
    metrics$SC_RICH_TRICHOP <- cf_score(metrics$RICH_TRICHOP, "DECREASE", 7, 1)
    metrics$SC_RICH_TRICHOP <- ifelse(metrics$SC_RICH_TRICHOP > 100, 100, metrics$SC_RICH_TRICHOP)
    metrics$SC_RICH_CLING <- cf_score(metrics$RICH_CLING, "DECREASE", 20, 4)
    metrics$SC_RICH_CLING <- ifelse(metrics$SC_RICH_CLING > 100, 100, metrics$SC_RICH_CLING)
    metrics$SC_HBI <- cf_score(metrics$HBI, "INCREASE", 6.18, 2.23)
    metrics$SC_HBI <- ifelse(metrics$SC_HBI > 100, 100, metrics$SC_HBI)
    metrics$SC_PCT_DOM5 <- cf_score(metrics$PCT_DOM5, "INCREASE", 92, 48)
    metrics$SC_PCT_DOM5 <- ifelse(metrics$SC_PCT_DOM5 > 100, 100, metrics$SC_PCT_DOM5)
    metrics$SC_PCT_EPHEM <- cf_score(metrics$PCT_EPHEM, "DECREASE", 59.7, 0.5)
    metrics$SC_PCT_EPHEM <- ifelse(metrics$SC_PCT_EPHEM > 100, 100, metrics$SC_PCT_EPHEM)
    metrics$SC_PCT_ORTHOCLADIINAE <- cf_score(metrics$PCT_ORTHOCLADIINAE, "INCREASE", 52.7, 0.5)
    metrics$SC_PCT_ORTHOCLADIINAE <- ifelse(metrics$SC_PCT_ORTHOCLADIINAE > 100, 100, metrics$SC_PCT_ORTHOCLADIINAE)
    metrics$SC_RICH_SCRAPE <- cf_score(metrics$RICH_SCRAPE, "DECREASE", 20, 4)
    metrics$SC_RICH_SCRAPE <- ifelse(metrics$SC_RICH_SCRAPE > 100, 100, metrics$SC_RICH_SCRAPE)
    
    sc_names <- c("SC_RICH_INTOL", "SC_RICH_EPHEM", "SC_RICH_PLECOP",
                  "SC_RICH_TRICHOP", "SC_RICH_CLING", "SC_HBI",
                  "SC_PCT_DOM5", "SC_PCT_EPHEM", "SC_PCT_ORTHOCLADIINAE",
                  "SC_RICH_SCRAPE")
  }
 
  if(index %in% "MT_SU"){
    metrics$SC_RICH <- cf_score(metrics$RICH, "DECREASE", 38, 14)
    metrics$SC_RICH_INTOL <- cf_score(metrics$RICH_INTOL, "DECREASE", 15, 0)
    metrics$SC_RICH_EPHEM <- cf_score(metrics$RICH_EPHEM, "DECREASE", 9, 0)
    metrics$SC_RICH_PLECOP <- cf_score(metrics$RICH_PLECOP, "DECREASE", 7, 0)
    metrics$SC_RICH_CLING <- cf_score(metrics$RICH_CLING, "DECREASE", 19, 5)
    metrics$SC_HBI <- cf_score(metrics$HBI, "INCREASE", 6.27, 2.89)
    metrics$SC_PCT_DOM5 <- cf_score(metrics$PCT_DOM5, "INCREASE", 91.7, 51)
    metrics$SC_PCT_EPT_CHEUMATO <- cf_score(metrics$PCT_EPT_CHEUMATO, "DECREASE", 85.4, 5.3)
    metrics$SC_PCT_ORTHOCLADIINAE <- cf_score(metrics$PCT_ORTHOCLADIINAE, "INCREASE", 36.8, 0.4)
    metrics$SC_RICH_SHRED <- cf_score(metrics$RICH_SHRED, "DECREASE", 5, 0)
    sc_names <- c("SC_RICH", "SC_RICH_INTOL", "SC_RICH_EPHEM", "SC_RICH_PLECOP",
                  "SC_RICH_CLING", "SC_HBI", "SC_PCT_DOM5", "SC_PCT_EPT_CHEUMATO",
                  "SC_PCT_ORTHOCLADIINAE", "SC_RICH_SHRED")
  }
  
  if(index %in% "MT_SU_60"){
    metrics$SC_RICH_INTOL <- cf_score(metrics$RICH_INTOL, "DECREASE", 11, 1)
    metrics$SC_RICH_INTOL <- ifelse(metrics$SC_RICH_INTOL > 100, 100, metrics$SC_RICH_INTOL)
    metrics$SC_RICH_EPT <- cf_score(metrics$RICH_EPT, "DECREASE", 18, 5)
    metrics$SC_RICH_EPT <- ifelse(metrics$SC_RICH_EPT > 100, 100, metrics$SC_RICH_EPT)
    metrics$SC_RICH_CLING <- cf_score(metrics$RICH_CLING, "DECREASE", 19, 8)
    metrics$SC_RICH_CLING <- ifelse(metrics$SC_RICH_CLING > 100, 100, metrics$SC_RICH_CLING)
    metrics$SC_HBI <- cf_score(metrics$HBI, "INCREASE", 5.87, 4.06)
    metrics$SC_HBI <- ifelse(metrics$SC_HBI > 100, 100, metrics$SC_HBI)
    metrics$SC_PCT_DOM5 <- cf_score(metrics$PCT_DOM5, "INCREASE", 86.2, 49)
    metrics$SC_PCT_DOM5 <- ifelse(metrics$SC_PCT_DOM5 > 100, 100, metrics$SC_PCT_DOM5)
    metrics$SC_PCT_EPT_CHEUMATO <- cf_score(metrics$PCT_EPT_CHEUMATO, "DECREASE", 76.9, 13.8)
    metrics$SC_PCT_EPT_CHEUMATO <- ifelse(metrics$SC_PCT_EPT_CHEUMATO > 100, 100, metrics$SC_PCT_EPT_CHEUMATO)
    metrics$SC_PCT_CHIRO <- cf_score(metrics$PCT_CHIRO, "INCREASE", 46.1, 1.5)
    metrics$SC_PCT_CHIRO <- ifelse(metrics$SC_PCT_CHIRO > 100, 100, metrics$SC_PCT_CHIRO)
    sc_names <- c("SC_RICH_INTOL", "SC_RICH_EPT",
                  "SC_RICH_CLING", "SC_HBI", "SC_PCT_DOM5", "SC_PCT_EPT_CHEUMATO",
                  "SC_PCT_CHIRO")
  }
  
  if(index %in% "PL_SP"){
    metrics$SC_RICH_INTOL <- cf_score(metrics$RICH_INTOL, "DECREASE", 15, 1)
    metrics$SC_RICH_EPHEM <- cf_score(metrics$RICH_EPHEM, "DECREASE", 10, 1)
    metrics$SC_RICH_PLECOP <- cf_score(metrics$RICH_PLECOP, "DECREASE", 7, 0)
    metrics$SC_RICH_CLING <- cf_score(metrics$RICH_CLING, "DECREASE", 17, 3)
    metrics$SC_HBI <- cf_score(metrics$HBI, "INCREASE", 6.54, 2.52)
    metrics$SC_RICH_TOL <- cf_score(metrics$RICH_TOL, "INCREASE", 65, 0)
    metrics$SC_PCT_EPT_CHEUMATO <- cf_score(metrics$PCT_EPT_CHEUMATO, "DECREASE", 91.2, 2.6)
    metrics$SC_PCT_CHIRO_ANNELID <- cf_score(metrics$PCT_CHIRO_ANNELID, "INCREASE", 83.5, 1.8)
    sc_names <- c("SC_RICH_INTOL", "SC_RICH_EPHEM", "SC_RICH_PLECOP",
                  "SC_RICH_CLING", "SC_HBI", "SC_RICH_TOL",
                  "SC_PCT_EPT_CHEUMATO", "SC_PCT_CHIRO_ANNELID")
    }
  
  if(index %in% "PL_SU"){
    metrics$SC_RICH <- cf_score(metrics$RICH, "DECREASE", 34, 14)
    metrics$SC_RICH_INTOL <- cf_score(metrics$RICH_INTOL, "DECREASE", 7, 0)
    metrics$SC_RICH_EPHEM <- cf_score(metrics$RICH_EPHEM, "DECREASE", 7, 0)
    metrics$SC_RICH_CLING <- cf_score(metrics$RICH_CLING, "DECREASE", 15, 0)
    metrics$SC_HBI <- cf_score(metrics$HBI, "INCREASE", 6.24, 3.73)
    metrics$SC_PCT_DOM5 <- cf_score(metrics$PCT_DOM5, "INCREASE", 91.5, 53.3)
    metrics$SC_PCT_EPT_CHEUMATO <- cf_score(metrics$PCT_EPT_CHEUMATO, "DECREASE", 67.1, 1.4)
    metrics$SC_PCT_CHIRO <- cf_score(metrics$PCT_CHIRO, "INCREASE", 69.1, 4)
    metrics$SC_RICH_SCRAPE <- cf_score(metrics$RICH_SCRAPE, "DECREASE", 7, 4)
    sc_names <- c("SC_RICH", "SC_RICH_INTOL", "SC_RICH_EPHEM",
                  "SC_RICH_CLING", "SC_HBI", "SC_PCT_DOM5",
                  "SC_PCT_EPT_CHEUMATO", "SC_PCT_CHIRO",
                  "SC_RICH_SCRAPE")
  }
  

  metrics$glimpss <- apply(metrics[, sc_names], 1, FUN = mean)
  
  return(metrics)
}

#==============================================================================
#'Calculate the metrics from GLIMPSS
#'
#'@param metrics = a data frame containing metric values
#'@param disturbance = does the metric increase or decrease with disturbance? 
#'Requires: "INCREASE" or "DECREASE".
#'@param sc_ceiling = the upper threshold typically defined by the
#' 5th or 95th percentile.
#'@param sc_floor = the lower threshold typically defined by the
#' 5th or 95th percentile.
#'@return A gradient scoring approach.
#'@export

cf_score <- function(metrics, disturbance, sc_ceiling, sc_floor){
  
  if(disturbance %in% "DECREASE"){
    final.df <- ((metrics - sc_floor) / (sc_ceiling - sc_floor)) * 100
  }
  
  if(disturbance %in% "INCREASE"){
    final.df <- ((sc_ceiling - metrics) / (sc_ceiling - sc_floor)) * 100
  }
  
  return(final.df)
}


