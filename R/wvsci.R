#==============================================================================
#'Calculate the metrics from WVSCI
#'
#'@param Long = a long data frame of taxonomic data
#'@return Calculates the metrics from the West Virginia Stream Condition Index
#' (WVSCI) (Tetra Tech, Inc. 2000).
#'@export

wvsci_metrics <- function(Long){
  ord <- wide(Long, "ORDER")
  fam <- wide(Long, "FAMILY")
  
  metrics <- fam[, 1:5]
  metrics$RICH <- vegan::specnumber(fam[, 6:ncol(fam)])
  metrics$EPT_RICH <- BIBI::rich_ept(Long, "GENUS")
  metrics$PCT_EPT <- BIBI::pct_ept(ord)
  metrics$PCT_CHIRONOMIDAE <- BIBI::pct_chironomidae(fam)
  metrics$PCT_DOM2 <- BIBI::pct_dom2(fam)
  metrics$HBI <- BIBI::tol_index(Long, BIBI::master)
  
  return(metrics)
}

#==============================================================================
#'Calculate and score the metrics from WVSCI
#'
#'@param Long = a long data frame of taxonomic data
#'@return Calculates and scores the metrics from the West Virginia Stream
#' Condition Index (WVSCI) (Tetra Tech, Inc. 2000).
#'@export

wvsci <- function(Long){
  metrics <- wvsci_metrics(Long)
  metrics$SC_RICH <- (metrics$RICH / 21) * 100
  metrics$SC_RICH <- ifelse(metrics$SC_RICH > 100, 100, metrics$SC_RICH)
  metrics$SC_EPT_RICH <- (metrics$EPT_RICH / 13) * 100
  metrics$SC_EPT_RICH <- ifelse(metrics$SC_EPT_RICH > 100, 100, metrics$SC_EPT_RICH)
  metrics$SC_PCT_EPT <- (metrics$PCT_EPT / 91.9) * 100
  metrics$SC_PCT_EPT <- ifelse(metrics$SC_PCT_EPT > 100, 100, metrics$SC_PCT_EPT)
  metrics$SC_PCT_CHIRONOMIDAE <- ((100 - metrics$PCT_CHIRONOMIDAE) / (100 - 0.98)) * 100
  metrics$SC_PCT_CHIRONOMIDAE<- ifelse(metrics$SC_PCT_CHIRONOMIDAE > 100, 100, metrics$SC_PCT_CHIRONOMIDAE)
  metrics$SC_DOM2 <- ((100 - metrics$PCT_DOM2) / (100 - 36)) * 100
  metrics$SC_DOM2 <- ifelse(metrics$SC_DOM2 > 100, 100, metrics$SC_DOM2)
  metrics$SC_HBI <- ((10 - metrics$HBI) / (10 - 2.9)) * 100
  metrics$SC_HBI <- ifelse(metrics$SC_HBI> 100, 100, metrics$SC_HBI)
  
  sc_names <- c("SC_RICH", "SC_EPT_RICH", "SC_PCT_EPT",
                "SC_PCT_CHIRONOMIDAE", "SC_EPT_RICH", "SC_HBI")
  metrics$wvsci <- apply(metrics[, sc_names], 1, FUN = mean)
  
  return(metrics)
}