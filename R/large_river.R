#==============================================================================
# Large River's
#==============================================================================
#'Large River's PREP
#'
#'@param Long = Taxonomic data in long format
#'@return  Subset the data to only sites in the Chesapeake Bay  watershed with
#' a Strahler Stream order >= 5
#'@export
prep_large <- function(Long) {
  sub.large <- subset(Long, Long$STRAHLER_STREAM_ORDER >= 5)
  
  agg.large <- bioregion_agg(sub.large)
  return(agg.large)
}

#==============================================================================
#'Piedmont Large River's PREP
#'
#'@param Long = Taxonomic data in long format
#'@return Subset the data to only sites in the Piedmont bioregion with a
#'Strahler Stream order >= 5
#'@export
prep_pied_large <- function(Long) {
  sub.L_pied <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                         c("45c", "45C", "45e", "45E", "45f", "45F", "45g",
                           "45G", "58h", "58H", "64d", "64D", "64c", "64C",
                           "64b", "64B", "64a", "64A") &
                         Long$STRAHLER_STREAM_ORDER >= 5)
  agg.L_pied <- bioregion_agg(sub.L_pied)
  return(agg.L_pied)
}

#==============================================================================
#'Ridges Large River's PREP
#'
#'@param Long = Taxonomic data in long format
#'@return Subset the data to only sites in the Ridges bioregion with a
#'Strahler Stream order >= 5
#'@export
prep_ridges_large <- function(Long) {
  sub.ridges <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                         c("66a", "66A", "66b", "66B", "67d", "67D", "67c", "67C",
                           "67i", "67I", "67h", "67H", "69a", "69A", "69b", "69B",
                           "70c", "70C") & Long$STRAHLER_STREAM_ORDER >= 5)
  agg.L_ridges <- bioregion_agg(sub.L_ridges)
  return(agg.L_ridges)
}

#==============================================================================
#'Valleys-All Large River's PREP
#'
#'@param Long = Taxonomic data in long format
#'@return Subset the data to only sites in the Valleys bioregion with a
#'Strahler Stream order >= 5
#'@export
prep_valleys_large <- function(Long) {
  sub.L_valleys <- subset(Long, Long$ECOREGION_LEVEL_4 %in%
                            c("67e", "67E", "67a", "67A", "67b", "67B", "67f", "67F",
                              "67g", "67G") & Long$STRAHLER_STREAM_ORDER >= 5)
  agg.L_valleys <- bioregion_agg(sub.L_valleys)
  return(agg.L_valleys)
}


