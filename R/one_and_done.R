#'One and done
#'
#'@param Info = Taxonomic information
#'@param Master = Taxa count data
#'@param Taxa = Taxonomic counts in long data format
#'@param Taxonomy = TSN #'s with necessary taxonomic levels (CLASS, ORDER, FAMILY)
#'@param WQ = Water quality data in a long data format
#'@param Habitat = Habitat data in a long data format
#'@param TAB_EVENT = EVENT_ID
#'@param TAB_STATIONS = Station information (STATION_ID, ECOREGION, STRAHLER)
#'@param TAB_PROJECT = Project information
#'@return Merge taxa water quality, and habitat data into one large data frame.
#'@export
one_and_done <- function(Info, Master, Taxa, WQ, Habitat,
                         TAB_EVENT, TAB_STATIONS, TAB_PROJECT) {
  #see prep_data function
  Data_PREP <- prep_data(Master, Taxa, WQ, Habitat, TAB_EVENT,
                         TAB_STATIONS, TAB_PROJECT)

  #Middle Atlantic Coastal Plain Bioregion
  #MA_Coast <- mac(Info, Data_PREP)
  #Southeastern Plains Bioregion
  #SE_Plains <- se_plains(Info, Data_PREP)
  #Piedmont Bioregion
  Piedmont <- pied(Info, Data_PREP)
  #North Central Appalachians Bioregion
  NCA <- nca(Info, Data_PREP)
  # Northern Appalachian Plateau and Uplands Bioregion
  NAPU <- napu(Info, Data_PREP)
  #Ridges Bioregion
  Ridges <- ridges(Info, Data_PREP)
  #All of the data for the Valley Bioregion
  Valleys_All <- valleys(Info, Data_PREP)
  #Valley Bioregion data excluding limestone sites
  #Valleys_NLV <- valleys_nlv(Info, Data_PREP)
  #Valley Bioregion data exclusively limestone sites
  #Valleys_LV <- valleys_lv(Info, Data_PREP)

  #Append the output from each bioregion to form one data frame
  new <- rbind(Piedmont, NCA, NAPU, Ridges, Valleys_All)
  #, Valleys_NLV, Valleys_LV, MA_Plain, SE_Plains)

  #merge site related information to the data frame
  #merg1 <- merge(TAB_STATIONS[, c("STATION_ID", "ECOREGION_LEVEL_4",
  #                                "WATERBODY_NAME", "STRAHLER_STREAM_ORDER",
  #                               "LATITUDE", "LONGITUDE")], new,
  #               by = "STATION_ID", all.y = TRUE)

  #merg2 <- merge(EVENT[, c("EVENT_ID", "STATION_ID")], merg1,
  #               by = c("EVENT_ID", "STATION_ID"), all.y = TRUE)

  #Reorder the data frame for final output
  #merg2 <- merg2[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
  #                  "AGENCY_CODE", "WATERBODY_NAME", "LATITUDE", "LONGITUDE",
  #                   "ECOREGION_LEVEL_4", "STRAHLER_STREAM_ORDER",
  #                   "METRIC", "METRIC_VALUE", "METRIC_SCORE", "MEAN",
  #                   "RATING")]
  merg2 <- new[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                   "AGENCY_CODE", "METRIC", "METRIC_VALUE", "METRIC_SCORE", "MEAN",
                   "RATING")]

  #Rename the columns for final output
  #colnames(merg2) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
  #                     "AGENCY_CODE", "WATERBODY_NAME", "LATITUDE", "LONGITUDE",
  #                     "ECOREGION_LEVEL_4", "STRAHLER_STREAM_ORDER", "METRIC",
  #                     "METRIC_VALUE", "METRIC_SCORE", "BIBI_SCORE",
  #                     "BIBI_RATING")
  colnames(merg2) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                       "AGENCY_CODE", "METRIC", "METRIC_VALUE", "METRIC_SCORE",
                       "BIBI_SCORE", "BIBI_RATING")

  return(merg2)
}
