#==============================================================================
#'Import MS Access Data Sheets
#'
#'@param odbc.name = A character string that refers to the odbc connection
#'name.
#'@return Assigns Chessie BIBI related ecoregion and bioregion distinctions.
#'@export
import_msa_data <- function(odbc.name){
  # Create an empty list to store the database sheets.
  final.list <- list()
  #============================================================================
  # Establish the connection with the Chessie BIBI Microsoft Access database.
  channel <- RODBC::odbcConnect(odbc.name)
  #============================================================================
  # Import the necessary database sheets (tabs) into a list.
  final.list[["tab.event"]] <- RODBC::sqlFetch(channel, "TAB_EVENT",
                                               stringsAsFactors = FALSE)
  final.list[["tab.stations"]] <- RODBC::sqlFetch(channel, "TAB_STATIONS",
                                                  stringsAsFactors = FALSE)
  final.list[["tab.project"]] <- RODBC::sqlFetch(channel, "TAB_PROJECT",
                                                 stringsAsFactors = FALSE)
  final.list[["tab.wq"]] <- RODBC::sqlFetch(channel, "TAB_WQ_DATA",
                                            stringsAsFactors = FALSE)
  final.list[["tab.wq_param"]] <- RODBC::sqlFetch(channel, "TAB_PARAMETER_WQ",
                                                  stringsAsFactors = FALSE)
  final.list[["tab.habitat"]] <- RODBC::sqlFetch(channel, "TAB_HABITAT_ASSESSMENT",
                                                 stringsAsFactors = FALSE)
  final.list[["tab.hab_param"]] <- RODBC::sqlFetch(channel, "TAB_HABITAT_PARAMETERS",
                                                   stringsAsFactors = FALSE)
  final.list[["tab.taxa"]] <- RODBC::sqlFetch(channel, "TAB_TAXONOMIC_COUNT",
                                              stringsAsFactors = FALSE)
  final.list[["tab.huc"]] <- RODBC::sqlFetch(channel, "TAB_HUC_12",
                                             stringsAsFactors = FALSE)
  final.list[["tab.program"]] <- RODBC::sqlFetch(channel, "TAB_PROGRAM",
                                                 stringsAsFactors = FALSE)
  #============================================================================
  # Close the ODBC connection
  RODBC::odbcCloseAll()
  #============================================================================
  # End import_msa_data function.
  return(final.list)
}

#==============================================================================
#'Ecoregion and Bioregion Modifications
#'
#'@param Long = Long data frame.
#'@return Assigns Chessie BIBI related ecoregion and bioregion distinctions.
#'@export

prep_ecoregion <- function(long){
  long$ECO3 <- substring(long$ECOREGION_LEVEL_4, 1, 2)
  long$ECO3 <- ifelse(!(long$SUBREGION_DESCRIPTION %in% "SUSQUEHANNA") & long$ECO3 == 67, "67S",
                      ifelse((long$SUBREGION_DESCRIPTION %in% "SUSQUEHANNA") & long$ECO3 == 67, "67N", long$ECO3))
  long$ECO3 <- ifelse(!(long$SUBREGION_DESCRIPTION %in% c("SUSQUEHANNA", "UPPER CHESAPEAKE")) & long$ECO3 == 64, "64S",
                      ifelse((long$SUBREGION_DESCRIPTION %in% c("SUSQUEHANNA", "UPPER CHESAPEAKE")) & long$ECO3 == 64, "64N", long$ECO3))
  long$ECO3 <- ifelse(!(long$SUBREGION_DESCRIPTION %in% c("SUSQUEHANNA", "UPPER CHESAPEAKE")) & long$ECO3 == 69, "69S",
                      ifelse((long$SUBREGION_DESCRIPTION %in% c("SUSQUEHANNA", "UPPER CHESAPEAKE")) & long$ECO3 == 69, "69N", long$ECO3))
  #long$ECO3 <- ifelse((long$ECOREGION_LEVEL_4 %in% c("67c", "67d", "67h", "67i")) & long$ECO3 %in% "67N", "67NR",
  #                      ifelse((long$ECOREGION_LEVEL_4 %in% c("67c", "67d", "67h", "67i")) & long$ECO3 %in% "67S", "67SR",
  #                             ifelse((long$ECOREGION_LEVEL_4 %in% c("67a", "67b", "67e", "67f", "67g")) & long$ECO3 %in% "67N", "67NV",
  #                                    ifelse((long$ECOREGION_LEVEL_4 %in% c("67a", "67b", "67e", "67f", "67g")) & long$ECO3 %in% "67S", "67SV", long$ECO3))))
  long$ECO3 <- ifelse(long$ICPRB_BIOREGION_ID %in% "SGV", "67SGV",
                      ifelse(long$ICPRB_BIOREGION_ID %in% "NGV", "67NGV",
                             ifelse(long$ICPRB_BIOREGION_ID %in% "BLUE", "67BLUE", long$ECO3)))
  
  long$BIOREGION <- ifelse(long$ECO3 %in% c("58", "60", "83"), "NAPU",
                           ifelse(long$ECO3 %in% c("62"), "NCA",
                                  #ifelse(long$ECO3 %in% c("69N"), "UCA",
                                  #ifelse(long$ECO3 %in% c("69S"), "LCA",
                                  ifelse(long$ECO3 %in% c("67N"), "NRV",
                                         ifelse(long$ECO3 %in% c("64N", "67NGV"), "UNP",  
                                                #ifelse(long$ECO3 %in% c("63", "65", "COAST", 
                                                ifelse(long$ECO3 %in% c("69S", "69N"), "CA",   
                                                       ifelse(long$ECO3 %in% c("67S"), "SRV",
                                                              #ifelse(long$ECO3 %in% c("64N"), "UNP",
                                                              ifelse(long$ECO3 %in% c("64S"), "LNP",
                                                                     ifelse(long$ECO3 %in% "45", "SPIED", 
                                                                            ifelse(long$ECO3 %in% c("66", "67BLUE"), "BLUE",
                                                                                   ifelse(long$ECO3 %in% "65", "SEP",
                                                                                          #ifelse(long$ECO3 %in% "67NGV", "NGV",
                                                                                          ifelse(long$ECO3 %in% "67SGV", "SGV",  
                                                                                                 ifelse(long$ECO3 %in% "63", "MAC", "ERROR"))))))))))))
  
  return(long)  
}

#==============================================================================
#'Check for Differences in the Taxa Identified by Agency
#'
#'@param Long = Long data frame, typically representing a single bioregion.
#'@param Level = The taxonomic rank to perform the comparison.
#'@return Compares the taxa identified by each agency.
#'@export

agency_taxa_diff <- function(Long, Level){
  wide.df <- wide(Long, Level)
  agg<- aggregate(wide.df[, 6:ncol(wide.df)], by = list(wide.df$AGENCY_CODE), FUN = sum)
  agg[, 2:ncol(agg)] <- ifelse(agg[, 2:ncol(agg)] > 0, 1, 0)
  cols_to_drop <- c(rep(TRUE, 1), colSums(agg[, 2:ncol(agg)]) < nrow(agg) &
                      colSums(agg[, 2:ncol(agg)]) > 0)
  final.df <- agg[, cols_to_drop]
  return(final.df)
}

#==============================================================================
#'Clean Up Character Vectors
#'
#'@param x = A vector.
#'@return Cleans up a character vector by verifying that it is class character,
#'all characters are uppercase, and all leading/trailing white space has been
#'removed.
#'@importFrom magrittr "%>%"
#'@export

clean_char <- function(x){
  final.vec <- as.character(x) %>% toupper() %>% trimws()
  return(final.vec)
}
#==============================================================================
#'Prepare Taxonomic Counts
#'
#'@param Master = Taxa count data
#'@param Taxa = Taxonomic counts in long data format
#'@return Merges a taxonomic counts data frame with a master taxa list
#'containing the following information for each taxon when applicable:
#'TSN, Phylum, Subphylum, Class, Subclass, Order, Suborder, Family,
#'Subfamily, Tribe, Genus, and Species.
#'@import data.table
#'@export

prep_taxa <- function(master.df, taxa.df, clean.taxa){
  # Taxonomic columns.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
                 "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS",
                 "SPECIES")
  #============================================================================
  # Remove any leading or trailing white space from the final TSN prior to
  # matching or merging.
  master.df$TSN_FINAL <- trimws(master.df$TSN_FINAL)
  taxa.df$TSN_FINAL <- trimws(taxa.df$TSN_FINAL)
  #============================================================================
  # Identify any TSN values that are present in the taxonomic counts data but
  # not in the master taxa list.
  miss.tsn <- unique(taxa.df$TSN_FINAL[!taxa.df$TSN_FINAL %in% master.df$TSN_FINAL])
  # If there are TSN values missing from the master taxa list provide a warning.
  if(length(miss.tsn) > 0){
    warning(paste0("Missing TSN: ", paste(miss.tsn, collapse = ", ")))
  } 
  # The warning indicates that there are issues but this script allows us to
  # move forward before we get a chance to update the master taxa list.
  # The missing values must be fixed or we risk misrepresenting the data.
  taxa.df <- taxa.df[!taxa.df$TSN_FINAL %in% miss.tsn, ]
  #============================================================================
  # Convert data.frames to class data.table to speed up merging below.
  # Subset the master.df to only include taxonomic hierarchy columns.
  master.df <- data.table::data.table(master.df[, c("TSN_FINAL", "TSN_R",
                                                    taxa.cols)])
  taxa.df <- data.table::data.table(taxa.df)
  #============================================================================
  # Join the master taxa list with the taxonomic counts by the final TSN.
  merged <- data.frame(merge(master.df, taxa.df, by = "TSN_FINAL",
                             all.y = TRUE))
  # Change the column name "TSN_FINAL" to just "TSN".
  colnames(merged)[colnames(merged) %in% "TSN_FINAL"] <- "TSN"
  # Make all column names are uppercase.
  names(merged) <- toupper(colnames(merged))
  # 2-26-16 removed "TSN" from subset below
  # At this point TSN is represented by TSN_DB which is only needed
  # to link to the BIBI data base.  TSN_R has the most up to date TSN
  # or a unique negative integer assigned
  # Subset of the merged data frames to represent only the columns specific
  # to sample information, taxonomic identification, and taxonomic counts.
  final.df <- merged[, c("EVENT_ID", "SAMPLE_NUMBER", "G_METHOD", "TSN_R",
                         taxa.cols, "REPORTING_VALUE")]
  #============================================================================
  
  final.df[, c("EVENT_ID", taxa.cols)] <- lapply(c("EVENT_ID", taxa.cols),
                                                 function(x){
                                                   clean_char(final.df[, x])
                                                 })
  #============================================================================
  # Convert data.table to class data.frame.
  final.df <- data.frame(final.df)
  # Replace NA's with "UNIDENTIFIED." This will prevent the wide function
  # from return different lengths.
  final.df[,c(taxa.cols)][is.na(final.df[, c(taxa.cols)])] <- "UNIDENTIFIED"
  # Rename the columns.
  names(final.df) <- c("EVENT_ID", "SAMPLE_NUMBER", "G_METHOD", "TSN",
                       taxa.cols, "REPORTING_VALUE")
  # All EVENT_ID's to upper case.  This will make sure all tabs merge properly
  # merged.
  final.df$EVENT_ID <- clean_char(final.df$EVENT_ID)
  #============================================================================
  # Remove and roll taxa up to the appropriate taxonomic level for the 
  # Chessie BIBI development.
  if(clean.taxa == TRUE) final.df <- Chessie::clean_taxa(final.df)
  #============================================================================
  # End prep_taxa function.
  return(final.df)
}


#==============================================================================
#'Prepare Water Quality Data
#'
#'@param WQ = Water quality data in long format
#'@return Transform water quality data from long to wide format
#'@import data.table
#'@export

prep_wq <- function(tab.wq){
  # Convert to data.table to speed up aggregation.
  wq <- data.table::data.table(tab.wq)
  # Aggregate by EVENT_ID and REPORTING_PARAMETER and find the mean 
  # REPORTED_VALUE (data.table style). 
  water.q <- wq[, mean(REPORTED_VALUE, na.rm = TRUE),
                by = list(EVENT_ID, REPORTING_PARAMETER)]
  # Rename columns.
  names(water.q) <- c("EVENT_ID", "REPORTING_PARAMETER",
                      "REPORTED_VALUE")
  # Long to wide data format.
  wide.wq <- tidyr::spread(water.q, REPORTING_PARAMETER, REPORTED_VALUE)
  # Make sure EVENT_ID is of class character, is all uppercase, and has no
  # leading/trailing white space.
  wide.wq$EVENT_ID <- clean_char(wide.wq$EVENT_ID)
  # Convert back to a class data.frame.
  final.df <- data.frame(wide.wq)
  # End prep_wq function.
  return(final.df)
}

#==============================================================================
#'Prepare Habitat Data
#'
#'@param Habitat = Habitat data in long format
#'@param agg.sample.num = Should the sample number be treated as seperate
#'samples (FALSE) or aggregated into a single sample (TRUE)?
#'@return Transform habitat data from long to wide format
#'@import data.table
#'@export

prep_habitat <- function(habitat.df, agg.sample.num = TRUE, bibi.version){
  # Change class to data.table to speed up aggreations below.
  habitat.df <- data.table::data.table(habitat.df)
  #============================================================================
  if (agg.sample.num == FALSE) {
    # If agg.sample.num is set to FALSE.
    # Find the mean reported value by EVENT_ID, SAMPLE_NUMBER, and parameter.
    hab.agg <- habitat.df[, mean(REPORTING_PARAMETER_VALUE, na.rm = TRUE),
                          by = list(EVENT_ID, SAMPLE_NUMBER, 
                                    HABITAT_REPORTING_PARAMETER)]
    # Update column names.
    names(hab.agg) <- c("EVENT_ID", "SAMPLE_NUMBER", "REPORTING_PARAMETER",
                        "REPORTED_VALUE")
    # change format from long to wide.
    hab.wide <- tidyr::spread(hab.agg, REPORTING_PARAMETER, REPORTED_VALUE)
    
  }else{
    # If agg.sample.num is NOT set to FALSE.
    # Find the mean reported value by EVENT_ID and parameter.
    hab.agg <- habitat.df[, mean(REPORTING_PARAMETER_VALUE, na.rm = TRUE),
                          by = list(EVENT_ID, HABITAT_REPORTING_PARAMETER)]
    # Update column names.
    names(hab.agg) <- c("EVENT_ID", "REPORTING_PARAMETER", "REPORTED_VALUE")
    # Create a new SAMPLE_NUMBER column and fill it with ones.
    hab.agg$SAMPLE_NUMBER <- 1
    # change format from long to wide.
    hab.wide <- tidyr::spread(hab.agg, REPORTING_PARAMETER, REPORTED_VALUE)
  }
  #============================================================================
  # Convert data.table to class data.frame.
  final.df <- data.frame(hab.wide)
  # Make sure EVNET_ID is class character, is all caps, and has no 
  # leading/trainling white space.
  final.df$EVENT_ID <- clean_char(final.df$EVENT_ID)
  # Make sure that SAMPLE_NUMBER is class numeric.
  final.df$SAMPLE_NUMBER <- as.numeric(final.df$SAMPLE_NUMBER)
  #============================================================================
  # End prep_habitat function.
  return(final.df)
}

#==============================================================================
#'Merge taxa, water quality, and habitat data into one
#'
#'@param TAB_EVENT = EVENT_ID
#'@param TAXA_PREP = Taxa data in a long data format
#'@param TAB_STATIONS = Station information (STATION_ID, ECOREGION, STRAHLER)
#'@param TAB_WQ = Water quality data in a wide format
#'@param TAB_HABITAT = Habitat data in a wide format
#'@param TAB_PROJECT = Project information
#'@param TAB_HUC = HUC information merged on HUC-12
#'@return Merge taxa water quality, and habitat data into one large data frame.
#'@import data.table
#'@export
prep_merge <- function(tab.event, taxa.prep, tab.stations, wq.wide,
                       habitat.wide, tab.project, tab.huc) {
  #============================================================================
  # Convert all data.frames to data.table to speed up the merge function.
  tab.event <- data.table::data.table(tab.event)
  taxa.prep <- data.table::data.table(taxa.prep)
  tab.stations <- data.table::data.table(tab.stations)
  wq.wide <- data.table::data.table(wq.wide)
  habitat.wide <- data.table::data.table(habitat.wide)
  tab.project <- data.table::data.table(tab.project)
  tab.huc <- data.table::data.table(tab.huc)
  #============================================================================
  # Sampling event lat/long and datum may differ slightly from the station
  # lat/long datum; therfore, a prefix is added to reduce confusion.
  # Additionally, the R_DATE refers to the date the data was recorded into the
  # Chessie BIBI database. Many of the stations were entered into the database
  # during the first Chessie BIBI iteration (2011) and are reused in the latest
  # iteration (2015-2017).  Therefore, R_DATEs will not always match between
  # the event tab and the stations tab. A prefix has been added to reduce 
  # confusion.
  add.prefix <- c("LL_DATUM", "R_DATE")
  names(tab.event)[names(tab.event) %in% add.prefix] <- c("EVENT_LL_DATUM", "EVENT_R_DATE")
  names(tab.stations)[names(tab.stations) %in% add.prefix] <- c("STATION_LL_DATUM", "STATION_R_DATE")
  #============================================================================
  # Join the event info with the taxanomic counts by EVENT_ID.
  merge.1 <- merge(tab.event, taxa.prep, by = "EVENT_ID")
  # Join the previous merge with the station information by STATION_ID.
  # As mentioned above, multiple sampling events have been collected at many
  # sampling locations because replicates were collected or the station was
  # revisted at a differnt time.
  merge.2 <- merge(tab.stations, merge.1, by = "STATION_ID", all.y = TRUE)
  # Join the previous merge with the water quality data in a wide format.
  merge.3 <- merge(wq.wide, merge.2, by = "EVENT_ID", all.y = TRUE)
  #============================================================================
  # HABITAT
  #============================================================================
  # Joining the habitat data requires extra attention.
  # The many of the habitat parameters are specific to sampling event and
  # sample number but not all. In some cases it appears that one set of habitat
  # data was used for multiple samples collected at the same time.
  # (same EVENT_ID but different SAMPLE_NUMBER)
  #----------------------------------------------------------------------------
  # First, all of the habitat data with SAMPLE_NUMBER == 1 are subset.
  tab.hab.1 <- habitat.wide[!habitat.wide$SAMPLE_NUMBER > 1, ]
  # The previous merge(merge.3) is subset to represent only the EVENT_IDs and
  # SAMPLE_NUMBERs that meet these requirements.
  sub.merge.3.1 <- merge.3[merge.3$EVENT_ID %in% tab.hab.1$EVENT_ID & 
                             merge.3$SAMPLE_NUMBER %in% tab.hab.1$SAMPLE_NUMBER, ]
  #----------------------------------------------------------------------------
  # Second, all of the habitat data with a SAMPLE_NUMBER > 1 are subset.
  tab.hab.2 <- habitat.wide[habitat.wide$SAMPLE_NUMBER > 1, ]
  # The previous merge(merge.3) is subset to represent only the EVENT_IDs and
  # SAMPLE_NUMBERs that meet these requirements.
  sub.merge.3.2 <- merge.3[merge.3$EVENT_ID %in% tab.hab.2$EVENT_ID & 
                             merge.3$SAMPLE_NUMBER %in% tab.hab.2$SAMPLE_NUMBER, ]
  #----------------------------------------------------------------------------
  # Finally, the previous merge (merge.3) is subset to represent all of the
  # remaining EVENT_ID and SAMPLE_NUMBER combinations that are not represented
  # by tab.hab.1 or tab.hab.2.
  sub.merge.3.3 <- merge.3[!(merge.3$EVENT_ID %in% tab.hab.2$EVENT_ID & 
                               merge.3$SAMPLE_NUMBER %in% tab.hab.2$SAMPLE_NUMBER) &
                             !(merge.3$EVENT_ID %in% tab.hab.1$EVENT_ID & 
                                 merge.3$SAMPLE_NUMBER %in% tab.hab.1$SAMPLE_NUMBER), ]
  #----------------------------------------------------------------------------
  # Join the habitat data to the respective subseted merge.3 data frames.
  # Two columns frequently used to merge data frames.
  merge.cols <- c("EVENT_ID", "SAMPLE_NUMBER")
  # Join the habitat data and previous merge for sample_number == 1.
  merge.4.1 <- merge(tab.hab.1, sub.merge.3.1, by = merge.cols, all.y = TRUE)
  # Join the habitat data and previous merge for sample_number > 1.
  merge.4.2 <- merge(tab.hab.2, sub.merge.3.2, by = merge.cols, all.y = TRUE)
  # Data.table method for removing column names.
  tab.hab.1[, SAMPLE_NUMBER:=NULL]
  # Join the habitat data with the subset of the previously merged data without
  # a a corresponding set of habitat variables.
  merge.4.3 <- merge(tab.hab.1, sub.merge.3.3, by = c("EVENT_ID"), all.y = TRUE)
  # Join all of the subsets together.
  merge.4 <- rbind(merge.4.1, merge.4.2, merge.4.3)
  #============================================================================
  # PROJECT_IDs to all caps and remove leading/trailing white space.
  merge.4$PROJECT_ID <- clean_char(merge.4$PROJECT_ID)
  tab.project$PROJECT_ID <- clean_char(tab.project$PROJECT_ID)
  # Join the previous merge with project information.
  merge.5 <- merge(tab.project, merge.4, by = "PROJECT_ID", all.y = TRUE)
  #============================================================================
  # HUC_12s to all caps and remove leading/trailing white space.
  merge.5$HUC_12 <- paste0(0,clean_char(merge.5$HUC_12))
  tab.huc$HUC_12 <- paste0(0, clean_char(tab.huc$HUC_12))
  # Join the previous merge with the HUC information.
  final.df <- merge(merge.5, tab.huc, by = "HUC_12", all.x = TRUE)
  #============================================================================
  # End prep_merge function.
  return(final.df)
}

#==============================================================================
#'Merge taxa, water quality, and habitat data into one
#'
#'@param Master = Taxa count data
#'@param odbc.name = The name of the ODBC connection.
#'@param agg_sample_num = Should the sample number be treated as seperate
#'samples (FALSE) or aggregated into a single sample (TRUE)?
#'@param development = during index development, "development" should be
#'set to TRUE.  When TRUE the sample_number (replicate) with the largest reporting
#'value (# of individual observed) is selected for.  During final scoring development
#'should be set to FALSE.  When FALSE all replicates are included in the final
#'output.
#'@param month = If true data from Nov, DEC, and Jan eliminated.
#'@param bibi.version = If set to 2011 the water quality and habitat parameters will be
#'handled as defined in the 2011 Chessie BIBI.  If set to 2016 the water quality
#'and habitat parameters will be handled as in the 2016 Chessie BIBI.
#'@return Merge taxa water quality, and habitat data into one large data frame
#'@import data.table
#'@export

prep_data <- function(master.df, odbc.name = "BIBI_2017",
                      sample.date = "SAMPLE_DATE_TIME",
                      agg.samp.num = TRUE, clean.taxa = TRUE,
                      development = TRUE, bibi.version = 2016, month = TRUE){
  #============================================================================
  print("[1/6] Preparing Data Tabs")
  # Function to make sure key columns are of class character, are all caps, and
  # all leading/trailing white space has been removed.
  clean_tab <- function(x){
    if("EVENT_ID" %in% names(x)){
      x$EVENT_ID <- clean_char(x$EVENT_ID)
    }
    if("STATION_ID" %in% names(x)){
      x$STATION_ID <- clean_char(x$STATION_ID)
    }
    if("AGENCY_CODE" %in% names(x)){
      x$AGENCY_CODE <- clean_char(x$AGENCY_CODE)
    }
    if("PROGRAM_CODE" %in% names(x)){
      x$PROGRAM_CODE <- clean_char(x$PROGRAM_CODE)
    }
    names(x) <- clean_char(names(x))
    return(x)
  }
  #============================================================================
  # Import the data sheets from the MS Access database.
  data.base <- import_msa_data(odbc.name)
  #============================================================================
  # Apply the clean_tab function to make sure key columns are of class
  # character, are all caps, and all leading/trailing white space has been
  # removed.
  tab.taxa <- clean_tab(data.base$tab.taxa)
  tab.wq <- clean_tab(data.base$tab.wq)
  tab.habitat <- clean_tab(data.base$tab.habitat)
  tab.event <- clean_tab(data.base$tab.event)
  tab.stations <- clean_tab(data.base$tab.stations)
  tab.project <- clean_tab(data.base$tab.project)
  tab.huc <- clean_tab(data.base$tab.huc)
  #============================================================================
  # Prepare the taxonomic counts to be joined with the other data sheets.
  taxa.prep <- prep_taxa(master.df, tab.taxa, clean.taxa)
  #============================================================================
  print("[3/6] Preparing Water Quality Data")
  # Prepare the water quality data to be joined witht the other data sheets.
  # The if statements handle the data according to the Chessie BIBI version of
  # interest.
  if(bibi.version == 2016){
    wq.wide <- prep_wq(tab.wq)
    wq.wide <- wq.wide[, c("EVENT_ID", "PH", "SPCOND", "DO", "WTEMP")]
  }
  if(bibi.version == 2011){
    wq.wide <- prep_wq_old(tab.wq)
    wq.wide <- wq.wide[, c("EVENT_ID", "PH", "SPCOND")]
  }
  #============================================================================
  print("[4/6] Preparing tab.habitat Data")
  # Prepare the habitat data to be joined witht the other data sheets.
  # The if statements handle the data according to the Chessie BIBI version of
  # interest.
  if(bibi.version == 2016){
    habitat.wide <- prep_habitat(tab.habitat, agg.samp.num, bibi.version)
  }
  if(bibi.version == 2011){
    habitat.wide <- prep_habitat_old(tab.habitat, agg.samp.num)
  }
  #============================================================================
  # Join all of the data sheets necessary for developing the Chessie BIBI.
  prep.df <- prep_merge(tab.event, taxa.prep, tab.stations, wq.wide, habitat.wide,
                        tab.project, tab.huc)
  # Make sure output is class data.frame.
  prep.df <- data.frame(prep.df)
  # Create column with the just the sampling date.
  prep.df$DATE <- format(as.Date(prep.df[, sample.date]), "%m/%d/%Y")
  #============================================================================
  print("[5/6] Checking for SAMPLE_NUMBER column")
  # Update the SAMPLE_NUMBER column.
  if(agg.samp.num  == FALSE){
    if("SAMPLE_NUMBER" %in% colnames(prep.df)){
      # If there are any NAs in the SAMPLE_NUMBER column, then fill these rows
      # with ones.
      if(sum(is.na(prep.df$SAMPLE_NUMBER)) > 0){
        warning("Warning: Missing SAMPLE_NUMBER values replaced with a 1.\
                Please review the values in the SAMPLE_NUMBER column before proceeding.")
        rep <- prep.df$SAMPLE_NUMBER
        rep[is.na(rep)] <- 1
        prep.df$SAMPLE_NUMBER <- rep
        
      }else{
        # If there are no NAs in the SAMPLE_NUMBER column, then no change occurs.
        prep.df$SAMPLE_NUMBER <- prep.df$SAMPLE_NUMBER
      }
    }else{
      # If the SAMPLE_NUMBER column does not exist, then create one and fill it
      # with ones.
      warning("Warning: Missing a SAMPLE_NUMBER column.\
              A SAMPLE_NUMBER column was created and filled with 1's as a placeholder.\
              Please review the SAMPLE_NUMBER values before proceeding.")
      prep.df$SAMPLE_NUMBER <- 1
    }
  }else{
    # If agg.samp.num is TRUE, then fill the SAMPLE_NUMBER column with ones.
    # When the samples are aggreated during the metric calculation stage the
    # samples that originally would have been kept seperate due to differing 
    # SAMPLE_NUMBERs will now be grouped.
    print("All sample numbers changed to 1, as specified by agg.samp.num = FALSE.
          All samples collected by the same agency at the same site and on the same date
          will be summed to form a composite sample.")
    prep.df$SAMPLE_NUMBER <- 1
  }
  
  #============================================================================
  print("[6/6] Checking for AGENCY_CODE column")
  # Update AGENCY_CODE column.
  if("AGENCY_CODE" %in% colnames(prep.df)){
    # If there are any NAs in the AGENCY_CODE column, then fill these rows with
    # NAs.
    if(sum(is.na(prep.df$AGENCY_CODE)) > 0){
      warning("Warning: Missing Agency Code replaced with a 1.\
              Please review the values in the AGENCY_CODE column before proceeding.")
      AC <- prep.df$AGENCY_CODE
      AC[is.na(AC)] <- 1
      prep.df$AGENCY_CODE <- AC
    }else{
      # If no NAs are present, then concatenate AGENCY_CODE and PROGRAM_CODE.
      prep.df$AGENCY_CODE <- paste(prep.df$AGENCY_CODE, prep.df$PROGRAM_CODE, sep = "_")
    }
  }else{
    warning("Warning: Missing a AGENCY_CODE column.\
            A AGENCY_CODE column was created and filled with 1's as a placeholder.\
            Please review the AGENCY_CODE column before proceeding.")
    # If the column AGENCY_CODE does not exist then, fill with ones.
    prep.df$AGENCY_CODE <- 1
  }
  #============================================================================
  # Final Preperation
  #============================================================================
  # Keep only samples collected with a kick-net or a similar sampling device.
  final.df <- prep_subset(prep.df)
  # Create a column representing the month the sample was collected in.
  final.df$MONTH <- as.numeric(format(as.Date(final.df$DATE, "%m/%d/%Y"), "%m"))
  # Remove Jan, Feb, and Dec becuase these months were infrequently sampled.                                
  if(month == TRUE) final.df <- final.df[!(final.df$MONTH %in% c(1, 2, 12)), ]
  # Create a column representing the julian day the sample was collected on.
  final.df$JULIAN <- lubridate::yday(lubridate::mdy(final.df$DATE))
  #============================================================================
  # If development is equal to TRUE then keep only the sample number with the
  # largest taxonomic count. Averaging taxa counts does not seem appropriate
  # and randomly selecting a sample number adds variability with each run of 
  # this function. The sample with the greatest taxonomic count arguably 
  # provides the best representation of the sampling event. Additionally,
  # there are many instances when sample count in one sample number is very low
  # (< 70), while another sample number for the same event is greater than 100.
  # We argue that samples < 70 could skew percent metrics and should be omitted
  # all together (we also found no significant relationship between degradation
  # and taxonomic counts < 70). Some of these replicates differ greatly, 
  # suggesting that one sample may have been collected to serve a purpose besides
  # reporting a standard taxonomic count (possibly some QA/QC measure?).
  if(development == TRUE){
    # Subset the data frame to include only necessary columns.
    prep.dev <- final.df[, c("EVENT_ID", "SAMPLE_NUMBER", "REPORTING_VALUE")]
    # Sum all of the REPORTING_VALUEs when the values in the remaining 
    # columns (.) match.
    agg.dev.1 <- aggregate(REPORTING_VALUE ~ ., data = prep.dev, sum)
    # Change to class data.table to speed up remaining aggregation and merge.
    agg.dev.1 <- data.table::data.table(agg.dev.1)
    # Find the max value of all of the replicates (sample_number) associated
    # with a single sampling event.
    agg.dev.2 <- agg.dev.1[, max(REPORTING_VALUE, na.rm = TRUE), by = EVENT_ID]
    # Rename the columns.
    names(agg.dev.2) <- c("EVENT_ID", "REPORTING_VALUE")
    # Join the two aggregated data tables by EVENT_ID and REPORTING_VALUE to
    # indicate which SAMPLE_NUMBER should be kept to represent each EVENT_ID.
    agg.final <- merge(agg.dev.1, agg.dev.2,  by = c("EVENT_ID",
                                                     "REPORTING_VALUE"))
    # Remove the REPORTING_VALUE column (data.table style).
    agg.final[, REPORTING_VALUE:=NULL]
    # Keep the selected EVENT_ID/SAMPLE number by joining the representative
    # table with the orgininal final.df.
    final.df <- merge(agg.final, final.df, by = c("EVENT_ID", "SAMPLE_NUMBER"))
  }
  #============================================================================
  # Make sure final.df is class data.frame.
  final.df <- data.frame(final.df)
  #============================================================================
  # Identify duplicated samples.
  dups.df <- identify_duplicates(final.df)
  # If no duplicates are found using the identify_duplicates function, then
  # the function returns NULL.  If the output does not equal NULL, then
  # preform the following modification to the final.df.
  if (!is.null(dups.df)) {
    # Join the identified duplicates with the final data.frame.
    merged.df <- merge(dups.df, final.df,
                       by = c("EVENT_ID", "SAMPLE_NUMBER", "AGENCY_CODE"),
                       all.y = TRUE)
    # Remove any rows that the identify_duplicates functions suggested to be
    # removed. (NOTE = "REMOVE")
    merged.df <- merged.df[!merged.df$NOTE %in% "REMOVE", ]
    # If the NEW_EVENT_ID column from the identify_duplicates function
    # is not NA then updated the EVENT_ID column with the value from the 
    # NEW_EVENT_ID column.
    merged.df$EVENT_ID <- ifelse(!is.na(merged.df$NEW_EVENT_ID), 
                                 merged.df$NEW_EVENT_ID, merged.df$EVENT_ID)
    # If the NEW_SAMPLE_NUMBER column from the identify_duplicates function
    # is not NA then updated the SAMPLE_NUMBER column with the value from the 
    # NEW_SAMPLE_NUMBER column.
    merged.df$SAMPLE_NUMBER <- ifelse(!is.na(merged.df$NEW_SAMPLE_NUMBER), 
                                      merged.df$NEW_SAMPLE_NUMBER,
                                      merged.df$SAMPLE_NUMBER)
    # Remove any traces of the identify_duplicates function.
    final.df <- merged.df[, !names(merged.df) %in% c("DUPLICATE",
                                                     "PCT_DUPLICATED", "NOTE",
                                                     "NEW_EVENT_ID",
                                                     "NEW_SAMPLE_NUMBER")]
  }
  #============================================================================
  # End prep_data function.
  return(final.df)
  }

#==============================================================================
#'Clean Taxa
#'
#'@param taxa.long = a data frame, in a long data format, containing taxonomic
#'hierarchy columns.
#'@return Standardizes taxa according to the Chessie BIBI protocol.
#'@export

clean_taxa <- function(taxa.long){
  # Taxonomic hierarchy column names.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                 "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  #============================================================================
  # If any NAs exist in the taxonomic columns, replace the NAs with
  # "UNIDENTIFIED."
  taxa.long[, taxa.cols] <- apply(taxa.long[, taxa.cols], 2, function(x){
    ifelse(is.na(x), "UNIDENTIFIED", as.character(x))
  } )
  #============================================================================
  # Keep only taxa from the phyla Annelida, Arthropoda, Mollosca,
  # and Platyhelminthes.
  phy.keep <- c("ANNELIDA", "ARTHROPODA", "MOLLUSCA", "PLATYHELMINTHES")
  taxa.long <- taxa.long[taxa.long$PHYLUM %in% phy.keep, ]
  #============================================================================
  # Keep only taxa from the subphyla Clitellata, Crustacea, Hexapoda, 
  # Rhabditophora, and taxa unidentified at this level.
  subphy.keep <- c("CLITELLATA", "CRUSTACEA", "HEXAPODA", "RHABDITOPHORA",
                   "UNIDENTIFIED")
  taxa.long <- taxa.long[taxa.long$SUBPHYLUM %in% subphy.keep, ]
  #============================================================================
  # Remove any taxa from the class Branchiopoda, Maxillopoda, and Ostracoda.
  class.exc <- c("BRANCHIOPODA", "MAXILLOPODA", "OSTRACODA")
  taxa.long <- taxa.long[!(taxa.long$CLASS %in% class.exc), ]
  #============================================================================
  # Remove any taxa from the order Hymenoptera. Aquatic Hymenoptera are 
  # generally small, parasitic organisms. Therefore, they may go easily 
  # unnoticed during sorting. We decided it was best to remove these organisms
  # from the analysis.
  order.exc <- c("HYMENOPTERA")
  taxa.long <- taxa.long[!(taxa.long$ORDER  %in% order.exc), ]
  #============================================================================
  # Remove any taxa from the families Gerridae, Hebridae, Veliidae, 
  # Hydrometridae, and Saldidae.  These taxa are generally considerd 
  # semi-aquatic because they live on the surface of the water. Therefore,
  # it is not appropriate to include these organisims in a benthic IBI.
  family.exc <- c("GERRIDAE", "HEBRIDAE", "VELIIDAE", "HYDROMETRIDAE",
                  "SALDIDAE")
  taxa.long <- taxa.long[!(taxa.long$FAMILY  %in% family.exc), ]
  #============================================================================
  # Remove any taxa from the genus Stenus. These taxa are considered 
  # semi-aquatic, and thus, were removed from the analysis.
  genus.exc <- c("STENUS")
  taxa.long <- taxa.long[!(taxa.long$GENUS  %in% genus.exc), ]
  #============================================================================
  # These taxa were not consitently identified to the same taxonomic rank
  # by the agencies that contributed data. Therefore, the samples were rolled
  # up to the lowest common denominator.
  # These columns will be influenced by the common denominator taxa.
  class.spp <- c("CLASS", "SUBCLASS", "ORDER", "SUBORDER",
                 "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  taxa.long[taxa.long$CLASS %in% "BIVALVIA", class.spp] <- "BIVALVIA"
  taxa.long[taxa.long$CLASS %in% "GASTROPODA", class.spp] <- "GASTROPODA"
  taxa.long[taxa.long$CLASS %in% "OLIGOCHAETA", class.spp] <- "OLIGOCHAETA"
  taxa.long[taxa.long$CLASS %in% "TREPAXONEMATA", class.spp] <- "TREPAXONEMATA"
  #============================================================================
  # These taxa were not consitently identified to the same taxonomic rank
  # by the agencies that contributed data. Therefore, the samples were rolled
  # up to the lowest common denominator.
  # These columns will be influenced by the common denominator taxa.
  order.spp <- c("ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  taxa.long[taxa.long$ORDER %in% "COLLEMBOLA", order.spp] <- "COLLEMBOLA"
  taxa.long[taxa.long$ORDER %in% "LEPIDOPTERA", order.spp] <- "LEPIDOPTERA"
  taxa.long[taxa.long$ORDER %in% "NEUROPTERA", order.spp] <- "NEUROPTERA"
  taxa.long[taxa.long$ORDER %in% "NEOOPHORA", order.spp] <- "NEOOPHORA"
  #============================================================================
  # End clean_taxa function.
  return(taxa.long)
}

#==============================================================================
#'Bioregion aggregation
#'
#'@param bioregion.df = Bioregion data frame to be aggregated
#'@return Prepares the bioregion data for metric calculations.
#'@export

bioregion_agg <- function(bioregion.df){
  # Columns to keep.
  keep.cols <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                 "AGENCY_CODE", "TSN", "PHYLUM", "SUBPHYLUM", "CLASS",
                 "SUBCLASS", "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES", "REPORTING_VALUE")
  # Subset the data frame to only include the specified columns.
  bioregion.taxa <- bioregion.df[, keep.cols]
  
  bioregion.taxa$FAMILY <- toupper(bioregion.taxa$FAMILY)
  
  final.df <- aggregate(REPORTING_VALUE ~ EVENT_ID + STATION_ID + DATE +
                          SAMPLE_NUMBER + AGENCY_CODE +TSN +
                          PHYLUM + SUBPHYLUM + CLASS + SUBCLASS + ORDER +
                          SUBORDER + FAMILY + SUBFAMILY + TRIBE + GENUS +
                          SPECIES, data = bioregion.taxa, FUN = sum,
                        na.rm = TRUE)
  colnames(final.df) <- keep.cols
  final.df[final.df == ""] <- NA
  # End bioregion_agg function.
  return(final.df)
}

#==============================================================================
#'Prepare subset specific to BIBI
#'
#'@param long.df = Data in a long data format.
#'@return A subset of the data containing sampling events with strahler stream
#'order less than or equal to 4 and samples collected with a kick net or 
#'another method considered comparable.
#'@export

prep_subset <- function(long.df){
  # Remove samples collected in streams with a strahler order greater than 4.
  # Remove samples that were not collected with a kick net or similar 
  # sampling method.
  # 57 = D-Frame Net (500 micron mesh; 12 inch diameter)
  # 58 = Rectangular Dip Net (0.5 x 0.5 m)
  # 86 = Kick Net (0.8 x 0.9mm mesh; 23 x 46cm size)
  # 87 = Kick Net (unspecified)
  # 89 = D-Frame Net (unspecified)
  # 92 = Kick Seine
  # 93 = D-Frame Net (600 micron; 12 inch diameter)
  # 94 = Kick Net (600 micron; 1 square meter kick screen)
  # 101 = D-frame Net (800-900 micron; 12 inch diameter)
  # 103 = Slack Sampler (500 micron; 1.25 square meter sampling area)
  # 104 = Rectangular Aquatic Net (0.8 x 0.9mm mesh; 9 x 18 inch size)
  final.df <- subset(long.df,  long.df$STRAHLER_STREAM_ORDER <= 4 &
                       long.df$G_METHOD %in% c(57, 58, 86, 87, 89, 92,
                                               93, 94, 101, 103, 104))
  # End prep_subset function.
  return(final.df)
}

#==============================================================================
#'Identify Duplicate Samples
#'
#'@param long.df = Data in a long data format.
#'@return Identify duplicate samples and suggest how the samples should be
#'handled.
#'@export


identify_duplicates <- function(long.df){
  # Extract only site related info.
  site.info <- unique(long.df[, c("EVENT_ID", "STATION_ID", "AGENCY_CODE",
                                  "SAMPLE_NUMBER", "DATE")])
  # Identify rows with duplicated values for all of the specified columns.
  dups.df <- site.info[duplicated(site.info[, c("STATION_ID", "AGENCY_CODE",
                                                "SAMPLE_NUMBER", "DATE")]), ]
  # Remove the EVENT_ID column.
  dups.df <- unique(dups.df[, !names(dups.df) %in% "EVENT_ID"])
  
  if (nrow(dups.df) == 0) return(NULL)
  #dups.df <- dups.df[dups.df$STATION_ID %in% "2-SMN002.19",]
  # Create an empty list to store each iteration from the loop.
  #event.list <- list()
  #============================================================================
  # Loop through dups.df by each row, using the STATION_ID, AGENCY_CODE,
  # SAMPLE_NUMBER, and DATE columns to subset the long.df. The subset is first
  # tested for multiple entry dates (2011 vs 2015). Second, if the sample is
  # 100% match to another sample. If the samples are different but collected
  # from the same station on the same day, then the samples are given unique
  # SAMPLE_NUMBERs.
  event.list <- lapply(1:nrow(dups.df), function(i){
    # Print out the current row being used in the loop.
    # Helps to track progress and identify which row is the issue if an
    # error occurs.
    print(paste(i, "/", nrow(dups.df)))
    # Join the duplicated row with the long.df to subset the long.df to 
    # represent the rows that have the potential of being duplicates.
    sub.df <- merge(dups.df[i, ], long.df, by = c("STATION_ID", "AGENCY_CODE",
                                                  "SAMPLE_NUMBER", "DATE"))
    #============================================================================
    # Identify the unique EVENT_R_DATE (the date the data was recorded into
    # the Chessie BIBI database).
    uniq.rec.date <- unique(sub.df$EVENT_R_DATE)
    # If there is more than one unique EVENT_R_DATE, then select the
    # sample that was entered most recently. It appears that several
    # agencies have modified/updated taxonomic counts since the firest
    # iteration of the Chessie BIBI. These samples do not appear as true 
    # duplicates because they do not produce a 100% match. However, we can be
    # confident that if a sample recorded into the database in 2011 was
    # represented by one replicate (SAMPLE_NUMBER) and again in 2015-2017
    # the sample has the same STATION_ID, AGENCY_CODE, SAMPLE_NUMBER, and DATE,
    # then the newest submission represents the most up-to-date version of the
    # sample. Therefore, the older submission is elminated.
    if (length(uniq.rec.date) > 1) {
      # Keep only the necessary columns.
      final.df <- unique(sub.df[, c("EVENT_ID", "SAMPLE_NUMBER",
                                    "AGENCY_CODE", "EVENT_R_DATE")])
      # The majority of these columns are not applicable to this if statement
      # but they are used in other sections of the script if there is only
      # one uniq.rec.date. Therefore, these columns are used as place holders
      # and filled with NA.
      final.df[, c("DUPLICATE", "PCT_DUPLICATED", "NOTE",
                   "NEW_EVENT_ID", "NEW_SAMPLE_NUMBER")] <- NA
      # Identify the newest submission to the Chessie BIBI database.
      newest.sub <- max(final.df$EVENT_R_DATE)
      # If the EVENT_R_DATE matches the newst.sub, then keep that event_id
      # otherwise remove the event_id.
      final.df$NOTE <- ifelse(final.df$EVENT_R_DATE %in% newest.sub,
                              "KEEP", "REMOVE")
      # Remove the EVENT_R_DATE column.
      final.df <- final.df[, !names(final.df) %in% "EVENT_R_DATE"]
      # Identify the unique EVENT_IDs.
      events.vec <- as.character(unique(final.df$EVENT_ID))
      # Specify which EVENT_ID(s) are duplicated with the EVENT_ID of interest.
      final.df$DUPLICATE <- unlist(lapply(final.df$EVENT_ID, function(x){
        paste(events.vec[!events.vec %in% x], collapse = ", ")
      }))
      # Add the data.frame to the list.
      #event.list[[i]] <- final.df
      return(final.df)
      # Stop the current iteration and move to the next iteration of the loop.
      #next
    }
    #============================================================================
    # Taxonomic columns.
    taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
                   "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS",
                   "SPECIES") 
    # If any of the taxonomic columns contain "UNIDENTIFIED", replace those
    # strings with NA.
    sub.df[, taxa.cols] <- apply(sub.df[, taxa.cols], 2, function(x){
      x <- ifelse(x %in% "UNIDENTIFIED", NA, x)
    })
    # Fill the NAs with the previous taxonomic rank.
    fill.df <- BIBI::fill_taxa(sub.df)
    #============================================================================
    # KEEP???????????????????????????????????????????????????????????????????????
    # For each EVENT_ID identify if there are 
    for (j in unique(sub.df$EVENT_ID)) {
      sub.event <- fill.df[fill.df$EVENT_ID %in% j, ]
      dups.1 <- sub.event[duplicated(sub.event[, c("SPECIES", "REPORTING_VALUE")]), ]
      check.1 <- merge(dups.1, sub.event, by = c("SPECIES", "REPORTING_VALUE"))
      pct.dups.1 <- (nrow(check.1) / nrow(sub.event)) * 100
    }
    #============================================================================
    # Sum the taxonomic counts by sampling event, replicate number, and final
    # taxonomic id.
    agg.df <- aggregate(REPORTING_VALUE ~ EVENT_ID + SAMPLE_NUMBER + SPECIES,
                        data = fill.df, sum)
    # Vector of habitat and water quality variable column names.
    # Environmental parameters.
    env.cols <- c('AESTH', 'BANKS', 'BANKV', 'CH_ALT', 'COVER', 'EMBED',
                  'EPI_SUB', 'FLOW', 'GRAZE', 'INSTR', 'P_SUB', 'POOL',
                  'REMOT', 'RIFF', 'ROOT', 'SED', 'SHAD', 'SINU', 'THAL',
                  'VEL_D', 'WOOD', 'WWID',
                  #"HAB_HETERO", 'INSTR_COND', 'RIP_ZONE',
                  'PH', 'SPCOND', 'DO', 'WTEMP')
    # Remove unnecessary columns.
    non.taxa.df <- unique(fill.df[, c("EVENT_ID", "STATION_ID", "SAMPLE_NUMBER",
                                      "AGENCY_CODE", "DATE", env.cols)])
    # Join the sample information with aggregated taxonomic counts.
    fill.df <- merge(non.taxa.df, agg.df,
                     by = c("EVENT_ID", "SAMPLE_NUMBER"))
    # Identify duplicated final taxonomic ids and associated taxonomic counts.
    dups.2 <- fill.df[duplicated(fill.df[, c("SPECIES", "REPORTING_VALUE")]), ]
    # Join the duplicated columns with the taxonomic counts to subset the
    # taxonomic counts (fill.df) to only represent duplicated rows.
    check.2 <- merge(dups.2, fill.df, by = c("SPECIES", "REPORTING_VALUE"))
    # Identify the percentage of rows that represent duplicates.
    pct.dups <- (nrow(check.2) / nrow(fill.df)) * 100
    #============================================================================
    # Create a new data.frame to store the output.
    final.df <- unique(fill.df[, c("EVENT_ID", "SAMPLE_NUMBER", "AGENCY_CODE")])
    # Identify the unique EVENT_IDs.
    events.vec <- as.character(unique(final.df$EVENT_ID))
    # Specify which EVENT_ID(s) are duplicated with the EVENT_ID of interest.
    final.df$DUPLICATE <- unlist(lapply(final.df$EVENT_ID, function(x){
      paste(events.vec[!events.vec %in% x], collapse = ", ")
    }))
    #============================================================================
    # Specify the percentage of rows that are duplicated between the two EVENT_IDs.
    final.df$PCT_DUPLICATED <- pct.dups
    # If the %duplicated is greater than or equal to fifty but less than 
    # one hundred, indicate that the sample should be revied. If the samples are
    # one hundred percent duplicated, then indicate that both samples should be
    # removed. Later in the script one of the samples will be selected to keep.
    # If the samples have less than fifty percent of the rows representing
    # duplicated taxa/taxonomic counts, then keep both samples.
    final.df$NOTE <- ifelse(final.df$PCT_DUPLICATED >= 50 & final.df$PCT_DUPLICATED < 100,
                            "REVIEW", 
                            ifelse(final.df$PCT_DUPLICATED == 100, "REMOVE",
                                   ifelse(final.df$PCT_DUPLICATED < 50, "KEEP",
                                          "ERROR")))
    # Create place holders for two new columns
    final.df[, c("NEW_EVENT_ID", "NEW_SAMPLE_NUMBER")] <- NA
    #============================================================================
    if(pct.dups == 100){
      #------------------------------------------------------------------------
      # If the percent of duplicated samples is equal to 100, then identify if
      # one EVENT_ID has more habitat and water quality variables. Samples with
      # more reported habitat and water quality variables should be favored
      # because these variables are used during the disturbance gradient 
      # classification process during IBI development.
      #------------------------------------------------------------------------
      # Keep only unique rows with EVENT_ID and the environmental parameters.
      # env.cols specified above.
      env.df <- unique(fill.df[, c("EVENT_ID", env.cols)])
      # Count the number of columns filled with NAs for each EVENT_ID.
      env.df$na.count <- apply(env.df[, env.cols], 1, function(x) length(x[is.na(x)]))
      # Identify which EVENT_ID to keep.
      if (min(env.df$na.count) == max(env.df$na.count)){
        # If both samples have the same environmental parameters reported,
        # then choose the EVENT_ID that was entered into the database first
        # (i.e., minimum EVENT_ID).
        env.df$NOTE <- ifelse(env.df$EVENT_ID %in% min(env.df$EVENT_ID), "KEEP",
                              "REMOVE")
      } else {
        # If the number of environmental parameters columns filled with NAs
        # differ between the EVENT_IDs, then select the EVENT_ID with the
        # fewest NAs (i.e., minimum number of NA counts).
        env.df$NOTE <- ifelse(env.df$na.count == min(env.df$na.count), "KEEP", 
                              ifelse(env.df$na.count > min(env.df$na.count), "REMOVE",
                                     "ERROR"))
        # Remove the "NOTE" column from the final.df so that it is not
        # duplicated by the following merge.
        final.df <- final.df[, !names(final.df) %in% "NOTE"]
        # Join the final.df with env.df, which specifies which EVENT_IDs to
        # keep.
        final.df <- merge(final.df, env.df[, c("EVENT_ID", "NOTE")], by = "EVENT_ID")
      }
    } else {
      #------------------------------------------------------------------------
      # If the EVENT_IDs are not a 100% match, then assign the first EVENT_ID
      # entered into the database (minimum EVENT_ID) to both samples and
      # assign unique SAMPLE_NUMBERs to both samples.
      #------------------------------------------------------------------------
      # Identify the first EVENT_ID created and entered into the database.
      final.df$NEW_EVENT_ID <- min(final.df$EVENT_ID)
      # Make sure the SAMPLE_NUMBER column is class numeric.
      final.df$SAMPLE_NUMBER <- as.numeric(final.df$SAMPLE_NUMBER)
      samp.num.i <- unique(final.df$SAMPLE_NUMBER)
      samp.num.df <- unique(long.df[long.df$EVENT_ID %in% final.df$EVENT_ID,
                                    c("EVENT_ID", "SAMPLE_NUMBER")])
      max.samp.num <- max(samp.num.df$SAMPLE_NUMBER)
      
      if (samp.num.i == max.samp.num) {
        # Identify the SAMPLE_NUMBER reported for both samples.
        # Create a vector that assigns a unique integer to each sample by
        # specifying a range from the reported SAMPLE_NUMBER to the number of 
        # additional samples being represented.
        samp.num.fill <- min(final.df$SAMPLE_NUMBER):(min(final.df$SAMPLE_NUMBER) + (nrow(final.df) - 1)) 
        # Add the new SAMPLE_NUMBERs.
        final.df$NEW_SAMPLE_NUMBER <- samp.num.fill
      } else {
        print("LOOK AT ME")
        # If the sample number for the samples is not equal to the max sample
        # number for any samples with the EVENT_IDs being evaluated, then
        # we do not want to create more issue by duplicating the sample number.
        # For example if our duplicates in question have a sample number of 1,
        # then the code above would assign sample numbers 1:2. However, if an
        # additional sample (e.g., SAMPLE_NUMBER = 2) with the same EVENT_ID
        # as the samples being evaluated already exits, then we would duplicate
        # the SAMPLE_NUMBER and these two samples would be incorrectly aggregated
        # during later development stages of the Chessie BIBI. To avoid this
        # I added one to the max sample number of any of the EVENT_IDs being 
        # ivestigated during a given iteration of the loop as the first new
        # SAMPLE_NUMBER. I then assign subsquent integers depending on how many
        # unique EVENT_IDs exist.
        samp.num.fill <- max.samp.num + 1:(max.samp.num + (nrow(final.df) - 1)) 
        # Add the new SAMPLE_NUMBERs.
        final.df$NEW_SAMPLE_NUMBER <- samp.num.fill
      }
      
    }
    # Store the data.frame in the list.
    #event.list[[i]] <- final.df
    return(final.df)
  }
  )
  #============================================================================
  # Join all of the data.frames in the list by appending them to one another.
  final.df <- do.call(rbind, event.list)
  #============================================================================
  # End identify_duplicates function.
  return(final.df)
}

