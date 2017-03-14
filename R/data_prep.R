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
#'Prepare Taxonomic Counts
#'
#'@param Master = Taxa count data
#'@param Taxa = Taxonomic counts in long data format
#'@return Merges a taxonomic counts data frame with a master taxa list
#'containing the following information for each taxon when applicable:
#'TSN, Phylum, Subphylum, Class, Subclass, Order, Suborder, Family,
#'Subfamily, Tribe, Genus, and Species.
#'@export

prep_taxa <- function(Master, Taxa){
  #Taxa$TSN <- as.numeric(Taxa$TSN)
  #merged <- merge(Master, Taxa, by.x = "TSN_DB", by.y = "TSN", all.y = TRUE)
  names(Taxa) <- toupper(names(Taxa))

  merged <- merge(Master, Taxa, by = "TSN_FINAL", all.y = TRUE)
  colnames(merged)[colnames(merged) == "TSN_FINAL"] <- "TSN"
  
  names(merged) <- toupper(colnames(merged))
  # 2-26-16 removed "TSN" from subset below
  # At this point TSN is represented by TSN_DB which is only needed
  # to link to the BIBI data base.  TSN_R has the most up to date TSN
  # or a unique negative integer assigned
  p_taxa <- merged[, c("EVENT_ID", "SAMPLE_NUMBER", "G_METHOD", "TSN_R",
                       "PHYLUM", "SUBPHYLUM", "CLASS","SUBCLASS", "ORDER",
                       "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS",
                       "SPECIES", "REPORTING_VALUE")]

  p_taxa[,"PHYLUM"] <- toupper(p_taxa[,"PHYLUM"])
  p_taxa[,"SUBPHYLUM"] <- toupper(p_taxa[,"SUBPHYLUM"])
  p_taxa[,"CLASS"] <- toupper(p_taxa[,"CLASS"])
  p_taxa[,"SUBCLASS"] <- toupper(p_taxa[,"SUBCLASS"])
  p_taxa[,"ORDER"] <- toupper(p_taxa[,"ORDER"])
  p_taxa[,"SUBORDER"] <- toupper(p_taxa[,"SUBORDER"])
  p_taxa[,"FAMILY"] <- toupper(p_taxa[,"FAMILY"])
  p_taxa[,"SUBFAMILY"] <- toupper(p_taxa[,"SUBFAMILY"])
  p_taxa[,"TRIBE"] <- toupper(p_taxa[,"TRIBE"])
  p_taxa[,"GENUS"] <- toupper(p_taxa[,"GENUS"])
  p_taxa[,"SPECIES"] <- toupper(p_taxa[,"SPECIES"])

  p_taxa[p_taxa == ""] <- NA
  # Replace NA's with "UNIDENTIFIED." This will prevent the wide function
  # from return different lengths.
  p_taxa[,c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")][is.na(p_taxa[,
          c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
          "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")])] <- "UNIDENTIFIED"

  # agg_taxa <- aggregate(REPORTING_VALUE ~ addNA(EVENT_ID) +
  #                      addNA(SAMPLE_NUMBER) + addNA(G_METHOD) +
  #                      addNA(TSN_R) + addNA(PHYLUM) + addNA(SUBPHYLUM) +
  #                      addNA(CLASS) + addNA(SUBCLASS) + addNA(ORDER) +
  #                      addNA(SUBORDER) + addNA(FAMILY) + addNA(SUBFAMILY) +
  #                      addNA(TRIBE) + addNA(GENUS) + addNA(SPECIES),
  #                      FUN = sum, data = p_taxa)
  names(p_taxa) <- c("EVENT_ID", "SAMPLE_NUMBER", "G_METHOD", "TSN",
                     "PHYLUM", "SUBPHYLUM", "CLASS","SUBCLASS", "ORDER",
                     "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS",
                     "SPECIES", "REPORTING_VALUE")
  # All EVENT_ID's to upper case.  This will make sure all tabs merge properly
  p_taxa$EVENT_ID <- toupper(p_taxa$EVENT_ID)
  
  return(p_taxa)
}


#==============================================================================
#'Prepare Water Quality Data
#'
#'@param WQ = Water quality data in long format
#'@return Transform water quality data from long to wide format
#'@export

prep_wq <- function(WQ){
  
  water.q <- aggregate(WQ$REPORTED_VALUE ~ WQ$EVENT_ID + WQ$REPORTING_PARAMETER,
                       FUN = mean, na.rm = TRUE)
  colnames(water.q) <- c("EVENT_ID", "REPORTING_PARAMETER", "REPORTED_VALUE")

  wide.wq <- reshape2::dcast(water.q, EVENT_ID ~ REPORTING_PARAMETER,
                             value.var = "REPORTED_VALUE")
  wide.wq$EVENT_ID <- toupper(wide.wq$EVENT_ID)
  return(wide.wq)
}

#==============================================================================
#'Prepare Habitat Data
#'
#'@param Habitat = Habitat data in long format
#'@param agg_sample_num = Should the sample number be treated as seperate
#'samples (FALSE) or aggregated into a single sample (TRUE)?
#'@return Transform habitat data from long to wide format
#'@export
prep_habitat <- function(Habitat, agg_sample_num = TRUE){
  wide.habitat <- if(agg_sample_num == FALSE){
    hab <- aggregate(data = Habitat, REPORTING_PARAMETER_VALUE ~ EVENT_ID +
                       SAMPLE_NUMBER +
                       HABITAT_REPORTING_PARAMETER,
                     FUN = mean, na.rm = TRUE)
    names(hab) <- c("EVENT_ID", "SAMPLE_NUMBER", "REPORTING_PARAMETER", "REPORTED_VALUE")
    wide.habitat <- reshape2::dcast(hab, EVENT_ID + SAMPLE_NUMBER ~ REPORTING_PARAMETER,
                                    value.var = "REPORTED_VALUE")
  }else{
    hab <- aggregate(data = Habitat, REPORTING_PARAMETER_VALUE ~ EVENT_ID +
                       HABITAT_REPORTING_PARAMETER, FUN = mean, na.rm = TRUE)
    names(hab) <- c("EVENT_ID", "REPORTING_PARAMETER", "REPORTED_VALUE")
    hab$SAMPLE_NUMBER <- 1
    wide.habitat <- data.frame(reshape2::dcast(hab, EVENT_ID + SAMPLE_NUMBER
                                               ~ REPORTING_PARAMETER,
                                    value.var = "REPORTED_VALUE"))
  }

  #wide2_Habitat<-wide_Habitat[rowSums(wide_Habitat[,2:ncol(wide_Habitat)]) != 0,]
  habitat.data <- data.frame(wide.habitat)
  habitat.data$HAB_HETERO <- ifelse (! is.na(habitat.data$RIFF) & ! is.na(habitat.data$POOL),
                                     (habitat.data$RIFF + habitat.data$POOL) / 2,
                                     ifelse (is.na(habitat.data$RIFF) & ! is.na(habitat.data$POOL),
                                             habitat.data$POOL,
                                             ifelse (! is.na(habitat.data$RIFF) & is.na(habitat.data$POOL),
                                                     habitat.data$RIFF, NA)))
  habitat.data$INSTR_COND <- ifelse (! is.na(habitat.data$EPI_SUB) & ! is.na(habitat.data$COVER),
                                     (habitat.data$EPI_SUB + habitat.data$COVER) / 2,
                                     ifelse (is.na(habitat.data$EPI_SUB) & ! is.na(habitat.data$COVER),
                                             habitat.data$COVER,
                                             ifelse (! is.na(habitat.data$EPI_SUB) & is.na(habitat.data$COVER),
                                                     habitat.data$EPI_SUB, NA)))
  #habitat.data$RIP_ZONE <- habitat.data$RIP_SC

  habitat.data$RIP_ZONE <- ifelse(!is.na(habitat.data$RIP_SC) & !is.na(habitat.data$RIP_W),
                                (habitat.data$RIP_SC + habitat.data$RIP_W) / 2,
                                       ifelse(!is.na(habitat.data$RIP_SC)& is.na(habitat.data$RIP_W),
                                              habitat.data$RIP_SC,
                                              ifelse(is.na(habitat.data$RIP_SC)& !is.na(habitat.data$RIP_W),
                                                     habitat.data$RIP_W, NA)))
  habitat.data <- habitat.data[, !(names(habitat.data) %in% c("RIP_SC", "RIP_W"))]


  #env.data$SUM <- rowSums(env.data[, 2:ncol(env.data)], na.rm = TRUE)
  
  habitat.data$EVENT_ID <- toupper(habitat.data$EVENT_ID)
  habitat.data$SAMPLE_NUMBER <- toupper(habitat.data$SAMPLE_NUMBER)
  
  #final.df <- aggregate(.~ EVENT_ID + SAMPLE_NUMBER, data = habitat.data, mean)
  #test <- habitat.data[duplicated(habitat.data[, c("EVENT_ID", "SAMPLE_NUMBER")]),]
  
  return(habitat.data)
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
#'@export
prep_merge <- function(TAB_EVENT, TAXA_PREP, TAB_STATIONS, TAB_WQ,
                       TAB_HABITAT, TAB_PROJECT, TAB_HUC) {
  merg3 <- merge(TAB_EVENT, TAXA_PREP, by = "EVENT_ID")
  merg4 <- merge(TAB_STATIONS, merg3, by = "STATION_ID", all.y = TRUE)
  merg5 <- merge(TAB_WQ, merg4, by = "EVENT_ID", all.y = TRUE)
  #names(TAB_HABITAT)[names(TAB_HABITAT) %in% "SAMPLE_NUMBER"] <- "HAB_SAMPLE_NUMBER"
  #============================================================================
  
  TAB_HABITAT_1 <- TAB_HABITAT[!TAB_HABITAT$SAMPLE_NUMBER > 1, ]
  merg5_1 <- merg5[merg5$EVENT_ID %in% TAB_HABITAT_1$EVENT_ID & 
                     merg5$SAMPLE_NUMBER %in% TAB_HABITAT_1$SAMPLE_NUMBER, ]
  TAB_HABITAT_2 <- TAB_HABITAT[TAB_HABITAT$SAMPLE_NUMBER > 1, ]
  merg5_2 <- merg5[merg5$EVENT_ID %in% TAB_HABITAT_2$EVENT_ID & 
                       merg5$SAMPLE_NUMBER %in% TAB_HABITAT_2$SAMPLE_NUMBER, ]
  merg5_3 <- merg5[!(merg5$EVENT_ID %in% TAB_HABITAT_2$EVENT_ID & 
                     merg5$SAMPLE_NUMBER %in% TAB_HABITAT_2$SAMPLE_NUMBER) &
                     !(merg5$EVENT_ID %in% TAB_HABITAT_1$EVENT_ID & 
                         merg5$SAMPLE_NUMBER %in% TAB_HABITAT_1$SAMPLE_NUMBER), ]
 
  
  merg6.1 <- merge(TAB_HABITAT_1, merg5_1, by = c("EVENT_ID"), all.y = TRUE)
  merg6.1 <- merg6.1[, !names(merg6.1) %in% "SAMPLE_NUMBER.x"]
  names(merg6.1)[names(merg6.1) %in% "SAMPLE_NUMBER.y"] <- "SAMPLE_NUMBER"
  merg6.2 <- merge(TAB_HABITAT_2, merg5_2, by = c("EVENT_ID", "SAMPLE_NUMBER"), all.y = TRUE)
  merg6.3 <- merge(TAB_HABITAT_1, merg5_3, by = c("EVENT_ID"), all.y = TRUE)
  merg6.3 <- merg6.3[, !names(merg6.3) %in% "SAMPLE_NUMBER.x"]
  names(merg6.3)[names(merg6.3) %in% "SAMPLE_NUMBER.y"] <- "SAMPLE_NUMBER"
  merg6 <- rbind(merg6.1, merg6.2, merg6.3)
  
  #nrow(unique(tab_taxa[, c("EVENT_ID", "SAMPLE_NUMBER")]))
  #nrow(unique(merg6[, c("EVENT_ID", "SAMPLE_NUMBER")]))
  
  #merg6 <- merge(TAB_HABITAT, merg5, by = c("EVENT_ID"), all.y = TRUE)
  #============================================================================
  
  merg6$PROJECT_ID <- toupper(merg6$PROJECT_ID)
  TAB_PROJECT$PROJECT_ID <- toupper(TAB_PROJECT$PROJECT_ID)
  merg7 <- merge(TAB_PROJECT, merg6, by = "PROJECT_ID", all.y = TRUE)
  merg7$HUC_12 <- toupper(merg7$HUC_12)
  TAB_HUC$HUC_12 <- toupper(TAB_HUC$HUC_12)
  merg8 <- merge(merg7, TAB_HUC, by = "HUC_12", all.x = TRUE)
  return(merg8)
}

#==============================================================================
#'Merge taxa, water quality, and habitat data into one
#'
#'@param Master = Taxa count data
#'@param Taxa = Taxonomic counts in long data format
#'@param Taxonomy = TSN #'s with necessary taxonomic levels (CLASS, ORDER, FAMILY)
#'@param WQ = Water quality data in a long data format
#'@param Habitat = Habitat data in a long data format
#'@param TAB_EVENT = EVENT_ID
#'@param TAB_STATIONS = Station information (STATION_ID, ECOREGION, STRAHLER)
#'@param agg_sample_num = Should the sample number be treated as seperate
#'samples (FALSE) or aggregated into a single sample (TRUE)?
#'@param development = during index development, "development" should be
#'set to TRUE.  When TRUE the sample_number (replicate) with the largest reporting
#'value (# of individual observed) is selected for.  During final scoring development
#'should be set to FALSE.  When FALSE all replicates are included in the final
#'output.
#'@param month = If true data from Nov, DEC, and Jan eliminated.
#'@param bibi = If set to 2011 the water quality and habitat parameters will be
#'handled as defined in the 2011 Chessie BIBI.  If set to 2016 the water quality
#'and habitat parameters will be handled as in the 2016 Chessie BIBI.
#'@return Merge taxa water quality, and habitat data into one large data frame
#'@export

prep_data <- function(Master, Taxa, WQ, Habitat, TAB_EVENT, TAB_STATIONS,
                      TAB_PROJECT, TAB_HUC,
                      Sample_Date = "SAMPLE_DATE_TIME",
                      agg_sample_num = TRUE, clean_taxa = TRUE,
                      development = TRUE, bibi = 2016, month = TRUE){
  
  print("[1/6] Preparing Data Tabs")
  upper.case <- function(x){
    if("EVENT_ID" %in% names(x)){
      x$EVENT_ID <- lapply(x$EVENT_ID, as.character)
      x$EVENT_ID <- toupper(x$EVENT_ID)
      x$EVENT_ID <- gsub("^\\s+|\\s+$", "", x$EVENT_ID)
    }
    if("STATION_ID" %in% names(x)){
      x$STATION_ID <- lapply(x$STATION_ID, as.character)
      x$STATION_ID <- toupper(x$STATION_ID)
      x$STATION_ID <- gsub("^\\s+|\\s+$", "", x$STATION_ID)
    }
    if("AGENCY_CODE" %in% names(x)){
      x$AGENCY_CODE <- lapply(x$AGENCY_CODE, as.character)
      x$AGENCY_CODE <- toupper(x$AGENCY_CODE)
      x$AGENCY_CODE <- gsub("^\\s+|\\s+$", "", x$AGENCY_CODE)
    }
    if("PROGRAM_CODE" %in% names(x)){
      x$PROGRAM_CODE <- lapply(x$PROGRAM_CODE, as.character)
      x$PROGRAM_CODE <- toupper(x$PROGRAM_CODE)
      x$PROGRAM_CODE <- gsub("^\\s+|\\s+$", "", x$PROGRAM_CODE)
    }
    return(x)
  } 
  
  Taxa <- upper.case(Taxa)
  WQ <- upper.case(WQ)
  Habitat <- upper.case(Habitat)
  TAB_EVENT <- upper.case(TAB_EVENT)
  TAB_STATIONS <- upper.case(TAB_STATIONS)
  TAB_PROJECT <- upper.case(TAB_PROJECT)
  
  print("[2/6] Preparing Taxanomic Counts")
  PT <- prep_taxa(Master, Taxa)

  print("[3/6] Preparing Water Quality Data")
  if(bibi == 2016){
    wide.wq <- prep_wq(WQ)
    wide.wq <- wide.wq[, c("EVENT_ID", "PH", "SPCOND", "DO", "WTEMP")]
  }
  if(bibi == 2011){
    wide.wq <- prep_wq_old(WQ)
    wide.wq <- wide.wq[, c("EVENT_ID", "PH", "SPCOND")]
  }


  print("[4/6] Preparing Habitat Data")
  if(bibi == 2016){
    H <- prep_habitat(Habitat, agg_sample_num)
  }
  if(bibi == 2011){
    H <- prep_habitat_old(Habitat, agg_sample_num)
  }

  Prep <- prep_merge(TAB_EVENT, PT, TAB_STATIONS, wide.wq, H, TAB_PROJECT, TAB_HUC)
  just_date <- as.Date(Prep[, Sample_Date], format = "%m/%d/%Y")
  Prep$DATE <- format(just_date, "%m/%d/%Y")
  
  print("[5/6] Checking for SAMPLE_NUMBER column")
  Prep$SAMPLE_NUMBER <- if(agg_sample_num  == FALSE){
    if("SAMPLE_NUMBER" %in% colnames(Prep)){
      if(sum(is.na(Prep$SAMPLE_NUMBER)) > 0){
        warning("Warning: Missing SAMPLE_NUMBER values replaced with a 1.\
                Please review the values in the SAMPLE_NUMBER column before proceeding.")
        rep <- Prep$SAMPLE_NUMBER
        rep[is.na(rep)] <- 1
        Prep$SAMPLE_NUMBER <- rep

      }else{
        Prep$SAMPLE_NUMBER
      }
    }else{
      warning("Warning: Missing a SAMPLE_NUMBER column.\
              A SAMPLE_NUMBER column was created and filled with 1's as a placeholder.\
              Please review the SAMPLE_NUMBER values before proceeding.")
      1
    }
  }else{
    print("All sample numbers changed to 1, as specified by agg_sample_num = FALSE.
          All samples collected by the same agency at the same site and on the same date
          will be summed to form a composite sample.")
    1
  }
  
  
  print("[6/6] Checking for AGENCY_CODE column")
  Prep$AGENCY_CODE <- if("AGENCY_CODE" %in% colnames(Prep)){
    if(sum(is.na(Prep$AGENCY_CODE)) > 0){
      warning("Warning: Missing Agency Code replaced with a 1.\
              Please review the values in the AGENCY_CODE column before proceeding.")
      AC <- Prep$AGENCY_CODE
      AC[is.na(AC)] <- 1
      Prep$SAMPLE_NUMBER <- AC
    }else{
      paste(Prep$AGENCY_CODE, Prep$PROGRAM_CODE, sep = "_")
    }
  }else{
    warning("Warning: Missing a AGENCY_CODE column.\
            A AGENCY_CODE column was created and filled with 1's as a placeholder.\
            Please review the AGENCY_CODE column before proceeding.")
    1
  }
  cnames <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  #rm_blanks <- function(Long){
   # gsub("^\\s+|\\s+$", "", Long)
  #}
  #Prep[, cnames] <- apply(Prep[, cnames], 2, rm_blanks)
  Prep$HUC_12 <- factor(Prep$HUC_12)
  
  
  #exclude.taxa <- c("CATENULIDA", "SARCOMASTIGOPHORA","BRYOZOA", "CNIDARIA", "NEMATOMORPHA",
  #                  "PORIFERA", "MYRIAPODA", "POLYCHAETA", "BRANCHIOPODA",
   #                 "MAXILLOPODA", "OSTRACODA")
 # phy_spp <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
  #             "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  #Prep[Prep$PHYLUM %in% exclude.taxa, phy_spp] <- "UNIDENTIFIED"
  
  #subphy_spp <- c("SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "SUBORDER",
   #            "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  #Prep[Prep$SUBPHYLUM %in% exclude.taxa, subphy_spp] <- "UNIDENTIFIED"
  
  if(clean_taxa == TRUE) Prep <- Chessie::clean_taxa(Prep)
  
  # Only Kick-net samples
  final.df <- prep_subset(Prep)
  
  # Remove Jan, Feb, and Dec becuase these months were infrequently sampled
  final.df$MONTH <- lubridate::month(lubridate::mdy(final.df$DATE))
                                     
  if(month %in% TRUE){
    final.df <- final.df[!(final.df$MONTH %in% c(1, 2, 12)), ]
  }
  # Add Julian Day
  final.df$JULIAN <- lubridate::yday(lubridate::mdy(final.df$DATE))
  # Strahler Order <= 4
  final.df <- final.df[final.df$STRAHLER_STREAM_ORDER <= 4, ]
  
  #final.df <- prep_ecoregion(final.df)
  if(development == TRUE){
    prep.test <- final.df[, c("EVENT_ID", "SAMPLE_NUMBER", "REPORTING_VALUE")]
    agg.test <- aggregate(REPORTING_VALUE ~ ., data = prep.test, sum)
    agg.test2 <- aggregate(REPORTING_VALUE ~ EVENT_ID, data = agg.test, max)
    agg.final <- merge(agg.test, agg.test2, by = c("EVENT_ID", "REPORTING_VALUE"))
    agg.final <- agg.final[, !names(agg.final) %in% "REPORTING_VALUE"]
    
    final.df <- merge(agg.final, final.df, by = c("EVENT_ID", "SAMPLE_NUMBER"))
  }
  
  return(final.df)
}

#==============================================================================
#'Clean Taxa
#'
#'@param Bioregion.df = Bioregion data frame to be aggregated
#'@return Standardizes taxa according to the Chessie BIBI standards.
#'@export

clean_taxa <- function(long){
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                 "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  
  long[, taxa.cols] <- apply(long[, taxa.cols], 2, function(x){
    ifelse(is.na(x), "UNIDENTIFIED", as.character(x))
  } )
  
  phy.keep <- c("ANNELIDA", "ARTHROPODA", "MOLLUSCA", "PLATYHELMINTHES")
  long <- long[long$PHYLUM %in% phy.keep, ]
  
  subphy.keep <- c("CLITELLATA", "CRUSTACEA", "HEXAPODA", "RHABDITOPHORA", "UNIDENTIFIED")
  long <- long[long$SUBPHYLUM %in% subphy.keep, ]
  
  class.exc <- c("BRANCHIOPODA", "MAXILLOPODA", "OSTRACODA")
  long <- long[!(long$CLASS %in% class.exc), ]
  
  order.exc <- c("HYMENOPTERA")
  long <- long[!(long$ORDER  %in% order.exc), ]
  
  family.exc <- c("GERRIDAE", "HEBRIDAE", "VELIIDAE", "HYDROMETRIDAE", "SALDIDAE")
  long <- long[!(long$FAMILY  %in% family.exc), ]
  
  genus.exc <- c("STENUS")
  long <- long[!(long$GENUS  %in% genus.exc), ]
  
  
  class_spp <- c("CLASS", "SUBCLASS", "ORDER", "SUBORDER",
                 "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  long[long$CLASS %in% "BIVALVIA", class_spp] <- "BIVALVIA"
  long[long$CLASS %in% "GASTROPODA", class_spp] <- "GASTROPODA"
  long[long$CLASS %in% "OLIGOCHAETA", class_spp] <- "OLIGOCHAETA"
  long[long$CLASS %in% "TREPAXONEMATA", class_spp] <- "TREPAXONEMATA"
  
  
  order_spp <- c("ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  long[long$ORDER %in% "COLLEMBOLA", order_spp] <- "COLLEMBOLA"
  long[long$ORDER %in% "LEPIDOPTERA", order_spp] <- "LEPIDOPTERA"
  long[long$ORDER %in% "NEUROPTERA", order_spp] <- "NEUROPTERA"
  long[long$ORDER %in% "NEOOPHORA", order_spp] <- "NEOOPHORA"
  
  return(long)
}
#==============================================================================
#'Bioregion aggregation
#'
#'@param Bioregion.df = Bioregion data frame to be aggregated
#'@return Prepares the bioregion data for metric calculations.
#'@export

bioregion_agg <- function(Bioregion.df){
  bioregion.taxa <- Bioregion.df[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                     "AGENCY_CODE", "TSN", "PHYLUM",
                                     "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
                                     "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE",
                                     "GENUS", "SPECIES", "REPORTING_VALUE")]
  bioregion.taxa$FAMILY <- toupper(bioregion.taxa$FAMILY)

  agg.bioregion <- aggregate(REPORTING_VALUE ~ EVENT_ID + STATION_ID + DATE +
                               SAMPLE_NUMBER + AGENCY_CODE +TSN +
                               PHYLUM + SUBPHYLUM + CLASS + SUBCLASS + ORDER +
                               SUBORDER + FAMILY + SUBFAMILY + TRIBE + GENUS +
                               SPECIES, data = bioregion.taxa, FUN = sum,
                             na.rm = TRUE)
  colnames(agg.bioregion) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                               "AGENCY_CODE", "TSN", "PHYLUM",
                               "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
                               "SUBORDER", "FAMILY", "SUBFAMILY", "TRIBE",
                               "GENUS", "SPECIES", "REPORTING_VALUE")
  agg.bioregion[agg.bioregion == ""] <- NA
  return(agg.bioregion)
}

#==============================================================================
#'Prepare subset specific to BIBI
#'
#'@param Long = Data in a long format.
#'@return A subset of the data containing sampling events with strahler stream
#'order less than or equal to 4 and methods comparable to kick net samples.
#'@export

prep_subset <- function(Long){
  sub.set <- subset(Long,  Long$STRAHLER_STREAM_ORDER <= 4 &
                      Long$G_METHOD %in% c(49, 57, 58, 59, 86,
                                           87, 89, 92, 93, 94, 101, 103, 104))
  # added 93 and 104
  return(sub.set)
}
