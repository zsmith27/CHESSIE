#==============================================================================
#Site Classification
#==============================================================================
#'Site Classification
#'
#'@param env.param = a data frame of environmental parameters representitive
#'of each sampling event.  Needs to include water quality parameters
#'(Specific Conductivity, pH, and Dissolved Oxygen), habitat parameters
#' (Rapid Habitat Assessment), and (watershed parameters()).
#'@return Each sampling event is classified as Reference, Mixed, or Degraded
#'based on environmental parameters. If environmental parameters are missing
#'the site cannot be accurately classified as Reference or Degraded.
#'Therefore, the site is classified as Mixed.  The Mixed class does not
#'influence the selection of metrics but will be scored and assigned an
#'IBI rating.
#'@export

site_classification <- function(env.param){
  env.param$new_sc <- ifelse(env.param$SPCOND <= 100, 6,
                             ifelse(env.param$SPCOND > 100 & env.param$SPCOND <= 200, 5,
                                    ifelse(env.param$SPCOND > 200 & env.param$SPCOND <= 300, 4,
                                           ifelse(env.param$SPCOND > 300 & env.param$SPCOND <= 400, 3,
                                                  ifelse(env.param$SPCOND > 400 & env.param$SPCOND <= 500, 2,
                                                         ifelse(env.param$SPCOND > 500 & env.param$SPCOND <= 600, 1,
                                                                ifelse(env.param$SPCOND > 600, 0, 100)))))))

  env.param$new_ph <- ifelse(env.param$PH >= 6 & env.param$PH <= 8, 6,
                             ifelse(env.param$PH < 6 & env.param$PH > 5 , 3.5,
                                    ifelse(env.param$PH > 8 & env.param$PH < 9 , 3.5,
                                           ifelse(env.param$PH <= 5 | env.param$PH >= 9, 0, 100))))

  env.param$new_do <- ifelse(env.param$DO >= 5 & env.param$DO <= 10, 6,
                             ifelse(env.param$DO < 5 & env.param$DO > 3 , 3.5,
                                    ifelse(env.param$DO > 10 & env.param$DO < 12, 3.5,
                                           ifelse(env.param$DO <= 3 | env.param$DO >= 12, 0, 100))))


  env.param$ave <- apply(env.param[,c("new_sc", "new_ph", "new_do")], 1, mean)


  #Habitat Params================================================================
  env.param$new_banks <- ifelse(env.param$BANKS <= 2, 0,
                                ifelse(env.param$BANKS > 2 & env.param$BANKS <= 5, 1,
                                       ifelse(env.param$BANKS > 5 & env.param$BANKS <= 8, 2,
                                              ifelse(env.param$BANKS > 8 & env.param$BANKS <= 11, 3,
                                                     ifelse(env.param$BANKS > 11 & env.param$BANKS <= 14, 4,
                                                            ifelse(env.param$BANKS > 14 & env.param$BANKS <= 17, 5,
                                                                   ifelse(env.param$BANKS > 17 & env.param$BANKS <= 20, 6, 100)))))))

  env.param$new_ch_alt <- ifelse(env.param$CH_ALT <= 2, 0,
                                 ifelse(env.param$CH_ALT > 2 & env.param$CH_ALT <= 5, 1,
                                        ifelse(env.param$CH_ALT > 5 & env.param$CH_ALT <= 8, 2,
                                               ifelse(env.param$CH_ALT > 8 & env.param$CH_ALT <= 11, 3,
                                                      ifelse(env.param$CH_ALT > 11 & env.param$CH_ALT <= 14, 4,
                                                             ifelse(env.param$CH_ALT > 14 & env.param$CH_ALT <= 17, 5,
                                                                    ifelse(env.param$CH_ALT > 17 & env.param$CH_ALT <= 20, 6, 100)))))))

  env.param$new_hab_hetero <- ifelse(env.param$HAB_HETERO <= 2, 0,
                                     ifelse(env.param$HAB_HETERO > 2 & env.param$HAB_HETERO <= 5, 1,
                                            ifelse(env.param$HAB_HETERO > 5 & env.param$HAB_HETERO <= 8, 2,
                                                   ifelse(env.param$HAB_HETERO > 8 & env.param$HAB_HETERO <= 11, 3,
                                                          ifelse(env.param$HAB_HETERO > 11 & env.param$HAB_HETERO <= 14, 4,
                                                                 ifelse(env.param$HAB_HETERO > 14 & env.param$HAB_HETERO <= 17, 5,
                                                                        ifelse(env.param$HAB_HETERO > 17 & env.param$HAB_HETERO <= 20, 6, 100)))))))

  env.param$new_instr_cond <- ifelse(env.param$INSTR_COND <= 2, 0,
                                     ifelse(env.param$INSTR_COND > 2 & env.param$INSTR_COND <= 5, 1,
                                            ifelse(env.param$INSTR_COND > 5 & env.param$INSTR_COND <= 8, 2,
                                                   ifelse(env.param$INSTR_COND > 8 & env.param$INSTR_COND <= 11, 3,
                                                          ifelse(env.param$INSTR_COND > 11 & env.param$INSTR_COND <= 14, 4,
                                                                 ifelse(env.param$INSTR_COND > 14 & env.param$INSTR_COND <= 17, 5,
                                                                        ifelse(env.param$INSTR_COND > 17 & env.param$INSTR_COND <= 20, 6, 100)))))))

  env.param$new_rip_zone <- ifelse(env.param$RIP_ZONE <= 2, 0,
                                   ifelse(env.param$RIP_ZONE > 2 & env.param$RIP_ZONE <= 5, 1,
                                          ifelse(env.param$RIP_ZONE > 5 & env.param$RIP_ZONE <= 8, 2,
                                                 ifelse(env.param$RIP_ZONE > 8 & env.param$RIP_ZONE <= 11, 3,
                                                        ifelse(env.param$RIP_ZONE > 11 & env.param$RIP_ZONE <= 14, 4,
                                                               ifelse(env.param$RIP_ZONE > 14 & env.param$RIP_ZONE <= 17, 5,
                                                                      ifelse(env.param$RIP_ZONE > 17 & env.param$RIP_ZONE <= 20, 6, 100)))))))

  env.param$new_embed <- ifelse(env.param$EMBED <= 2, 0,
                                ifelse(env.param$EMBED > 2 & env.param$EMBED <= 5, 1,
                                       ifelse(env.param$EMBED > 5 & env.param$EMBED <= 8, 2,
                                              ifelse(env.param$EMBED > 8 & env.param$EMBED <= 11, 3,
                                                     ifelse(env.param$EMBED > 11 & env.param$EMBED <= 14, 4,
                                                            ifelse(env.param$EMBED > 14 & env.param$EMBED <= 17, 5,
                                                                   ifelse(env.param$EMBED > 17 & env.param$EMBED <= 20, 6, 100)))))))


  env.param$new_habitat <- ifelse(env.param$SUM >= 100 , 6,
                                  ifelse(env.param$SUM < 100 & env.param$SUM > 90, 5,
                                         ifelse(env.param$SUM <= 90 & env.param$SUM > 80, 4,
                                                ifelse(env.param$SUM <= 80 & env.param$SUM > 70, 3,
                                                       ifelse(env.param$SUM <= 70 & env.param$SUM > 60, 2,
                                                              ifelse(env.param$SUM <= 60 & env.param$SUM > 50, 1,
                                                                     ifelse(env.param$SUM <= 50, 0, 100)))))))


  #env.param <- na.omit(env.param)

  env.param$CATEGORY <- ifelse(env.param$ave >= 5 & env.param$new_habitat >= 5, "REF",
                               ifelse(env.param$ave >= 5 & env.param$new_habitat < 5, "MIXED",
                                      ifelse(env.param$ave <= 5 & env.param$new_habitat > 3, "MIXED",
                                             ifelse(env.param$ave <= 3 | env.param$new_habitat <= 3, "DEG", "ERROR"))))

  env.param$CATEGORY <- ifelse(is.na(env.param$ave) | is.na(env.param$new_habitat),
                               "MIXED", env.param$CATEGORY)

  env.class <- env.param[, c("EVENT_ID", "CATEGORY")]
  return(env.class)
}

#==============================================================================
#'Just Important Environmental Parameters
#'
#'@param Prep.Data = the output of the prep_data function.
#'@return A data frame of environmental parameters of interest.
#'@export

shrink_env <- function(Prep.Data){
  #shrink <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER",
  #            "PROJECT_ID", "ECOREGION_LEVEL_4", "LATITUDE", "LONGITUDE", "KARST_TYPE",
  #            "STRAHLER_STREAM_ORDER","HUC_12", "PCT_FOREST_2006", "BANKS",
  #            "CH_ALT", "HAB_HETERO", "INSTR_COND", "RIP_ZONE", "EMBED", "SUM",
  #            "PH", "SPCOND", "DO")
  #shrink <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER",
  #            "PROJECT_ID", "ECOREGION_LEVEL_4", "LATITUDE", "LONGITUDE",
  #            "KARST_TYPE", "STRAHLER_STREAM_ORDER","HUC_12", "PCT_FOREST_2006",
  #            "BANKS", "BANKV", "CH_ALT", "EMBED", "EPI_SUB", "FLOW", "RIFF",
  #            "SED", "RIP_ZONE", "VEL_D", "AESTH", "POOL", "SHAD", "INSTR",
  #            "HAB_HETERO", "INSTR_COND",
  #            "PH", "SPCOND", "DO", "WTEMP")
  if ("ICPRB_BIOREGION_ID" %in% names(Prep.Data)) {
    shrink <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER",
                "PROJECT_ID", "ECOREGION_LEVEL_4", "LATITUDE", "LONGITUDE",
                "STRAHLER_STREAM_ORDER", "HUC_12", "SUBREGION_DESCRIPTION",
                "ICPRB_BIOREGION_ID",
                "BANKS", "BANKV", "CH_ALT", "EMBED", "EPI_SUB", "FLOW", "RIFF",
                "SED", "RIP_ZONE", 
                "PH", "SPCOND", "DO", "WTEMP")
  } else {
    shrink <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER",
                "PROJECT_ID", "ECOREGION_LEVEL_4", "LATITUDE", "LONGITUDE",
                "STRAHLER_STREAM_ORDER", "HUC_12", "SUBREGION_DESCRIPTION",
                #"ICPRB_BIOREGION_ID",
                "BANKS", "BANKV", "CH_ALT", "EMBED", "EPI_SUB", "FLOW", "RIFF",
                "SED", "RIP_ZONE", 
                "PH", "SPCOND", "DO", "WTEMP")
  }
  
  env.final <- unique(Prep.Data[, shrink])
  return(env.final)
}

#==============================================================================
#'Add a Bioregion Column
#'
#'@param Long = a long data frame with Ecoregion level 4 classes.
#'@return A column classifying each site into a bioregion based on
#'the Ecoregion level 4 class.
#'@export

prep_bioregion <- function(Long){
  
  Long$ECOREGION_LEVEL_4 <- tolower(Long$ECOREGION_LEVEL_4)
  
  # returns string w/o leading or trailing whitespace
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  Long$ECOREGION_LEVEL_4 <- trim(Long$ECOREGION_LEVEL_4)
  
  NAPU.list <- c("60a", "60b", "60d", "60e", "83f")
  NCA.list <- c("62a", "62b", "62c",  "62d")
  PIED.list <- c("45c", "45e", "45f","45g", "58h",  "64d",  "64c",
                 "64b",  "64a")
  VAL.list <- c("67e", "67a", "67b", "67f", "67g")
  RIDGE.list <- c("66a", "66b",  "67d",  "67c", "67i",  "67h",  "69a",
                  "69b", "70c", "70C")
  SEP.list <- c("65m", "65n")
  MAC.list <- c("63b", "63c", "63d", "63e", "63f")
  
  
  Long$BIOREGION <- ifelse(Long$ECOREGION_LEVEL_4 %in% NAPU.list, "NAPU", 
                           ifelse(Long$ECOREGION_LEVEL_4 %in% NCA.list, "NCA",
                                  ifelse(Long$ECOREGION_LEVEL_4 %in% PIED.list, "PIED",
                                         ifelse(Long$ECOREGION_LEVEL_4 %in% VAL.list, "VAL",
                                                ifelse(Long$ECOREGION_LEVEL_4 %in% RIDGE.list, "RIDGE",
                                                       ifelse(Long$ECOREGION_LEVEL_4 %in% SEP.list, "SEP",
                                                              ifelse(Long$ECOREGION_LEVEL_4 %in% MAC.list, "MAC",
                                                                     "ERROR")))))))
  return(Long)
  
}

#==============================================================================
#'Remove Environmental Paramater Errors
#'
#'@param env.df = a data frame contatining environmental parameters.
#'@return Removes rows with environmental parameter errors.
#'@export

remove_errors <- function(env.df){
  env.df <- env.df[!is.na(env.df$STRAHLER_STREAM_ORDER), ]
  env.df <- env.df[env.df$STRAHLER_STREAM_ORDER <= 4, ]
  env.df$DO <- ifelse(env.df$DO < 22, env.df$DO, NA)
  env.df$SPCOND <- ifelse(env.df$SPCOND > 25, env.df$SPCOND, NA)
  env.df$WTEMP <- ifelse(env.df$WTEMP < 35, env.df$WTEMP, NA)
  env.df$PH <- ifelse(env.df$PH >= 0 & env.df$PH <= 14, env.df$PH, NA)
  
  #env.df$AESTH <- ifelse(env.df$AESTH >= 0 & env.df$AESTH <= 20, env.df$AESTH, NA)
  env.df$BANKS <- ifelse(env.df$BANKS >= 0 & env.df$BANKS <= 20, env.df$BANKS, NA)
  env.df$BANKV <- ifelse(env.df$BANKV >= 0 & env.df$BANKV <= 20, env.df$BANKV, NA)
  env.df$CH_ALT <- ifelse(env.df$CH_ALT >= 0 & env.df$CH_ALT <= 20, env.df$CH_ALT, NA)
  env.df$EMBED <- ifelse(env.df$EMBED >= 0 & env.df$EMBED <= 20, env.df$EMBED, NA)
  env.df$EPI_SUB <- ifelse(env.df$EPI_SUB >= 0 & env.df$EPI_SUB <= 20, env.df$EPI_SUB, NA)
  env.df$FLOW <- ifelse(env.df$FLOW >= 0 & env.df$FLOW <= 20, env.df$FLOW, NA)
  env.df$RIFF <- ifelse(env.df$RIFF >= 0 & env.df$RIFF <= 20, env.df$RIFF, NA)
  env.df$SED <- ifelse(env.df$SED >= 0 & env.df$SED <= 20, env.df$SED, NA)
  #env.df$RIP_ZONE <- ifelse(env.df$RIP_ZONE >= 0 & env.df$RIP_ZONE <= 20, env.df$RIP_ZONE, NA)
  #env.df$INSTR <- ifelse(env.df$INSTR >= 0 & env.df$INSTR <= 20, env.df$INSTR, NA)
  #env.df$POOL <- ifelse(env.df$POOL >= 0 & env.df$POOL <= 20, env.df$POOL, NA)
  #env.df$SHAD <- ifelse(env.df$SHAD >= 0 & env.df$SHAD <= 20, env.df$SHAD, NA)
  
  #env.df$HAB_HETERO <- ifelse(env.df$HAB_HETERO >= 0 & env.df$HAB_HETERO <= 20, env.df$HAB_HETERO, NA)
  #env.df$INSTR_COND <- ifelse(env.df$INSTR_COND >= 0 & env.df$INSTR_COND <= 20, env.df$INSTR_COND, NA)


  #Env.df <- Env.df[Env.df$SUM >= 0 & Env.df$SUM <= 120 | is.na(Env.df$SUM), ]
  return(env.df)
}

#==============================================================================
#'Prepare Environmental Parameters
#'
#'@param Prep.Data = the output of the prep_data function.
#'@return A data frame containing environmental parameters of interest,
#'bioregion classes, and no errors.
#'@export

prep_env <- function(Prep.data){
  prep.env <- shrink_env(Prep.data)
  #prep.bio <- prep_bioregion(prep.env)
  prep.bio <- prep_ecoregion(prep.env)
  final.df <- remove_errors(prep.bio)
  return(final.df)
}

#==============================================================================
#'Prepare Class Multi
#'
#'@param Prep_ENV = the output of the prep_env function.
#'@return A data frame containing environmental data to classify sites.
#'@export

prep_class_multi <- function(Prep_Env){
  
  h_params <- c("BANKS", "BANKV", "CH_ALT", "EMBED", "EPI_SUB", "FLOW", "RIFF", "SED")
  env.df <- Prep_Env
  env.df$MEAN_HAB <- apply(env.df[, h_params], 1, function(x) mean(x, na.rm = TRUE))
  env.df$MEDIAN_HAB <- apply(env.df[, h_params], 1, function(x) median(x, na.rm = TRUE))
  # Number of Missing habitat values
  env.df$NA_HAB <- apply(is.na(env.df[, h_params]), 1, sum)
  #==============================================================================
  # Number of Missing Water quality values
  wq_params <- c("SPCOND", "PH", "DO")
  env.df$NA_WQ <- apply(is.na(env.df[, wq_params]), 1, sum)
  
  env.df$DEG_SPCOND <- ifelse(env.df$SPCOND >= 1000, 3, 
                              ifelse(env.df$SPCOND >= 750, 2,
                                     ifelse(env.df$SPCOND > 500, 1,
                                            ifelse(is.na(env.df$SPCOND), NA, 0))))
  
  env.df$DEG_PH <- ifelse(env.df$PH > 9.5 | env.df$PH < 4, 3,
                          ifelse(env.df$PH > 9 | env.df$PH < 5, 2,
                                 ifelse(env.df$PH > 8.5 | env.df$PH < 6, 1,
                                        ifelse(is.na(env.df$PH), NA, 0))))
  env.df$DEG_DO <- ifelse(env.df$DO <= 5, 1,
                          ifelse(is.na(env.df$DO), NA, 0))
  
  env.df[c("MEAN_HAB", "DEG_SPCOND", "DEG_PH", "DEG_DO")][is.na(env.df[c("MEAN_HAB","DEG_SPCOND", "DEG_PH", "DEG_DO")])] <- 0
  
  env.df$DEG_HAB <- ifelse(env.df$MEAN_HAB >= 16 & env.df$NA_HAB <= 3, 0,
                           ifelse(env.df$MEAN_HAB >= 14 & env.df$NA_HAB <= 3, 1,
                                  ifelse(env.df$MEAN_HAB > 12 & env.df$NA_HAB <= 3, 2,
                                         ifelse(env.df$MEAN_HAB <= 12 & env.df$NA_HAB <= 3, 3, NA))))
  env.df$DEG_HAB2 <- ifelse(env.df$MEDIAN_HAB >= 16 & env.df$NA_HAB <= 3, 0,
                            ifelse(env.df$MEDIAN_HAB >= 14 & env.df$NA_HAB <= 3, 2,
                                   ifelse(env.df$MEDIAN_HAB > 12 & env.df$NA_HAB <= 3, 3,
                                          #ifelse(env.df$MEDIAN_HAB > 12 & env.df$NA_HAB <= 3, 3,   
                                          ifelse(env.df$MEDIAN_HAB <= 12 & env.df$NA_HAB <= 3, 6, NA))))
  env.df$DEG_6 <- rowSums(ifelse(env.df[, h_params] <= 6, 1, 0), na.rm = TRUE)
  env.df$DEG_10 <- rowSums(ifelse(env.df[, h_params] <= 10, 1, 0), na.rm = TRUE)
  env.df$DEG_12 <- rowSums(ifelse(env.df[, h_params] <= 12, 1, 0), na.rm = TRUE)
  #env.df$DEG_12G <- rowSums(ifelse(env.df[, h_params] > 12, 1, 0), na.rm = TRUE)
  env.df$DEG_16 <- rowSums(ifelse(env.df[, h_params] >= 16, 1, 0), na.rm = TRUE)
  env.df$DEG_HAB3 <- ifelse(env.df$DEG_6 >= 2, 3,
                            ifelse(env.df$DEG_10 >= 2, 2,
                                   ifelse(env.df$DEG_16 >= 2, 1, 0)))
  
  env.df$DEG_SUM <- rowSums(env.df[, c("DEG_HAB", "DEG_SPCOND", "DEG_PH", "DEG_DO")], na.rm = TRUE)
  
  env.df$CATEGORY <- ifelse(rowSums(env.df[, c("DEG_SPCOND", "DEG_PH", "DEG_DO")], na.rm = TRUE) == 0 &
                              (env.df$DEG_16 / (8 - env.df$NA_HAB)) * 100 >= 66,
                            "MIN", "MIX")
  
  env.df$CATEGORY <- ifelse(env.df$DEG_SPCOND == 0 & env.df$DEG_PH == 0 & env.df$DEG_DO == 0 &
                              (((env.df$DEG_16) / ((8 - env.df$NA_HAB))) * 100 >= 75) &
                              env.df$DEG_12 == 0, "REF", env.df$CATEGORY)
  
  env.df$CATEGORY <- ifelse(env.df$DEG_SPCOND == 0 & env.df$DEG_PH == 0 & env.df$DEG_DO == 0 & 
                              (env.df$DEG_16 == (8 - env.df$NA_HAB)) &
                              env.df$DEG_12 == 0, "REF+", env.df$CATEGORY)
  
  
  env.df$CATEGORY <- ifelse(rowSums(env.df[, c("DEG_SPCOND", "DEG_PH", "DEG_DO")], na.rm = TRUE) == 0 & 
                              (env.df$DEG_12 / (8 - env.df$NA_HAB)) * 100 >= 50,
                            "MOD", env.df$CATEGORY)
  
  env.df$CATEGORY <- ifelse(((env.df$DEG_6 / (8 - env.df$NA_HAB)) * 100 >= 50) |
                              env.df$DEG_SPCOND >= 1 | env.df$DEG_PH >= 1 | env.df$DEG_DO >= 1,
                            "DEG", env.df$CATEGORY)
  
  env.df$CATEGORY <- ifelse((env.df$DEG_SPCOND >= 1 | env.df$DEG_PH >= 1 | env.df$DEG_DO >= 1) & 
                              (env.df$DEG_6 / (8 - env.df$NA_HAB)) * 100 >= 50,
                            "DEG+", env.df$CATEGORY)
  env.df$CATEGORY <- ifelse((env.df$DEG_SPCOND  + env.df$DEG_PH + env.df$DEG_DO >= 2) |
                              (env.df$DEG_6 / (8 - env.df$NA_HAB)) * 100 >= 50,
                            "DEG2", env.df$CATEGORY)
  
  table(env.df$CATEGORY)
  #env.df$CATEGORY <- ifelse(env.df$DEG_SUM == 0 , "REF", NA)
  #env.df$CATEGORY <- ifelse(env.df$DEG_SUM > 0, "NEAR", env.df$CATEGORY)
  #env.df$CATEGORY <- ifelse(env.df$DEG_SUM > 0, "MIN", env.df$CATEGORY)
  #env.df$CATEGORY <- ifelse(env.df$DEG_SUM >= 2, "MOD", env.df$CATEGORY)
  #env.df$CATEGORY <- ifelse(env.df$DEG_SUM >= 3, "SEV", env.df$CATEGORY)
  #env.df$CATEGORY <- ifelse(env.df$EVENT_ID == 8258, "MIX", env.df$CATEGORY)
  env.df$CATEGORY <- ifelse((8 - env.df$NA_HAB) < 3, #| (8 - env.df$NA_WQ) >= 0,
                            "MIX", env.df$CATEGORY)
  
  #env.df$CATEGORY <- ifelse(env.df$DEG_6 >= 3 | env.df$SPCOND  >= 1000 | (env.df$PH < 6  | env.df$PH > 9), "SEV",
  #                          ifelse(env.df$DEG_16 == 0 | env.df$SPCOND  <= 500 | (env.df$PH >= 6 & env.df$PH <= 9), "REF", "MIX"))
  
  return(env.df)
}


