

#==============================================================================
#'Prepare Water Quality Data
#'
#'@param WQ = Water quality data in long format
#'@param agg_sample_num = Should the sample number be treated as seperate
#'samples (FALSE) or aggregated into a single sample (TRUE)?
#'@return Transform water quality data from long to wide format
#'@export

prep_wq_old <- function(WQ){
  water.q <- aggregate(WQ$REPORTED_VALUE ~ WQ$EVENT_ID + WQ$REPORTING_PARAMETER,
                       FUN = mean, na.rm = TRUE)
  colnames(water.q) <- c("EVENT_ID", "REPORTING_PARAMETER", "REPORTED_VALUE")
  
  wide.wq <- reshape(water.q, v.names = "REPORTED_VALUE", idvar = "EVENT_ID",
                     timevar = "REPORTING_PARAMETER", direction = "wide")
  #wide_WQ[is.na(wide_WQ)] <- 0 #NA = zero
  ColnamesRemovingPrefix <- function(df, prefix) {
    names <- colnames(df)
    indices <- (substr(names, 1, nchar(prefix)) == prefix)
    names[indices] <- substr(names[indices], nchar(prefix) + 1,
                             nchar(names[indices]))
    return(names)
  }
  
  wide.df <- reshape2::dcast(water.q, EVENT_ID ~ REPORTING_PARAMETER,
                             value.var = "REPORTED_VALUE")
  
  colnames(wide.wq) <- ColnamesRemovingPrefix(wide.wq, "REPORTED_VALUE.")
  #wide2_WQ <- wide_WQ[rowSums(wide_WQ[,2:ncol(wide_WQ)]) != 0,]
  
  return(wide.wq)
}


#==============================================================================
#'Prepare Habitat Data
#'
#'@param Habitat = Habitat data in long format
#'@return Transform habitat data from long to wide format
#'@export
prep_habitat_old <- function(Habitat, agg_sample_num){
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

  #wide_Habitat[is.na(wide_Habitat)]<-0 #NA = zero
  ColnamesRemovingPrefix <- function(df, prefix) {
    names <- colnames(df)
    indices <- (substr(names, 1, nchar(prefix)) == prefix)
    names[indices] <- substr(names[indices], nchar(prefix) + 1,
                             nchar(names[indices]))
    return(names)
  }
  
  #wide<-subset(wide, select=-REPORTED_VALUE.)
  colnames(wide.habitat) <- ColnamesRemovingPrefix(wide.habitat,
                                                   "REPORTED_VALUE.")
  #wide2_Habitat<-wide_Habitat[rowSums(wide_Habitat[,2:ncol(wide_Habitat)]) != 0,]
  habitat.data <- wide.habitat[ , c("EVENT_ID", "SAMPLE_NUMBER", "BANKS", "CH_ALT",
                                    "RIFF", "POOL", "EPI_SUB", "COVER",
                                    "INSTR", "RIP_SC", "RIP_W", "EMBED")]
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
  habitat.data$RIP_ZONE <- habitat.data$RIP_SC
  #Habitat_Data$RIP_ZONE<-ifelse( !is.na(Habitat_Data$RIP_Sc) & !is.na(Habitat_Data$RIP_SC) & !is.na(Habitat_Data$RIP_W), (Habitat_Data$RIP_Sc+Habitat_Data$RIP_SC+Habitat_Data$RIP_W)/3,
  #                               ifelse(is.na(Habitat_Data$RIP_Sc) & !is.na(Habitat_Data$RIP_SC)& !is.na(Habitat_Data$RIP_W), (Habitat_Data$RIP_SC+Habitat_Data$RIP_W)/2,
  #                                      ifelse(!is.na(Habitat_Data$RIP_Sc) & !is.na(Habitat_Data$RIP_SC)& is.na(Habitat_Data$RIP_W), (Habitat_Data$RIP_Sc+Habitat_Data$RIP_SC)/2,
  #                                             ifelse(!is.na(Habitat_Data$RIP_Sc) & is.na(Habitat_Data$RIP_SC)& !is.na(Habitat_Data$RIP_W), (Habitat_Data$RIP_Sc+Habitat_Data$RIP_W)/2,
  #                                                    ifelse(!is.na(Habitat_Data$RIP_Sc) & is.na(Habitat_Data$RIP_SC)& is.na(Habitat_Data$RIP_W), Habitat_Data$RIP_Sc,
  #                                                           ifelse(is.na(Habitat_Data$RIP_Sc) & !is.na(Habitat_Data$RIP_SC)& is.na(Habitat_Data$RIP_W), Habitat_Data$RIP_SC,
  #                                                                  ifelse(is.na(Habitat_Data$RIP_Sc) & is.na(Habitat_Data$RIP_SC)& !is.na(Habitat_Data$RIP_W), Habitat_Data$RIP_W, NA)))))))
  
  env.data <- habitat.data[, c("EVENT_ID", "SAMPLE_NUMBER", "BANKS", "CH_ALT", "HAB_HETERO",
                               "INSTR_COND", "RIP_ZONE", "EMBED")]
  env.data$SUM <- rowSums(env.data[, 2:ncol(env.data)])
  env <- env.data[! is.na(env.data$SUM), ]
  return(env)
}

#==============================================================================
#'Site Classification (2011)
#'
#'@param env.parm = a data frame containing the necessary habitat and water
#'quality variables to classify sites based on the 2011 Chessie BIBI
#'protocol.
#'@return A dataframe containg site classes (CATEGORY).
#'@export

site_classification_2011 <- function(env.param){
  env.param <- unique(env.param)
  env.param <- prep_bioregion(env.param)
  #============================================================================== 
  pied_val.class <- function(env.df){
    hab.params <- c("BANKS", "CH_ALT", "EMBED", "HAB_HETERO",
                    "INSTR_COND", "RIP_ZONE")
    env.df$HAB_REF <- rowSums(env.df[, hab.params] >= 16, na.rm = T)
    env.df$HAB_DEG <- rowSums(env.df[, hab.params] <= 5, na.rm = T)
    env.df$CATEGORY <- ifelse(env.df$HAB_REF == 6 & 
                                   env.df$SPCOND <= 500 &
                                   env.df$PH <= 9 & env.df$PH >= 6, "REF",
                                 ifelse(env.df$HAB_DEG >= 3 |
                                          env.df$SPCOND > 1000 |
                                          env.df$PH > 9 |
                                          env.df$PH < 6, "SEV", "MIX"))
    return(env.df)
  }
  #==============================================================================
  ridges.class <- function(env.df){
    hab.params <- c("BANKS", "CH_ALT", "EMBED", "HAB_HETERO",
                    "INSTR_COND", "RIP_ZONE")
    env.df$HAB_REF <- rowSums(env.df[, hab.params] >= 16)
    env.df$HAB_DEG <- rowSums(env.df[, hab.params] <= 5)
    env.df$CATEGORY <- ifelse(env.df$HAB_REF == 6 & 
                                   env.df$SPCOND <= 500 &
                                   env.df$PH <= 9 & env.df$PH >= 6, "REF",
                                 ifelse(env.df$HAB_DEG >= 3 |
                                          env.df$SPCOND > 1000 |
                                          env.df$PH > 9.5 |
                                          env.df$PH < 5, "SEV", "MIX"))
    return(env.df)
  }
  #==============================================================================
  nca.class <- function(env.df){
    hab.params <- c("BANKS", "CH_ALT", "EMBED", "HAB_HETERO",
                    "INSTR_COND", "RIP_ZONE")
    env.df$HAB_REF <- rowSums(env.df[, hab.params] >= 15)
    env.df$HAB_DEG <- rowSums(env.df[, hab.params] <= 5)
    env.df$CATEGORY <- ifelse(env.df$HAB_REF == 6 & 
                                   env.df$SPCOND <= 500 &
                                   env.df$PH <= 9 & env.df$PH >= 6, "REF",
                                 ifelse(env.df$HAB_DEG >= 3 |
                                          env.df$SPCOND > 500 |
                                          env.df$PH > 9 |
                                          env.df$PH < 5.5, "SEV", "MIX"))
    return(env.df)
  }
  #==============================================================================
  napu.class <- function(env.df){
    hab.params <- c("BANKS", "CH_ALT", "EMBED", "HAB_HETERO",
                    "INSTR_COND", "RIP_ZONE")
    env.df$HAB_REF <- rowSums(env.df[, hab.params] >= 15)
    env.df$HAB_DEG <- rowSums(env.df[, hab.params] <= 5)
    env.df$CATEGORY <- ifelse(env.df$EMBED >= 6 & 
                                   env.df$SPCOND <= 500 &
                                   env.df$PH <= 9 & env.df$PH >= 6, "REF",
                                 ifelse(env.df$EMBED <=15 &
                                          env.df$SPCOND > 500 |
                                          env.df$PH > 9 |
                                          env.df$PH < 5.5, "SEV", "MIX"))
    return(env.df)
  }
  #==============================================================================
  if("PIED" %in% env.param$BIOREGION | "VAL" %in% env.param$BIOREGION){
    env.pied_val <- env.param[env.param$BIOREGION %in% c("PIED", "VAL"), ]
    env.class.pied_val <- pied_val.class(env.pied_val)
  } 
  if("RIDGE" %in% env.param$BIOREGION){
    env.ridges <- env.param[env.param$BIOREGION %in% c("RIDGE"), ]
    env.class.ridges <- ridges.class(env.ridges)
  }

  if("NCA" %in% env.param$BIOREGION){
    env.nca <- env.param[env.param$BIOREGION %in% c("NCA"), ]
    env.class.nca <- nca.class(env.nca)
  } 
  if("NAPU" %in% env.param$BIOREGION){
    env.napu <- env.param[env.param$BIOREGION %in% c("NAPU"), ]
    env.class.napu <- napu.class(env.napu)
  } 
  #==============================================================================
  
  env.class <- rbind(if(exists("env.class.pied_val")){env.class.pied_val},
                     if(exists("env.class.ridges")){env.class.ridges},
                     if(exists("env.class.nca")){env.class.nca},
                     if(exists("env.class.napu")){env.class.napu})

  #==============================================================================
  env.class$CATEGORY <- ifelse(is.na(env.class$CATEGORY), "MIX", as.character(env.class$CATEGORY))
  
  final.df <- unique(env.class)
  return(final.df)
  
}


