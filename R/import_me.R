#==============================================================================
#'Import Output
#'
#'@param job = When equal to NULL the data is not modified beyond the file imported, just 
#'aggregated with the other bioregion. When equal to "METRICS" the data is imported, transformed
#'to a long data type, and aggregated with the other bioregions.
#'@param high.res.front = specify the common prefix to the high resolution bioregion output files.
#'@param low.res.front = specify the common prefix to the low resolution bioregion output files.
#'@param high.res.proc_date = specify the date associated with the high resolution bioregion output files.
#'@param low.res.proc_date = specify the date associated with the low resolution bioregion output files.
#'@param high.res.data_set = specify the data set of interest associated with the high resolution bioregion output files.
#'@param low.res.data_set = specify the data set of interest associated with the low resolution bioregion output files.
#'@return Import and aggregate BIBI output files.
#'@export
read.me <- function(high.res.front, high.res.proc_date, high.res.data_set,
                    low.res.front, low.res.proc_date, low.res.data_set, job = NULL){
  #==============================================================================
  tidy.me <- function(bioregion) tidyr::gather(bioregion, METRIC, SCORE, 8:ncol(bioregion))
  #==============================================================================
  high.res.bioregions <- c("BLUE", "CA", "LNP", "MAC", "NAPU", "NCA", "NRV", "PIED", "SRV", "UNP")
  
  datalist <- list()
  for (i in high.res.bioregions){
    bio <- read.csv(paste(high.res.front, i, high.res.data_set, high.res.proc_date, ".csv", sep = ""))
    bio <- bio[, !names(bio) %in% "X"]
    bio$BIOREGION <- i
    bio <- bio[, c("BIOREGION", (1:ncol(bio))[-"BIOREGION"])]
    if(job %in% "METRICS") bio <- tidy.me(bio)
    datalist[[i]] <- bio # add it to your list
  }
  high.res <- unique(do.call(rbind, datalist))
  #==============================================================================
  low.res.bioregions <- c("COAST", "NON_COAST", "BASIN")
  
  datalist <- list()
  for (i in low.res.bioregions){
    bio <- read.csv(paste(low.res.front, i, low.res.data_set, low.res.proc_date, ".csv", sep = ""))
    bio <- bio[, !names(bio) %in% "X"]
    bio$BIOREGION <- i
    bio <- bio[, c("BIOREGION", (1:ncol(bio))[-"BIOREGION"])]
    if(job %in% "METRICS") bio <- tidy.me(bio)
    datalist[[i]] <- bio # add it to your list
  }
  low.res <- unique(do.call(rbind, datalist))
  #==============================================================================
  bound <- rbind(high.res, low.res)
  
  return(bound)
}


#==============================================================================
#==============================================================================
#'Group Output Files
#'
#'@param job = When equal to NULL the data is not modified beyond the file imported, just 
#'aggregated with the other bioregion. When equal to "METRICS" the data is imported, transformed
#'to a long data type, and aggregated with the other bioregions.
#'@param high.res.front = specify the common prefix to the high resolution bioregion output files.
#'@param low.res.front = specify the common prefix to the low resolution bioregion output files.
#'@param high.res.proc_date = specify the date associated with the high resolution bioregion output files.
#'@param low.res.proc_date = specify the date associated with the low resolution bioregion output files.
#'@param high.res.data_set = specify the data set of interest associated with the high resolution bioregion output files.
#'@param low.res.data_set = specify the data set of interest associated with the low resolution bioregion output files.
#'@return Import, aggregate, and modify the BIBI output into a standard table.
#'@export
group.me <- function(job, high.res.front, high.res.proc_date, high.res.data_set,
                     low.res.front, low.res.proc_date, low.res.data_set){
  
  bound <- read.me(high.res.front, high.res.proc_date, high.res.data_set,
                   low.res.front, low.res.proc_date, low.res.data_set, job)
  #==============================================================================
  if("ORGINAL_VALUE" %in% names(bound) | "ORGINAL_THRESHOLD" %in% names(bound)){
    names(bound)[names(bound) %in% "ORGINAL_VALUE"|
                   names(bound) %in% "ORGINAL_THRESHOLD"] <- "ORIGINAL_VALUE"
  } 
  
  if("MEAN_SIM_THRESHOLD" %in% names(bound)){
    names(bound)[names(bound) %in% "MEAN_SIM_THRESHOLD"] <- "MEAN_SIM_VALUE"
  } 
  return(bound)
}


#==============================================================================
