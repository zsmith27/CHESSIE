#==============================================================================
#==============================================================================
# Author: Zachary M. Smith
# Maintainer: Zachary M. Smith
# Organization: ICPRB
# Created: December 2016
# Updated: 3-15-2017
# Purpose: Create a subdirectory to help organize BIBI output.
#==============================================================================
#==============================================================================
#'Create a Subdirectory for BIBI Outpt
#'
#'@param main.dir <- Specificy the main directory (file path).
#'@param spatial.res <- Specify the spatial resolution of the data.
#'@param taxon.res <- Specify the taxonomic resolution of the data.
#'@return 
#'@export

create_subdir <- function(main.dir, spatial.res, taxon.res){
  month.year <- format(Sys.Date(), format="%B_%Y")
  todays.date <- format(Sys.Date(), format="%m_%d_%Y")
  #========================================================
  # Create a folder for the Month and Year.
  if (!dir.exists(file.path(main.dir, month.year))){
    dir.create(file.path(main.dir, month.year))
  } 
  
  main.dir <- paste0(main.dir, "/", month.year)
  #========================================================
  # Create a folder for todays date.
  if (!dir.exists(file.path(main.dir, todays.date))){
    dir.create(file.path(main.dir, todays.date))
  } 
  
  main.dir <- paste0(main.dir, "/", todays.date)
  #========================================================
  # Create a folder for the specified spatial resolution.
  if (!dir.exists(file.path(main.dir, spatial.res))){
    dir.create(file.path(main.dir, spatial.res))
  } 
  
  main.dir <- paste0(main.dir, "/", spatial.res)
  #========================================================
  # Create a folder for the specified taxonomic resolution.
  if (!dir.exists(file.path(main.dir, taxon.res))){
    dir.create(file.path(main.dir, taxon.res))
  } 
  main.dir <- paste0(main.dir, "/", taxon.res)
  # create a folder for each region or bioregion.
  #========================================================
  return(main.dir)
}