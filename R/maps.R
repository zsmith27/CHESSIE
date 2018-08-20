#==============================================================================
#'Prep event ratings
#'
#'@export

organize_all_event_ratings <- function(all.scores, spatial, taxon.rank, todays.date){
  # Keep only the necessary columns.
  all.scores <- all.scores[, c("BIOREGION", "EVENT_ID", "STATION_ID", "HUC_8", "HUC_10", "HUC_12",
                               "DATE", "SAMPLE_NUMBER", "AGENCY_CODE", "CATEGORY",
                               "RANDOM_SAMPLING", "FINAL_SCORE", "DEG_50",
                               "REF_10", "REF_25", "REF_50", "RATING")]
  # DEG_50 is no longer used, so replace with "HALF_REF_10", which refers to the
  # value that is half of the 10th Reference percentile.
  names(all.scores)[names(all.scores) %in% "DEG_50"] <- "HALF_REF_10"
  # Multiply the scores and Reference thresholds by 100 and round to the second 
  # decimal place.
  all.scores[, 12:16] <- round(all.scores[, 12:16] * 100, 2)
  # Convert the Disturbance Gradient categories to class factor.
  # This allows for grouping by the CATEGORY and presenting the results in the
  # specified order (REF to SEV).
  all.scores$CATEGORY <- factor(all.scores$CATEGORY, levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
  # Sort the data frame by Bioregion (alphabetical), Category (by factor),
  # and Score (decending).
  all.scores <- all.scores[order(all.scores$BIOREGION, all.scores$CATEGORY,
                                 -all.scores$FINAL_SCORE), ]
  # Save the rating and scores data frame as a .csv table.
  write.csv(all.scores, paste(spatial, taxon.rank, todays.date,
                              "All_Event_Ratings.csv", sep = "_"), row.names = F)
  # End organize_all_event_ratings function.
  return(all.scores)
}
#==============================================================================
#'Caclulate Reference and Severely Degraded percentiles
#'
#'@export

ref_sev_percentiles <- function(all.scores, spatial, taxon.rank, todays.date) {
  # Create a new data frame with each row representing a unique Bioregion.
  # Bioregion also represents Basin and Regions.
  new.df <- data.frame(BIOREGION = unique(all.scores$BIOREGION))
  # Create a blank list to store the for loop output.
  data.list <- list()
  # The for loop
  for (i in unique(new.df$BIOREGION)) {
    # Subset the data frame to only include the bioregion represented by i.
    sub.scores <- all.scores[all.scores$BIOREGION %in% i, ]
    
    # Identify outliers, remove outliers, calculate percentiles without outliers.
    quick_percentiles <- function(scores.df, bioregion, category) {
      outlier <- function(scores.df, mult.val, job = "REF"){
        # Create a blank list to store for loop output.
        datalist = list()
        for (i in unique(scores.df$BIOREGION)) {
          bioregion.df <- scores.df[scores.df$BIOREGION %in% i, ]
          # Find the IQR.
          quantile.vec <- quantile(bioregion.df$FINAL_SCORE, c(0.25, 0.75))
          # 25th - (IQR * 1.5)
          if (job %in% "REF"){
            outlier.thresh <- quantile.vec[1] - (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
          }
          # 75th + (IQR * 1.5)
          if (job %in% "SEV"){
            outlier.thresh <- quantile.vec[2] + (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
          }
          # Create a new data frame and fill the BIOREGION column with i.
          outlier.df <- data.frame(BIOREGION = i)
          # Identify the outlier threshold for the bioregion.
          outlier.df$THRESH <- outlier.thresh
          # add the outlier.df data frame to the list.
          datalist[[i]] <- outlier.df
        }
        # Append the list of data frames together.
        final.df <- do.call(rbind, datalist)
        # End the outlier function.
        return(final.df)
      }
      
      # Subset the data frame to only include rows that represent the specified
      # disturbance gradient category (i.e., REF or SEV).
      sub.df <- scores.df[scores.df$CATEGORY %in% category, ]
      # Changed job = "REF" to job = category (3/20/17)
      out.df <- outlier(sub.df, 1.5, job = category)
      # If category is REF, then we are only worried about the lower outlier threshold.
      if (category %in% "REF") sub.df <- sub.df[sub.df$FINAL_SCORE >= out.df$THRESH, ]
      # If category is SEV, then we are only worried about the upper outlier threshold.
      if (category %in% "SEV") sub.df <- sub.df[sub.df$FINAL_SCORE <= out.df$THRESH, ]
      # Calculate the specified percentiles of the remaining scores.
      sub.df <- data.frame(t(quantile(sub.df$FINAL_SCORE, c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95))))
      # Create a "BIOREGION" column and fill with the specified bioregion.
      sub.df$BIOREGION <- bioregion
      # Create a "CATEGORY" column and fill with the specified category.
      sub.df$CATEGORY <- category
      # Rearrange the columns.
      sub.df <- sub.df[, c(7, 8, 1:6)]
      # Re-name the columns.
      names(sub.df) <- c("Bioregion", "Category", "5%", "10%", "25%", "50%", "75%", "95%")
      # End quick_percentiles function.
      return(sub.df)
    }
    
    # Reference percentiles.
    sub.ref <- quick_percentiles(sub.scores, i, "REF")
    # Severely Degraded percentiles.
    sub.sev <- quick_percentiles(sub.scores, i, "SEV")
    # Append the two data frames together.
    sub.final <- rbind(sub.ref, sub.sev)
    # Add the data frame to the list.
    data.list[[i]] <- sub.final
  }
  # Append the list of data frames together.
  ref_sev.pct <- do.call(rbind, data.list)
  # Round the percentile values to the second decimal place.
  ref_sev.pct[, 3:8] <- round(ref_sev.pct[, 3:8], 2)
  # Export the Reference and Severely Degraded percentile values as a .csv table.
  write.csv(ref_sev.pct, paste(spatial, taxon.rank,
                               todays.date, "REF_SEV_Score_Percentiles.csv",
                               sep = "_"), row.names = F)
  # End ref_sev_percentiles function.
  return(ref_sev.pct)
}

#==============================================================================
#'Prepare HUC data to be mapped
#'
#'@export

agg_hucs_map <- function(huc.data, final.scores){
  map.me <- unique(merge(huc.data, final.scores, by = c("EVENT_ID", "STATION_ID"), all.y = T))
  #map.me$HUC12_BIO <- paste("0", map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  map.me$HUC12_BIO <- paste(map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  huc12 <- data.frame(unique(map.me$HUC12_BIO))
  # Find the mean score per station
  map.me <- aggregate(FINAL_SCORE ~ STATION_ID + HUC_8 + HUC_10 + HUC_12 +
                        HUC12_BIO + ACRES + SQ_MILE, data = map.me, mean)
  #test <- data.frame(huc12[!huc12$unique.map.me.HUC_12. %in% map.me$HUC_12,])
  #==============================================================================
  map.me <- unique(merge(huc.data, final.scores, by = c("EVENT_ID", "STATION_ID"), all = TRUE))
  # NAs are present because the sample had <= 70 total taxa reported.
  # We decided that <= 70 could skew the metric value and thus it was safer to 
  # exclude the sample.
  map.me <- map.me[!is.na(map.me$BIOREGION), ]
  #map.me$HUC12_BIO <- paste("0", map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  map.me$HUC12_BIO <- paste(map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  return(map.me)
}
#==============================================================================
#'Aggregate scores by station
#'
#'@export
agg_station_map <- function(map.me, pct_sum, spatial, taxon.rank, todays.date){
  # Find the mean of the scores by station ID and bioregion.
  map.station <- aggregate(FINAL_SCORE ~ STATION_ID + BIOREGION, data = map.me, mean)
  # Find the total number of scores reported for a Station.
  map.station.length <- aggregate(FINAL_SCORE ~ STATION_ID + BIOREGION, data = map.me, length)
  # Change the column names.
  names(map.station.length) <- c("STATION_ID", "COUNT")
  # Join the data frame that contains the mean score by station ID and bioregion with
  # the data frame that contains the total count of scores for a particular station.
  map.station <- merge(map.station, map.station.length, by = "STATION_ID")
  # Apply the rating scheme.
  map.station <- rate.me2(map.station, pct_sum, taxon.rank)
  # Export the data frame as a .csv table.
  write.csv(map.station, paste(spatial, taxon.rank, todays.date,
                               "STATION.csv", sep = "_"), row.names = F)
  # End agg_station_map function/
  return(map.station)
}
#==============================================================================
#'Aggregate scores by HUC12
#'
#'@export
agg_huc12_map <- function(map.station, map.me, pct_sum, spatial, taxon.rank, todays.date,
                          agg.col = "STATION_ID", catch.count){
  # Append the HUC12_BIO column to the data frame containing
  # the mean BIBI scores for each station.
  if (agg.col %in% "STATION_ID") {
    merged.huc12 <- merge(map.station,
                          unique(map.me[, c(agg.col, "HUC12_BIO")]),
                          by = agg.col)
  } else {
    merged.huc12 <- map.station
  }

  # Find the mean score for each HUC12 Bioregion.
  map.huc12 <- aggregate(FINAL_SCORE ~ HUC12_BIO + BIOREGION,
                         data = merged.huc12,  mean)
  # Find the total number of stations found within the HUC12 Bioregion.
  map.huc12.length <- aggregate(FINAL_SCORE ~ HUC12_BIO + BIOREGION,
                                data = merged.huc12,  length)
  # Change the column names.
  names(map.huc12.length) <- c("HUC12_BIO", "BIOREGION", "COUNT")
  # Join the average HUC12 scores with the count of stations within each HUC12.
  map.huc12 <- merge(map.huc12, map.huc12.length[, c("HUC12_BIO", "COUNT")],
                     by = c("HUC12_BIO"))
  # Apply rating scheme.
  map.huc12 <- rate.me2(map.huc12, pct_sum, taxon.rank)
  # If the number of stations within a HUC12 is less than five, then
  # specify that there is an insufficient number of samples to
  # classify the given HUC12.
  map.huc12$RATING <- ifelse(map.huc12$COUNT < catch.count, "Insufficient Samples",
                             as.character(map.huc12$RATING))
  # Make sure the "HUC12_BIO" column is a character class.
  map.huc12$HUC12_BIO <- as.character(map.huc12$HUC12_BIO)
  # Append the HUC_12 column on to the data frame.
  final.df <- merge(map.huc12, unique(map.me[, c("HUC12_BIO", "HUC_12",
                                                 "ACRES", "SQ_MILE")]),
                    by = "HUC12_BIO")
  # Export the data frame as a .csv table.
  write.csv(final.df, paste(spatial, taxon.rank, todays.date,
                             "NEW_HUC12.csv", sep = "_"), row.names = F)
  # End agg_huc12_map function.
  return(final.df)
}






#==============================================================================
#'Aggregate scores by catchment
#'
#'@export
agg_catchment_map <- function(map.station, map.me, all.data, pct_sum, spatial,
                              taxon.rank, todays.date, agg.col = "FEATUREID",
                              catch.count = 3, sub.dir){
  #==============================================================================
  # The location of the Catchment Shapefiles.
  setwd("D:/ZSmith/Projects/Chessie_BIBI/GIS/spatial_02")
  # Import the Catchment Shapefiles.
  catchments <- rgdal::readOGR(".", "Catchments02", stringsAsFactors = FALSE)
  # Transform the HUC12 Shapefiles to prepare for leaflet.
  catchments <- sp::spTransform(catchments, sp::CRS("+init=epsg:4326"))
  #==============================================================================
  stations.df <- merge(map.station, 
                       unique(all.data[, c("STATION_ID", "LATITUDE",
                                           "LONGITUDE", "HUC_12")]),
                       by = "STATION_ID")
  # Specify the coordinates of each station.
  sp::coordinates(stations.df) <- ~ LONGITUDE + LATITUDE
  # Set the projection of the SpatialPointsDataFrame using the projection of 
  # the shapefile.
  sp::proj4string(stations.df) <- sp::proj4string(catchments)
  #==============================================================================
  # Identify the catchments that contain at least one sample.
  #keep.catchments <- sp::over(stations.df, catchments)
  # Keep only the catchments that were found to contain one or more samples.
  #catchments <- catchments[catchments@data$FEATUREID %in% keep.catchments$FEATUREID, ]
  # Identify the Catchment that the station is located within.
  stations.df$FEATUREID <- sp::over(stations.df, catchments)$FEATUREID
  
  stations.df$CATCH_AREA <- sp::over(stations.df, catchments)$AreaSqKM
  #==============================================================================
  stations.df <- as.data.frame(stations.df)
  #==============================================================================
  catch.huc.list <- lapply(unique(stations.df$FEATUREID), function(x) {
    sub.catch <- stations.df[stations.df$FEATUREID %in% x, ]
    sub.catch$HUC_12 <- as.character(sub.catch$HUC_12)
    final.df <- as.data.frame(table(sub.catch$HUC_12))
    final.df$Var1 <- as.character(final.df$Var1)
    final.df$FEATUREID <- x
    max.freq <- max(final.df$Freq)
    final.df <- final.df[final.df$Freq == max.freq, c("FEATUREID", "Var1")]
    names(final.df) <- c("FEATUREID", "STATION_HUC12")
    if (nrow(final.df) > 1) final.df <- final.df[1, ]
    return(final.df)
  })
  catch.huc.df <- do.call(rbind, catch.huc.list)
  stations.df <- merge(stations.df, catch.huc.df, by = "FEATUREID")
  #==============================================================================
  # Append the HUC12_BIO column to the data frame containing
  # the mean BIBI scores for each station.
  merged.catchment <- merge(map.station,
                            unique(stations.df[, c("STATION_ID", "FEATUREID",
                                                   "CATCH_AREA", "STATION_HUC12")]),
                            by = "STATION_ID")
  # Find the mean score for each catchment Bioregion.
  map.catchment <- aggregate(FINAL_SCORE ~ FEATUREID + BIOREGION + STATION_HUC12,
                             data = merged.catchment,  mean)
  rate.catchment <- rate.me2(map.catchment, pct_sum, taxon.rank)
  # Export the data frame as a .csv table.
  setwd(sub.dir)
  write.csv(rate.catchment, paste(spatial, taxon.rank, todays.date,
                                  "CATCHMENT.csv", sep = "_"), row.names = F)
  # Find the total number of stations found within the catchment Bioregion.
  #map.catchment.length <- aggregate(FINAL_SCORE ~ FEATUREID + BIOREGION,
  #                                  data = merged.catchment,  length)
  # Change the column names.
  #names(map.catchment.length) <- c("FEATUREID", "BIOREGION", "COUNT")
  # Join the average catchment scores with the count of stations within each catchment.
  #map.catchment <- merge(map.catchment, map.catchment.length[, c("FEATUREID", "COUNT")],
  #                       by = c("FEATUREID"))
  map.catchment$HUC12_BIO <- paste(map.catchment$STATION_HUC12, map.catchment$BIOREGION, sep = "_")
  #==============================================================================
  map.catchment <- agg_huc12_map(map.catchment, map.me, pct_sum, spatial, taxon.rank, todays.date,
                                 agg.col = "FEATUREID", catch.count)
  #==============================================================================
  # Apply rating scheme.
  #map.catchment <- rate.me2(map.catchment, pct_sum, taxon.rank)
  # If the number of stations within a catchment is less than five, then
  # specify that there is an insufficient number of samples to
  # classify the given catchment.
  map.catchment$RATING <- ifelse(map.catchment$COUNT < catch.count, "Insufficient Samples",
                                 as.character(map.catchment$RATING))
  # Make sure the "FEATUREID" column is a character class.
  #map.catchment$FEATUREID <- as.character(map.catchment$FEATUREID)
  # Append the HUC_12 column on to the data frame.
  #final.df <- merge(map.catchment, unique(map.me[, c("FEATUREID", "HUC_12",
  #                                                   "Acres", "SQ_MILE")]),
  #                  by = "FEATUREID")
  # Export the data frame as a .csv table.
  final.df <- map.catchment
  setwd(sub.dir)
  write.csv(final.df, paste(spatial, taxon.rank, todays.date,
                            "HUC12.csv", sep = "_"), row.names = F)
  # End agg_catchment_map function.
  return(final.df)
}
#==============================================================================
#'Aggregate scores by catchment
#'
#'@export
agg_catchment_map_old <- function(map.station, map.me, all.data, pct_sum, spatial, taxon.rank, todays.date,
                              agg.col = "FEATUREID", catch.count = 3, sub.dir){
  #==============================================================================
  # The location of the Catchment Shapefiles.
  setwd("D:/ZSmith/Projects/Chessie_BIBI/GIS/spatial_02")
  # Import the Catchment Shapefiles.
  catchments <- rgdal::readOGR(".", "Catchments02", stringsAsFactors = FALSE)
  # Transform the HUC12 Shapefiles to prepare for leaflet.
  catchments <- sp::spTransform(catchments, sp::CRS("+init=epsg:4326"))
  #==============================================================================
  stations.df <- merge(map.station, 
                       unique(all.data[, c("STATION_ID", "LATITUDE",
                                           "LONGITUDE", "HUC_12")]),
                       by = "STATION_ID")
  # Specify the coordinates of each station.
  sp::coordinates(stations.df) <- ~ LONGITUDE + LATITUDE
  # Set the projection of the SpatialPointsDataFrame using the projection of 
  # the shapefile.
  sp::proj4string(stations.df) <- sp::proj4string(catchments)
  #==============================================================================
  # Identify the catchments that contain at least one sample.
  keep.catchments <- sp::over(stations.df, catchments)
  # Keep only the catchments that were found to contain one or more samples.
  catchments <- catchments[catchments@data$FEATUREID %in% keep.catchments$FEATUREID, ]
  # Identify the Catchment that the station is located within.
  stations.df$FEATUREID <- sp::over(stations.df, catchments)$FEATUREID
  #==============================================================================
  stations.df <- as.data.frame(stations.df)
  #==============================================================================
  # Append the HUC12_BIO column to the data frame containing
  # the mean BIBI scores for each station.
  merged.catchment <- merge(map.station,
                            unique(stations.df[, c("STATION_ID", "FEATUREID", "HUC_12")]),
                            by = "STATION_ID")
  # Find the mean score for each catchment Bioregion.
  map.catchment <- aggregate(FINAL_SCORE ~ FEATUREID + BIOREGION + HUC_12,
                             data = merged.catchment,  mean)
  rate.catchment <- rate.me2(map.catchment, pct_sum, taxon.rank)
  # Export the data frame as a .csv table.
  setwd(sub.dir)
  write.csv(rate.catchment, paste(spatial, taxon.rank, todays.date,
                               "CATCHMENT.csv", sep = "_"), row.names = F)
  # Find the total number of stations found within the catchment Bioregion.
  #map.catchment.length <- aggregate(FINAL_SCORE ~ FEATUREID + BIOREGION,
  #                                  data = merged.catchment,  length)
  # Change the column names.
  #names(map.catchment.length) <- c("FEATUREID", "BIOREGION", "COUNT")
  # Join the average catchment scores with the count of stations within each catchment.
  #map.catchment <- merge(map.catchment, map.catchment.length[, c("FEATUREID", "COUNT")],
  #                       by = c("FEATUREID"))
  map.catchment$HUC12_BIO <- paste(map.catchment$HUC_12, map.catchment$BIOREGION, sep = "_")
  #==============================================================================
  map.catchment <- agg_huc12_map(map.catchment, map.me, pct_sum, spatial, taxon.rank, todays.date,
                                 agg.col = "FEATUREID", catch.count)
  #==============================================================================
  # Apply rating scheme.
  #map.catchment <- rate.me2(map.catchment, pct_sum, taxon.rank)
  # If the number of stations within a catchment is less than five, then
  # specify that there is an insufficient number of samples to
  # classify the given catchment.
  map.catchment$RATING <- ifelse(map.catchment$COUNT < catch.count, "Insufficient Samples",
                                 as.character(map.catchment$RATING))
  # Make sure the "FEATUREID" column is a character class.
  #map.catchment$FEATUREID <- as.character(map.catchment$FEATUREID)
  # Append the HUC_12 column on to the data frame.
  #final.df <- merge(map.catchment, unique(map.me[, c("FEATUREID", "HUC_12",
  #                                                   "Acres", "SQ_MILE")]),
  #                  by = "FEATUREID")
  # Export the data frame as a .csv table.
  final.df <- map.catchment
  setwd(sub.dir)
  write.csv(final.df, paste(spatial, taxon.rank, todays.date,
                            "HUC12.csv", sep = "_"), row.names = F)
  # End agg_catchment_map function.
  return(final.df)
}
#==============================================================================
#'Generate BIBI Map
#'
#'@export
create_map <- function(huc.df, sub.dir, title.me, spatial, taxon.rank, todays.date){
  # The location of the HUC12 shapefile.
  setwd("//Pike/data/Projects/Chessie_BIBI/CBIBI_Package_July2016/GIS_July2016/ShapeFiles_Environmental")
  # Import the HUC12 shapefile.
  bioregions <- rgdal::readOGR(".", "ChesBay_WBD12_032015", stringsAsFactors = FALSE)
  # Reset the working directory back to the appropriate folder.
  setwd(sub.dir)
  # Remove the Bioregion from the HUC12 Bioregion string to report
  # only the HUC12 Code.
  huc.df$HUC_12 <- gsub( "_.*$", "", huc.df$HUC12_BIO)
  # Joing the HUC12 Shapefile with the data frame containing scores and ratings.
  bio <- sp::merge(bioregions, huc.df, by.x = "HUC12", by.y = "HUC_12", duplicateGeoms = TRUE)
  
  # Assign colors based on the Rating.
  bio$COLORS <- ifelse(bio$HUC12 %in% c("020801010000", "020600010000", "020700111001"), "steelblue1",
                       ifelse(bio$RATING %in% "Excellent", "darkgreen",
                              ifelse(bio$RATING %in% "Good", "green3",
                                     ifelse(bio$RATING %in% "Fair", "yellow2",
                                            ifelse(bio$RATING %in% "Poor", "orange2",
                                                   ifelse(bio$RATING %in% "VeryPoor", "red3",
                                                          ifelse(bio$RATING %in% c("Insufficient Samples", "TBD") |
                                                                   is.na(bio$RATING) , "white", "deeppink")))))))
  
  # Create a .png image of the mapped ratings.
  png(paste(spatial, taxon.rank, todays.date, "rating_map.png", sep = "_"),
      units = "in", res = 720,  width = 6.5, height = 9)
  
  sp::plot(bio, main = title.me, 
       #axes = TRUE,
       #xaxt = 'n', yaxt = 'n',
       #border =" black",
       col = bio@data$COLORS,  #xlim = c(0, 0),
       xlim=c(200000, 400000),
       ylim = c(4020000, 4750000))
  #ylim = c(4050000, 4750000))
  #ylim = c(0, 0))
  
  dev.off()
  #graphics.off()
  # End create_map function.
}

#==============================================================================
#'Map the BIBI Ratings
#'
#'@param my.data =
#'@param index_res = Spatial Resolution (i.e., "BASIN", "COAST", "INLAND", or
#' "BIOREGION")
#'@param bioregion = If index_res is set to "BIOREGION," then a list of 
#'bioregions to keep must be specified.
#'@return Multiple tables related to aggregated scores and ratings, as wel as,
#'a .png image of the map representing the HUC12 ratings.
#'@export

map_this <- function(spatial, calc.date, taxon.rank, all.data, tab.events, todays.date,
                     title.me, very.poor.thresh, random_event, sub.dir, agg.col){
  # Import the appropriate data based on the specified spatial resolution (spatial),
  # taxonomic resolution (taxon.rank), and the date the scores were calculated (calc.date).
  final.scores <- group.me2(spatial, taxon.rank, "scored_metrics", calc.date, final_score = TRUE)
  #==============================================================================
  # Establish the rating thresholds based on the Reference distribution.
  # Outliers in this case reflect ONLY the Refernce distribution.
  pct_sum <- percentile_summary(final.scores, very.poor.thresh, outlier = TRUE)
  # Save the Reference percentile summary as a .csv.
  write.csv(pct_sum, paste(spatial, taxon.rank, todays.date, "Thresholds.csv", sep = "_"), row.names = F)
  #==============================================================================
  # random_events was built specifically to select only randomly sampled events
  # according to the SITE_TYPE_CODE in the Chessie BIBI database.
  # It appears that there are currently (3/20/17) inaccuracies associated with
  # the SITE_TYPE_CODE; therefore, random_events should alwasy be set to FALSE
  # until the randomly selected sampling events can be verified.
  if (random_event == TRUE) rand.events <- random_events(tab.events)
  #==============================================================================
  # Apply the rating scheme.
  all.scores <- rate.me2(final.scores, pct_sum, taxon.rank)
  #==============================================================================
  # Create a new data frame with only event, station, and huc information.
  huc_data <- unique(all.data[, c("EVENT_ID", "STATION_ID", "HUC_8", "HUC_10",
                                  "HUC_12", "ACRES", "SQ_MILE")])
  # Join the scores with the HUC information.
  all.scores <- unique(merge(all.scores, huc_data, by = c("EVENT_ID", "STATION_ID")))
  # Add the HUC12 zero back on to the beginning of the HUC12 number.
  #all.scores$HUC_12 <- paste("0", all.scores$HUC_12, sep = "")
  #==============================================================================
  # Identify which events were considered randomly sampled based on the 
  # SITE_TYPE_CODE. If random_event is NOT True, then fill column with NAs.
  if (random_event == TRUE) {
    all.scores$RANDOM_SAMPLING <- ifelse(all.scores$EVENT_ID %in% rand.events$EVENT_ID, "YES", "NO") 
  } else {
    all.scores$RANDOM_SAMPLING <- NA
  }
  #==============================================================================
  # Organize and export the sampling event ratings.
  all.scores <- organize_all_event_ratings(all.scores, spatial, taxon.rank,
                                           todays.date)
  #==============================================================================
  # Identify outliers, remove outliers, calculate percentiles without outliers for
  # each Reference and Severely Degraded distributions within each bioregion.
  ref_sev.pct <- ref_sev_percentiles(all.scores, spatial, taxon.rank, todays.date)
  #==============================================================================
  # Keep only randomly sampled sampling events.
  if (random_event == TRUE){
    final.scores <- final.scores[final.scores$EVENT_ID %in% rand.events$EVENT_ID, ] 
  } 
  #==============================================================================
  # Create a new data frame containing just event, station, and HUC information.
  huc.data <- unique(all.data[, c("EVENT_ID", "STATION_ID", "HUC_8",
                                  "HUC_10", "HUC_12", "ACRES", "SQ_MILE")])
  #==============================================================================
  # Find the mean score by station, bioregion, and huc.
  map.me <- agg_hucs_map(huc.data, final.scores)
  #==============================================================================
  # Find the mean score by station.
  map.station <- agg_station_map(map.me, pct_sum, spatial, taxon.rank, todays.date)
  #==============================================================================
  if (agg.col %in% "STATION_ID") {
    # Find the mean score of the mean station scores for each HUC12.
    # Exclude HUC12s with less than five unique stations.
    map.huc12 <- agg_huc12_map(map.station, map.me, pct_sum, spatial, taxon.rank,
                               todays.date, agg.col,
                               catch.count = 5)
  }
  if (agg.col %in% "FEATUREID") {
    # Find the mean score of the mean catchment scores for each HUC12.
    # Exclude HUC12s with less than three unique catchments.
    map.huc12 <- agg_catchment_map(map.station, map.me, all.data, pct_sum, spatial, taxon.rank,
                               todays.date, agg.col,
                               catch.count = 3, sub.dir)
  }

  #==============================================================================
  # Create and export the map.
  create_map(map.huc12, sub.dir, title.me, spatial, taxon.rank, todays.date)
  
  # End map_this function.
}