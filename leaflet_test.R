


map.this.leaflet <- function(spatial, calc.date, taxon.rank, all.data, tab.events, todays.date,
                     title.me, very.poor.thresh, random_event, sub.dir){
  #==============================================================================
  # Create/specify the appropriate subdirectory
  new.dir <- create_subdir(main.dir, spatial, i)
  setwd(new.dir)
  #==============================================================================
  
  final.scores <- group.me2(spatial, taxon.rank, "scored_metrics", calc.date, final_score = TRUE)
  
  
  
  #rand.events <- random_events(tab.events)
  #final.scores <- final.scores[final.scores$EVENT_ID %in% rand.events$EVENT_ID, ]
  pct_sum <- percentile_summary(final.scores, very.poor.thresh, outlier = TRUE)
  write.csv(pct_sum, paste(spatial, taxon.rank, todays.date, "Thresholds.csv", sep = "_"), row.names = F)
  #==============================================================================
  if (random_event == TRUE) rand.events <- random_events(tab.events)
  all.scores <- rate.me2(final.scores, pct_sum, taxon.rank)
  
  huc_data <- unique(all.data[, c("EVENT_ID", "STATION_ID", "HUC_8", "HUC_10", "HUC_12")])
  all.scores <- unique(merge(all.scores, huc_data, by = c("EVENT_ID", "STATION_ID")))
  all.scores$HUC_12 <- paste("0", all.scores$HUC_12, sep = "")
  
  if (random_event == TRUE) {
    all.scores$RANDOM_SAMPLING <- ifelse(all.scores$EVENT_ID %in% rand.events$EVENT_ID, "YES", "NO") 
  } else {
    all.scores$RANDOM_SAMPLING <- NA
  }
  all.scores <- all.scores[, c("BIOREGION", "EVENT_ID", "STATION_ID", "HUC_8", "HUC_10", "HUC_12",
                               "DATE", "SAMPLE_NUMBER", "AGENCY_CODE", "CATEGORY",
                               "RANDOM_SAMPLING", "FINAL_SCORE", "DEG_50",
                               "REF_10", "REF_25", "REF_50", "RATING")]
  names(all.scores)[names(all.scores) %in% "DEG_50"] <- "HALF_REF_10"
  #all.scores[, 12:16] <- round(all.scores[, 12:16] * 100, 2)
  all.scores[, 12:16] <- all.scores[, 12:16] * 100
  all.scores$CATEGORY <- factor(all.scores$CATEGORY, levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
  all.scores <- all.scores[order(all.scores$BIOREGION, all.scores$CATEGORY, -all.scores$FINAL_SCORE), ]
  write.csv(all.scores, paste(spatial, taxon.rank, todays.date, "All_Event_Ratings.csv", sep = "_"), row.names = F)
  
  new.df <- data.frame(BIOREGION = unique(all.scores$BIOREGION))
  data.list <- list()
  for (i in unique(new.df$BIOREGION)) {
    sub.scores <- all.scores[all.scores$BIOREGION %in% i, ]
    quick.percentiles <- function(scores.df, bioregion, category) {
      outlier <- function(scores.df, mult.val, job = "REF"){
        datalist = list()
        for (i in unique(scores.df$BIOREGION)) {
          bioregion.df <- scores.df[scores.df$BIOREGION %in% i, ]
          
          quantile.vec <- quantile(bioregion.df$FINAL_SCORE, c(0.25, 0.75))
          if (job %in% "REF"){
            outlier.thresh <- quantile.vec[1] - (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
          }
          if (job %in% "SEV"){
            outlier.thresh <- quantile.vec[2] + (abs(quantile.vec[1] - quantile.vec[2]) * mult.val)
          }
          
          outlier.df <- data.frame(BIOREGION = i)
          outlier.df$THRESH <- outlier.thresh
          datalist[[i]] <- outlier.df # add it to your list
        }
        final.df <- do.call(rbind, datalist)
        return(final.df)
      }
      
      sub.df <- scores.df[scores.df$CATEGORY %in% category, ]
      out.df <- outlier(sub.df, 1.5, job = category)
      if (category %in% "REF") sub.df <- sub.df[sub.df$FINAL_SCORE >= out.df$THRESH, ]
      if (category %in% "SEV") sub.df <- sub.df[sub.df$FINAL_SCORE <= out.df$THRESH, ]
      sub.df <- data.frame(t(quantile(sub.df$FINAL_SCORE, c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95))))
      sub.df$BIOREGION <- bioregion
      sub.df$CATEGORY <- category
      sub.df <- sub.df[, c(7, 8, 1:6)]
      names(sub.df) <- c("Bioregion", "Category", "5%", "10%", "25%", "50%", "75%", "95%")
      return(sub.df)
    }
    
    sub.ref <- quick.percentiles(sub.scores, i, "REF")
    sub.sev <- quick.percentiles(sub.scores, i, "SEV")
    sub.final <- rbind(sub.ref, sub.sev)
    data.list[[i]] <- sub.final
  }
  ref_sev.pct <- do.call(rbind, data.list)
  ref_sev.pct[, 3:8] <- round(ref_sev.pct[, 3:8], 2)
  write.csv(ref_sev.pct, paste(spatial, taxon.rank, todays.date, "REF_SEV_Score_Percentiles.csv", sep = "_"), row.names = F)
  #==============================================================================
  if (random_event == TRUE) {
    rand.events <- random_events(tab.events)
    final.scores <- final.scores[final.scores$EVENT_ID %in% rand.events$EVENT_ID, ]
  }
  
  #==============================================================================
  #bio_huc <- read.csv("HUC12_BIOREGIONS.csv")
  huc_data <- unique(all.data[, c("EVENT_ID", "STATION_ID", "HUC_8", "HUC_10", "HUC_12")])
  #huc_data <- merge(huc_data, bio_huc, by = c("STATION_ID"))
  #huc_data <- unique(huc_data[, c("EVENT_ID", "STATION_ID", "HUC_8", "HUC_10", "HUC_12", "HUC12BIOREGION_ID")])
  #==============================================================================
  map.me <- unique(merge(huc_data, final.scores, by = c("EVENT_ID", "STATION_ID"), all.y = T))
  map.me$HUC12_BIO <- paste("0", map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  huc12 <- data.frame(unique(map.me$HUC12_BIO))
  # Find the mean score per station
  map.me <- aggregate(FINAL_SCORE ~ STATION_ID + HUC_8 + HUC_10 + HUC_12 + HUC12_BIO, data = map.me, mean)
  #test <- data.frame(huc12[!huc12$unique.map.me.HUC_12. %in% map.me$HUC_12,])
  #==============================================================================
  map.me <- unique(merge(huc_data, final.scores, by = c("EVENT_ID", "STATION_ID"), all.x = T,  all.y = T))
  map.me <- map.me[!is.na(map.me$BIOREGION), ]
  map.me$HUC12_BIO <- paste("0", map.me$HUC_12, "_", map.me$BIOREGION, sep = "")
  #==============================================================================
  map.station <- aggregate(FINAL_SCORE ~ STATION_ID + BIOREGION, data = map.me, mean)
  map.station.length <- aggregate(FINAL_SCORE ~ STATION_ID, data = map.me, length)
  names(map.station.length) <- c("STATION_ID", "COUNT")
  map.station <- merge(map.station, map.station.length, by = "STATION_ID")
  map.station <- rate.me2(map.station, pct_sum, taxon.rank)
  #map.station <- rate.me(map.station, index = 2016, "ORDER")
  #map.station2 <- rate.me(map.station, index = 2011)
  write.csv(map.station, paste(spatial, taxon.rank, todays.date, "STATION.csv", sep = "_"), row.names = F)
  #==============================================================================
  merged.huc12 <- merge(map.station, unique(map.me[, c("STATION_ID", "HUC12_BIO")]), by = "STATION_ID")
  
  map.huc12 <- aggregate(FINAL_SCORE ~ HUC12_BIO + BIOREGION, data = merged.huc12,  mean)
  map.huc12.length <- aggregate(FINAL_SCORE ~ HUC12_BIO, data = merged.huc12,  length)
  names(map.huc12.length) <- c("HUC12_BIO", "COUNT")
  map.huc12 <- merge(map.huc12, map.huc12.length, by = c("HUC12_BIO"))
  map.huc12 <- rate.me2(map.huc12, pct_sum, taxon.rank)
  #map.huc12 <- rate.me(map.huc12, index = 2016, "GENUS")
  map.huc12$RATING <- ifelse(map.huc12$COUNT < 5, "Insufficient Samples", as.character(map.huc12$RATING))
  #map.huc12$HUC12_BIO<- as.character(paste("0", map.huc12$HUC12_BIO, sep = ""))
  write.csv(map.huc12, paste(spatial, taxon.rank, todays.date, "NEW_HUC12.csv", sep = "_"), row.names = F)
  #==============================================================================
  #map.huc12 <- aggregate(FINAL_SCORE ~ HUC12_BIO + BIOREGION, data = map.me,  mean)
  #map.huc12.length <- aggregate(FINAL_SCORE ~ HUC12_BIO, data = map.me,  length)
  #names(map.huc12.length) <- c("HUC12_BIO", "COUNT")
  #map.huc12 <- merge(map.huc12, map.huc12.length, by = c("HUC12_BIO"))
  #map.huc12 <- rate.me2(map.huc12, pct_sum, taxon.rank)
  #map.huc12 <- rate.me(map.huc12, index = 2016, "GENUS")
  #map.huc12$RATING <- ifelse(map.huc12$COUNT < 5, "Insufficient Samples", as.character(map.huc12$RATING))
  #map.huc12$HUC12_BIO<- as.character(paste("0", map.huc12$HUC12_BIO, sep = ""))
  #write.csv(map.huc12, paste(spatial, taxon.rank, todays.date, "HUC12.csv", sep = "_"), row.names = F)
  #==============================================================================
  
  
  fam <- map.huc12
  fam$HUC12_BIO <- as.character(fam$HUC12_BIO)
  fam <- merge(fam, unique(map.me[, c("HUC12_BIO", "HUC_12")]), by = "HUC12_BIO")
  
  
  
  library(maps)
  library(ggplot2)
  library(mapdata)
  library(rgdal)
  print("readOGR")
  setwd("//Pike/data/Projects/Chessie_BIBI/CBIBI_Package_July2016/GIS_July2016/ShapeFiles_Environmental")
  bioregions <- rgdal::readOGR(".", "ChesBay_WBD12_032015", stringsAsFactors = FALSE)
  bioregions <- spTransform(bioregions, CRS("+init=epsg:4326"))
  setwd("D:/ZSmith/Projects/Chessie_BIBI/GIS/spatial_02")
  catchments <- rgdal::readOGR(".", "Catchments02", stringsAsFactors = FALSE)
  catchments <- spTransform(catchments, CRS("+init=epsg:4326"))
  #Set Working directory (The folder files are imported and exported to)
  #setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
  test <- point.in.polygon(stations.df$LATITUDE, stations.df$LONGITUDE, catchments@polygons)
  
  test <- over(stations.df, catchments@polygons)
  # Assignment modified according
  coordinates(stations.df) <- ~ LONGITUDE + LATITUDE
  # Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
  proj4string(stations.df) <- proj4string(catchments)
  
  test <- over(stations.df, catchments)
  stations.df$TEST <- over(stations.df, catchments)$FEATUREID
  test2 <- catchments[catchments$FEATUREID %in% test$FEATUREID, ]
  
  setwd(sub.dir)
  print("bioregions.arcgis as dataframe")
  #bioregions <- data.frame(bioregions.arcgis)
  #bioregions$HUC12 <- as.numeric(as.character(bioregions$HUC12))
  print("fix HUC12")
  fam$HUC_12 <- gsub( "_.*$", "", fam$HUC12_BIO)
  print("Merge")
  bio <- sp::merge(bioregions, fam, by.x = "HUC12", by.y = "HUC_12", duplicateGeoms=TRUE)
  print("Colors")
  bio$COLORS <- ifelse(bio$HUC12 %in% c(20801010000, 20600010000, 20700111001), "steelblue1",
                       ifelse(bio$RATING %in% "Excellent", "darkgreen",
                              ifelse(bio$RATING %in% "Good", "green3",
                                     ifelse(bio$RATING %in% "Fair", "yellow2",
                                            ifelse(bio$RATING %in% "Poor", "orange2",
                                                   ifelse(bio$RATING %in% "VeryPoor", "red3",
                                                          ifelse(bio$RATING %in% c("Insufficient Samples", "TBD") |
                                                                   is.na(bio$RATING) , "white", "deeppink")))))))
  print("Colors2")
  bio$COLORS2 <- ifelse(bio$HUC12 %in% c(20801010000, 20600010000, 20700111001), "transparent",
                        ifelse(bio$RATING %in% "Excellent", "#007505",
                               ifelse(bio$RATING %in% "Good", "#24FF0A",
                                      ifelse(bio$RATING %in% "Fair", "#FFEC06",
                                             ifelse(bio$RATING %in% "Poor", "#FF8E11",
                                                    ifelse(bio$RATING %in% "VeryPoor", "#FF0000",
                                                           ifelse(bio$RATING %in% c("Insufficient Samples", "TBD") |
                                                                    is.na(bio$RATING) , "transparent", "#000000")))))))
  
  
  
  
  
  bio$COLORS <- as.factor(bio$COLORS2)
  
  stations.df <- merge(map.station, unique(all.data[, c("STATION_ID", "LATITUDE", "LONGITUDE", "HUC_12")]), by = "STATION_ID")
  stations.df$COLORS2 <- ifelse(stations.df$RATING %in% "Excellent", "#007505",
                               ifelse(stations.df$RATING %in% "Good", "#24FF0A",
                                      ifelse(stations.df$RATING %in% "Fair", "#FFEC06",
                                             ifelse(stations.df$RATING %in% "Poor", "#FF8E11",
                                                    ifelse(stations.df$RATING %in% "VeryPoor", "#FF0000",
                                                           ifelse(stations.df$RATING %in% c("Insufficient Samples", "TBD") |
                                                                    is.na(stations.df$RATING) , "transparent", "#000000"))))))
  
  huc.vec <- c("020501030605", "020600020202", "020700010301", "020700020205", "020700020302", "020700020403",
               "020700030102", "020700030501", "020700030801", "020700100805")
  
  highlight.huc <- bio[bio$HUC12 %in% huc.vec, ]
  
  library(leaflet)
  #bioregions <- readOGR("ChesBay_WBD12_032015.shp")
  leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery) %>%
    #addTiles() %>%
    #addPolygons(data = bio, weight = 0.5, color = "#000000", fillColor = ~bio$COLORS2, fillOpacity = 0.3,
    #            popup = ~as.character(bio$HUC12)) %>%
    #addPolygons(data = highlight.huc, weight = 4, color = "#000000", fillColor = ~highlight.huc$COLORS2, fillOpacity = 0.3,
    #            popup = ~as.character(highlight.huc$HUC12)) %>%
    addPolygons(data = test2, weight = 0.5, color = "#000000", #fillColor = ~bio$COLORS2, fillOpacity = 0.3,
                popup = ~as.character(test2$FEATUREID)) %>%
    addCircleMarkers(data = stations.df, lng = stations.df$LONGITUDE, lat = stations.df$LATITUDE,
                     stroke = TRUE, weight = 2, color = "#000000",  radius = 5,
                     fillOpacity = 1.0, fillColor = ~stations.df$COLORS2,
                     popup = ~as.character(stations.df$STATION_ID))
  
  #addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
  #            opacity = 1.0, fillOpacity = 0.5,
  #            #fillColor = ~colorQuantile("YlOrRd", ALAND)(ALAND),
  #            highlightOptions = highlightOptions(color = "white", weight = 2,
  #                                                bringToFront = TRUE))
  
  
  
}






library(maps)
library(ggplot2)
library(mapdata)
library(rgdal)
print("readOGR")
setwd("//Pike/data/Projects/Chessie_BIBI/CBIBI_Package_July2016/GIS_July2016/ShapeFiles_Environmental")
bioregions <- rgdal::readOGR(".", "ChesBay_WBD12_032015", stringsAsFactors = FALSE)
bioregions <- spTransform(bioregions, CRS("+init=epsg:4326"))
setwd("D:/ZSmith/Projects/Chessie_BIBI/GIS/spatial_02")
catchments <- rgdal::readOGR(".", "Catchments02", stringsAsFactors = FALSE)
catchments <- spTransform(catchments, CRS("+init=epsg:4326"))
#Set Working directory (The folder files are imported and exported to)
#setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
stations.df <- merge(map.station, unique(all.data[, c("STATION_ID", "LATITUDE", "LONGITUDE", "HUC_12")]), by = "STATION_ID")
stations.df$COLORS2 <- ifelse(stations.df$RATING %in% "Excellent", "#007505",
                              ifelse(stations.df$RATING %in% "Good", "#24FF0A",
                                     ifelse(stations.df$RATING %in% "Fair", "#FFEC06",
                                            ifelse(stations.df$RATING %in% "Poor", "#FF8E11",
                                                   ifelse(stations.df$RATING %in% "VeryPoor", "#FF0000",
                                                          ifelse(stations.df$RATING %in% c("Insufficient Samples", "TBD") |
                                                                   is.na(stations.df$RATING) , "transparent", "#000000"))))))

huc.vec <- c("020501030605", "020600020202", "020700010301", "020700020205", "020700020302", "020700020403",
             "020700030102", "020700030501", "020700030801", "020700100805")

highlight.huc <- bio[bio$HUC12 %in% huc.vec, ]

coordinates(stations.df) <- ~ LONGITUDE + LATITUDE
# Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
proj4string(stations.df) <- proj4string(catchments)

test <- over(stations.df, catchments)
catchments <- catchments[catchments$FEATUREID %in% test$FEATUREID, ]
catchments@data$COLORS2 <- over(stations.df, catchments)$COLORS2

stations.df@data$FEATUREID <- over(stations.df, catchments)$FEATUREID

test <- aggregate(FINAL_SCORE ~ FEATUREID + BIOREGION, data = stations.df, FUN = mean)
map.catch <- rate.me2(test, pct_sum, taxon.rank)
map.catch$COLORS2 <- ifelse(map.catch$RATING %in% "Excellent", "#007505",
                              ifelse(map.catch$RATING %in% "Good", "#24FF0A",
                                     ifelse(map.catch$RATING %in% "Fair", "#FFEC06",
                                            ifelse(map.catch$RATING %in% "Poor", "#FF8E11",
                                                   ifelse(map.catch$RATING %in% "VeryPoor", "#FF0000",
                                                          ifelse(map.catch$RATING %in% c("Insufficient Samples", "TBD") |
                                                                   is.na(map.catch$RATING) , "transparent", "#000000"))))))

test2 <- sp::merge(catchments, map.catch[, c("COLORS2", "FEATUREID")], by = "FEATUREID", duplicateGeoms=TRUE)

setwd(sub.dir)
print("bioregions.arcgis as dataframe")
#bioregions <- data.frame(bioregions.arcgis)
#bioregions$HUC12 <- as.numeric(as.character(bioregions$HUC12))
print("fix HUC12")
fam$HUC_12 <- gsub( "_.*$", "", fam$HUC12_BIO)
print("Merge")
bio <- sp::merge(bioregions, fam, by.x = "HUC12", by.y = "HUC_12", duplicateGeoms=TRUE)
print("Colors")
bio$COLORS <- ifelse(bio$HUC12 %in% c(20801010000, 20600010000, 20700111001), "steelblue1",
                     ifelse(bio$RATING %in% "Excellent", "darkgreen",
                            ifelse(bio$RATING %in% "Good", "green3",
                                   ifelse(bio$RATING %in% "Fair", "yellow2",
                                          ifelse(bio$RATING %in% "Poor", "orange2",
                                                 ifelse(bio$RATING %in% "VeryPoor", "red3",
                                                        ifelse(bio$RATING %in% c("Insufficient Samples", "TBD") |
                                                                 is.na(bio$RATING) , "white", "deeppink")))))))
print("Colors2")
bio$COLORS2 <- ifelse(bio$HUC12 %in% c(20801010000, 20600010000, 20700111001), "transparent",
                      ifelse(bio$RATING %in% "Excellent", "#007505",
                             ifelse(bio$RATING %in% "Good", "#24FF0A",
                                    ifelse(bio$RATING %in% "Fair", "#FFEC06",
                                           ifelse(bio$RATING %in% "Poor", "#FF8E11",
                                                  ifelse(bio$RATING %in% "VeryPoor", "#FF0000",
                                                         ifelse(bio$RATING %in% c("Insufficient Samples", "TBD") |
                                                                  is.na(bio$RATING) , "transparent", "#000000")))))))





bio$COLORS <- as.factor(bio$COLORS2)



library(leaflet)
#bioregions <- readOGR("ChesBay_WBD12_032015.shp")
leaflet() %>%
  #addProviderTiles(providers$Esri.WorldImagery) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
  #addTiles() %>%
  addPolygons(data = bio, weight = 0.5, color = "#000000", fillColor = "transparent", #fillOpacity = 0.3,
              popup = ~as.character(bio$HUC12)) %>%
  #addPolygons(data = highlight.huc, weight = 4, color = "#000000", fillColor = ~highlight.huc$COLORS2, fillOpacity = 0.3,
  #            popup = ~as.character(highlight.huc$HUC12)) %>%
  addPolygons(data = test2, weight = 0.5, color = "#000000", fillColor = test2$COLORS2, fillOpacity = 0.5,
              popup = ~as.character(test2$FEATUREID)) %>%
  addCircleMarkers(data = stations.df, lng = stations.df$LONGITUDE, lat = stations.df$LATITUDE,
                   stroke = TRUE, weight = 2, color = "#000000",  radius = 5,
                   fillOpacity = 1.0, fillColor = ~stations.df$COLORS2,
                   popup = ~as.character(stations.df$STATION_ID))

#addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
#            opacity = 1.0, fillOpacity = 0.5,
#            #fillColor = ~colorQuantile("YlOrRd", ALAND)(ALAND),
#            highlightOptions = highlightOptions(color = "white", weight = 2,
#                                                bringToFront = TRUE))




library(leaflet)
#bioregions <- readOGR("ChesBay_WBD12_032015.shp")
leaflet() %>%
  #addProviderTiles(providers$Esri.WorldImagery) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
  #addTiles() %>%
  addPolygons(data = bio, weight = 0.5, color = "#000000", fillColor = bio$COLORS2, fillOpacity = 1,
              popup = ~as.character(bio$HUC12)) %>%
  #addPolygons(data = highlight.huc, weight = 4, color = "#000000", fillColor = ~highlight.huc$COLORS2, fillOpacity = 0.3,
  #            popup = ~as.character(highlight.huc$HUC12)) %>%
  #addPolygons(data = test2, weight = 0.5, color = "#000000", fillColor = test2$COLORS2, fillOpacity = 0.5,
  #            popup = ~as.character(test2$FEATUREID)) %>%
  addCircleMarkers(data = stations.df, lng = stations.df$LONGITUDE, lat = stations.df$LATITUDE,
                   stroke = TRUE, weight = 2, color = "#000000",  radius = 5,
                   fillOpacity = 1.0, fillColor = ~stations.df$COLORS2,
                   popup = ~as.character(stations.df$STATION_ID))



















test <- stations.df[stations.df$HUC_12 %in% "20700010304", ]
mean(test$FINAL_SCORE)



spatial <- "REGION"
calc.date <- "3_17_17_500i"
todays_date = format(Sys.time(), "%m_%d_%y")
todays.date <- todays_date
random_event = FALSE 
main.dir <- "D:/ZSmith/Projects/Chessie_BIBI/Output"
taxon.rank <- "FAMILY"
very.poor.thresh <- "half_10"
all.data <- prep.data
i <- "FAMILY"
sub.dir <- create_subdir(main.dir, "REGION", "FAMILY")
title.me <- "TEST"
m.c
new.dir <- create_subdir(main.dir, spatial, i)
setwd(new.dir)
tab.event <- tab_event
setwd("D:/ZSmith/Projects/Chessie_BIBI/Output/March_2017/03_17_2017/REGION/FAMILY")

map.this.leaflet(spatial, calc.date, taxon.rank = i, all.data, tab.event, todays_date,
         title.me = paste(i, "-level ", spatial, " Ratings (Random)", sep = ""),
         very.poor.thresh = "half_10", random_event, sub.dir = new.dir)
