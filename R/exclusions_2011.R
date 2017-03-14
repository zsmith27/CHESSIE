#setwd("C:\\Users\\zsmith\\Desktop\\BIBI\\R\\chessie\\Tables")
#exclusions_2011 <- read.csv("BIBI_2011_Taxa_Exclusions.csv")
#exclusions_2011$TSN <- gsub("(?<![0-9])0+", "", exclusions_2011$TSN, perl = TRUE)
#exclusions_2011 <- exclusions_2011[, !names(exclusions_2011) %in% "X"]

#setwd("C:\\Users\\zsmith\\Desktop\\BIBI\\R\\chessie")
#devtools::use_data(exclusions_2011, overwrite = TRUE)