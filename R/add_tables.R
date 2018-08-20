#Set Working directory (The folder files are imported and exported to)
#setwd("//pike/data/Projects/Chessie_BIBI/BIBI_June_2016")
#metric.class <- read.csv("Metrics_List_9_19_16.csv")
#m.c <- metric.class[, c(1, 3)]
#names(m.c) <- c("METRICS", "METRIC_CLASS")
setwd("C:/Users/zsmith/Desktop/BIBI/R/GitHub/CHESSIE")
m.c <- read.csv("metric_class.csv", stringsAsFactors = FALSE)
#==============================================================================
#Set Working directory (The folder files are imported and exported to)
setwd("C:/Users/zsmith/Desktop/BIBI/R/GitHub/CHESSIE")
devtools::use_data(m.c, overwrite = TRUE)