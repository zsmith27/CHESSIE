#==============================================================================
# Check for Taxonomic Resolution Discrepancies Between Agencies/Programs
#==============================================================================
#'Prepare data for checking taxonomic resolution discrepancies between agencies/programs
#'
#'@param long = long data frame containing taxonomic counts.
#'@return Concatonates AGENCY_CODE and PROGRAM_CODE, and replaces AGENCY_CODE
#'in the long data frame with the new concatonation.  This allows the data
#'to be aggregated by unique a unique agency/program.
#'@export

prep_agency_taxa_res <- function(long){
  keep <- c("EVENT_ID", "AGENCY_CODE", "PROGRAM_CODE")
  keep.df <- long[, keep]
  # Concatonate AGENCY_CODE and PROGRAM_CODE
  keep.df$AP <- paste(keep.df$AGENCY_CODE, keep.df$PROGRAM_CODE, sep ="_")
  # Keep only unique instances of EVENT_ID and the concatonation
  keep.df <- unique(keep.df[, c("EVENT_ID", "AP")])
  # Merge the orignal data set with new concatonated strings
  final.df <- merge(keep.df, long, by = "EVENT_ID", all.y = T)
  # Remove the old Agency_Code and replace it with the concatonated strings
  final.df <- final.df[, !(names(final.df) %in% "AGENCY_CODE")]
  colnames(final.df)[names(final.df) == "AP"] <- "AGENCY_CODE"
  # Remove leading and trailing blanks
  rm_blanks <- function(Long){
    gsub("^\\s+|\\s+$", "", Long)
  }
  final.df$AGENCY_CODE <- rm_blanks(final.df$AGENCY_CODE)
  return(final.df)
}


#==============================================================================
#'Taxa Absent from Agency/Program
#'
#'@param long = long data frame containing taxonomic counts.
#'@param rank = Taxonomic rank used to aggregate the data.
#'@param col.ord = a list of names used to sort the final output.
#'@return Returns a data frame containing taxa that are observed by >= 1 
#'agency/program and not observed by >= 1 agency/program. 
#'@export

taxa_absence <- function(long, rank, col.ord = NULL){
  wide.df <- BIBI::wide(long, rank)
  agg <- aggregate(wide.df[, 6:ncol(wide.df)], by = list(wide.df$AGENCY_CODE), FUN = sum)
  agg2 <- agg
  agg[, 2:ncol(agg)] <- ifelse(agg[, 2:ncol(agg)] > 0, 1, 0)
  cols_to_drop <- c(rep(TRUE, 1), colSums(agg[, 2:ncol(agg)]) < nrow(agg) &
                      colSums(agg[, 2:ncol(agg)]) > 0)
  tricky <- names(agg[, cols_to_drop])
  final.df <- agg2[, tricky]
  if(!is.null(col.ord)){
    final.df <- final.df[match(col.ord, final.df[, 1]), ]
    names(final.df)[1] <- "Agency Code and Program"
    final.df
  } 

  return(final.df)
}

#==============================================================================
#'Agency Taxonomic Differences
#'
#'@param long = long data frame containing taxonomic counts.
#'@param rank = Taxonomic rank used to aggregate the data.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Returns a data frame containing taxa that are observed by >= 1 
#'agency/program and not observed by >= 1 agency/program. This data frame
#'also contains the percentage of data that each agency/program composes
#'in the bioregion being investigated.
#'@export

agency_taxa_diff2 <- function(long, bioregion, rank){
  agency.count <- agency_count(long, bioregion)
  prep.bio <- prep_bioregion(long)
  bioregion.df <- prep.bio[prep.bio$BIOREGION %in% bioregion, ]
  absent.taxa <- taxa_absence(bioregion.df, rank, agency.count$`Agency Code and Program`)
  names(absent.taxa)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, absent.taxa, by = "Agency Code and Program")
  final.name <- names(final.df[, 4:ncol(final.df)])
  corrected.names <- paste0(toupper(substr(final.name, 1, 1)), 
                            tolower(substr(final.name, 2, nchar(final.name))))
  names(final.df)[4:ncol(final.df)] <- corrected.names
  return(final.df)
}

#==============================================================================
#'Write csv files containing agency/program taxonomic discrepancies
#'
#'@param long = long data frame containing taxonomic counts.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Creates multiple data frames containing taxa that are observed by >= 1 
#'agency/program and not observed by >= 1 agency/program at the phylum,
#' subphylum, class, subclass, and order levels. These data frames
#'also contains the percentage of data that each agency/program composes
#'in the bioregion being investigated.
#'@export

output_agency_taxa_diff2 <- function(long, bioregion){
  atd <- agency_taxa_diff2(long, bioregion, "PHYLUM")
  write.csv(atd, "phylum_ridges.csv")
  write.csv(atd, file = paste(c(bioregion, "_phylum.csv"), collapse = ""))
  atd <- agency_taxa_diff2(long, bioregion, "SUBPHYLUM")
  write.csv(atd, file = paste(c(bioregion, "_subphylum.csv"), collapse = ""))
  atd <- agency_taxa_diff2(long, bioregion, "CLASS")
  write.csv(atd, file = paste(c(bioregion, "_class.csv"), collapse = ""))
  atd <- agency_taxa_diff2(long, bioregion, "SUBCLASS")
  write.csv(atd, file = paste(c(bioregion, "_subclass.csv"), collapse = ""))
  atd <- agency_taxa_diff2(long, bioregion, "ORDER")
  write.csv(atd, file = paste(c(bioregion, "_order.csv"), collapse = ""))
}

#==============================================================================
#'Data frame of Taxa of Interest
#'
#'@param long = long data frame containing taxonomic counts.
#'@param low.rank = The lower of taxonomic rank of interest.
#'@param upper.rank = The upper of taxonomic rank of interest; corresponds with
#'the taxon of interest.
#'@param taxon = The taxon of interest.
#'@return Returns a data frame at the taxonomic level specified by low.rank,
#' containing only lower level taxa of the taxon of interest.
#'@export

group.names <- function(long, low.rank, upper.rank, taxon){
  taxa.group <- split(long[, low.rank], long[, upper.rank])
  new.df <- unique(data.frame(taxa.group[taxon]))
  new.df <- new.df[, !(names(new.df) %in% "UNIDENTIFIED")]
  return(new.df)
}

#==============================================================================
#'Test Agency/Program Standard Taxonomic Resolution of Specific Taxa
#'
#'@param long = long data frame containing taxonomic counts.
#'@param low.rank = The lower of taxonomic rank of interest.
#'@param upper.rank = The upper of taxonomic rank of interest; corresponds with
#'the taxon of interest.
#'@param taxon = The taxon of interest.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Returns a list of names at the taxonomic level specified by low.rank,
#' aggregated by the taxon of interest. The data frame
#'also contains the percentage of data that each agency/program composes
#'in the bioregion being investigated.
#'@export

test_id <- function(long, low.rank, upper.rank, taxon, bioregion){
  agency.count <- agency_count(long, bioregion)
  #test <- agency_taxa_diff2(long, low.rank)
  taxa.names <- group.names(long, low.rank, upper.rank, taxon)
  if(length(taxa.names) > 0) taxa.names <- droplevels.factor(taxa.names[taxa.names != "UNIDENTIFIED"])
  wide.df <- BIBI::wide(long, low.rank)
  agg <- aggregate(wide.df[, 6:ncol(wide.df)], by = list(wide.df$AGENCY_CODE), FUN = sum)
  taxa.df <- agg[, c(names(agg) %in% taxa.names)]
  agency.taxa <- data.frame(cbind(agg[, 1], taxa.df), stringsAsFactors = F)
  names(agency.taxa)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, agency.taxa, by = "Agency Code and Program")
  if(ncol(final.df) == 4) return(as.numeric(as.character(final.df[, 4])))
  if(ncol(final.df) > 4) return(rowSums(final.df[, 4:ncol(final.df)]))
  if(ncol(final.df) < 4) return(as.numeric(as.character(0)))
}


#==============================================================================
#'Test Agency/Program Standard Taxonomic Resolution
#'
#'@param long = long data frame containing taxonomic counts.
#'@param col.ord = a list of names used to sort the final output.
#'@param taxon = The taxon of interest.
#'@return Returns the percentage of taxa identified by each agency/program
#'at the specified taxonomic rank.
#'@export


agency_resolution <- function(long, col.ord, rank){
  wide.df <- BIBI::wide(long, rank)
  t2 <- aggregate(UNIDENTIFIED ~ AGENCY_CODE, data = wide.df, sum)
  wide.df$SUM <- rowSums(wide.df[, 6:ncol(wide.df)])
  t3 <- aggregate(SUM ~ AGENCY_CODE, data = wide.df, sum)
  merged <- merge(t2, t3, by = "AGENCY_CODE")
  merged$'% Identified' <- (1 - (merged$UNIDENTIFIED / merged$SUM)) * 100
  final.df <- merged[match(col.ord$`Agency Code and Program`, merged[, 1]), ]
  names(final.df) <- c("Agency Code and Program", "# Unidentified",
                       "Total # of Taxa", "% Identified")
  return(final.df)
}

#==============================================================================
#'Write csv files of Agency/Program Standard Taxonomic Resolution
#'
#'@param long = long data frame containing taxonomic counts.
#'@param col.ord = a list of names used to sort the final output.
#'@param taxon = The taxon of interest.
#'@return Writes multiple csv files containing the percentage of taxa
#' identified by each agency/program at the order, family and genus
#' level.
#'@export

ag_res_output <- function(long, bioregion){
  bioregion.df <- long[long$BIOREGION %in% bioregion, ]
  ord <- wide(long , "ORDER")
  new <- data.frame(table(ord$AGENCY_CODE))
  names(new) <- c("Agency Code and Program", "Total Count")
  new$'Percentage of Data' <- round((new$`Total Count` / sum(new$`Total Count`)) * 100, 2)
  new <- new[order(new$`Total Count`), ]
  ord.df <- agency_resolution(long, new, "ORDER")
  write.csv(ord.df, file = paste(c(bioregion, "_ord_res.csv"), collapse = ""))
  fam.df <- agency_resolution(long, new, "FAMILY")
  write.csv(fam.df, file = paste(c(bioregion, "_fam_res.csv"), collapse = ""))
  gen.df <- agency_resolution(long, new, "GENUS")
  write.csv(gen.df, file = paste(c(bioregion, "_gen_res.csv"), collapse = ""))
}


#==============================================================================
#'Agency/Program count and percentage per bioregion
#'
#'@param long = long data frame containing taxonomic counts.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return The number and percentage of samples that each agency/bioregion
#' contributes to a bioregion of interest.
#'@export

agency_count <- function(long, bioregion){
  prep.bio <- prep_bioregion(long)
  unique(prep.bio$BIOREGION)
  bioregion.df <- prep.bio[prep.bio$BIOREGION %in% bioregion, ]
  phy <- wide(bioregion.df, "PHYLUM")
  agency.count <- data.frame(table(phy$AGENCY_CODE))
  names(agency.count) <- c("Agency Code and Program", "Total Count (N)")
  agency.count$'% of Data' <- round((agency.count$`Total Count` / sum(agency.count$`Total Count`)) * 100, 2)
  agency.count <- agency.count[order(agency.count$`Total Count`), ]
  return(agency.count)
}


#==============================================================================
#'Checks for common Agency/Program Discrepancies
#'
#'@param long = long data frame containing taxonomic counts.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Checks for common agency/program discrepancies at the order, family,
#'and genus level.
#'@export



common_agency_discrepancies <- function(long, bioregion){
  agency.count <- agency_count(long, bioregion)
  prep.bio <- prep_bioregion(long)
  bioregion.df <- prep.bio[prep.bio$BIOREGION %in% bioregion, ]
  agency.df <- data.frame(AGENCY_CODE = sort(unique(bioregion.df$AGENCY_CODE)))
  agency.df$O_Chelicerata <- test_id(bioregion.df, "ORDER", "SUBPHYLUM", "CHELICERATA", bioregion)
  agency.df$O_Gastropoda <- test_id(bioregion.df, "ORDER", "CLASS", "GASTROPODA", bioregion)
  agency.df$O_Hirudinea <- test_id(bioregion.df, "ORDER", "SUBCLASS", "HIRUDINEA", bioregion)
  agency.df$O_Oligochaeta <- test_id(bioregion.df, "ORDER", "CLASS", "OLIGOCHAETA", bioregion)
  agency.df$O_Nemertea <- test_id(bioregion.df, "ORDER", "PHYLUM", "NEMERTEA", bioregion)
  
  agency.df$F_Chelicerata <- test_id(bioregion.df, "FAMILY", "SUBPHYLUM", "CHELICERATA", bioregion)
  agency.df$F_Gastropoda <- test_id(bioregion.df, "FAMILY", "CLASS", "GASTROPODA", bioregion)
  agency.df$F_Hirudinea <- test_id(bioregion.df, "FAMILY", "SUBCLASS", "HIRUDINEA", bioregion)
  agency.df$F_Oligochaeta <- test_id(bioregion.df, "FAMILY", "CLASS", "OLIGOCHAETA", bioregion)
  agency.df$F_Nemertea <- test_id(bioregion.df, "FAMILY", "PHYLUM", "NEMERTEA", bioregion)
  
  agency.df$G_Chelicerata <- test_id(bioregion.df, "GENUS", "SUBPHYLUM", "CHELICERATA", bioregion)
  agency.df$G_Gastropoda <- test_id(bioregion.df, "GENUS", "CLASS", "GASTROPODA", bioregion)
  agency.df$G_Hirudinea <- test_id(bioregion.df, "GENUS", "SUBCLASS", "HIRUDINEA", bioregion)
  agency.df$G_Oligochaeta <- test_id(bioregion.df, "GENUS", "CLASS", "OLIGOCHAETA", bioregion)
  agency.df$G_Nemertea <- test_id(bioregion.df, "GENUS", "PHYLUM", "NEMERTEA", bioregion)
  agency.df$G_Chironomidae <- test_id(bioregion.df, "GENUS", "FAMILY", "CHIRONOMIDAE", bioregion)
  
  names(agency.df)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, agency.df, by = "Agency Code and Program")
  #return(final.df)
  write.csv(final.df, file = paste(c(bioregion, "_agency_res_exclusions.csv"), collapse = ""))
}

#==============================================================================
#'Agency/Program count and percentage per bioregion
#'
#'@param long = long data frame containing taxonomic counts.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return The number and percentage of samples that each agency/bioregion
#' contributes to a bioregion of interest.
#'@export

agency_count2 <- function(long){
  phy <- wide(long, "PHYLUM")
  agency.count <- data.frame(table(phy$AGENCY_CODE))
  names(agency.count) <- c("Agency Code and Program", "Total Count (N)")
  agency.count$'% of Data' <- round((agency.count$`Total Count` / 
                                     sum(agency.count$`Total Count`)) * 100, 2)
  agency.count <- agency.count[order(agency.count$`Total Count`), ]
  return(agency.count)
}

#==============================================================================
#'Test Agency/Program Standard Taxonomic Resolution of Specific Taxa
#'
#'@param long = long data frame containing taxonomic counts.
#'@param low.rank = The lower of taxonomic rank of interest.
#'@param upper.rank = The upper of taxonomic rank of interest; corresponds with
#'the taxon of interest.
#'@param taxon = The taxon of interest.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Returns a list of names at the taxonomic level specified by low.rank,
#' aggregated by the taxon of interest. The data frame
#'also contains the percentage of data that each agency/program composes
#'in the bioregion being investigated.
#'@export

test_id2 <- function(long, low.rank, upper.rank, taxon){
  agency.count <- agency_count2(long)
  #test <- agency_taxa_diff2(long, low.rank)
  taxa.names <- group.names(long, low.rank, upper.rank, taxon)
  if(length(taxa.names) > 0) taxa.names <- droplevels.factor(taxa.names[taxa.names != "UNIDENTIFIED"])
  wide.df <- BIBI::wide(long, low.rank)
  agg <- aggregate(wide.df[, 6:ncol(wide.df)], by = list(wide.df$AGENCY_CODE), FUN = sum)
  taxa.df <- agg[, c(names(agg) %in% taxa.names)]
  agency.taxa <- data.frame(cbind(agg[, 1], taxa.df), stringsAsFactors = F)
  names(agency.taxa)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, agency.taxa, by = "Agency Code and Program")
  if(ncol(final.df) == 4) return(as.numeric(as.character(final.df[, 4])))
  if(ncol(final.df) > 4) return(rowSums(final.df[, 4:ncol(final.df)]))
  if(ncol(final.df) < 4) return(as.numeric(0))
}



#==============================================================================
#'Checks for common Agency/Program Discrepancies
#'
#'@param long = long data frame containing taxonomic counts.
#'@param bioregion = The bioregion code used to aggregate the data.
#'@return Checks for common agency/program discrepancies at the order, family,
#'and genus level.
#'@export



common_agency_discrepancies2 <- function(long){
  agency.count <- agency_count2(long)
  agency.df <- data.frame(AGENCY_CODE = sort(unique(long$AGENCY_CODE)))

  agency.df$O_Gastropoda <- test_id2(long, "ORDER", "CLASS", "GASTROPODA")
  agency.df$O_Hirudinea <- test_id2(long, "ORDER", "SUBCLASS", "HIRUDINEA")
  agency.df$O_Oligochaeta <- test_id2(long, "ORDER", "CLASS", "OLIGOCHAETA")


  agency.df$F_Gastropoda <- test_id2(long, "FAMILY", "CLASS", "GASTROPODA")
  agency.df$F_Hirudinea <- test_id2(long, "FAMILY", "SUBCLASS", "HIRUDINEA")
  agency.df$F_Oligochaeta <- test_id2(long, "FAMILY", "CLASS", "OLIGOCHAETA")

  

  agency.df$G_Gastropoda <- test_id2(long, "GENUS", "CLASS", "GASTROPODA")
  agency.df$G_Hirudinea <- test_id2(long, "GENUS", "SUBCLASS", "HIRUDINEA")
  agency.df$G_Oligochaeta <- test_id2(long, "GENUS", "CLASS", "OLIGOCHAETA")

  agency.df$G_Chironomidae <- test_id2(long, "GENUS", "FAMILY", "CHIRONOMIDAE")
  
  names(agency.df)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, agency.df, by = "Agency Code and Program")
  #return(final.df)
  #write.csv(final.df, file = paste(c(bioregion, "_agency_res_exclusions.csv"), collapse = ""))
  return(final.df)
}

common_agency_discrepancies3 <- function(long, level){
  agency.count <- agency_count2(long)
  agency.df <- data.frame(AGENCY_CODE = sort(unique(long$AGENCY_CODE)))
  
  if(level %in% "CLASS"){
    class <- wide(long, "CLASS")
    class2 <- data.frame(aggregate(class[, 6:ncol(class)] , by =list(class$AGENCY_CODE),  sum))
    write.csv(class2, "class_counts.csv")
    
    agency.order <- data.frame(sapply(unique(long$CLASS), function(x) test_id2(long, "ORDER", "CLASS", x)))
    names(agency.order) <- paste("ORD", names(agency.order), sep = "_")
    agency.family <- data.frame(sapply(unique(long$CLASS), function(x) test_id2(long, "FAMILY", "CLASS", x)))
    names(agency.family) <- paste("FAM", names(agency.family), sep = "_")
    agency.genus <- data.frame(sapply(unique(long$CLASS), function(x) test_id2(long, "GENUS", "CLASS", x)))
    names(agency.genus) <- paste("GEN", names(agency.genus), sep = "_")

    grouped <- cbind(agency.df, agency.order, agency.family, agency.genus)
  }

  if(level %in% "ORDER"){
    ord <- wide(long, "ORDER")
    ord2 <- data.frame(aggregate(ord[, 6:ncol(ord)] , by =list(ord$AGENCY_CODE),  sum))
    write.csv(ord2, "Ordinal_counts.csv")
    
    agency.family <- data.frame(sapply(unique(long$ORDER), function(x) test_id2(long, "FAMILY", "ORDER", x)))
    names(agency.family) <- paste("FAM", names(agency.family), sep = "_")
    agency.genus <- data.frame(sapply(unique(long$ORDER), function(x) test_id2(long, "GENUS", "ORDER", x)))
    names(agency.genus) <- paste("GEN", names(agency.genus), sep = "_")
    
    grouped <- cbind(agency.df, agency.family, agency.genus)
  }
  
  if(level %in% "FAMILY"){
    fam <- wide(long, "FAMILY")
    agency.genus <- data.frame(sapply(unique(long$FAMILY), function(x) test_id2(long, "GENUS", "FAMILY", x)))
    #names(agency.genus) <- paste("GEN", names(agency.genus), sep = "_")
    
    grouped <- cbind(agency.df, agency.genus)
  }
  
  
  names(grouped)[1] <- "Agency Code and Program"
  final.df <- merge(agency.count, grouped, by = "Agency Code and Program")
  #return(final.df)
  #write.csv(final.df, file =  "class_agency_res_exclusions.csv")
  return(final.df)
}