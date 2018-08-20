#==============================================================================
#'Index Plot
#'
#'@param index.df
#'@param condition.colname
#'@return 
#'@import ggplot2
#'@export
index_plot2 <- function(index.df, condition.colname, condition.order = NULL,
                        index.colname = "INDEX", threshold, plot.title) {
  
  if (!is.null(condition.order)) {
    index.df[, condition.colname] <- factor(index.df[, condition.colname],
                                            levels = condition.order)
  }
  index.df <- index.df[!index.df$CATEGORY %in% "MIX", ]
  #index.df$CATEGORY <- ifelse(index.df$CATEGORY %in% "REF", "Reference",
  #                            ifelse(index.df$CATEGORY %in% "MIN", "Minimally Degraded",
  #                                   ifelse(index.df$CATEGORY %in% "MOD", "Moderately Degraded",
  #                                          ifelse(index.df$CATEGORY %in% "DEG", "Degraded", "ERROR"))))
  #----------------------------------------------------------------------------
  agg.df <- aggregate(index.df[, index.colname] ~ index.df[, condition.colname],
                      data = index.df, FUN = max)
  names(agg.df) <- c(condition.colname, "MAX")
  len.df <- aggregate(index.df[, index.colname] ~ index.df[, condition.colname],
                      data = index.df, FUN = length)
  names(len.df) <- c(condition.colname, "N")
  n.size <- merge(agg.df, len.df, by = condition.colname)
  #----------------------------------------------------------------------------
  ref.50 <- unique(index.df$REF_50)
  ref.25 <- unique(index.df$REF_25)
  ref.10 <- unique(index.df$REF_10)
  ref.half <- unique(index.df$HALF_REF_10)
  #----------------------------------------------------------------------------
  final.plot <- ggplot(index.df, aes_string(condition.colname, index.colname)) +
    geom_rect(xmin = 0, xmax = 5, ymin = ref.50, ymax = 100,
              fill = "darkgreen", alpha = 0) +
    annotate("rect", xmin = 0, xmax = 5, ymin = ref.50, ymax = 100, alpha = 0.3,
             fill = "darkgreen") +
    geom_rect(xmin = 0, xmax = 5, ymin = ref.25, ymax = ref.50 - 0.001,
              fill = "green3", alpha = 0) +
    annotate("rect", xmin = 0, xmax = 5, ymin = ref.25, ymax = ref.50 - 0.001, alpha = 0.3,
             fill = "green3") +
    geom_rect(xmin = 0, xmax = 5, ymin = ref.10, ymax = ref.25 - 0.001,
              fill = "yellow2", alpha = 0) +
    annotate("rect", xmin = 0, xmax = 5, ymin = ref.10, ymax = ref.25 - 0.001, alpha = 0.3,
             fill = "yellow2") +
    geom_rect(xmin = 0, xmax = 5, ymin = ref.half, ymax = ref.10 - 0.001,
              fill = "orange2", alpha = 0) +
    annotate("rect", xmin = 0, xmax = 5, ymin = ref.half, ymax = ref.10 - 0.001, alpha = 0.3,
             fill = "orange2") +
    geom_rect(xmin = 0, xmax = 5, ymin = 0, ymax = ref.half - 0.001,
              fill = "red3", alpha = 0) +
    annotate("rect", xmin = 0, xmax = 5, ymin = 0, ymax = ref.half - 0.001, alpha = 0.3,
             fill = "red3") +
    stat_boxplot(geom = 'errorbar', width = 0.3, size = 0.8) + 
    geom_boxplot(notch = FALSE, size = 0.8) +
    scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          # text = element_text(size = 20),
          #axis.text.x = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
          axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
          plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm")) +
    labs(title = plot.title) +
    ylab("Index") +
    xlab("")+
    geom_text(data = n.size, aes(y = 101, label = N), vjust = 0)
  
  
  
  if (!is.null(threshold)) {
    final.plot <- final.plot + 
      geom_hline(aes(yintercept = threshold),
                 size = 1.25, color = "red2", linetype = "dotted")
  }
  
  return(final.plot)
}
