#==============================================================================
#'Import Output
#'
#'@param scores.df = a data frame containing final index scores.
#'@param bioregion = specify the bioregion of interest.
#'@return Box-and-Whisker plots of final index scores aggregated by stream condition class.
#'@export

plot.me <- function(scores.df, bioregion){
  library(ggplot2)
  library(stringr)
  scores.df <- fam.plot[fam.plot$BIOREGION %in% bioregion, ]
  
  scores.df$CATEGORY <- factor(scores.df$CATEGORY,
                               levels = c("REF", "MIN", "MIX", "MOD", "SEV"))
  
  ggplot(data = scores.df, aes(CATEGORY, SCORE))+
    stat_boxplot(geom ='errorbar', width = 0.5) + 
    geom_boxplot(notch = FALSE) + #, aes(fill = factor(CATEGORY))) +
    geom_hline(aes(yintercept = THRESHOLD), size = 1.25, color = "red2") +
    #scale_fill_manual( values = c("forestgreen", "gold", "goldenrod2", "red", "red4")) +
    scale_x_discrete(labels = c("Reference", "Minimally\nDegraded",
                                "Mixed", "Moderately\nDegraded", "Degraded")) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          # text = element_text(size = 20),
          #axis.text.x = element_text(size = 15),
          axis.title.y = element_text(margin= margin(0, 10, 0, 0)),
          axis.title.x = element_text(margin= margin(10, 0, 0, 0)),
          plot.margin  =unit(c(0.5, 1, 0.5, 0.5), "cm")) +
    ylim(0, 100) +
    ylab("Final Score (%)") +
    xlab("Class") +
    labs(title = scores.df$BIOREGION)
}