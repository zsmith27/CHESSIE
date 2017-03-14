#==============================================================================
#' 1-3-5 Scoring
#'
#'@param score.df = a dataframe containing raw metric values.
#'@param thresh.df = a matrix containing thresholds for scoring.
#'@param metric = a string identifying the metric of interest.
#'@param distrubance = A string identifying if the metric increases
#' ("INCREASE") or decreases ("DECREASE") with distrubance.
#'@return Scores metrics on a 1-3-5 scale based on the thresholds from the 2011
#'Chessie BIBI
#'@export

score_1_3_5 <- function(score.df, thresh.df, metric, disturbance){
  XT <- thresh.df$XT
  XM <- thresh.df$XM
  if(disturbance %in% "INCREASE"){
    score.df[, metric] <- ifelse (score.df[ , metric] < XT & score.df[ , metric] > XM, 3,
                                  ifelse (score.df[ , metric] >= XT, 1,
                                          ifelse (score.df[ , metric] <= XM, 5, "ERROR")))
  }else{
    score.df[, metric] <- ifelse (score.df[ , metric] > XT & score.df[ , metric] < XM, 3,
                                  ifelse (score.df[ , metric] >= XT, 5,
                                          ifelse (score.df[ , metric] <= XM, 1, "ERROR")))
  }
  
}