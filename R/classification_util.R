#' @title Returns Accuracy/Precision/Recall/F1-Score table
#' @name get_prediction_results
#' @param x vector of values computed by classification model
#' @param y vector of actual values to compare with \code{x}
#' @param positive character value (class) to be considered as positive in \code{x} and \code{y}
#' @return values of prediction metrics
#' @examples {
#' get_prediction_results(c(1, 1, 0, 0, 1, 1, 1), c(1, 1, 0, 0, 1, 0, 1), '1')
#' }
#' @import caret
#' @export
get_prediction_results <- function(x, y, positive) {
  x <- as.factor(x)
  y <- as.factor(y)
  levels(y) <- levels(x)
  matrix <- confusionMatrix(x, y, positive)
  accuracy <- matrix$overall[1] # Correctness of model
  precision <- matrix$byClass[3] # Positive prediction value
  neg_precision <- matrix$byClass[4] # Negative prediction value
  sensitivity <- matrix$byClass[1] # True positive recognition rate (aka recall)
  specificity <- matrix$byClass[2] # True negative recognition rate
  f1_score <- 2 * ((precision * sensitivity) / (precision + sensitivity))
  names(f1_score) <- 'F1-Score'
  results <- c(accuracy, precision, sensitivity, specificity, f1_score)
  results
}

#' @title Discretizes a vector naturally
#' @name discretize_naturally
#' @description Discretizes a numeric vector into 3, 5 or 7 and returns a factor vector with values as: LOW, MEDIUM, HIGH; VERY_LOW, LOW, MEDIUM, HIGH, VERY_HIGH; VERY_LOW, LOW, MEDIUM_LOW, MEDIUM, MEDIUM_HIGH, HIGH, VERY_HIGH
#' @param x vector to discretize
#' @param breaks is the number of classes. This can only be 3, 5 or 7
#' @return descritized vector
#' @import arules
#' @export
discretize_naturally <- function(x, breaks=7) {
  if (mode(x) != 'numeric') {
    x
  } else {
    if (breaks == 3) {
      discretize(x, method="interval", breaks, labels=c('L','M','H'))
    } else if (breaks == 5) {
      discretize(x, method="interval", breaks, labels=c('VL','L','M','H','VH'))
    } else {
      discretize(x, method="interval", breaks, labels=c('VL','L','ML','M','MH','H','VH'))
    }
  }
}
