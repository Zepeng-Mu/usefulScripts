#' @param pvector a vector for observed P-values
#' @param col a string for point color in QQ-plot
#' @param add a Boolean for whether to overlay on an existing plot
#' @param ylim a numeric vector specifying the range of y axis
qqplot <- function(pvector, col = "black", add = F, ylim = NULL){
  expectedP <- -log10(ppoints(length(pvector)))
  observedP <- -log10(sort(pvector, decreasing = F))
  if (add == F) {
    plot(x = expectedP, y = observedP, col = col,
         xlab = expression(Expected ~ ~-log[10](italic(p))),
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         pch = 16, cex = 0.7, ylim = NULL)
    abline(0, 1, col = "red")
  }else{
    points(x = expectedP, y = observedP, col = col,
           xlab = expression(Expected ~ ~-log[10](italic(p))),
           ylab = expression(Observed ~ ~-log[10](italic(p))),
           pch = 16, cex = 0.7, ylim = NULL) 
  }
}