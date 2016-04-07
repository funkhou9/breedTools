#' dist_plot
#' 
#' Plotting function to inspect the distribution of results from a genome-wide
#' scan of breed composition. Plots the joint distribution (scatterplot) of breed composition
#' results of breed specified and R squared values. Marginal distributions (histograms)
#' on the axes.
#'
#'  @param dat data.frame of genome-wide breed composition results
#'  @param breed character naming the desired breed composition result in dat to plot
#'  @param xlab character giving label for x-axis
#'  @return NONE
#'  @export
dist_plot <- function(dat,
                      breed,
                      xlab = "Estimated Breed Composition") {
  
  # Marginal distributions as histograms
  yhist <- hist(dat[, "R2"], prob = T, breaks = 100)
  xhist <- hist(dat[, breed], prob = T, breaks = 200)
  
  # Reset plotting window
  plot.new()
  
  # 1st layer - joint scatterplot
  par(fig = c(0, 0.8, 0, 0.8), new = TRUE)
  plot(dat[, breed], dat[, "R2"],
       xlab = xlab, ylab = "R Squared")
  abline(lm(dat[, "R2"] ~ dat[, breed]))
  
  # 2nd layer - x marginal distribution (composition results)
  par(fig = c(0, 0.8, 0.60, 1), new = TRUE)
  barplot(xhist$counts, axes = FALSE, space = 0)
  
  # 3rd layer - y marginal distribution (R squared results)
  par(fig = c(0.67, 1, 0, 0.8), new = TRUE)
  barplot(yhist$counts, axes = FALSE, space = 0, horiz = TRUE)
}