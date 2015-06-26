#' sim_plot
#' 
#' Generates a scatterplot to inspect the results of repeatedly simulated 
#' QP_SimCross output.
#' 
#' @param sim data.frame returned from _____, containing breed composition
#'  estimates from simulated offspring.
#' @param y character naming the column of sim to be plotted on y. By default,
#'  "XBreed1.qp" refers to the estimated breed composition of breed 1.
#' @return NONE
#' @export
sim_plot <- function(sim, y = "XBreed1.qp") {
  
  # Randomly selected animals for use in simulation may not be
  #   completely purebred according to our methods. Therefore
  #   we need to calculate the actual breed of the simulated offspring
  act_b <- sim[, "ActBreed1"] * sim[, "Breed1.qp"]
  
  # Plot actual breed compos on x, simulated on y
  plot(act_b, sim[, "XBreed1.qp"], xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Actual Breed Composition", ylab = "Predicted Breed Composition",
       col = (c("black","green3","blue","magenta")[sim[, "idx"]]))
  
  legend(0.8, 0.2, c("Duroc", "Hampshire", "Landrace", "Yorkshire"),
         text.col = c("black","green3","blue","magenta"), bty = "n")
  
  # To show correlation between actual and simulated - draw a (1, 1) line
  abline(0, 1, col = "red")
}