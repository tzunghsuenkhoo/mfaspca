
#' Title
#'
#' @param W
#' @param XS
#' @param CM
#' @param years
#' @param k
#' @param l
#'
#' @return
#' @export
#'
#' @examples
F_moran_plots <- function(W, XS, CM, years, k, l){

  #Transform the Spatial weight matrix
  W2=W/apply(W,1,sum)


  FMMoranall <- FMORAN_I(XS, W2, CM$mylistX, CM$Function_proj, k, l)
  FMMoran <- FMMoranall$FMMoran
  BMMoran <- FMMoranall$BVMoran

  ### Bivariate Functional Moran
  # Extract functional data values for bivariate functional Moran
  BMMoran_values <- BMMoran@X
  # Extract functional data values for multivariate functional Moran
  FMMoran_values <- FMMoran@X

  # Set up plotting parameters
  par(cex.axis = 0.8, cex.lab=0.8, cex=1.5,  # Reduce size of axis tick labels
      tcl = -0.2)       # Set tick mark size (negative values make ticks smaller)

  # Define the y-axis tick marks
  y_ticks <- seq(round(min(BMMoran_values)-0.01, digits=2), round(max(BMMoran_values)+0.01,digits=2), by=0.01)

  # Create the plot with years on the x-axis
  matplot(years, t(BMMoran_values), type = "l", col = 1:ncol(FMMoran_values),
          ylab = "Bivariate Functional Moran's I statistic",
          xlab = "Year",
          xlim = c(min(years), max(years)),  # Set x-axis limits to cover the range of years
          ylim = c(min(BMMoran_values)-0.01, max(BMMoran_values)+0.005),  # Set y-axis limits
          yaxt = "n",            # Suppress default y-axis
          xaxt = "n")            # Suppress default x-axis

  # Add custom y-axis with specified tick marks
  axis(side = 2, at = y_ticks, labels = y_ticks, col.axis = "black")

  # Add custom x-axis with years
  axis(side = 1, at = years, labels = years, cex.axis = 0.6, col.axis = "black")

  ### Multivariate Functional Moran's I
  # Set up plotting parameters
  par(cex.axis = 0.8, cex.lab=0.8, cex=1.5,  # Reduce size of axis tick labels
      tcl = -0.2)       # Set tick mark size (negative values make ticks smaller)

  # Define the y-axis tick marks
  y_ticks <- seq(round(min(FMMoran_values)-0.01, digits=2), round(max(FMMoran_values)+0.01,digits=2), by=0.01)

  # Create the plot with years on the x-axis
  matplot(years, t(FMMoran_values), type = "l", col = 1:ncol(FMMoran_values),
          ylab = "Multivariate Functional Moran's I statistic",
          xlab = "Year",
          xlim = c(min(years), max(years)),  # Set x-axis limits to cover the range of years
          ylim = c(min(FMMoran_values)-0.01, max(FMMoran_values)+0.005),  # Set y-axis limits
          yaxt = "n",            # Suppress default y-axis
          xaxt = "n")            # Suppress default x-axis

  # Add custom y-axis with specified tick marks
  axis(side = 2, at = y_ticks, labels = y_ticks, col.axis = "black")

  # Add custom x-axis with years
  axis(side = 1, at = years, labels = years, cex.axis = 0.6, col.axis = "black")
}

