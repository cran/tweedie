#' @title Plot Tweedie Models
#' @name tweedie_plot
#' @description This function produced a plot of the specified Tweedie distribution.
#'
#' @usage tweedie_plot(y, xi = NULL, mu, phi, type = "pdf", power = NULL, add = FALSE, ...)
#' 
#' @details If \eqn{1 < p < 2}{1 < power < 2}, the mass at \eqn{Y=0}{Y = 0} is automatically added.
#'
#' @param y the values for \eqn{y}{y} in the plot.
#' @param power the variance power \eqn{p}{power}.
#' @param mu the mean of the distribution \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param xi a synonym for \code{power}.
#' @param type the type of plot, either \code{PDF} (the default) or \code{CDF}.
#' @param add logical; if \code{TRUE}, the plot is added to the current plot; if \code{FALSE} (the default) the plot is produced on a fresh plot.
#' @param ... plotting parameters passed to \code{plot()}.
#' 
#' @examples
#' y <- seq(0, 4, length = 50)
#' tweedie_plot(y, power = 1.1, mu = 1, phi = 1)
#' 
#' 
#' @importFrom graphics lines rug par mtext abline axis  points plot
#'
#' @export
tweedie_plot <- function(y, xi = NULL, mu, phi, type = "pdf", power = NULL, 
                         add =FALSE, ...) {
  
  # Sort out the xi/power notation
  if ( is.null(power) & is.null(xi) ) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if ( is.null(power) ) {
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if ( is.null(xi) ) {
    xi.notation <- FALSE
    xi <- power
  }
  if ( xi != power ) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  index.par       <- ifelse( xi.notation, "xi", "p")
  index.par.long  <- ifelse( xi.notation, "xi", "power")
  
  if ( ( power < 0 ) | ( ( power > 0 ) & ( power < 1 ) ) ) {
    stop( paste("Plots cannot be produced for", index.par.long, "=", power, "\n") )
  }   
  
  is.pg <- ( power > 1 ) & ( power < 2 )
  
  if ( type == "pdf") {
    fy <- dtweedie( y = y, 
                    power = power, 
                    mu = mu,
                    phi = phi)
  } else {
    fy <- ptweedie( q = y, 
                    power = power, 
                    mu = mu,
                    phi = phi)
  }
  
  if ( !add ) {
    if ( is.pg ) {
      plot( range(fy) ~ range( y ), 
            type = "n", ...)
      if ( any( y == 0 ) ) { # The exact zero
        points( fy[y == 0] ~ y[y == 0], 
                pch = 19, ... )
      }
      if ( any( y > 0 ) ) { # The exact zero
        lines( fy[y > 0]   ~ y[y > 0], 
               pch = 19, ... )
      }
    } else {  # Not a Poison-gamma dist
      plot( range(fy) ~ range( y ), 
            type = "n", ...)
      lines( fy ~ y, 
             pch = 19, ... )
    }
  } else {# Add; no new plot
    if ( is.pg ) {
      if ( any( y == 0 ) ) { # The exact zero
        points( fy[y == 0] ~ y[y == 0], 
                pch = 19, ... )
      }
      if ( any( y > 0 ) ) { # The exact zero
        lines( fy[y > 0]   ~ y[y > 0], 
               pch = 19, ... )
      }
    } else {  # Not a Poison-gamma dist
      lines( fy ~ y, 
             pch = 19, ... )
    }
    
  }
  return(invisible(list(y = fy, 
                        x = y) ))
  
}





#' @rdname tweedie_plot
#' @export
tweedie.plot <- function(y, xi = NULL, mu, phi, type = "pdf", power = NULL, 
                         add =FALSE, ...){
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "tweedie.plot()", 
                            with = "tweedie_plot()")
  if ( is.null(power)) power <- xi
  tweedie_plot(y = y, 
              power = power, 
              mu = mu, 
              phi = phi, 
              type = type, 
              add = add, 
              ... = ...)
}

