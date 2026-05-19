
#' @title Plot Tweedie Models
#' @name tweedie_plot
#' @description This function produced a plot of the specified Tweedie distribution.
#'
#' @usage tweedie_plot(y, xi = NULL, mu, phi, type = "pdf", power = NULL, 
#'                     add = FALSE, 
#'                     plot_args = list(), point_args = list(), line_args = list())
#'
#' @details If \eqn{1 < p < 2}{1 < power < 2}, the mass at \eqn{Y=0}{Y = 0} is automatically added.
#'
#' @param y the values for \eqn{y}{y} in the plot.
#' @param power the variance power \eqn{p}{power}.
#' @param mu the mean of the distribution \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param xi a synonym for \code{power}.
#' @param type the type of plot, either \code{pdf} (the default) or \code{cdf}.
#' @param add logical; if \code{TRUE}, the plot is added to the current plot; if \code{FALSE} (the default) the plot is produced on a fresh plot.
#' @param plot_args A named list of arguments controlling the main plot.
#'   These are passed to \code{\link[graphics]{plot}}.
#' @param point_args A named list of graphical parameters for plotted points.
#'   These are passed to \code{\link[graphics]{points}}.
#' @param line_args A named list of graphical parameters for plotted lines.
#'   These are passed to \code{\link[graphics]{lines}}.
#'
#' @examples
#' y <- seq(0, 4, length = 50)
#' tweedie_plot(y, power = 1.1, mu = 1, phi = 1,
#'              line_args = list(lwd = 2))
#'
#'
#' @importFrom graphics lines rug par mtext abline axis  points plot
#'
#' @export
tweedie_plot <- function(y,
                         xi = NULL,
                         mu,
                         phi,
                         type = "pdf",
                         power = NULL,
                         add = FALSE,
                         plot_args  = list(),
                         point_args = list(),
                         line_args  = list()) {
  # Sort out the xi/power notation
  if (is.null(power) &
      is.null(xi))
    stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if (is.null(power)) {
    power <- xi
  } else {
    xi.notation <- FALSE
  }
  if (is.null(xi)) {
    xi.notation <- FALSE
    xi <- power
  }
  if (xi != power) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  index.par       <- ifelse(xi.notation, "xi", "p")
  index.par.long  <- ifelse(xi.notation, "xi", "power")
  
  if ((power < 0) | ((power > 0) & (power < 1))) {
    stop(paste(
      "Plots cannot be produced for",
      index.par.long,
      "=",
      power,
      "\n"
    ))
  }
  
  is.pg <- (power > 1) & (power < 2)
  
  if (type == "pdf") {
    fy <- dtweedie(
      y = y,
      power = power,
      mu = mu,
      phi = phi
    )
  } else {
    fy <- ptweedie(
      q = y,
      power = power,
      mu = mu,
      phi = phi
    )
  }
  
  # BEGIN: Check for some given parameters supplied via args
  plot_defaults <- list(
    col = "black",
    main = "Tweedie distribution",
    xlab = "Y",
    ylab = "Prob fn"
  )
  point_defaults <- list(pch = 19,
                         col = "black",
                         cex = 1)
  line_defaults  <- list(col = "black", 
                         lwd = 1)
  
  plot_args  <- utils::modifyList(plot_defaults, 
                                  plot_args)
  point_args <- utils::modifyList(point_defaults, 
                                  point_args)
  line_args  <- utils::modifyList(line_defaults, 
                                  line_args)
  # END: Check for some given parameters supplied via args
  
  
  if (!add) {
    if (is.pg) {
      do.call(plot, 
              c(list(
                x = range(y),
                y = range(fy),
                type = "n"
                ), 
                plot_args))
      if (any(y == 0)) {
        # The exact zero
        do.call(points, 
                c(list(x = y[y == 0], 
                       y = fy[y == 0]), 
                  point_args))
      }
      if (any(y > 0)) {
        # The exact zero
        do.call(lines, 
                c(list(x = y[y > 0], 
                       y = fy[y > 0]), 
                  line_args))
      }
    } else {
      # Not a Poison-gamma dist
      do.call(plot, 
              c(list(
                x = range(y),
                y = range(fy),
                type = "n"
                ), 
                plot_args))
      do.call(lines, 
              c(list(x = y, 
                     y = fy), 
                line_args))
    }
  } else {
    # Add; no new plot
    if (is.pg) {
      if (any(y == 0)) {
        # The exact zero
        do.call(points, 
                c(list(x = y[y == 0], 
                       y = fy[y == 0]), 
                  point_args))
      }
      if (any(y > 0)) {
        # The exact zero
        do.call(lines, 
                c(list(x = y[y > 0], 
                       y = fy[y > 0]), 
                  line_args))
      }
    } else {
      # Not a Poison-gamma dist
      do.call(lines, 
              c(list(x = y, 
                     y = fy), 
                line_args))
    }
    
  }
  return(invisible(list(y = fy, x = y)))
  
}





#' @rdname tweedie_plot
#' @param ... Additional graphical arguments, passed to \code{plot_args} for 
#'   backward compatibility. Deprecated; use \code{plot_args} in 
#'   \code{\link{tweedie_plot}} instead.
#' @export
tweedie.plot <- function(y,
                         xi = NULL,
                         mu,
                         phi,
                         type = "pdf",
                         power = NULL,
                         add = FALSE,
                         ...) {
  lifecycle::deprecate_warn(when = "3.0.5",
                            what = "tweedie.plot()",
                            with = "tweedie_plot()")
  if (is.null(power))
    power <- xi
  
  # Capture ... and pass to plot_args for backward compatibility
  dots <- list(...)
  
  tweedie_plot(
    y = y,
    power = power,
    mu = mu,
    phi = phi,
    type = type,
    add = add,
    plot_args = dots
  )
}
