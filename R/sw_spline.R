#' Return a linear spline function with specified knots and slopes
#'
#' @description Spline starts at the point (0,0) and is constrained to end with
#'     a slope of zero after the last knot
#' @param x Location at which spline should be evaluated
#' @param knots Num/int vector that specifies x-coordinates of knots (excludes
#'     knot at x=0)
#' @param slopes Num vector that specifies slopes between knots (including the
#'     knot at x=0)
#' @return
#' @export
sw_spline <- function(x, knots, slopes) {

  if (length(knots)!=length(slopes)) {
    stop("Length of `knots` argument must equal length of `slopes` argument")
  }

  y <- slopes[1] * pmax(0,x)

  for (i in 1:length(knots)) {

    if (i<length(knots)) {
      y <- y + (slopes[i+1]-slopes[i]) * pmax(0, x-knots[i])
    } else {
      y <- y - slopes[i] * pmax(0, x-knots[i])
    }

  }

  return(y)

}
