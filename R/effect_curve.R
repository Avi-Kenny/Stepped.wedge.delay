#' Return the effect curve value (% of Tx effect reached as a fn of time steps)
#'
#' @param x Value at which to evaluate the effect curve
#' @param type One of the following character strings:
#'     - "exp": Hussey & Hughes; ignores delay
#'     - "spline": "exposure treatment indicators" model
#' @param params List of parameters, which differs by "type"
#'     - For "exp", `d` is the rate parameter
#'     - For "spline", `knots` specifies the knot locations (the first of
#'       which must be zero) and `slopes` represents the slopes between
#'       knots. Spline starts at (0,0) and has slope zero after the last knot.
#'     - For "parabola", paramaters `a`, `b`, and `c` correspond to the parabola
#'       given by y = a(x^2) + bx + c
#' @return A number between zero and one, representing the % of the Tx effect
#'     reached after x time steps

effect_curve <- function(x, type, params) {

  p <- params

  if (type=="exp") {
    y <- ifelse(x>6, 1, (1-exp(-x/p$d)) / (1-exp(-6/p$d)))
  }

  if (type=="spline") {

    if ((length(p$knots)-1)!=length(p$slopes)) {
      stop("Length of `knots` must equal length of `slopes` plus one")
    }

    y <- p$slopes[1] * pmax(0,x)

    for (i in 2:length(p$knots)) {

      if (i<length(p$knots)) {
        y <- y + (p$slopes[i] - p$slopes[i-1]) * pmax(0,x-p$knots[i])
      } else {
        y <- y - p$slopes[i-1] * pmax(0,x-p$knots[i])
      }

    }

  }

  if (type=="parabola") {

    y <- ifelse(x>6, 1, (p$a * x^2) + (p$b * x) + p$c)

  }

  if (type=="non-monotonic") {

    y <- ifelse(x>6, 0.75, (-1/16*x^2)+(1/2*x))

  }

  if (!(type %in% c("exp", "spline", "parabola", "non-monotonic"))) {
    stop("type must be in c('exp', 'spline', 'parabola', 'non-monotonic')")
  }

  return (y)

}
