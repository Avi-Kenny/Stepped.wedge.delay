#' Return the log likelihood for the spline model
#'
#' @param data A dataset returned by generate_dataset()
#' @return TO DO
#' @export
# FN: generate_dataset
log_lik_spline <- function(sigma_v, sigma_e, alpha, beta, theta, p_x, p_y, g_x) {

  data <- data$data

  # Data: v_i, y

  loop_sum <- 0
  for (i in 1:I) {
    v_i <- 999
    loop_sum <- loop_sum - (v_i^2)
  }
  loop_sum <- loop_sum / (2*(sigma_v^2))
  for (i in 1:I) {
    c_i <- 999
    v_i <- 999
    for (j in 1:J) {
      beta_j <- 999
      r_ij <- r_ij(p_x, p_y, g_x, c_i, j)
      for (k in 1:K) {
        y <- 999
        term <- (y-(alpha + v_i + beta_j + (theta*r_ij)))^2 / (2*(sigma_e^2))
        loop_sum <- loop_sum - term
      }
    }
  }


}



# Helper function to calculate the spline log likelihood (above)
r_ij <- function(p_x, p_y, g_x, c_i, j) {

  I1 <- ifelse(c_i<j & j<=c_i+p_x, 1, 0)
  I2 <- ifelse(c_i+p_x<j & j<=c_i+g_x, 1, 0)
  I3 <- ifelse(c_i+g_x<j, 1, 0)
  (p_y/p_x)*(j-c_i)*I1 + (((c_i-j+1)*(p_y-1))/(g_x-p_x)+1)*I2 + I3

}
