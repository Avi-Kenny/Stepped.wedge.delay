#' Return the log likelihood for the spline model
#'
#' @param sigma_v !!!!! TO DO
#' @return !!!!! TO DO
#' @export
log_lik_spline <- function(
  sigma_v, sigma_e, alpha,
  beta_1, beta_2, beta_3, beta_4, beta_5,
  theta, p_x, p_y, g_x, data
) {

  I <- data$params$n_clusters
  J <- data$params$n_time_points
  K <- data$params$n_ind_per_cluster
  beta <- c(beta_1, beta_2, beta_3, beta_4, beta_5)

  # !!!!! Can speed this up by calculating this in advance
  df <- data$data %>% mutate(
    r_ij = r_ij(p_x, p_y, g_x, c_i, j),
    s = y - alpha - beta[j] - theta*r_ij,
    t = s^2
  )

  df_sum <- df %>% group_by(i) %>% summarize(
    s = sum(s),
    t = sum(t)
  )

  s_i <- df_sum$s
  t_i <- df_sum$t

  t1 <- (I/2) *
        log((2*pi*(sigma_e^2)*(sigma_v^2))/((sigma_e^2)+(J*K*(sigma_v^2))))
  t2 <- (sigma_v^2)/(2*(sigma_e^2)*((sigma_e^2)+J*K*(sigma_v^2)))
  t3 <- sum(s_i^2-t_i)

  return(t1 + t2*t3)

}



# Helper function to calculate the spline log likelihood (above)
r_ij <- function(p_x, p_y, g_x, c_i, j) {

  I1 <- ifelse(c_i<j & j<=c_i+p_x, 1, 0)
  I2 <- ifelse(c_i+p_x<j & j<=c_i+g_x, 1, 0)
  I3 <- ifelse(c_i+g_x<j, 1, 0)
  (p_y/p_x)*(j-c_i)*I1 + (((c_i-j+1)*(p_y-1))/(g_x-p_x)+1)*I2 + I3

}
