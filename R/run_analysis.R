#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @return TO DO
#' @export
# FN: generate_dataset
run_analysis <- function(data) {

  # # Step 1: Estimate the time-on-intervention fixed effects
  # model_normal_gee2 <- geeglm(
  #   y ~ factor(j) + factor(l),
  #   data = data$data,
  #   id = i,
  #   family = "gaussian",
  #   corstr = "independence"
  #   # corstr = "exchangeable"
  # )
  #
  # # Extract estimates and covariance matrix
  # coeff_names <- names(model_normal_gee2$coefficients)
  # theta_l_hat <- as.numeric(model_normal_gee2$coefficients)
  # sigma_l_hat <- model_normal_gee2$geese$vbeta
  # indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  # coeff_names <- coeff_names[indices]
  # theta_l_hat <- theta_l_hat[indices]
  # sigma_l_hat <- sigma_l_hat[indices,indices]

  # Step 2: Use nonlinear GLS to estimate long-term effect and lag duration








}
