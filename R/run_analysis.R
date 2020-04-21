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




  # # !!!!! Testing: start

  model_normal_lm2 <- lm(
    y ~ factor(j) + factor(l),
    data = data$data
  )

  coeff_names <- names(model_normal_lm2$coefficients)
  theta_l_hat <- as.numeric(model_normal_lm2$coefficients)
  sigma_l_hat <- vcov(model_normal_lm2)
  indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  coeff_names <- coeff_names[indices]
  theta_l_hat <- theta_l_hat[indices]
  sigma_l_hat <- sigma_l_hat[indices,indices]

  # # !!!!! Testing: end



  # Use numerical ML to estimate mu_hat
  neg_log_lik <- function(theta, d, J, theta_l_hat, sigma_l_hat) {

    l_times <- 1:(J-1)
    mu_d <- 1-exp(-l_times/d)
    log_lik <- -0.5 * t(theta_l_hat - theta*mu_d) %*%
      solve(sigma_l_hat) %*% (theta_l_hat - theta*mu_d)

    return(-1 * log_lik)

  }

  # Generate ML estimates of theta and d
  opt <- optim(
    par = c(theta=-0.1, d=1),
    fn = function(par) {
      return (
        neg_log_lik(
          theta = par[1],
          d = par[2],
          J = data$params$n_time_points,
          theta_l_hat = theta_l_hat,
          sigma_l_hat = sigma_l_hat
        )
      )
    }
  )
  theta_hat <- opt$par[["theta"]]
  d_hat <- opt$par[["d"]]

  # Use hessian to estimate SEs
  h <- hessian(
    func = function(par) {
      return (
        neg_log_lik(
          theta = par[1],
          d = par[2],
          J = data$params$n_time_points,
          theta_l_hat = theta_l_hat,
          sigma_l_hat = sigma_l_hat
        )
      )
    },
    c(theta_hat,d_hat)
  )
  se_theta_hat <- sqrt(solve(h)[1,1])
  se_d_hat <- sqrt(solve(h)[2,2])

  # Construct theta_hat2 estimator (corresponds to "approach 2")
  l_times <- 1:(data$params$n_time_points-1)
  mu_d_hat <- 1-exp(-l_times/d_hat)
  theta_hat2 <- (t(theta_l_hat) %*% solve(sigma_l_hat) %*% mu_d_hat) /
    (t(mu_d_hat) %*% solve(sigma_l_hat) %*% mu_d_hat)

  # # Construct var_theta_hat2 estimator (corresponds to "approach 2")
  var_d_hat <- se_d_hat^2
  grad_g <- -l_times * d_hat^(-2) * exp(-l_times/d_hat)
  m <- mu_d_hat
  si <- solve(sigma_l_hat)
  tl <- theta_l_hat
  l <- l_times
  grad_h <- t(sapply(l_times, function(j) {
    ( ((t(m)%*%si%*%m)*(t(tl)%*%si[,j]) -
         2*(t(m)%*%si[,j])*(t(tl)%*%si%*%m) ) /
        ((t(m)%*%si%*%m)^2) )[1,1]
  }))
  var_theta_hat2 <- (grad_h%*%grad_g%*%var_d_hat%*%t(grad_g)%*%t(grad_h))
  se_theta_hat2 <- sqrt(var_theta_hat2[1,1])

  # ggplot(theta_hats, aes(x=l, y=theta_l_est)) +
  #   geom_point(color="turquoise")

  return (list(
    theta_hat = theta_hat,
    d_hat = d_hat,
    se_d_hat = se_d_hat,
    se_theta_hat = se_theta_hat,
    se_theta_hat2 = se_theta_hat2
  ))

}
