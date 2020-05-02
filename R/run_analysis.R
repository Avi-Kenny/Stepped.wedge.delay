#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param type A character string corresponding to one of the following: "2S LM"
#'     (two-stage using a simple linear model in first stage), "2S GEE"
#'     (two-stage using GEE in first stage), "2S LMM" (two-stage using LMM in
#'     first stage) "PP SPL" (the "Purple point spline" method), "FX SPL" (the
#'     "fixed X-coordinate" spline method), "IG LM" (simple linear model that
#'     ignores the time-lag effect), "IG GEE" (GEE that ignores the time-lag
#'     effect)
#' @param L Passed via simba; list of simulation levels
#' @param C Passed via simba; list of simulation constants
#' @return TO DO
#' @export
# FN: generate_dataset
run_analysis <- function(data, type, L, C) {

  # # !!!!! Testing
  # data <- generate_dataset(
  #   alpha = log(0.1),
  #   tau = 0,
  #   theta = log(0.5),
  #   d = 0,
  #   n_clusters = 24,
  #   n_time_points = 5,
  #   n_ind_per_cluster = 100,
  #   type = "normal",
  #   sigma = 0.1
  # )

# print("check 2")
# print('exists("neg_log_lik")')
# print(exists("neg_log_lik"))

  if (type %in% c("2S LM", "2S GEE", "2S LMM")) {

    if (type=="2S LM") {

      # Run normal linear model
      model_normal_lm2 <- lm(
        y ~ factor(j) + factor(l),
        data = data$data
      )

      # Extract coefficients and SEs
      coeff_names <- names(model_normal_lm2$coefficients)
      theta_l_hat <- as.numeric(model_normal_lm2$coefficients)
      sigma_l_hat <- vcov(model_normal_lm2)

    }

    if (type=="2S GEE") {

      # Run GEE model with exchangeable working covariance matrix
      model_normal_gee2 <- geeglm(
        y ~ factor(j) + factor(l),
        data = data$data,
        id = i,
        family = "gaussian",
        corstr = "exchangeable"
      )

      # Extract estimates and covariance matrix
      coeff_names <- names(model_normal_gee2$coefficients)
      theta_l_hat <- as.numeric(model_normal_gee2$coefficients)
      sigma_l_hat <- model_normal_gee2$geese$vbeta

    }

    if (type=="2S LMM") {

      model_lmm <- lme(
        fixed = y ~ factor(j) + factor(l), # !!!!! Force beta_1=0
        random = ~1|i,
        data = data$data,
        method = "ML"
      )
      summary(model_lmm)

      # Extract estimates and covariance matrix
      coeff_names <- names(model_lmm$coefficients$fixed)
      theta_l_hat <- as.numeric(model_lmm$coefficients$fixed)
      sigma_l_hat <- model_lmm$varFix

    }

    # Truncate theta_l_hat and sigma_l_hat
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]

    # Generate ML estimates of theta and d
    # neg_log_lik <- sim$methods$neg_log_lik # !!!!! Temp fix
    # neg_log_lik <- sim_obj$methods$neg_log_lik # !!!!! Temp fix
    nll <- function(par) {
      return (
        neg_log_lik(
          theta = par[1],
          d = par[2],
          J = L$n_time_points,
          theta_l_hat = theta_l_hat,
          sigma_l_hat = sigma_l_hat
        )
      )
    }
    opt <- optim(par=c(theta=-0.1, d=1), fn=nll)
    theta_hat <- opt$par[["theta"]]
    d_hat <- opt$par[["d"]]
    h <- optimHess(par=c(theta_hat, d_hat), fn=nll)

    # Use hessian to estimate SEs and extract standard errors
    # Note: tryCatch block necessary because sometimes estimated variances are <0
    se_theta_hat <- NA
    se_d_hat <- NA
    tryCatch(
      expr = {
        se_theta_hat <- sqrt(solve(h)[1,1])
        se_d_hat <- sqrt(solve(h)[2,2])
      },
      error = function(cond) {}, # !!!!! log error
      warning = function(cond) {} # !!!!! log error
    )

    return (list(
      d_hat = d_hat,
      theta_hat = theta_hat,
      se_d_hat = se_d_hat,
      se_theta_hat = se_theta_hat
    ))

  }

  if (type == "2S SPL") {

    # !!!!! Next three blocks are identical with "2S LM" method above

    # Run normal linear model
    model_normal_lm2 <- lm(
      y ~ factor(j) + factor(l),
      data = data$data
    )

    # Extract coefficients and SEs
    coeff_names <- names(model_normal_lm2$coefficients)
    theta_l_hat <- as.numeric(model_normal_lm2$coefficients)
    sigma_l_hat <- vcov(model_normal_lm2)

    # Truncate theta_l_hat and sigma_l_hat
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]

    # Generate ML estimates of theta, p_x, and p_y
    neg_log_lik_spl <- sim$methods$neg_log_lik_spl # !!!!! Temp fix
    # neg_log_lik_spl <- sim_obj$methods$neg_log_lik_spl # !!!!! Temp fix
    nll <- function(par) {
      return (
        neg_log_lik_spl(
          theta = par[1],
          p_x = par[2],
          p_y = par[3],
          J = L$n_time_points,
          theta_l_hat = theta_l_hat,
          sigma_l_hat = sigma_l_hat
        )
      )
    }
    opt <- optim(par=c(theta=-0.1, p_x=2, p_y=0.5), fn=nll)
    theta_hat <- opt$par[["theta"]]
    p_x_hat <- opt$par[["p_x"]]
    p_y_hat <- opt$par[["p_y"]]
    h <- optimHess(par=c(theta_hat, p_x_hat, p_y_hat), fn=nll)

    # Use hessian to estimate SEs and extract standard errors
    # Note: tryCatch block necessary because sometimes estimated variances are <0
    se_theta_hat <- NA
    se_p_x_hat <- NA
    se_p_y_hat <- NA
    tryCatch(
      expr = {
        se_theta_hat <- sqrt(solve(h)[1,1])
        se_p_x_hat <- sqrt(solve(h)[2,2])
        se_p_y_hat <- sqrt(solve(h)[3,3])
      },
      error = function(cond) {}, # !!!!! log error
      warning = function(cond) {} # !!!!! log error
    )

    return (list(
      p_x_hat = p_x_hat,
      p_y_hat = p_y_hat,
      theta_hat = theta_hat,
      se_p_x_hat = se_p_x_hat,
      se_p_y_hat = se_p_y_hat,
      se_theta_hat = se_theta_hat
    ))

  }

  if (type %in% c("IG LM", "IG GEE")) {

    if (type=="IG LM") {

      # Run normal linear model
      model_normal_lm <- lm(
        y ~ factor(j) + x_ij,
        data = data$data
      )

      # Extract coefficients and SEs
      theta_hat <- summary(model_normal_lm)$coefficients["x_ij",1]
      se_theta_hat <- summary(model_normal_lm)$coefficients["x_ij",2]

    }

    if (type=="IG GEE") {

      # !!!!! TO DO

    }

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

  if (type=="PP SPL") {

    # !!!!! TO DO

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

  if (type=="FX SPL") {

    # !!!!! Unfinished

    # Generate ML estimates of theta and d
    nll <- function(par) {
      return (
        -1 * log_lik_spline(
          sigma_v = par[1],
          sigma_e = par[2],
          alpha = par[3],
          beta_1 = par[4],
          beta_2 = par[5],
          beta_3 = par[6],
          beta_4 = par[7],
          beta_5 = par[8],
          theta = par[9],
          p_x = par[10],
          p_y = par[11],
          g_x = 5,
          data = data
        )
      )
    }
    opt <- optim(
      par = c(
        sigma_v = 0.01,
        sigma_e = 0.01,
        alpha = 1,
        beta_1 = 0.2,
        beta_2 = 0.2,
        beta_3 = 0.2,
        beta_4 = 0.2,
        beta_5 = 0.2,
        theta = -0.1,
        p_x = 2,
        p_y = 0.5
      ),
      fn = nll
    )
    theta_hat <- opt$par[["theta"]]
    print(theta_hat)



    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

}
