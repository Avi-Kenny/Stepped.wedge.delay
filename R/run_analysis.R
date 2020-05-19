#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param analysis_type A character string corresponding to one of the
#'     following: "2S LM" (two-stage using a simple linear model in first
#'     stage), "2S GEE EX" (two-stage using GEE in first stage; exchangeable
#'     correlation), "2S GEE ID" (two-stage using GEE in first stage;
#'     independence correlation), "2S LMM ML" (two-stage using LMM in first
#'     stage; maximum likelihood), "2S LMM REML" (two-stage using LMM in first
#'     stage; REML), "2S SPL" (two-stage using spline in second stage) "SPL 1K"
#'     (single stage linear spline with a single knot), "SPL 2K" (single stage
#'     linear spline with two knots), "IG LM" (simple linear model that ignores
#'     the time-lag effect), "IG GEE" (GEE that ignores the time-lag effect),
#'     "2S HL" (two-stage using H-likelihood in first stage), "Staircase" (new
#'     one-stage method using inequality constraints)
#' @param L Passed via simba; list of simulation levels
#' @param C Passed via simba; list of simulation constants
#' @return TO DO
#' @export
# FN: generate_dataset
run_analysis <- function(data, analysis_type, data_type, L, C) {

  if (!(data_type %in% c("normal", "binomial"))) {
    stop ("`data_type` must be either 'normal' or 'binomial'")
  }

  if (analysis_type %in% c("2S LM", "2S GEE EX", "2S GEE ID", "2S LMM ML",
                           "2S LMM REML", "2S HL")) {

    if (analysis_type=="2S LM") {

      # Run linear model
      if (data_type=="normal") {
        model <- lm(
          y ~ factor(j) + factor(l),
          data = data$data
        )
      } else if (data_type=="binomial") {
        model <- glm(
          y ~ factor(j) + factor(l),
          data = data$data,
          family = binomial(link="log")
        )
      }

      # Extract coefficients and SEs
      coeff_names <- names(model$coefficients)
      theta_l_hat <- as.numeric(model$coefficients)
      sigma_l_hat <- vcov(model)

    }

    if (analysis_type %in% c("2S GEE EX","2S GEE ID")) {

      corstr <- ifelse(
        analysis_type=="2S GEE EX",
        "exchangeable",
        "independence"
      )

      # Run GEE model with exchangeable working covariance matrix
      if (data_type=="normal") {
        model <- geeglm(
          y ~ factor(j) + factor(l),
          data = data$data,
          id = i,
          family = "gaussian",
          corstr = corstr
        )
      } else if (data_type=="binomial") {
        model <- geeglm(
          y ~ factor(j) + factor(l),
          data = data$data,
          id = i,
          family = binomial(link="log"),
          corstr = corstr
        )
      }

      # Extract estimates and covariance matrix
      coeff_names <- names(model$coefficients)
      theta_l_hat <- as.numeric(model$coefficients)
      sigma_l_hat <- model$geese$vbeta

    }

    if (analysis_type %in% c("2S LMM ML","2S LMM REML")) {

      reml <- ifelse(analysis_type=="2S LMM REML", TRUE, FALSE)

      if (analysis_type=="2S LMM REML" && data_type=="binomial") {
        # warning(paste("ML being used (rather than REML) since REML is not",
        #               "well-defined for GLMMs"))
      }

      if (data_type=="normal") {
        model <- lmer(
          y ~ factor(j) + factor(l) + (1|i),
          data = data$data,
          REML = reml
        )
      } else if (data_type=="binomial") {
        model <- glmer(
          y ~ factor(j) + factor(l) + (1|i),
          data = data$data,
          family = binomial(link="log")
        )
      }

      # Extract estimates and covariance matrix
      coeff_names <- names(summary(model)$coefficients[,1])
      theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
      sigma_l_hat <- vcov(model)

    }

    if (analysis_type=="2S HL") {

      # !!!!!

      # Run linear model
      if (data_type=="normal") {

        model <- glmmTMB(
          y ~ factor(j) + factor(l) + (1|i),
          data = data$data
        )

      } else if (data_type=="binomial") {

        model <- glmmTMB(
          y ~ factor(j) + factor(l) + (1|i),
          data = data$data,
          family = binomial(link="log"),
          REML = TRUE,
          start = list(beta=rep(-1,2*L$n_time_points-1)) # Avoids a gradient error
        )
        # as.numeric(diag(sigma_l_hat))/as.numeric(diag(sigma_l_hat2)) # !!!!! Directly compare SE estimates with REML=FALSE

      }

      # Extract coefficients and SEs
      coeff_names <- dimnames(summary(model)$coefficients$cond)[[1]]
      theta_l_hat <- as.numeric(summary(model)$coefficients$cond[,1])
      sigma_l_hat <- unname(vcov(model)[[1]])

    }

    # Truncate theta_l_hat and sigma_l_hat
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]
    sigma_l_hat <- as.matrix(sigma_l_hat)

    # Generate ML estimates of theta and d
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

  if (analysis_type=="Staircase") {

    # Run linear model
    if (data_type=="normal") {

      model <- lm(
        y ~ factor(j) + factor(l),
        data = data$data
      )

    } else if (data_type=="binomial") {

      model <- glm(
        y ~ factor(j) + factor(l),
        data = data$data,
        family = binomial(link="log")
      )

    }

    # !!!!! Currently hard-coded for num_times=7
    res <- restriktor(
      object = model,
      constraints = rbind(
        c(0,0,0,0,0,0,0,1,-1,0,0,0,0),
        c(0,0,0,0,0,0,0,0,1,-1,0,0,0),
        c(0,0,0,0,0,0,0,0,0,1,-1,0,0),
        c(0,0,0,0,0,0,0,0,0,0,1,-1,0),
        c(0,0,0,0,0,0,0,0,0,0,0,1,-1)
      ),
      rhs = c(0,0,0,0,0)
    )
    s <- summary(res)$coefficients

    return (list(
      d_hat = NA,
      theta_hat = s[nrow(s),"Estimate"],
      se_d_hat = NA,
      se_theta_hat = s[nrow(s),"Std. Error"]
    ))

  }

  if (analysis_type == "2S SPL") {

    # !!!!! Next three blocks are identical with "2S LM" method above

    # Run linear model
    if (data_type=="normal") {
      model <- lm(
        y ~ factor(j) + factor(l),
        data = data$data
      )
    } else if (data_type=="binomial") {
      model <- glm(
        y ~ factor(j) + factor(l),
        data = data$data,
        family = binomial(link="log")
      )
    }

    # Extract coefficients and SEs
    coeff_names <- names(model$coefficients)
    theta_l_hat <- as.numeric(model$coefficients)
    sigma_l_hat <- vcov(model)

    # Truncate theta_l_hat and sigma_l_hat
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]

    # Generate ML estimates of theta, p_x, and p_y
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

  if (analysis_type %in% c("IG LM", "IG GEE")) {

    if (analysis_type=="IG LM") {

      # Run linear model
      if (data_type=="normal") {
        model <- lm(
          y ~ factor(j) + x_ij,
          data = data$data
        )
      } else if (data_type=="binomial") {
        model <- glm(
          y ~ factor(j) + x_ij,
          data = data$data,
          family = binomial(link="log")
        )
      }

      # Extract coefficients and SEs
      theta_hat <- summary(model)$coefficients["x_ij",1]
      se_theta_hat <- summary(model)$coefficients["x_ij",2]

    }

    if (analysis_type=="IG GEE") {

      # !!!!! TO DO

    }

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

  if (analysis_type=="FX SPL") {

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
