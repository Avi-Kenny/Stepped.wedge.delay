#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param analysis_type A character string corresponding to one of the
#'     following: "2S LM" (two-stage using a simple linear model in first
#'     stage), "2S GEE EX" (two-stage using GEE in first stage; exchangeable
#'     correlation), "2S GEE ID" (two-stage using GEE in first stage;
#'     independence correlation), "2S LMM ML" (two-stage using LMM in first
#'     stage; maximum likelihood), "2S LMM REML" (two-stage using LMM in first
#'     stage; REML), "SPL 1K" (single stage linear spline with a single knot),
#'     "SPL 2K" (single stage linear spline with two knots), "IG LM" (simple
#'     linear model that ignores the time-lag effect), "IG GEE" (GEE that
#'     ignores the time-lag effect), "2S HL" (two-stage using H-likelihood in
#'     first stage), "Staircase" (new one-stage method using inequality
#'     constraints)
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

  if (analysis_type=="SPL 1K") {

    # Add spline covariates to dataset
    # !!!!! Currently specific to the case where J=7
    df <- data$data %>% mutate(
      s1 = pmax(0,l), # Equal to l
      s2 = pmax(0,l-3)
    )

    # Run linear model
    # !!!!! Change this to GLMM
    if (data_type=="normal") {

      model <- lm(
        y ~ factor(j) + s1 + s2,
        data = df
      )

    } else if (data_type=="binomial") {

      model <- glm(
        y ~ factor(j) + s1 + s2,
        data = df,
        family = binomial(link="log")
      )

    }

    # Extract coefficients and SEs
    coeff_names <- names(model$coefficients)
    s_hat <- as.numeric(model$coefficients)
    sigma_s_hat <- vcov(model)

    # Truncate s_hat and sigma_s_hat
    indices <- c((length(coeff_names)-1):length(coeff_names))
    coeff_names <- coeff_names[indices]
    s_hat <- s_hat[indices]
    sigma_s_hat <- sigma_s_hat[indices,indices]
    sigma_s_hat <- as.matrix(sigma_s_hat)

    # Calculate estimators
    theta_hat <- (6*s_hat[1]) + (3*s_hat[2])
    se_theta_hat <- sqrt(
      36 * sigma_s_hat[1,1] +
      9 * sigma_s_hat[1,2] +
      36 * sigma_s_hat[2,2]
    )

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
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

}
