#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param analysis A list containing `type` (a character string) and
#'     optionally `params` (a list that differs in structure depending on type).
#'       Possible values of `type` include:
#'         - "2S LM" (two-stage with linear model in first stage)
#'         - "2S GEE" (two-stage with GEE in first stage)
#'         - "2S LMM" (two-stage with linear mixed model in first stage)
#'         - "2S HL" (two-stage with H-likelihood in first stage)
#'         - "SPL" (one-stage linear spline model)
#'         - "Last" (!!!!! Testing: use last time point only !!!!!)
#'       Possible values of `params` include:
#'         - For type="2S GEE", params should equal list(corr=c), where values
#'           of `c` include "exchangeable" and "independence"
#'         - For type="2S LMM", params should equal list(REML=r), where values
#'           of `r` include TRUE (fit using REML) and FALSE (fit using ML)
#'         - For type="SPL", params show equal list(knots=k, mono=TRUE), where
#'           `k` is a vector of knots (x-coordinates; excluding the knot at x=0)
#'           and mono is TRUE/FALSE (indicating whether to force the spline to
#'           be monotonic)
#'
#'         - !!!!! Need to add ignore=TRUE to params
#' @param data_type Type of data being analyzed ("binomial" or "normal")
#' @param L Passed via simba; list of simulation levels
#' @param C Passed via simba; list of simulation constants
#' @return TO DO
#' @export
# FN: generate_dataset
run_analysis <- function(data, analysis, data_type, L, C) {

  if (analysis$type %in% c("2S LM", "2S GEE", "2S LMM", "2S HL")) {

    if (analysis$type=="2S LM") {

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

    if (analysis$type=="2S GEE") {

      if (data_type=="normal") {
        model <- geeglm(
          y ~ factor(j) + factor(l),
          data = data$data,
          id = i,
          family = "gaussian",
          corstr = analysis$params$corr
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

    if (analysis$type=="2S LMM") {

      if (analysis$params$REML && data_type=="binomial") {
        stop("REML is not well-defined for binomial GLMMs")
      }

      if (data_type=="normal") {
        model <- lmer(
          y ~ factor(j) + factor(l) + (1|i),
          data = data$data,
          REML = analysis$params$REML
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

    if (analysis$type=="2S HL") {

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
          start = list(beta=rep(-1,2*L$n_time_points-1)) # !!!!! Avoids a gradient error; perhaps implement this in the 2S LMM as well
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
      error = function(cond) {}, # !!!!! catch/process error
      warning = function(cond) {} # !!!!! catch/process error
    )

    return (list(
      d_hat = d_hat,
      theta_hat = theta_hat,
      se_d_hat = se_d_hat,
      se_theta_hat = se_theta_hat
    ))

  }

  # !!!!! Testing: Last !!!!!
  if (analysis$type=="Last") {

    # Run GLMM
    if (data_type=="normal") {
      model <- lmer(
        y ~ factor(j) + factor(l) + (1|i),
        data = data$data
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

    # Truncate theta_l_hat and sigma_l_hat
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]
    sigma_l_hat <- as.matrix(sigma_l_hat)

    # Use last theta_l_hat as estimate
    theta_hat <- theta_l_hat[length(theta_l_hat)]
    se_theta_hat <- sqrt(sigma_l_hat[nrow(sigma_l_hat),nrow(sigma_l_hat)])

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

  # !!!!! Save and repurpose this code when imposing inequality constraints

  # if (analysis$type=="Staircase") {
  #
  #   # Run linear model
  #   if (data_type=="normal") {
  #
  #     model <- lm(
  #       y ~ factor(j) + factor(l),
  #       data = data$data
  #     )
  #
  #   } else if (data_type=="binomial") {
  #
  #     model <- glm(
  #       y ~ factor(j) + factor(l),
  #       data = data$data,
  #       family = binomial(link="log")
  #     )
  #
  #   }
  #
  #   # !!!!! Currently hard-coded for num_times=7
  #   res <- restriktor(
  #     object = model,
  #     constraints = rbind(
  #       c(0,0,0,0,0,0,0,1,-1,0,0,0,0),
  #       c(0,0,0,0,0,0,0,0,1,-1,0,0,0),
  #       c(0,0,0,0,0,0,0,0,0,1,-1,0,0),
  #       c(0,0,0,0,0,0,0,0,0,0,1,-1,0),
  #       c(0,0,0,0,0,0,0,0,0,0,0,1,-1)
  #     ),
  #     rhs = c(0,0,0,0,0)
  #   )
  #   s <- summary(res)$coefficients
  #
  #   return (list(
  #     d_hat = NA,
  #     theta_hat = s[nrow(s),"Estimate"],
  #     se_d_hat = NA,
  #     se_theta_hat = s[nrow(s),"Std. Error"]
  #   ))
  #
  # }

  if (analysis$type=="SPL") {

    # Dynamically build formula
    # !!!!! monotonic spline currently only works with lm/glm
    if (analysis$params$mono) {
      formula <- "y ~ factor(j) + t0"
    } else {
      formula <- "y ~ factor(j) + (1|i) + t0"
    }

    # Add spline covariates to dataset
    k <- analysis$params$knots
    df <- data$data
    df$t0 <- df$l
    if (length(k)>1) {
      for (i in 1:(length(k)-1)) {
        df[paste0("t",i)] <- pmax(0,df$l-k[i]) - pmax(0,df$l-k[length(k)])
        formula <- paste0(formula," + t",i)
      }
    }

    # !!!!! monotonic spline currently only works with lm/glm
    if (analysis$params$mono) {

      # Run LM with spline terms
      if (data_type=="normal") {
        model <- lm(
          formula = formula,
          data = df
        )
      } else if (data_type=="binomial") {
        model <- glm(
          formula = formula,
          data = df,
          family = binomial(link="log")
        )
      }

      # !!!!! Currently hard-coded for num_times=7
      res <- restriktor(
        object = model,
        constraints = rbind(
          c(0,0,0,0,0,0,0,-1,0,0,0,0,0),
          c(0,0,0,0,0,0,0,-1,-1,0,0,0,0),
          c(0,0,0,0,0,0,0,-1,-1,-1,0,0,0),
          c(0,0,0,0,0,0,0,-1,-1,-1,-1,0,0),
          c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0),
          c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1)
        ),
        rhs = c(0,0,0,0,0,0)
      )
      s <- summary(res)

      # Extract coefficients and SEs
      coeff_names <- names(summary(res)$coefficients[,1])
      coeffs <- as.numeric(summary(res)$coefficients[,1])
      sigma_hat <- summary(res)$V

    } else {

      # Run GLMM with spline terms
      if (data_type=="normal") {
        model <- lmer(
          formula = formula,
          data = df
        )
      } else if (data_type=="binomial") {
        model <- glmer(
          formula = formula,
          data = df,
          family = binomial(link="log")
        )
      }

      # Extract estimates and covariance matrix
      coeff_names <- names(summary(model)$coefficients[,1])
      coeffs <- as.numeric(summary(model)$coefficients[,1])
      sigma_hat <- vcov(model)

    }

    # Truncate s_hat and sigma_s_hat
    indices <- c((length(coeff_names)-length(k)+1):length(coeff_names))
    coeff_names <- coeff_names[indices]
    coeffs <- coeffs[indices]
    sigma_hat <- sigma_hat[indices,indices]
    sigma_hat <- as.matrix(sigma_hat)

    # Calculate estimators
    if (length(k)>1) {
      kdiffs <- k[length(k)] - c(0, k[1:(length(k)-1)])
      theta_hat <- (kdiffs %*% coeffs)[1,1]
      se_theta_hat <- sqrt((kdiffs %*% sigma_hat %*% kdiffs)[1,1])
    } else {
      theta_hat <- k * coeffs
      se_theta_hat <- k * sqrt(sigma_hat[1,1])
    }

    return (list(
      d_hat = NA,
      theta_hat = theta_hat,
      se_d_hat = NA,
      se_theta_hat = se_theta_hat
    ))

  }

  # !!!!! Update the "ignore" code

  # if (analysis$type %in% c("IG LM", "IG GEE")) {
  #
  #   if (analysis$type=="IG LM") {
  #
  #     # Run linear model
  #     if (data_type=="normal") {
  #       model <- lm(
  #         y ~ factor(j) + x_ij,
  #         data = data$data
  #       )
  #     } else if (data_type=="binomial") {
  #       model <- glm(
  #         y ~ factor(j) + x_ij,
  #         data = data$data,
  #         family = binomial(link="log")
  #       )
  #     }
  #
  #     # Extract coefficients and SEs
  #     theta_hat <- summary(model)$coefficients["x_ij",1]
  #     se_theta_hat <- summary(model)$coefficients["x_ij",2]
  #
  #   }
  #
  #   if (analysis$type=="IG GEE") {
  #
  #     # !!!!! TO DO
  #
  #   }
  #
  #   return (list(
  #     d_hat = NA,
  #     theta_hat = theta_hat,
  #     se_d_hat = NA,
  #     se_theta_hat = se_theta_hat
  #   ))
  #
  # }

}
