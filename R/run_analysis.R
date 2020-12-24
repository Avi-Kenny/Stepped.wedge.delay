#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param analysis A list containing `type` (a character string) and
#'     optionally `params` (a list that differs in structure depending on type).
#'       Possible values of `type` include:
#'         - "HH" (Hussey & Hughes; ignores delay)
#'         - "ETI" ("exposure treatment indicators" model)
#'         - "SS" (smoothing spline model)
#'       ARCHIVED values of `type` include:
#'         - "SPL" (linear spline model)
#'         - "MSS" (monotonic smoothing spline model)
#'         - "2S LM" (two-stage with linear model in first stage)
#'         - "2S GEE" (two-stage with GEE in first stage)
#'         - "2S LMM" (two-stage with linear mixed model in first stage)
#'         - "2S HL" (two-stage with H-likelihood in first stage)
#'         - "Last" (use last time point only)
#'         - "WASH" ("washout" model)
#'       Possible values of `params` include:
#'         - For type="2S GEE", params should equal list(corr=c), where values
#'           of `c` include "exchangeable" and "independence"
#'         - For type="WASH", params should equal list(length=k), where `k` is
#'           the length of the washout period (i.e. the number of time steps to
#'           discard).
#'         - For type="2S LMM", params should equal list(REML=r), where values
#'           of `r` include TRUE (fit using REML) and FALSE (fit using ML).
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

  if (analysis$type=="HH") {

    # Run GLMM
    if (data_type=="normal") {
      model <- lmer(
        y ~ factor(j) + x_ij + (1|i),
        data = data$data
      )
    } else if (data_type=="binomial") {
      model <- glmer(
        y ~ factor(j) + x_ij + (1|i),
        data = data$data,
        family = "binomial" # binomial(link="log")
      )
    }

    # Extract estimate and SE
    ate_hat <- summary(model)$coefficients["x_ij",1] # same as theta_hat
    se_ate_hat <- summary(model)$coefficients["x_ij",2] # same as se_theta_hat

    return (list(
      ate_hat = ate_hat,
      se_ate_hat = se_ate_hat
    ))

  }

  if (analysis$type=="ETI") {

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
        family = "binomial" # binomial(link="log")
      )
    }

    coeff_names <- names(summary(model)$coefficients[,1])
    theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_l_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]
    sigma_l_hat <- as.matrix(sigma_l_hat)

    # Trapezoid sum: sum of first len-1 values plus half the last value
    len <- length(theta_l_hat)
    A <- matrix(c(rep(1/len,len-1),(1/len)/2), nrow=1)
    ate_hat <- (A %*% theta_l_hat)[1,1]
    se_ate_hat <- sqrt(A %*% sigma_l_hat %*% t(A))[1,1]

    return (list(
      ate_hat = ate_hat,
      se_ate_hat = se_ate_hat
    ))

  }

  if (analysis$type=="SS") {

    J <- L$n_time_points

    if (data_type=="normal") {
      if (analysis$params$t==1) {
        model <- gamm(
          y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
          random = list(i=~1),
          data = data$data
        )
      } else if (analysis$params$t==2) {
        model <- gamm(
          y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
            s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
          random = list(i=~1),
          data = data$data
        )
      }
    } else if (data_type=="binomial") {

      if (analysis$params$t==1) {
        model <- gamm(
          y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
          random = list(i=~1),
          data = data$data,
          family = "binomial" # binomial(link="log")
        )
      } else if (analysis$params$t==2) {
        model <- gamm(
          y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
            s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
          random = list(i=~1),
          data = data$data,
          family = "binomial" # binomial(link="log")
        )
      }

    }

    theta_l_hat <- sapply(c(1:(J-1)), function(l) {
      predict(model$gam, newdata=list(j=1, l=l), type = "terms")[2]
    })
    sigma_l_hat <- vcov(model$gam, freq=FALSE) # freq=TRUE seems to make little difference
    coeff_names <- dimnames(sigma_l_hat)[[1]]
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,4)=="s(l)"]
    coeff_names <- coeff_names[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]

    # Trapezoid sum: sum of first len-1 values plus half the last value
    len <- length(theta_l_hat)
    A <- matrix(c(rep(1/len,len-1),(1/len)/2), nrow=1)
    ate_hat <- (A %*% theta_l_hat)[1,1]
    se_ate_hat <- sqrt(A %*% sigma_l_hat %*% t(A))[1,1]

    return (list(
      ate_hat = ate_hat,
      se_ate_hat = se_ate_hat
    ))

  }

  if (analysis$type=="MON") {

    # !!!!!
    summary(lmer(
      y ~ j + (1|i),
      # y ~ factor(j) + factor(l) + (1|i),
      data = data$data
    ))

    library(rjags)
    jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta1*j[n] + alpha[i[n]], 1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }
        beta1 ~ dnorm(0.0, 1.0E-4)
        beta0 ~ dnorm(0.0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data$data$i)),
        N = nrow(data$data),
        y = data$data$y,
        i = data$data$i,
        j = data$data$j
      ),
      quiet = FALSE,
      n.chains = 1,
      n.adapt = 1000
    )
    output <- coda.samples(
      model = jm,
      variable.names = c("beta0", "beta1", "sigma", "tau"),
      n.iter = 1000,
      thin = 1
    )
    summary(output)

  }

}
