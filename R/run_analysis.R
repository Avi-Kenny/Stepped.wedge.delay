#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param data_type Type of data being analyzed ("binomial" or "normal")
#' @param method A list of the form list(method="HH", enforce=NULL). Method
#'     is one of the following character strings:
#'       - "HH": Hussey & Hughes; ignores delay
#'       - "ETI": "exposure treatment indicators" model
#'       - "SS": smoothing spline model
#'       - "MCMC-SPL": linear spline (ETI) using JAGS
#'       - "MCMC-SPL-MON": monotonic linear spline (ETI) using JAGS
#'       - "MCMC-STEP-MON": monotonic step function (ETI) using JAGS
#'     The `enforce` argument is only used for "MCMC-STEP-MON" and controls how
#'     monotonicity is enforced. It is one of the following:
#'       - "prior; gamma prior": Monotonicity enforced via prior; gamma prior
#'       - "prior; unif prior": Monotonicity enforced via prior; uniform prior
#'       - "exp; exp prior": Monotonicity enforced via exponentiation;
#'         exponential prior
#'       - "exp; N(1,10) prior": Monotonicity enforced via exponentiation;
#'         normal prior (INLA default)
#'     The `effect_reached` integer argument is only used for "ETI" and "SS". If
#'     >0, this represents the maximum number of time points it takes to reach
#'     the maximum effect; if =0, no assumption is made about time to maximum
#'     effect
#'     The `re` argument is only used for "ETI" and "SS". It can be set to
#'     "height" or "height+shape".
#' @param return_extra An empty list, or a list of the form
#'     list(rte=TRUE, whole_curve=TRUE). Include rte=TRUE to return estimates of
#'     random effect parameters. Include whole_curve=TRUE to return estimates of
#'     the entire effect curve
#' @return A list of key-value pairs, containing the following:
#'     - ate_hat: ATE estimate
#'     - se_ate_hat: Standard error of ATE estimate
#'     - lte_hat: LTE estimate
#'     - se_lte_hat: Standard error of LTE estimate

run_analysis <- function(data, data_type, method, return_extra) {

  # Globally set MCMC tuning parameters
  if_null <- function(x,y) {ifelse(!is.null(x),x,y)}
  mcmc <- list(
    n.adapt = if_null(method$mcmc$n.adapt, 1000), # 5000
    n.iter = if_null(method$mcmc$n.iter, 1000), # 2000
    n.burn = if_null(method$mcmc$n.burn, 1000), # 2000
    thin = if_null(method$mcmc$thin, 1),
    n.chains = if_null(method$mcmc$n.chains, 2)
  )

  # Set defaults for effect_reached and re
  if (is.null(method$effect_reached)) {
    method$effect_reached <- 0
  }
  if (is.null(method$re)) {
    method$re <- NA
  }

  # Helper function to calculate ATE and LTE
  res <- function(theta_l_hat, sigma_l_hat, effect_reached) {

    J <- L$n_time_points
    R <- effect_reached
    len <- length(theta_l_hat)

    if (R>0 && method$method!="ETI") {
      stop("Error: R>0 && method$method!='ETI'")
    }

    # # Tradezoidal Riemann sum
    # if (R>0) {
    #   # Hard-coded for J=7
    #   a <- rep(0,R)
    #   for (r in c(1:R)) {
    #     a[r-1] <- a[r-1] + 0.5
    #     a[r] <- a[r] + 0.5
    #   }
    #   a[R] <- a[R] + (6-R)
    #   a <- (1/6)*a
    #   A <- t(as.matrix(a))
    # } else if (R==0) {
    #   A <- (1/(J-1)) * matrix(c(rep(1,len-1),0.5), nrow=1)
    # }

    # Right-hand Riemann sum
    if (R>0) {
      A <- (1/(J-1)) * matrix(c(rep(1,R-1),J-R))
    } else if (R==0) {
      A <- (1/(J-1)) * matrix(rep(1,len), nrow=1)
    }

    return (list(
      ate_hat = (A %*% theta_l_hat)[1,1],
      se_ate_hat = sqrt(A %*% sigma_l_hat %*% t(A))[1,1],
      lte_hat = theta_l_hat[len],
      se_lte_hat = sqrt(sigma_l_hat[len,len])
      # n_kept = NA,  # only relevant for "rejection" method
      # n_tossed = NA # only relevant for "rejection" method
    ))

  }

  if (method$method=="HH") {

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
        family = "binomial"
      )
    }

    # Extract estimates and SEs
    return (list(
      ate_hat = summary(model)$coefficients["x_ij",1],
      se_ate_hat = summary(model)$coefficients["x_ij",2],
      lte_hat = summary(model)$coefficients["x_ij",1],
      se_lte_hat = summary(model)$coefficients["x_ij",2]
    ))

  }

  if (method$method=="ETI") {

    J <- L$n_time_points

    # If n_extra_time_points>0, recode l terms
    if (data$params$n_extra_time_points>0) {
      data$data %<>% mutate(
        l = ifelse(j>J, J-1, l)
      )
    }

    # If effect_reached>0, recode l terms
    if (method$effect_reached>0) {
      data$data %<>% mutate(
        l = ifelse(l>method$effect_reached, method$effect_reached, l)
      )
    }

    # Run GLMM
    if (is.na(method$re)) {
      formula <- y ~ factor(j) + factor(l) + (1|i)
    } else if (method$re=="height") {
      formula <- y ~ factor(j) + factor(l) + (x_ij|i)
    }
    if (data_type=="normal") {
      model <- lmer(formula, data=data$data)
    } else if (data_type=="binomial") {
      model <- glmer(formula, data=data$data, family="binomial")
    }

    coeff_names <- names(summary(model)$coefficients[,1])
    theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_l_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
    coeff_names <- coeff_names[indices]
    theta_l_hat <- theta_l_hat[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]
    sigma_l_hat <- as.matrix(sigma_l_hat)

    res <- res(theta_l_hat,sigma_l_hat,method$effect_reached)

    if (is.null(return_extra$rte) && is.null(return_extra$whole_curve)) {
      return (res)
    }

    if (!is.null(return_extra$rte)) {

      if (is.na(method$re)) {
        s <- summary(model)
        return(c(res,list(
          sigma_hat = s$sigma,
          rho1_hat = NA,
          rho2_hat = NA,
          tau_hat = sqrt(s$varcor[[1]][1,1]),
          nu_hat = NA
        )))
      } else if (method$re=="height") {
        s <- summary(model)
        return(c(res,list(
          sigma_hat = s$sigma,
          rho1_hat = attr(s$varcor[[1]],"correlation")[1,2],
          rho2_hat = NA,
          tau_hat = sqrt(s$varcor[[1]][1,1]),
          nu_hat = sqrt(s$varcor[[1]][2,2])
        )))
      }

    }

    if (!is.null(return_extra$whole_curve)) {
      return(c(res,list(
        theta_1_hat = theta_l_hat[1],
        theta_2_hat = theta_l_hat[2],
        theta_3_hat = theta_l_hat[3],
        theta_4_hat = theta_l_hat[4],
        theta_5_hat = theta_l_hat[5],
        theta_6_hat = theta_l_hat[6]
      )))
    }

  }

  if (method$method=="CUBIC-4df") {

    # !!!!! INCOMPLETE; TESTING !!!!!
    J <- L$n_time_points
    data$data %<>% mutate(
      b1 = l,
      b2 = l^2,
      b3 = l^3,
      b4 = pmax(0,(l-3)^3)
    )
    formula <- y ~ factor(j) + b1+b2+b3+b4 + (1|i)
    model <- lmer(formula, data=data$data)
    coeff_names <- names(summary(model)$coefficients[,1])

    b_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_b_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,1)=="b"]
    coeff_names <- coeff_names[indices]
    b_hat <- b_hat[indices]
    sigma_b_hat <- sigma_b_hat[indices,indices]
    sigma_b_hat <- as.matrix(sigma_b_hat)

    B <- rbind(
      c(1,1^2,1^3,max(0,(1-3)^3)),
      c(2,2^2,2^3,max(0,(2-3)^3)),
      c(3,3^2,3^3,max(0,(3-3)^3)),
      c(4,4^2,4^3,max(0,(4-3)^3)),
      c(5,5^2,5^3,max(0,(5-3)^3)),
      c(6,6^2,6^3,max(0,(6-3)^3))
    )
    theta_l_hat <- as.numeric(B %*% b_hat)
    sigma_l_hat <- B %*% sigma_b_hat %*% t(B)
    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  # Notes for all MCMC methods:
  # !!!!! Only coded for Normal data with J=7
  # !!!!! Does not currently handle case when n_extra_time_points>0
  # !!!!! Does not currently handle case when effect_reached>0

  if (method$method=="MCMC-SPL") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = l,
      s_2 = pmax(0,l-1),
      s_3 = pmax(0,l-2),
      s_4 = pmax(0,l-3),
      s_5 = pmax(0,l-4),
      s_6 = pmax(0,l-5)
    )

    jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
          + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
          beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
          1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }
        beta_s_6 ~ dnorm(0, 1.0E-4)
        beta_s_5 ~ dnorm(0, 1.0E-4)
        beta_s_4 ~ dnorm(0, 1.0E-4)
        beta_s_3 ~ dnorm(0, 1.0E-4)
        beta_s_2 ~ dnorm(0, 1.0E-4)
        beta_s_1 ~ dnorm(0, 1.0E-4)
        beta_j_7 ~ dnorm(0, 1.0E-4)
        beta_j_6 ~ dnorm(0, 1.0E-4)
        beta_j_5 ~ dnorm(0, 1.0E-4)
        beta_j_4 ~ dnorm(0, 1.0E-4)
        beta_j_3 ~ dnorm(0, 1.0E-4)
        beta_j_2 ~ dnorm(0, 1.0E-4)
        beta0 ~ dnorm(0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        s_1 = data_mod$s_1,
        s_2 = data_mod$s_2,
        s_3 = data_mod$s_3,
        s_4 = data_mod$s_4,
        s_5 = data_mod$s_5,
        s_6 = data_mod$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    # Extract beta_s means
    beta_s_hat <- c()
    for (i in 1:6) {
      beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
    }

    # Construct covariance matrix of s terms
    sigma_s_hat <- matrix(NA, nrow=6, ncol=6)
    n_samp <- length(output[[1]][,1])
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_s_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",j)]}))
        )
      }
    }

    # Calculate theta_l_hat vector and sigma_l_hat matrix
    B = rbind(
      c(1,0,0,0,0,0),
      c(2,1,0,0,0,0),
      c(3,2,1,0,0,0),
      c(4,3,2,1,0,0),
      c(5,4,3,2,1,0),
      c(6,5,4,3,2,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  if (method$method=="MCMC-SPL-MON") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = l,
      s_2 = pmax(0,l-1),
      s_3 = pmax(0,l-2),
      s_4 = pmax(0,l-3),
      s_5 = pmax(0,l-4),
      s_6 = pmax(0,l-5)
    )

    jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
          + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
          beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
          1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }
        beta_s_6 <- exp(alpha_s_5) - exp(alpha_s_6)
        beta_s_5 <- exp(alpha_s_4) - exp(alpha_s_5)
        beta_s_4 <- exp(alpha_s_3) - exp(alpha_s_4)
        beta_s_3 <- exp(alpha_s_2) - exp(alpha_s_3)
        beta_s_2 <- exp(alpha_s_1) - exp(alpha_s_2)
        beta_s_1 <- - exp(alpha_s_1)
        alpha_s_6 <- log(10) - e_s_6
        alpha_s_5 <- log(10) - e_s_5
        alpha_s_4 <- log(10) - e_s_4
        alpha_s_3 <- log(10) - e_s_3
        alpha_s_2 <- log(10) - e_s_2
        alpha_s_1 <- log(10) - e_s_1
        e_s_6 ~ dexp(1)
        e_s_5 ~ dexp(1)
        e_s_4 ~ dexp(1)
        e_s_3 ~ dexp(1)
        e_s_2 ~ dexp(1)
        e_s_1 ~ dexp(1)
        beta_j_7 ~ dnorm(0, 1.0E-4)
        beta_j_6 ~ dnorm(0, 1.0E-4)
        beta_j_5 ~ dnorm(0, 1.0E-4)
        beta_j_4 ~ dnorm(0, 1.0E-4)
        beta_j_3 ~ dnorm(0, 1.0E-4)
        beta_j_2 ~ dnorm(0, 1.0E-4)
        beta0 ~ dnorm(0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        s_1 = data_mod$s_1,
        s_2 = data_mod$s_2,
        s_3 = data_mod$s_3,
        s_4 = data_mod$s_4,
        s_5 = data_mod$s_5,
        s_6 = data_mod$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    # MCMC diagnostics
    if (FALSE) {
      n_samp <- length(output[[1]][,1])
      var <- "beta_s_6" # beta0 tau
      c1 <- output[[1]][1:n_samp,var]
      c2 <- output[[2]][1:n_samp,var]
      c3 <- output[[3]][1:n_samp,var]
      c4 <- output[[4]][1:n_samp,var]
      ggplot(
        data.frame(
          x = rep(c(1:n_samp),4),
          y = c(c1,c2,c3,c4),
          chain = rep(c(1:4), each=n_samp)
        ),
        aes(x=x, y=y, color=factor(chain))) +
        geom_line() +
        labs(title=var)
    }

    # Extract beta_s means
    beta_s_hat <- c()
    for (i in 1:6) {
      beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
    }

    # Construct covariance matrix of s terms
    sigma_s_hat <- matrix(NA, nrow=6, ncol=6)
    n_samp <- length(output[[1]][,1])
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_s_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",j)]}))
        )
      }
    }

    # Calculate theta_l_hat vector and sigma_l_hat matrix
    B = rbind(
      c(1,0,0,0,0,0),
      c(2,1,0,0,0,0),
      c(3,2,1,0,0,0),
      c(4,3,2,1,0,0),
      c(5,4,3,2,1,0),
      c(6,5,4,3,2,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  if (method$method=="MCMC-RTE-height") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% dummy_cols(select_columns="l", remove_first_dummy=TRUE)

    jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] +
            beta_j_7*j_7[n] + (beta_l_1+re[i[n],2])*l_1[n] +
            (beta_l_2+re[i[n],2])*l_2[n] + (beta_l_3+re[i[n],2])*l_3[n] +
            (beta_l_4+re[i[n],2])*l_4[n] + (beta_l_5+re[i[n],2])*l_5[n] +
            (beta_l_6+re[i[n],2])*l_6[n] + re[i[n],1],
            1/(sigma^2))
          }
          for (n in 1:I) {
            re[n,1:2] ~ dmnorm.vcov(c(0,0),Sigma)
          }
          Sigma[1,1] <- tau^2
          Sigma[1,2] <- rho*tau*nu
          Sigma[2,1] <- rho*tau*nu
          Sigma[2,2] <- nu^2
          beta_l_6 ~ dnorm(0, 1.0E-4)
          beta_l_5 ~ dnorm(0, 1.0E-4)
          beta_l_4 ~ dnorm(0, 1.0E-4)
          beta_l_3 ~ dnorm(0, 1.0E-4)
          beta_l_2 ~ dnorm(0, 1.0E-4)
          beta_l_1 ~ dnorm(0, 1.0E-4)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          rho ~ dunif(-1,1)
          nu <- 1/sqrt(nu_prec)
          nu_prec ~ dgamma(1.0E-3, 1.0E-3)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        l_1 = data_mod$l_1,
        l_2 = data_mod$l_2,
        l_3 = data_mod$l_3,
        l_4 = data_mod$l_4,
        l_5 = data_mod$l_5,
        l_6 = data_mod$l_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_l_1", "beta_l_2", "beta_l_3", "beta_l_4",
                         "beta_l_5", "beta_l_6", "tau", "nu", "rho", "sigma"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    n_samp <- length(output[[1]][,1])

    if (return_extra$rte) {
      rho_hat <- summary(output)$statistics["rho","Mean"]
      tau_hat <- summary(output)$statistics["tau","Mean"]
      nu_hat <- summary(output)$statistics["nu","Mean"]
      sigma_hat <- summary(output)$statistics["sigma","Mean"]
    }

    # !!!!! MCMC diagnostics
    if (FALSE) {
      var <- "nu" # rho tau nu sigma
      c1 <- output[[1]][1:n_samp,var]
      c2 <- output[[2]][1:n_samp,var]
      c3 <- output[[3]][1:n_samp,var]
      ggplot(
        data.frame(
          x = rep(c(1:n_samp),3),
          y = c(c1,c2,c3),
          chain = rep(c(1:3), each=n_samp)
        ),
        aes(x=x, y=y, color=factor(chain))) +
        geom_line() +
        labs(title=var)
    }

    # Extract beta_s means
    theta_l_hat <- c()
    for (i in 1:6) {
      theta_l_hat[i] <- mean(
        unlist(lapply(output, function(l) {
          l[1:n_samp,paste0("beta_l_",i)]
        })),
        na.rm = TRUE
      )
    }

    # Construct covariance matrix of l terms
    sigma_l_hat <- matrix(NA, nrow=6, ncol=6)
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_l_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_l_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_l_",j)]})),
          use = "complete.obs"
        )
      }
    }

    res <- res(theta_l_hat,sigma_l_hat,method$effect_reached)
    if (return_extra$rte) {
      res$sigma_hat <- sigma_hat
      res$rho1_hat <- rho_hat
      res$rho2_hat <- NA
      res$tau_hat <- tau_hat
      res$nu_hat <- nu_hat
    }

    return (res)

  }

  if (method$method=="MCMC-RTE-height+shape") {

    print(paste("Check 3:",Sys.time())) # !!!!!

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% dummy_cols(select_columns="l", remove_first_dummy=TRUE)

    jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] +
            beta_j_7*j_7[n] + (beta_l_1+re[i[n],2])*l_1[n] +
            (beta_l_2+re[i[n],3])*l_2[n] + (beta_l_3+re[i[n],4])*l_3[n] +
            (beta_l_4+re[i[n],5])*l_4[n] + (beta_l_5+re[i[n],6])*l_5[n] +
            (beta_l_6+re[i[n],7])*l_6[n] + re[i[n],1],
            1/(sigma^2))
          }
          for (n in 1:I) {
            re[n,1:7] ~ dmnorm.vcov(rep(0,7),Sigma)
          }
          Sigma[1,1] <- tau^2
          for (i in 2:7) {
            Sigma[1,i] <- rho1*tau*nu
            Sigma[i,1] <- rho1*tau*nu
            for (j in 2:7) {
              Sigma[i,j] <- ifelse(i==j, nu^2, rho2*nu^2)
            }
          }
          beta_l_6 ~ dnorm(0, 1.0E-4)
          beta_l_5 ~ dnorm(0, 1.0E-4)
          beta_l_4 ~ dnorm(0, 1.0E-4)
          beta_l_3 ~ dnorm(0, 1.0E-4)
          beta_l_2 ~ dnorm(0, 1.0E-4)
          beta_l_1 ~ dnorm(0, 1.0E-4)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          rho1 ~ dunif(-1,1)
          rho2 ~ dunif(-1,1)
          nu <- 1/sqrt(nu_prec)
          nu_prec ~ dgamma(1.0E-3, 1.0E-3)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    print(paste("Check 4:",Sys.time())) # !!!!!

    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        l_1 = data_mod$l_1,
        l_2 = data_mod$l_2,
        l_3 = data_mod$l_3,
        l_4 = data_mod$l_4,
        l_5 = data_mod$l_5,
        l_6 = data_mod$l_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    print(paste("Check 5:",Sys.time())) # !!!!!
    update(jm, n.iter = mcmc$n.burn)
    print(paste("Check 6:",Sys.time())) # !!!!!
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_l_1", "beta_l_2", "beta_l_3", "beta_l_4",
                         "beta_l_5", "beta_l_6", "tau", "nu", "rho1", "rho2",
                         "sigma"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    n_samp <- length(output[[1]][,1])

    if (return_extra$rte) {
      rho1_hat <- summary(output)$statistics["rho1","Mean"]
      rho2_hat <- summary(output)$statistics["rho2","Mean"]
      tau_hat <- summary(output)$statistics["tau","Mean"]
      nu_hat <- summary(output)$statistics["nu","Mean"]
      sigma_hat <- summary(output)$statistics["sigma","Mean"]
    }

    # !!!!! MCMC diagnostics
    if (FALSE) {
      var <- "sigma" # rho1 rho2 tau nu sigma
      c1 <- output[[1]][1:n_samp,var]
      c2 <- output[[2]][1:n_samp,var]
      c3 <- output[[3]][1:n_samp,var]
      ggplot(
        data.frame(
          x = rep(c(1:n_samp),3),
          y = c(c1,c2,c3),
          chain = rep(c(1:3), each=n_samp)
        ),
        aes(x=x, y=y, color=factor(chain))) +
        geom_line() +
        labs(title=var)
    }

    # Extract beta_s means
    theta_l_hat <- c()
    for (i in 1:6) {
      theta_l_hat[i] <- mean(
        unlist(lapply(output, function(l) {
          l[1:n_samp,paste0("beta_l_",i)]
        })),
        na.rm = TRUE
      )
    }

    # Construct covariance matrix of l terms
    sigma_l_hat <- matrix(NA, nrow=6, ncol=6)
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_l_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_l_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_l_",j)]})),
          use = "complete.obs"
        )
      }
    }

    res <- res(theta_l_hat,sigma_l_hat,method$effect_reached)
    if (return_extra$rte) {
      res$sigma_hat <- sigma_hat
      res$rho1_hat <- rho1_hat
      res$rho2_hat <- rho2_hat
      res$tau_hat <- tau_hat
      res$nu_hat <- nu_hat
    }

    return (res)

  }

  if (method$method=="MCMC-STEP-MON") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = as.numeric(l>=1),
      s_2 = as.numeric(l>=2),
      s_3 = as.numeric(l>=3),
      s_4 = as.numeric(l>=4),
      s_5 = as.numeric(l>=5),
      s_6 = as.numeric(l>=6)
    )

    if (method$enforce %in% c("none","rejection")) {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 ~ dnorm(0, 1.0E-4)
          beta_s_5 ~ dnorm(0, 1.0E-4)
          beta_s_4 ~ dnorm(0, 1.0E-4)
          beta_s_3 ~ dnorm(0, 1.0E-4)
          beta_s_2 ~ dnorm(0, 1.0E-4)
          beta_s_1 ~ dnorm(0, 1.0E-4)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="prior; gamma prior") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- - alpha_s_6
          beta_s_5 <- - alpha_s_5
          beta_s_4 <- - alpha_s_4
          beta_s_3 <- - alpha_s_3
          beta_s_2 <- - alpha_s_2
          beta_s_1 <- - alpha_s_1
          alpha_s_6 ~ dgamma(1.0E-2, 1.0E-2)
          alpha_s_5 ~ dgamma(1.0E-2, 1.0E-2)
          alpha_s_4 ~ dgamma(1.0E-2, 1.0E-2)
          alpha_s_3 ~ dgamma(1.0E-2, 1.0E-2)
          alpha_s_2 ~ dgamma(1.0E-2, 1.0E-2)
          alpha_s_1 ~ dgamma(1.0E-2, 1.0E-2)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="prior; unif prior") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 ~ dunif(-10,0)
          beta_s_5 ~ dunif(-10,0)
          beta_s_4 ~ dunif(-10,0)
          beta_s_3 ~ dunif(-10,0)
          beta_s_2 ~ dunif(-10,0)
          beta_s_1 ~ dunif(-10,0)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; exp prior") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- exp(alpha_s_6)
          beta_s_5 <- exp(alpha_s_5)
          beta_s_4 <- exp(alpha_s_4)
          beta_s_3 <- exp(alpha_s_3)
          beta_s_2 <- exp(alpha_s_2)
          beta_s_1 <- exp(alpha_s_1)
          alpha_s_6 <- log(10) - e_s_6
          alpha_s_5 <- log(10) - e_s_5
          alpha_s_4 <- log(10) - e_s_4
          alpha_s_3 <- log(10) - e_s_3
          alpha_s_2 <- log(10) - e_s_2
          alpha_s_1 <- log(10) - e_s_1
          e_s_6 ~ dexp(1)
          e_s_5 ~ dexp(1)
          e_s_4 ~ dexp(1)
          e_s_3 ~ dexp(1)
          e_s_2 ~ dexp(1)
          e_s_1 ~ dexp(1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="new mix") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- ifelse(bern_6, e_l_6, -1*e_r_6)
          beta_s_5 <- ifelse(bern_5, e_l_5, -1*e_r_5)
          beta_s_4 <- ifelse(bern_4, e_l_4, -1*e_r_4)
          beta_s_3 <- ifelse(bern_3, e_l_3, -1*e_r_3)
          beta_s_2 <- ifelse(bern_2, e_l_2, -1*e_r_2)
          beta_s_1 <- ifelse(bern_1, e_l_1, -1*e_r_1)
          e_l_6 ~ dexp(10)
          e_l_5 ~ dexp(10)
          e_l_4 ~ dexp(10)
          e_l_3 ~ dexp(10)
          e_l_2 ~ dexp(10)
          e_l_1 ~ dexp(10)
          e_r_6 ~ dexp(0.1)
          e_r_5 ~ dexp(0.1)
          e_r_4 ~ dexp(0.1)
          e_r_3 ~ dexp(0.1)
          e_r_2 ~ dexp(0.1)
          e_r_1 ~ dexp(0.1)
          bern_6 ~ dbern(0.001)
          bern_5 ~ dbern(0.001)
          bern_4 ~ dbern(0.001)
          bern_3 ~ dbern(0.001)
          bern_2 ~ dbern(0.001)
          bern_1 ~ dbern(0.001)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; mix prior 0.1") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- ifelse(bern_6, 0, exp(alpha_s_6))
          beta_s_5 <- ifelse(bern_5, 0, exp(alpha_s_5))
          beta_s_4 <- ifelse(bern_4, 0, exp(alpha_s_4))
          beta_s_3 <- ifelse(bern_3, 0, exp(alpha_s_3))
          beta_s_2 <- ifelse(bern_2, 0, exp(alpha_s_2))
          beta_s_1 <- ifelse(bern_1, 0, exp(alpha_s_1))
          alpha_s_6 <- log(10) - e_s_6
          alpha_s_5 <- log(10) - e_s_5
          alpha_s_4 <- log(10) - e_s_4
          alpha_s_3 <- log(10) - e_s_3
          alpha_s_2 <- log(10) - e_s_2
          alpha_s_1 <- log(10) - e_s_1
          bern_6 ~ dbern(0.1)
          bern_5 ~ dbern(0.1)
          bern_4 ~ dbern(0.1)
          bern_3 ~ dbern(0.1)
          bern_2 ~ dbern(0.1)
          bern_1 ~ dbern(0.1)
          e_s_6 ~ dexp(1)
          e_s_5 ~ dexp(1)
          e_s_4 ~ dexp(1)
          e_s_3 ~ dexp(1)
          e_s_2 ~ dexp(1)
          e_s_1 ~ dexp(1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; mix prior 0.2") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- ifelse(bern_6, 0, exp(alpha_s_6))
          beta_s_5 <- ifelse(bern_5, 0, exp(alpha_s_5))
          beta_s_4 <- ifelse(bern_4, 0, exp(alpha_s_4))
          beta_s_3 <- ifelse(bern_3, 0, exp(alpha_s_3))
          beta_s_2 <- ifelse(bern_2, 0, exp(alpha_s_2))
          beta_s_1 <- ifelse(bern_1, 0, exp(alpha_s_1))
          alpha_s_6 <- log(10) - e_s_6
          alpha_s_5 <- log(10) - e_s_5
          alpha_s_4 <- log(10) - e_s_4
          alpha_s_3 <- log(10) - e_s_3
          alpha_s_2 <- log(10) - e_s_2
          alpha_s_1 <- log(10) - e_s_1
          bern_6 ~ dbern(0.2)
          bern_5 ~ dbern(0.2)
          bern_4 ~ dbern(0.2)
          bern_3 ~ dbern(0.2)
          bern_2 ~ dbern(0.2)
          bern_1 ~ dbern(0.2)
          e_s_6 ~ dexp(1)
          e_s_5 ~ dexp(1)
          e_s_4 ~ dexp(1)
          e_s_3 ~ dexp(1)
          e_s_2 ~ dexp(1)
          e_s_1 ~ dexp(1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; N(1,10) mix (0.1)") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- ifelse(bern_6, 0, -exp(e_s_6))
          beta_s_5 <- ifelse(bern_5, 0, -exp(e_s_5))
          beta_s_4 <- ifelse(bern_4, 0, -exp(e_s_4))
          beta_s_3 <- ifelse(bern_3, 0, -exp(e_s_3))
          beta_s_2 <- ifelse(bern_2, 0, -exp(e_s_2))
          beta_s_1 <- ifelse(bern_1, 0, -exp(e_s_1))
          bern_6 ~ dbern(0.1)
          bern_5 ~ dbern(0.1)
          bern_4 ~ dbern(0.1)
          bern_3 ~ dbern(0.1)
          bern_2 ~ dbern(0.1)
          bern_1 ~ dbern(0.1)
          e_s_6 ~ dnorm(1, 0.1)
          e_s_5 ~ dnorm(1, 0.1)
          e_s_4 ~ dnorm(1, 0.1)
          e_s_3 ~ dnorm(1, 0.1)
          e_s_2 ~ dnorm(1, 0.1)
          e_s_1 ~ dnorm(1, 0.1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; N(1,10) mix (0.2)") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- ifelse(bern_6, 0, -exp(e_s_6))
          beta_s_5 <- ifelse(bern_5, 0, -exp(e_s_5))
          beta_s_4 <- ifelse(bern_4, 0, -exp(e_s_4))
          beta_s_3 <- ifelse(bern_3, 0, -exp(e_s_3))
          beta_s_2 <- ifelse(bern_2, 0, -exp(e_s_2))
          beta_s_1 <- ifelse(bern_1, 0, -exp(e_s_1))
          bern_6 ~ dbern(0.2)
          bern_5 ~ dbern(0.2)
          bern_4 ~ dbern(0.2)
          bern_3 ~ dbern(0.2)
          bern_2 ~ dbern(0.2)
          bern_1 ~ dbern(0.2)
          e_s_6 ~ dnorm(1, 0.1)
          e_s_5 ~ dnorm(1, 0.1)
          e_s_4 ~ dnorm(1, 0.1)
          e_s_3 ~ dnorm(1, 0.1)
          e_s_2 ~ dnorm(1, 0.1)
          e_s_1 ~ dnorm(1, 0.1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; N(1,10) prior") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- - exp(e_s_6)
          beta_s_5 <- - exp(e_s_5)
          beta_s_4 <- - exp(e_s_4)
          beta_s_3 <- - exp(e_s_3)
          beta_s_2 <- - exp(e_s_2)
          beta_s_1 <- - exp(e_s_1)
          e_s_6 ~ dnorm(1, 0.1)
          e_s_5 ~ dnorm(1, 0.1)
          e_s_4 ~ dnorm(1, 0.1)
          e_s_3 ~ dnorm(1, 0.1)
          e_s_2 ~ dnorm(1, 0.1)
          e_s_1 ~ dnorm(1, 0.1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    if (method$enforce=="exp; N(0,10) prior") {

      jags_code <- quote("
        model {
          for (n in 1:N) {
            y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
            + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
            beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
            1/(sigma^2))
          }
          for (n in 1:I) {
            alpha[n] ~ dnorm(0, 1/(tau^2))
          }
          beta_s_6 <- - exp(e_s_6)
          beta_s_5 <- - exp(e_s_5)
          beta_s_4 <- - exp(e_s_4)
          beta_s_3 <- - exp(e_s_3)
          beta_s_2 <- - exp(e_s_2)
          beta_s_1 <- - exp(e_s_1)
          e_s_6 ~ dnorm(0, 0.1)
          e_s_5 ~ dnorm(0, 0.1)
          e_s_4 ~ dnorm(0, 0.1)
          e_s_3 ~ dnorm(0, 0.1)
          e_s_2 ~ dnorm(0, 0.1)
          e_s_1 ~ dnorm(0, 0.1)
          beta_j_7 ~ dnorm(0, 1.0E-4)
          beta_j_6 ~ dnorm(0, 1.0E-4)
          beta_j_5 ~ dnorm(0, 1.0E-4)
          beta_j_4 ~ dnorm(0, 1.0E-4)
          beta_j_3 ~ dnorm(0, 1.0E-4)
          beta_j_2 ~ dnorm(0, 1.0E-4)
          beta0 ~ dnorm(0, 1.0E-4)
          tau <- 1/sqrt(tau_prec)
          tau_prec ~ dgamma(1.0E-3, 1.0E-3)
          sigma <- 1/sqrt(sigma_prec)
          sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
        }
      ")

    }

    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        s_1 = data_mod$s_1,
        s_2 = data_mod$s_2,
        s_3 = data_mod$s_3,
        s_4 = data_mod$s_4,
        s_5 = data_mod$s_5,
        s_6 = data_mod$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4",
                         "beta_s_5", "beta_s_6"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    n_samp <- length(output[[1]][,1])

    if (method$enforce=="rejection") {

      # Throw out samples that don't satisfy parameter constraints
      # tol=0.2:  616/1,000 kept
      # tol=0.1:  105/1,000 kept
      # tol=0.05: 18/1,000 kept
      # output=o
      n_kept <- 0
      n_tossed <- 0
      for (i in 1:mcmc$n.chains) {
        for (j in 1:n_samp) {
          if (max(output[[i]][j,])>method$tolerance) {
            output[[i]][j,] <- NA
            n_tossed <- n_tossed + 1
          } else {
            n_kept <- n_kept + 1
          }
        }
      }
      if (n_kept==0) {
        stop("No posterior samples were kept")
      }

    }

    # Extract beta_s means
    beta_s_hat <- c()
    for (i in 1:6) {
      beta_s_hat[i] <- mean(
        unlist(lapply(output, function(l) {
          l[1:n_samp,paste0("beta_s_",i)]
        })),
        na.rm = TRUE
      )
    }

    # Construct covariance matrix of s terms
    sigma_s_hat <- matrix(NA, nrow=6, ncol=6)
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_s_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",j)]})),
          use = "complete.obs"
        )
      }
    }

    # Calculate theta_l_hat vector and sigma_l_hat matrix
    B = rbind(
      c(1,0,0,0,0,0),
      c(1,1,0,0,0,0),
      c(1,1,1,0,0,0),
      c(1,1,1,1,0,0),
      c(1,1,1,1,1,0),
      c(1,1,1,1,1,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    res <- res(theta_l_hat,sigma_l_hat,method$effect_reached)
    if (method$enforce=="rejection") {
      res$n_kept <- n_kept
      res$n_tossed <- n_tossed
    }

    if (is.null(return_extra$whole_curve)) {
      return (res)
    } else if (return_extra$whole_curve) {
      return(c(res,list(
        theta_1_hat = theta_l_hat[1],
        theta_2_hat = theta_l_hat[2],
        theta_3_hat = theta_l_hat[3],
        theta_4_hat = theta_l_hat[4],
        theta_5_hat = theta_l_hat[5],
        theta_6_hat = theta_l_hat[6]
      )))
    }

  }

  if (method$method=="MCMC-STEP-PAVA") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = as.numeric(l>=1),
      s_2 = as.numeric(l>=2),
      s_3 = as.numeric(l>=3),
      s_4 = as.numeric(l>=4),
      s_5 = as.numeric(l>=5),
      s_6 = as.numeric(l>=6)
    )

    jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
          + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
          beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
          1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }
        beta_s_6 ~ dnorm(0, 1.0E-4)
        beta_s_5 ~ dnorm(0, 1.0E-4)
        beta_s_4 ~ dnorm(0, 1.0E-4)
        beta_s_3 ~ dnorm(0, 1.0E-4)
        beta_s_2 ~ dnorm(0, 1.0E-4)
        beta_s_1 ~ dnorm(0, 1.0E-4)
        beta_j_7 ~ dnorm(0, 1.0E-4)
        beta_j_6 ~ dnorm(0, 1.0E-4)
        beta_j_5 ~ dnorm(0, 1.0E-4)
        beta_j_4 ~ dnorm(0, 1.0E-4)
        beta_j_3 ~ dnorm(0, 1.0E-4)
        beta_j_2 ~ dnorm(0, 1.0E-4)
        beta0 ~ dnorm(0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")

    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        s_1 = data_mod$s_1,
        s_2 = data_mod$s_2,
        s_3 = data_mod$s_3,
        s_4 = data_mod$s_4,
        s_5 = data_mod$s_5,
        s_6 = data_mod$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4",
                         "beta_s_5", "beta_s_6"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    # Run PAVA algorithm on each element of posterior sample
    n_samp <- length(output[[1]][,1])
    B = rbind(
      c(1,0,0,0,0,0),
      c(1,1,0,0,0,0),
      c(1,1,1,0,0,0),
      c(1,1,1,1,0,0),
      c(1,1,1,1,1,0),
      c(1,1,1,1,1,1)
    )
    if (method$wts=="equal") {
      wts <- rep(1,6)
    }
    if (method$wts=="samp_size") {
      wts <- c(sum(data_mod$s_1), sum(data_mod$s_2), sum(data_mod$s_3),
               sum(data_mod$s_4),sum(data_mod$s_5),sum(data_mod$s_6))
    }
    if (method$wts=="sqrt_samp_size") {
      wts <- sqrt(c(sum(data_mod$s_1), sum(data_mod$s_2), sum(data_mod$s_3),
                    sum(data_mod$s_4),sum(data_mod$s_5),sum(data_mod$s_6)))
    }
    theta_l_hat_sample <- lapply(output, function(l) {

      return(
        t(apply(l, 1, function(beta_s_hat) {

          theta_l_hat <- B %*% beta_s_hat

          return(pava(
            y = theta_l_hat,
            w = wts,
            decreasing = TRUE
          ))

        }))
      )

    })

    # Extract theta_l_hat means
    theta_l_hat <- c()
    for (i in 1:6) {
      theta_l_hat[i] <- mean(
        unlist(lapply(theta_l_hat_sample, function(l) {
          l[1:n_samp,i]
        }))
      )
    }

    # Construct covariance matrix
    sigma_l_hat <- matrix(NA, nrow=6, ncol=6)
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_l_hat[i,j] <- cov(
          unlist(lapply(theta_l_hat_sample, function(l) {l[1:n_samp,i]})),
          unlist(lapply(theta_l_hat_sample, function(l) {l[1:n_samp,j]}))
        )
      }
    }

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  if (method$method=="MCMC-Shively-Sager") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% dummy_cols(select_columns="l", remove_first_dummy=TRUE)

    jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
          + theta_l_1*l_1[n] + theta_l_2*l_2[n] + theta_l_3*l_3[n] +
          theta_l_4*l_4[n] + theta_l_5*l_5[n] + theta_l_6*l_6[n] + alpha[i[n]],
          1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }

        theta_l_6 <- theta_l_5 - (exp(ss_beta+6*w_5-5*w_6)/(w_6-w_5)) *
                     (exp(6*(w_6-w_5))-exp(5*(w_6-w_5)))
        theta_l_5 <- theta_l_4 - (exp(ss_beta+5*w_4-4*w_5)/(w_5-w_4)) *
                     (exp(5*(w_5-w_4))-exp(4*(w_5-w_4)))
        theta_l_4 <- theta_l_3 - (exp(ss_beta+4*w_3-3*w_4)/(w_4-w_3)) *
                     (exp(4*(w_4-w_3))-exp(3*(w_4-w_3)))
        theta_l_3 <- theta_l_2 - (exp(ss_beta+3*w_2-2*w_3)/(w_3-w_2)) *
                     (exp(3*(w_3-w_2))-exp(2*(w_3-w_2)))
        theta_l_2 <- theta_l_1 - (exp(ss_beta+2*w_1-1*w_2)/(w_2-w_1)) *
                     (exp(2*(w_2-w_1))-exp(1*(w_2-w_1)))
        theta_l_1 <- -1 * (exp(ss_beta)/w_1) * (exp(w_1)-1)
        w_6 <- w_5 + n_6
        w_5 <- w_4 + n_5
        w_4 <- w_3 + n_4
        w_3 <- w_2 + n_3
        w_2 <- w_1 + n_2
        w_1 ~ dnorm(1.0E-6, ss_tau2_prec)
        n_6 ~ dnorm(1.0E-6, ss_tau2_prec)
        n_5 ~ dnorm(1.0E-6, ss_tau2_prec)
        n_4 ~ dnorm(1.0E-6, ss_tau2_prec)
        n_3 ~ dnorm(1.0E-6, ss_tau2_prec)
        n_2 ~ dnorm(1.0E-6, ss_tau2_prec)
        ss_tau2_prec <- 1/ss_tau2
        ss_tau2 ~ dunif(0,1.0E3)
        ss_beta ~ dnorm(0, 1.0E-2)

        beta_j_7 ~ dnorm(0, 1.0E-4)
        beta_j_6 ~ dnorm(0, 1.0E-4)
        beta_j_5 ~ dnorm(0, 1.0E-4)
        beta_j_4 ~ dnorm(0, 1.0E-4)
        beta_j_3 ~ dnorm(0, 1.0E-4)
        beta_j_2 ~ dnorm(0, 1.0E-4)
        beta0 ~ dnorm(0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")

    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        I = length(unique(data_mod$i)),
        N = nrow(data_mod),
        y = data_mod$y,
        i = data_mod$i,
        j_2 = data_mod$j_2,
        j_3 = data_mod$j_3,
        j_4 = data_mod$j_4,
        j_5 = data_mod$j_5,
        j_6 = data_mod$j_6,
        j_7 = data_mod$j_7,
        l_1 = data_mod$l_1,
        l_2 = data_mod$l_2,
        l_3 = data_mod$l_3,
        l_4 = data_mod$l_4,
        l_5 = data_mod$l_5,
        l_6 = data_mod$l_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    update(jm, n.iter = mcmc$n.burn)
    output <- coda.samples(
      model = jm,
      # variable.names = c("theta_l_1", "theta_l_2", "theta_l_3", "theta_l_4",
      #                    "theta_l_5", "theta_l_6"),
      variable.names = c("theta_l_1", "theta_l_2", "theta_l_3", "theta_l_4",
                         "theta_l_5", "theta_l_6", "ss_tau2", "ss_beta"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    # !!!!! MCMC diagnostics
    if (FALSE) {
      n_samp <- length(output[[1]][,1])
      var <- "ss_tau2" # theta_l_6 ss_beta ss_tau2
      c1 <- output[[1]][1:n_samp,var]
      c2 <- output[[2]][1:n_samp,var]
      c3 <- output[[2]][1:n_samp,var]
      ggplot(
        data.frame(
          x = rep(c(1:n_samp),3),
          y = c(c1,c2,c3),
          chain = rep(c(1:3), each=n_samp)
        ),
        aes(x=x, y=y, color=factor(chain))) +
        geom_line() +
        labs(title=var)
    }

    # Extract theta_l_hat means
    n_samp <- length(output[[1]][,1])
    theta_l_hat <- c()
    for (i in 1:6) {
      theta_l_hat[i] <- mean(
        unlist(lapply(output, function(l) {
          l[1:n_samp,paste0("theta_l_",i)]
        }))
      )
    }

    # Construct covariance matrix of s terms
    sigma_l_hat <- matrix(NA, nrow=6, ncol=6)
    for (i in 1:6) {
      for (j in 1:6) {
        sigma_l_hat[i,j] <- cov(
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("theta_l_",i)]})),
          unlist(lapply(output, function(l) {l[1:n_samp,paste0("theta_l_",j)]}))
        )
      }
    }

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  if (method$method=="STEP-INLA") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = as.numeric(l>=1),
      s_2 = as.numeric(l>=2),
      s_3 = as.numeric(l>=3),
      s_4 = as.numeric(l>=4),
      s_5 = as.numeric(l>=5),
      s_6 = as.numeric(l>=6)
    )

    # Run INLA model
    model_inla <- inla(
      y ~ j_2 + j_3 + j_4 + j_5 + j_6 + j_7 + s_1 + s_2 + s_3 + s_4 + s_5 +
        s_6 + f(i, model="iid"),
      family = "gaussian",
      data = data_mod
    )

    # Extract beta_s means
    beta_s_hat <- c()
    for (i in 1:6) {
      beta_s_hat[i] <- model_inla$summary.fixed[paste0("s_",i),"mean"]
    }

    # Construct covariance matrix of s terms
    # !!!!! TO DO

    # Calculate theta_l_hat vector and sigma_l_hat matrix
    B = rbind(
      c(1,0,0,0,0,0),
      c(1,1,0,0,0,0),
      c(1,1,1,0,0,0),
      c(1,1,1,1,0,0),
      c(1,1,1,1,1,0),
      c(1,1,1,1,1,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

  if (method$method=="STEP-MON-INLA") {

    data_mod <- data$data
    data_mod %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_mod %<>% mutate(
      s_1 = as.numeric(l>=1),
      s_2 = as.numeric(l>=2),
      s_3 = as.numeric(l>=3),
      s_4 = as.numeric(l>=4),
      s_5 = as.numeric(l>=5),
      s_6 = as.numeric(l>=6)
    )

    # Run INLA model
    model_inla <- inla(
      y ~ j_2 + j_3 + j_4 + j_5 + j_6 + j_7 +
        f(s_1, model="clinear", range=c(-10,0)) +
        f(s_2, model="clinear", range=c(-10,0)) +
        f(s_3, model="clinear", range=c(-10,0)) +
        f(s_4, model="clinear", range=c(-10,0)) +
        f(s_5, model="clinear", range=c(-10,0)) +
        f(s_6, model="clinear", range=c(-10,0)) +
        f(i, model="iid"),
      family = "gaussian",
      data = data_mod
    )

    # Extract beta_s means
    beta_s_hat <- c()
    for (i in 1:6) {
      beta_s_hat[i] <- model_inla$summary.hyperpar[paste0("Beta for s_",i),
                                                   "mean"]
    }

    # Construct covariance matrix of s terms
    # !!!!! TO DO

    # Calculate theta_l_hat vector and sigma_l_hat matrix
    B = rbind(
      c(1,0,0,0,0,0),
      c(1,1,0,0,0,0),
      c(1,1,1,0,0,0),
      c(1,1,1,1,0,0),
      c(1,1,1,1,1,0),
      c(1,1,1,1,1,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    return (res(theta_l_hat,sigma_l_hat,method$effect_reached))

  }

}
