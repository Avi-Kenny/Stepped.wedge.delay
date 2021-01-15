#' Run the data analysis
#'
#' @param data A dataset returned by generate_dataset()
#' @param method One of the following character strings:
#'     - "HH": Hussey & Hughes; ignores delay
#'     - "ETI": "exposure treatment indicators" model
#'     - "SS": smoothing spline model
#'     - "MCMC-SPL": linear spline (ETI) using JAGS
#'     - "MCMC-SPL-MON": monotonic linear spline (ETI) using JAGS
#' @param data_type Type of data being analyzed ("binomial" or "normal")
#' @param L Passed via simba; list of simulation levels
#' @param C Passed via simba; list of simulation constants
#' @return A list of key-value pairs, containing the following:
#'     - ate_hat: ATE estimate
#'     - se_ate_hat: Standard error of ATE estimate
#'     - lte_hat: LTE estimate
#'     - se_lte_hat: Standard error of LTE estimate

run_analysis <- function(data, method, data_type, L, C) {

  # Globally set MCMC tuning parameters
  mcmc <- list(n.adapt=6000, n.iter=2000, thin=1, n.chains=2)

  # Helper function to calculate ATE and LTE
  res <- function(theta_l_hat, sigma_l_hat) {

    # Right-hand Riemann sum (i.e. average of theta_l_hats)
    len <- length(theta_l_hat)
    A <- matrix(rep(1/len,len), nrow=1)

    return (list(
      ate_hat = (A %*% theta_l_hat)[1,1],
      se_ate_hat = sqrt(A %*% sigma_l_hat %*% t(A))[1,1],
      lte_hat = theta_l_hat[len],
      se_lte_hat = sqrt(sigma_l_hat[len,len])
    ))

  }

  if (method=="HH") {

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

  if (method=="ETI") {

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
        family = "binomial"
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

    return (res(theta_l_hat,sigma_l_hat))

  }

  if (method=="SS") {

    J <- L$n_time_points

    if (data_type=="normal") {
      model <- gamm(
        y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
        random = list(i=~1),
        data = data$data
      )
    } else if (data_type=="binomial") {
      model <- gamm(
        y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
        random = list(i=~1),
        data = data$data,
        family = "binomial"
      )
    }

    theta_l_hat <- sapply(c(1:(J-1)), function(l) {
      predict(model$gam, newdata=list(j=1, l=l), type = "terms")[2]
    })
    sigma_l_hat <- vcov(model$gam, freq=FALSE) # freq=TRUE seems to make little difference
    coeff_names <- dimnames(sigma_l_hat)[[1]]
    indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,4)=="s(l)"]
    coeff_names <- coeff_names[indices]
    sigma_l_hat <- sigma_l_hat[indices,indices]

    return (res(theta_l_hat,sigma_l_hat))

  }

  if (method=="MCMC-SPL") {

    # !!!!! Only coded for Normal data with J=7

    data_jags <- data$data
    data_jags %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_jags %<>% mutate(
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
        I = length(unique(data_jags$i)),
        N = nrow(data_jags),
        y = data_jags$y,
        i = data_jags$i,
        j_2 = data_jags$j_2,
        j_3 = data_jags$j_3,
        j_4 = data_jags$j_4,
        j_5 = data_jags$j_5,
        j_6 = data_jags$j_6,
        j_7 = data_jags$j_7,
        s_1 = data_jags$s_1,
        s_2 = data_jags$s_2,
        s_3 = data_jags$s_3,
        s_4 = data_jags$s_4,
        s_5 = data_jags$s_5,
        s_6 = data_jags$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
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

    return (res(theta_l_hat,sigma_l_hat))

  }

  if (method=="MCMC-SPL-MON") {

    # !!!!! Only coded for Normal data with J=7

    data_jags <- data$data
    data_jags %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_jags %<>% mutate(
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
        I = length(unique(data_jags$i)),
        N = nrow(data_jags),
        y = data_jags$y,
        i = data_jags$i,
        j_2 = data_jags$j_2,
        j_3 = data_jags$j_3,
        j_4 = data_jags$j_4,
        j_5 = data_jags$j_5,
        j_6 = data_jags$j_6,
        j_7 = data_jags$j_7,
        s_1 = data_jags$s_1,
        s_2 = data_jags$s_2,
        s_3 = data_jags$s_3,
        s_4 = data_jags$s_4,
        s_5 = data_jags$s_5,
        s_6 = data_jags$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
    output <- coda.samples(
      model = jm,
      variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
      # variable.names = c("beta0", "beta_j_2", "beta_j_3", "beta_j_4", "beta_j_5", "beta_j_6", "beta_j_7", "beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6", "sigma", "tau"),
      n.iter = mcmc$n.iter,
      thin = mcmc$thin
    )

    # !!!!! MCMC diagnostics
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

    return (res(theta_l_hat,sigma_l_hat))

  }

  if (method=="MCMC-STEP-MON") {

    # !!!!! Only coded for Normal data with J=7

    data_jags <- data$data
    data_jags %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    data_jags %<>% mutate(
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
        beta_s_6 <- - exp(alpha_s_6)
        beta_s_5 <- - exp(alpha_s_5)
        beta_s_4 <- - exp(alpha_s_4)
        beta_s_3 <- - exp(alpha_s_3)
        beta_s_2 <- - exp(alpha_s_2)
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
        I = length(unique(data_jags$i)),
        N = nrow(data_jags),
        y = data_jags$y,
        i = data_jags$i,
        j_2 = data_jags$j_2,
        j_3 = data_jags$j_3,
        j_4 = data_jags$j_4,
        j_5 = data_jags$j_5,
        j_6 = data_jags$j_6,
        j_7 = data_jags$j_7,
        s_1 = data_jags$s_1,
        s_2 = data_jags$s_2,
        s_3 = data_jags$s_3,
        s_4 = data_jags$s_4,
        s_5 = data_jags$s_5,
        s_6 = data_jags$s_6
      ),
      n.chains = mcmc$n.chains,
      n.adapt = mcmc$n.adapt
    )
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
      c(1,1,0,0,0,0),
      c(1,1,1,0,0,0),
      c(1,1,1,1,0,0),
      c(1,1,1,1,1,0),
      c(1,1,1,1,1,1)
    )
    theta_l_hat <- B %*% beta_s_hat
    sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

    return (res(theta_l_hat,sigma_l_hat))

  }

}
