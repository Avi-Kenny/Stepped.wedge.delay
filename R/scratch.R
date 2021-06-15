
# 2 x 3 stepped wedge
{
  mu <- 2
  beta2 <- -5
  beta3 <- 10
  delta <- 1
  df <- data.frame(
    i = integer(),
    j_2 = integer(),
    j_3 = integer(),
    x_ij = integer(),
    y = double()
  )
  for (i in 1:2) {
    for (j in 1:3) {
      for (k in 1:10) {
        x_ij <- as.integer((i==1 && j>=2) || (i==2 && j>=3))
        j_2 <- as.integer(j==2)
        j_3 <- as.integer(j==3)
        m <- mu + j_2*beta2 + j_3*beta3 + x_ij*delta
        y <- rnorm(n=1, mean=m, sd=0.1)
        df[nrow(df)+1,] <- c(i, j_2, j_3, x_ij, y)
      }
    }
  }
  summary(lm(y~j_2+j_3+x_ij, data=df))
  mean(filter(df, i==1 & j_2==1)$y) - mean(filter(df, i==2 & j_2==1)$y)
}

# Dirichlet
x=rdirichlet(n=10000, alpha=c(4,4,1,1))
x=rdirichlet(n=10000, alpha=c(400,400,100,100))
x=rdirichlet(n=10000, alpha=c(.04,.04,.01,.01))
c(mean(x[,1]),mean(x[,2]),mean(x[,3]),mean(x[,4]))



# Test integration to reach theta_1
beta=0.3
w_1=0.9
w_2=-0.5
integrate( function(t) { exp(beta+w_1*t) }, lower=0, upper=1)$value
(exp(beta)/w_1) * (exp(w_1)-1)
integrate( function(t) {
  exp(beta + (w_2-w_1)*t + (2*w_1-1*w_2))
}, lower=1, upper=2)$value
(exp(beta+2*w_1-1*w_2)/(w_2-w_1)) * (exp(2*(w_2-w_1))-exp(1*(w_2-w_1)))

# theta_l_2 <- theta_l_1 + (exp(ss_beta+w_1)/(w_2-w_1)) *
#   (exp(2*(w_2-w_1))-exp(1*(w_2-w_1)))
# theta_l_1 <-



# Diagnose rejection method
{
  # Set up objects
  B = rbind(c(1,0,0,0,0,0),c(1,1,0,0,0,0),c(1,1,1,0,0,0),c(1,1,1,1,0,0),c(1,1,1,1,1,0),c(1,1,1,1,1,1))
  df <- data.frame("i"=integer(), "l"=integer(),"theta_l"=double(),"which"=character(),stringsAsFactors=FALSE)

  # Unrejected
  tol1 <- 0.1
  tol2 <- 0.2
  for (i in 1:(length(output[[1]][,1]))) {
    beta_s <- as.numeric(output[[1]][i,])
    theta_l <- round(as.numeric(B %*% beta_s),4)
    keep1 <- as.numeric(!(max(beta_s)>tol1))
    keep2 <- as.numeric(!(sum(pmax(beta_s,0))>tol2))
    for (j in 1:6) {
      df[nrow(df)+1,] <- list(i,j,theta_l[j],"all")
      if (keep1==1) {
        df[nrow(df)+1,] <- list(i,j,theta_l[j],"truncated 1")
      }
      if (keep2==1) {
        df[nrow(df)+1,] <- list(i,j,theta_l[j],"truncated 2")
      }
    }
  }

  # Group averages
  avgs <- df %>% group_by(l,which) %>% summarize(
    theta_l = mean(theta_l)
  )

  # True values of theta_l
  true_theta_ls <- -0.5 * effect_curve(x=c(1:6),type="exp",params=list(d=1.5))

  ggplot(df, aes(x=l, y=theta_l, group=i)) +
    facet_wrap(~which, ncol=3) +
    geom_line(alpha=0.1) +
    geom_line(aes(x=l,y=y), data.frame(l=c(1:6),y=true_theta_ls), color="purple", size=1) +
    geom_line(aes(x=l,y=theta_l), avgs, color="orange", size=1)

}







library(dplyr)
library(lme4)
library(rjags)

generate_dataset <- function(beta0, tau, theta, n_clusters, n_time_points,
                             n_ind_per_cluster, sigma, delay_model) {

  # Generate data frame
  data <- data.frame(
    "i" = integer(), # cluster
    "j" = integer(), # time; 1=baseline, J=endline
    "k" = integer(), # individual
    "l" = integer(), # time since intervention
    "v_i" = double(), # cluster random effect
    "y_ij" = double(), # cluster-level probability or mean
    "x_ij" = integer(), # treatment state indicator
    "c_i" = integer(), # the start time of the treatment
    "y" = integer() # binary outcome
  )

  # Generate crossover times (assumes a "balanced and complete" design)
  n_clust_per_time <- n_clusters/(n_time_points-1)
  crossover_times <- rep(2:n_time_points, each=n_clust_per_time)

  # Create beta_js (linear time trend from 0 down to -0.5)
  # Main constraint is that beta_1=0
  beta_js <- sapply(1:n_time_points, function(j){
    ((1-j)/(n_time_points-1)) * 0.5
  })

  # Create theta_ls (intervention effects) based on continuous fn "delay_model"
  theta_ls <- theta * effect_curve(
    x = 1:(n_time_points-1),
    type = delay_model$type,
    params = delay_model$params
  )

  # Loop through clusters, time, and individuals
  for (i in 1:n_clusters) {

    v_i <- rnorm(1, mean=0, sd=tau)
    c_i <- crossover_times[i]-1

    for (j in 1:(n_time_points)) {

      x_ij <- ifelse(j<crossover_times[i], 0, 1)
      l <- ifelse(j<crossover_times[i], 0, (j-crossover_times[i])+1)
      theta_l <- ifelse(l>0, theta_ls[l], 0)

      # y_ij <- beta0 + beta_js[j] + theta*x_ij + v_i
      y_ij <- beta0 + theta_l*x_ij
      k <- n_ind_per_cluster
      y <- rnorm(k, mean=y_ij, sd=sigma)
      data <- rbind(data, data.frame(cbind(
        i=rep(i,k), j=rep(j,k), k=c(1:k), l=rep(l,k), v_i=rep(v_i,k),
        y_ij=rep(y_ij,k), x_ij=rep(x_ij,k), c_i=c_i, y=y
      )))

    }

  }

  params <- as.list(match.call())
  params$crossover_times <- crossover_times

  return (list(
    "params" = params,
    "beta_js" = beta_js,
    "data" = data
  ))

}



# Generate data and run model
{
  data <- generate_dataset(
    beta0 = 0,
    tau = 0,
    theta = -0.5,
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    sigma = 1,
    delay_model = list(
      type = "exp",
      params = list(d=1.5)
    )
  )

  model <- lm(
    y ~ factor(l),
    # y ~ factor(j) + factor(l),
    data = data$data
  )
  coeff_names <- names(summary(model)$coefficients[,1])
  theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
  indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  coeff_names <- coeff_names[indices]
  theta_l_hat <- theta_l_hat[indices]
  print(theta_l_hat)

}
#





# Run JAGS model: exponentiation with uniform effective prior
{
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
          y[n] ~ dnorm(beta0 + beta_s_1*s_1[n] + beta_s_2*s_2[n] +
          beta_s_3*s_3[n] + beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n],
          1/(sigma^2))
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
        beta0 ~ dnorm(0, 1.0E-4)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
  jm <- jags.model(
    file = textConnection(jags_code),
    data = list(
      N = nrow(data_jags),
      y = data_jags$y,
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
    n.chains = 1,
    n.adapt = 1000
  )
  output <- coda.samples(
    model = jm,
    variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
    n.iter = 1000,
    thin = 1
  )

  # Extract beta_s means and calculate theta_l_hat
  beta_s_hat <- c()
  for (i in 1:6) {
    beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
  }
  B = rbind(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,1,1,0,0,0),
    c(1,1,1,1,0,0),
    c(1,1,1,1,1,0),
    c(1,1,1,1,1,1)
  )
  theta_l_hat_mon0 <- B %*% beta_s_hat
  print(as.numeric(theta_l_hat_mon1))
}



# Run JAGS model: exponentiation with mixture prior that puts point mass at zero
{
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
          y[n] ~ dnorm(beta0 + beta_s_1*s_1[n] + beta_s_2*s_2[n] +
          beta_s_3*s_3[n] + beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n],
          1/(sigma^2))
        }
        beta_s_6 <- ifelse(bern_6, 0, -exp(alpha_s_6))
        beta_s_5 <- ifelse(bern_5, 0, -exp(alpha_s_5))
        beta_s_4 <- ifelse(bern_4, 0, -exp(alpha_s_4))
        beta_s_3 <- ifelse(bern_3, 0, -exp(alpha_s_3))
        beta_s_2 <- ifelse(bern_2, 0, -exp(alpha_s_2))
        beta_s_1 <- ifelse(bern_1, 0, -exp(alpha_s_1))
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
        beta0 ~ dnorm(0, 1.0E-4)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
  jm <- jags.model(
    file = textConnection(jags_code),
    data = list(
      N = nrow(data_jags),
      y = data_jags$y,
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
    n.chains = 1,
    n.adapt = 1000
  )
  output <- coda.samples(
    model = jm,
    variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
    n.iter = 1000,
    thin = 1
  )

  # Extract beta_s means and calculate theta_l_hat
  beta_s_hat <- c()
  for (i in 1:6) {
    beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
  }
  B = rbind(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,1,1,0,0,0),
    c(1,1,1,1,0,0),
    c(1,1,1,1,1,0),
    c(1,1,1,1,1,1)
  )
  theta_l_hat_mon4 <- B %*% beta_s_hat
  print(as.numeric(theta_l_hat_mon1))
}



# Run JAGS model: dunif(0,10)
{
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
          y[n] ~ dnorm(beta0 + beta_s_1*s_1[n] + beta_s_2*s_2[n] +
          beta_s_3*s_3[n] + beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n],
          1/(sigma^2))
        }
        beta_s_6 <- - alpha_s_6
        beta_s_5 <- - alpha_s_5
        beta_s_4 <- - alpha_s_4
        beta_s_3 <- - alpha_s_3
        beta_s_2 <- - alpha_s_2
        beta_s_1 <- - alpha_s_1
        alpha_s_6 ~ dunif(0,10)
        alpha_s_5 ~ dunif(0,10)
        alpha_s_4 ~ dunif(0,10)
        alpha_s_3 ~ dunif(0,10)
        alpha_s_2 ~ dunif(0,10)
        alpha_s_1 ~ dunif(0,10)
        beta0 ~ dnorm(0, 1.0E-4)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
  jm <- jags.model(
    file = textConnection(jags_code),
    data = list(
      N = nrow(data_jags),
      y = data_jags$y,
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
    n.chains = 1,
    n.adapt = 1000
  )
  output <- coda.samples(
    model = jm,
    variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
    n.iter = 1000,
    thin = 1
  )

  # Extract beta_s means and calculate theta_l_hat
  beta_s_hat <- c()
  for (i in 1:6) {
    beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
  }
  B = rbind(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,1,1,0,0,0),
    c(1,1,1,1,0,0),
    c(1,1,1,1,1,0),
    c(1,1,1,1,1,1)
  )
  theta_l_hat_mon1 <- B %*% beta_s_hat
  print(as.numeric(theta_l_hat_mon1))
}



# Run JAGS model: gamma
{
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
          y[n] ~ dnorm(beta0 + beta_s_1*s_1[n] + beta_s_2*s_2[n] +
          beta_s_3*s_3[n] + beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n],
          1/(sigma^2))
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
        beta0 ~ dnorm(0, 1.0E-4)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
  jm <- jags.model(
    file = textConnection(jags_code),
    data = list(
      N = nrow(data_jags),
      y = data_jags$y,
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
    n.chains = 1,
    n.adapt = 1000
  )
  output <- coda.samples(
    model = jm,
    variable.names = c("beta_s_1", "beta_s_2", "beta_s_3", "beta_s_4", "beta_s_5", "beta_s_6"),
    n.iter = 1000,
    thin = 1
  )

  # Extract beta_s means and calculate theta_l_hat
  beta_s_hat <- c()
  for (i in 1:6) {
    beta_s_hat[i] <- summary(output)$statistics[paste0("beta_s_",i),"Mean"]
  }
  B = rbind(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,1,1,0,0,0),
    c(1,1,1,1,0,0),
    c(1,1,1,1,1,0),
    c(1,1,1,1,1,1)
  )
  theta_l_hat_mon2 <- B %*% beta_s_hat
  print(as.numeric(theta_l_hat_mon2))
}



# Plot curves against data
ggplot(data.frame(x=data_jags$l, y=data_jags$y), aes(x=x,y=y)) +
  geom_point(alpha=0.01) +
  coord_cartesian(ylim=c(-1,1)) +
  geom_point(
    aes(x=x,y=y),
    data.frame(
      x = c(0:6),
      y = c(mean(data_jags$y[data_jags$l==0]),
            mean(data_jags$y[data_jags$l==1]),
            mean(data_jags$y[data_jags$l==2]),
            mean(data_jags$y[data_jags$l==3]),
            mean(data_jags$y[data_jags$l==4]),
            mean(data_jags$y[data_jags$l==5]),
            mean(data_jags$y[data_jags$l==6]))
    ),
    color = "purple",
    size = 3
  ) +
  geom_line(
    aes(x=x,y=y, color=type),
    data.frame(
      x = rep(c(0:6), 3),
      y = c(0,as.numeric(theta_l_hat),
            0,as.numeric(theta_l_hat_mon1),
            0,as.numeric(theta_l_hat_mon2)
      ),
      type = rep(c("unconstr", "monotone U(0,1)", "mono (gamma)"), each=7)
    )
  )
# geom_line(
#     aes(x=x,y=y, color=type),
#     data.frame(
#       x = rep(c(0:6), 2),
#       y = c(0,as.numeric(theta_l_hat),0,as.numeric(theta_l_hat_mon)),
#       type = rep(c("unconstrained", "monotone"), each=7)
#     )
#   )

print(theta_l_hat)
print(as.numeric(theta_l_hat_mon1))
print(as.numeric(theta_l_hat_mon2))
print(mean(theta_l_hat))
print(mean(as.numeric(theta_l_hat_mon1)))
print(mean(as.numeric(theta_l_hat_mon2)))



# !!!!!
# Plot curves against data
# mon1 <- theta_l_hat_mon1
ggplot(data.frame(x=data_jags$l, y=data_jags$y), aes(x=x,y=y)) +
  geom_point(alpha=0.01) +
  coord_cartesian(ylim=c(-1,1)) +
  geom_point(
    aes(x=x,y=y),
    data.frame(
      x = c(0:6),
      y = c(mean(data_jags$y[data_jags$l==0]),
            mean(data_jags$y[data_jags$l==1]),
            mean(data_jags$y[data_jags$l==2]),
            mean(data_jags$y[data_jags$l==3]),
            mean(data_jags$y[data_jags$l==4]),
            mean(data_jags$y[data_jags$l==5]),
            mean(data_jags$y[data_jags$l==6]))
    ),
    color = "purple",
    size = 3
  ) +
  geom_line(
    aes(x=x,y=y, color=type),
    data.frame(
      x = rep(c(0:6), 2),
      y = c(0,as.numeric(theta_l_hat),
            0,as.numeric(theta_l_hat_mon4)
      ),
      type = rep(c("unconstrained", "constrained"), each=7)
    )
  )



# # !!!!!
# shape = 0.01
# rate = 0.01
# ggplot(data.frame(x=c(0,100)), aes(x=x)) + stat_function(fun = function(x) {
#   dgamma(x, shape=shape, rate=rate)
# }) + labs(title=paste0("Gamma density (shape: ",shape,", rate: ",rate,")"))

# jags_code <- quote("
#   model {
#     for (n in 1:N) {
#       y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
#       beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
#       + beta_s_1*s_1[n] + beta_s_2*s_2[n] + beta_s_3*s_3[n] +
#       beta_s_4*s_4[n] + beta_s_5*s_5[n] + beta_s_6*s_6[n] + alpha[i[n]],
#       1/(sigma^2))
#     }
#     for (n in 1:I) {
#       alpha[n] ~ dnorm(0, 1/(tau^2))
#     }
#     beta_s_6 <- - exp(alpha_s_6)
#     beta_s_5 <- - exp(alpha_s_5)
#     beta_s_4 <- - exp(alpha_s_4)
#     beta_s_3 <- - exp(alpha_s_3)
#     beta_s_2 <- - exp(alpha_s_2)
#     beta_s_1 <- - exp(alpha_s_1)
#     alpha_s_6 <- log(10) - e_s_6
#     alpha_s_5 <- log(10) - e_s_5
#     alpha_s_4 <- log(10) - e_s_4
#     alpha_s_3 <- log(10) - e_s_3
#     alpha_s_2 <- log(10) - e_s_2
#     alpha_s_1 <- log(10) - e_s_1
#     e_s_6 ~ dexp(1)
#     e_s_5 ~ dexp(1)
#     e_s_4 ~ dexp(1)
#     e_s_3 ~ dexp(1)
#     e_s_2 ~ dexp(1)
#     e_s_1 ~ dexp(1)
#     beta_j_7 ~ dnorm(0, 1.0E-4)
#     beta_j_6 ~ dnorm(0, 1.0E-4)
#     beta_j_5 ~ dnorm(0, 1.0E-4)
#     beta_j_4 ~ dnorm(0, 1.0E-4)
#     beta_j_3 ~ dnorm(0, 1.0E-4)
#     beta_j_2 ~ dnorm(0, 1.0E-4)
#     beta0 ~ dnorm(0, 1.0E-4)
#     tau <- 1/sqrt(tau_prec)
#     tau_prec ~ dgamma(1.0E-3, 1.0E-3)
#     sigma <- 1/sqrt(sigma_prec)
#     sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
#   }
# ")
