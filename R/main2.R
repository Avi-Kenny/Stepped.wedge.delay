# Title: "Stepped wedge lagged effect simulation"
# Author: Avi Kenny
# Date: 2020-05-20



#################.
##### Setup #####
#################.

# devtools::install_github(repo="Avi-Kenny/simba", dependencies=TRUE)

# Load packages
library(z.stepped.wedge)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lme4)
library(geepack)
library(stringr)
library(simba)
library(parallel)
library(glmmTMB, lib.loc=)
library(restriktor)

# Load functions
source(generate_dataset.R)
source(log_lik_spline.R)
source(one_simulation.R)
source(plot_outcome.R)
source(plot_sw_design.R)
source(run_analysis.R)

# Set code blocks to run
run_setup <- TRUE
run_main <- FALSE
run_coverage <- FALSE
run_misc <- FALSE
run_testing_staircase <- FALSE
run_testing_1KSPL <- TRUE
run_testing <- FALSE



###################################.
##### SETUP: Simulation setup #####
###################################.

if ( run_setup ) {

  # Set up and configure simba object
  sim <- new_sim()
  sim %<>% set_config(
    # num_sim = 100,
    num_sim = 1000,
    parallel = "outer",
    # parallel = "none",
    packages = c("dplyr", "magrittr", "stringr", "geepack", "lme4",
                 "z.stepped.wedge", "glmmTMB", "restriktor")
  )
  sim %<>% add_constants(
    alpha = log(0.1)
  )

  # Negative log lik corresponding to two-stage dissertation approach
  sim %<>% add_method(
    "neg_log_lik",
    function(theta, d, J, theta_l_hat, sigma_l_hat) {

      l_times <- 1:(J-1)
      mu_d <- 1-exp(-l_times/d)
      log_lik <- -0.5 * t(theta_l_hat - theta*mu_d) %*%
        solve(sigma_l_hat) %*% (theta_l_hat - theta*mu_d)

      return (-1*log_lik)

    }
  )

  # Negative log lik corresponding to two-stage spline approach
  sim %<>% add_method(
    "neg_log_lik_spl",
    function(theta, p_x, p_y, J, theta_l_hat, sigma_l_hat) {

      # !!!!! g_x hard-coded for now
      g_x <- J

      l_times <- 1:(J-1)
      mu_spl <- sapply(l_times, function(l) {
        I1 <- ifelse(0<l & l<=p_x, 1, 0)
        I2 <- ifelse(p_x<l & l<=g_x, 1, 0)
        I3 <- ifelse(g_x<l, 1, 0)
        (p_y/p_x)*l*I1 + ((1-p_y)*l+g_x+p_y-p_x-1)/(g_x-p_x)*I2 + I3
      })

      log_lik <- -0.5 * t(theta_l_hat - theta*mu_spl) %*%
        solve(sigma_l_hat) %*% (theta_l_hat - theta*mu_spl)

      return (-1*log_lik)

    }
  )

  # Add functions to simba object
  sim %<>% add_creator(generate_dataset)
  sim %<>% add_script(one_simulation)

}



#####################################.
##### MAIN: Reproduce table 3.1 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    d = c(0, 0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM", "IG LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### MAIN: Reproduce table 3.2 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    d = c(0, 0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM", "IG LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### MAIN: Reproduce table 3.3 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    d = c(0, 0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM", "IG LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### MAIN: Reproduce table 3.4 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    d = c(0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_d", truth="d", estimate="d_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### MAIN: Reproduce table 3.5 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    d = c(0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_d", truth="d", estimate="d_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### MAIN: Reproduce table 3.6 #####
#####################################.

if ( run_main ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    d = c(0.5, 1.4),
    tau = 0,
    sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = c("2S LM"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_d", truth="d", estimate="d_hat")
  ))

  print(round(difftime(Sys.time(), start_time),2))

}



#########################################################.
##### COVERAGE: Investigate coverage issue (no GEE) #####
#########################################################.

if ( FALSE ) {
# if ( run_coverage ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = seq(12, 120, 24),
    # n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = seq(20, 100, 20),
    # n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    d = 1.4,
    # d = c(0,1.4),
    tau = 0,
    # tau = c(0,0.25),
    sigma = 3,
    # sigma = 0.3,
    data_type = "binomial",
    # data_type = c("normal", "binomial"),
    analysis_type = "2S LM",
    # analysis_type = c("2S LM", "2S LMM REML"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-06 02_46_27")
  print(round(difftime(Sys.time(), start_time),2))

}



######################################################.
##### COVERAGE: Investigate coverage issue (GEE) #####
######################################################.

if ( FALSE ) {
# if ( run_coverage ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    theta = log(0.5),
    d = 1.4,
    tau = c(0,0.25),
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis_type = c("2S GEE EX", "2S GEE ID"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-06 03_41_40")
  print(round(difftime(Sys.time(), start_time),2))

}



##############################################################.
##### COVERAGE: Investigate coverage issue (REML vs. ML) #####
##############################################################.

if ( FALSE ) {
# if ( run_coverage ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    d = 1.4,
    tau = 0.25,
    # tau = c(0,0.25),
    sigma = 0.3,
    data_type = "normal",
    analysis_type = "2S LMM REML",
    # analysis_type = c("2S LMM REML", "2S LMM ML"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-06 10_15_49")
  print(round(difftime(Sys.time(), start_time),2))

}



###############################################################.
##### COVERAGE: Investigate coverage issue (H-likelihood) #####
###############################################################.

if ( run_coverage ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 12,
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    d = 1.4,
    tau = 0,
    sigma = 3,
    data_type = "binomial",
    analysis_type = c("2S LM", "2S HL"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-18 10_21_46")
  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### TESTING: Staircase method #####
#####################################.

if ( run_testing_staircase ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    d = 1.4,
    tau = 0,
    sigma = 3,
    data_type = c("normal", "binomial"),
    analysis_type = c("2S LM", "Staircase"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat"),
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-18 10_21_46")
  print(round(difftime(Sys.time(), start_time),2))

}



####################################.
##### TESTING: Spline (1 knot) #####
####################################.

if ( run_testing_1KSPL ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    d = c(0,1),
    tau = 0,
    sigma = 3,
    data_type = c("normal", "binomial"),
    analysis_type = "SPL 1K",
    delay_model = "spline"
  )

  sim %<>% run("one_simulation")

  print(summary(
    sim_obj = sim,
    means = list(all=TRUE, na.rm=TRUE),
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat"),
    coverage = list(
      name = "cov_theta",
      truth = "theta",
      estimate = "theta_hat",
      se = "se_theta_hat",
      na.rm = TRUE
    )
  ))

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-19 17_20_20")
  print(round(difftime(Sys.time(), start_time),2))

}



#####################################.
##### TESTING: Two-stage spline #####
#####################################.

if ( run_testing ) {

  start_time <- Sys.time()

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    d = c(0, 0.5, 1.4),
    tau = 0,
    analysis_type = c("2S SPL"),
    delay_model = "s-curve"
  )

  sim %<>% run("one_simulation")

  summary(
    sim_obj = sim,
    bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
  )

  print("")
  print(round(difftime(Sys.time(), start_time),2))

}



#######################################.
##### MISC: Check MLE calculation #####
#######################################.

if ( run_misc ) {

  J <- 5
  theta_hat_l <- matrix(1:J,ncol=1)
  mu_d <- matrix(0.3*c(3:(J+2)),ncol=1)
  A <- matrix(runif(J^2)*2-1, ncol=J)
  sigma_inv <- t(A) %*% A

  neg_log_lik <- function(theta) {
    return(
      (1/2) * t(theta_hat_l-(theta*mu_d)) %*%
        sigma_inv %*%
        (theta_hat_l-(theta*mu_d))
    )
  }

  optim(
    par = 6,
    fn = function(x) {
      return ( neg_log_lik(x) )
    },
    method = "BFGS"
  )

  mle_ank <- (
    (t(theta_hat_l) %*% sigma_inv %*% mu_d) +
      t((t(theta_hat_l) %*% sigma_inv %*% mu_d))) /
      (2* t(mu_d) %*% sigma_inv %*% mu_d)
  mle_tsg <- (t(theta_hat_l) %*% sigma_inv %*% mu_d) /
             (t(mu_d) %*% sigma_inv %*% mu_d)

  print(mle_ank)
  print(mle_tsg)

}



###########################################.
##### MISC: Linear spline alternative #####
###########################################.

if ( run_misc ) {

  ggplot(data.frame(x=c(0,6)), aes(x=x)) +
    stat_function(fun = function(x) {
      return ( 0.7*x + (0.3/5 - 0.7)*pmax(0,x-1) )
    }) +
    geom_point(aes(x=0, y=0)) +
    geom_point(aes(x=6, y=1), colour="green") +
    geom_point(aes(x=1, y=0.7), colour="purple") +
    labs(
      title = "Linear spline model for R_il",
      y = "R_il (% of treatment effect achieved)",
      x = "Time since implementation (l_i)"
    )

}



################################################.
##### ARCHIVE: Old code (to recycle later) #####
################################################.

if ( FALSE ) {

  print(paste("Number of available cores:", parallel::detectCores()))
  print(paste("SLURM_ARRAY_JOB_ID:", Sys.getenv("SLURM_ARRAY_JOB_ID")))
  print(paste("SLURM_CPUS_ON_NODE:", Sys.getenv("SLURM_CPUS_ON_NODE")))
  print(paste("SLURM_NODELIST:", Sys.getenv("SLURM_NODELIST")))
  print(paste("SLURM_NNODES:", Sys.getenv("SLURM_NNODES")))
  print(paste("SLURM_NTASKS:", Sys.getenv("SLURM_NTASKS")))

  summary(
    sim_obj = sim,
    sd = list(
      list(name="sd_theta_hat", x="theta_hat"),
      list(name="sd_d_hat", x="d_hat")
    ),
    bias = list(
      list(name="bias_theta", truth="theta", estimate="theta_hat"),
      list(name="bias_d", truth="d", estimate="d_hat")
    ),
    coverage = list(
      list(name="cov_theta", truth="theta",
           estimate="theta_hat", se="se_theta_hat"),
      list(name="cov_d", truth="d", estimate="d_hat", se="se_d_hat")
    )
  )

  # Plots
  plot_sw_design(data_1)
  plot_outcome(data_1, type="no error")
  plot_outcome(data_1, type="realized")

  # Binomial GLM
  model_binomial_gee1 <- geeglm(
    y ~ factor(j) + factor(x_ij),
    data = data$data,
    id = i,
    family = binomial(link = "log"),
    corstr = "exchangeable"
  )
  summary(model_binomial_gee1)
  system.time(
    model_binomial_gee2 <- geeglm(
      y ~ factor(j) + factor(l),
      data = data$data,
      id = i,
      family = binomial(link = "log"),
      corstr = "exchangeable"
    )
  )
  summary(model_gee2)

  # Only estimate theta (step function)
  model_normal_gee1 <- geeglm(
    y ~ factor(j) + factor(x_ij),
    data = data$data,
    id = i,
    family = "gaussian",
    corstr = "exchangeable"
  )
  summary(model_normal_gee1)

}
