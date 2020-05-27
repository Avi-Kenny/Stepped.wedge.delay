# Title: "Stepped wedge lagged effect simulation"
# Author: Avi Kenny
# Date: 2020-05-20



#################.
##### Setup #####
#################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jim Hughes/Project - Stepped wedge/z.stepped.wedge/R")
} else {
  setwd("z.stepped.wedge/R")
}

# Load packages
{
  library(z.stepped.wedge)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(lme4)
  library(geepack)
  library(stringr)
  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
  library(parallel)
  library(glmmTMB, lib.loc=)
  library(restriktor)
}

# Load functions
{
  source("generate_dataset.R")
  source("log_lik_spline.R")
  source("one_simulation.R")
  source("plot_outcome.R")
  source("plot_sw_design.R")
  source("run_analysis.R")
  source("sw_spline.R")
}

# Set code blocks to run
{
  run_setup <- TRUE
  run_results <- FALSE
  run_main_526 <- TRUE
  run_tables3 <- FALSE
  run_cov_nogee <- FALSE
  run_cov_gee <- FALSE
  run_cov_reml <- FALSE
  run_cov_hlik <- FALSE
  run_misc <- FALSE
  run_testing_staircase <- FALSE
  run_testing_1Kspl <- FALSE
  run_testing_2Sspl <- FALSE
}



###################################.
##### SETUP: Simulation setup #####
###################################.

if ( run_setup ) {

  # Set up and configure simba object
  sim <- new_sim()
  sim %<>% set_config(
    num_sim = 100,
    # num_sim = 1000,
    parallel = "none",
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

  # Create directory to store results
  if (!file.exists("../simba.out")) { dir.create("../simba.out") }

}



###########################################################.
##### Compile simulation results into a single object #####
###########################################################.

if (FALSE) {

  # Merge *.simba files
  sims <- list.files(
    path = "../simba.out/simba.out",
    pattern = "*.simba",
    full.names = TRUE,
    recursive = FALSE
  )
  sim <- NULL
  for (s in sims) {
    s <- readRDS(s)
    if (is.null(sim)) { sim <- s } else { sim <- merge(sim, s) }
  }
  saveRDS(sim, file="../simba.out/sim_tab3.1.simba")

}



#########################################.
##### MAIN: 5/26 set of simulations #####
#########################################.

if ( run_main_526 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 48,
    # n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    # n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    tau = 0,
    # tau = c(0,0.25), # !!!!! Is 0.25 the best level given the new sigma? Make this 50% of the variance
    sigma = 0.3,
    data_type = "normal",
    analysis = list(
      "SPL" = list(type="SPL", params=list(knots=c(1,6)))
    ),
    # analysis = list(
    #   "SPL 1" = list(...), # !!!!!
    #   "SPL 2" = list(...), # !!!!!
    #   "2S LMM" = list(...) # !!!!!
    # ),
    delay_model = list(
      "Spline model" = list(
        type = "spline",
        params = list(knots=c(1,6),slopes=c(0.8,0.04))
      )
    )
    # delay_model = list(
    #   "Exp model" = list(type="exp", params=list(d=0)),
    #   "Spline model" = list(
    #     type = "spline",
    #     params = list(knots=c(1,6),slopes=c(0.8,0.04))
    #   )
    # )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_gee.simba")
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
  }

}



#####################################.
##### MAIN: Reproduce table 3.1 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.1.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#####################################.
##### MAIN: Reproduce table 3.2 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.2.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#####################################.
##### MAIN: Reproduce table 3.3 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.3.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#####################################.
##### MAIN: Reproduce table 3.4 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.4.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#####################################.
##### MAIN: Reproduce table 3.5 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.5.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#####################################.
##### MAIN: Reproduce table 3.6 #####
#####################################.

if ( run_tables3 ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.6.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#########################################################.
##### COVERAGE: Investigate coverage issue (no GEE) #####
#########################################################.

if ( run_cov_nogee ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = seq(20, 100, 20),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 3,
    data_type = c("normal", "binomial"),
    analysis = list(
      "2S LM" = list(type="2S LM"),
      "2S LMM REML" = list(type="2S LMM", params=list(REML=TRUE))
    ),
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_nogee.simba")
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
  }

}



######################################################.
##### COVERAGE: Investigate coverage issue (GEE) #####
######################################################.

if ( run_cov_gee ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list(
      "2S GEE EXC" = list(type="2S GEE", params=list(corr="exchangeable")),
      "2S GEE IND" = list(type="2S GEE", params=list(corr="independence"))
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_gee.simba")
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
  }

}



##############################################################.
##### COVERAGE: Investigate coverage issue (REML vs. ML) #####
##############################################################.

if ( run_cov_reml ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 0.3,
    data_type = "normal",
    analysis = list(
      "2S LMM REML" = list(type="2S LMM", params=list(REML=TRUE)),
      "2S LMM ML" = list(type="2S LMM", params=list(REML=FALSE))
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_reml.simba")
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
  }

}



###############################################################.
##### COVERAGE: Investigate coverage issue (H-likelihood) #####
###############################################################.

if ( run_cov_hlik ) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 12,
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 3,
    data_type = "binomial",
    analysis = list(
      "2S LM" = list(type="2S LM"),
      "2S HL" = list(type="2S HL")
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_hlik.simba")
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
  }

}



#####################################.
##### TESTING: Staircase method #####
#####################################.

if ( run_testing_staircase ) {

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

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_staircase.simba")
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
  }

}



####################################.
##### TESTING: Spline (1 knot) #####
####################################.

if ( run_testing_1Kspl ) {

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
    delay_model = "spline K1"
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_1Kspl.simba")
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
  }

}



#####################################.
##### TESTING: Two-stage spline #####
#####################################.

if ( run_testing_2Sspl ) {

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

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_2Sspl.simba")
    summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    )
  }

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



################################.
##### MISC: Graphs for PPT #####
################################.

if ( run_misc ) {

  # Plot 1
  # Exported as 600w X 400h
  ggplot(
    data.frame(
      x = c(1:6),
      y = c(0,0,3,3,3,3),
      state = c( rep("Control",2), rep("Treatment",4) )
    ),
    aes(x=x, y=y, color=state)
  ) +
    stat_function(
      fun = function(x) { return(
        ifelse(x<=2, 0, ifelse(x>=3, 3, 3*(x-2)^2))
      )},
      color="#333333"
    ) +
    geom_point(size=3) +
    labs(
      title = "Immediate effect",
      y = "Treatment effect",
      x = "Time (steps)",
      color = "State"
    )

  # Plot 2
  # Exported as 600w X 400h
  ggplot(
    data.frame(
      x = c(1:6),
      y = c(0,0,0.693,2.135,3,3),
      state = c( rep("Control",2), rep("Treatment",4) )
    ),
    aes(x=x, y=y, color=state)
  ) +
    stat_function(
      fun = function(x) { return(
        ifelse(x<=2, 0, ifelse(x>=5, 3, 1.507543*(1+sin((x-2-(pi/2))))))
      )},
      color="#333333"
    ) +
    geom_point(size=3) +
    labs(
      title = "Lagged effect",
      y = "Treatment effect",
      x = "Time (steps)",
      color = "State"
    )

  # Plot 3
  # Exported as 600w X 400h
  set.seed(17)
  x = c(rep(3,32),rep(4,24),rep(5,16),rep(6,8))
  y = c(
    rnorm(n=32, mean=1.8, sd=0.3),
    rnorm(n=24, mean=2.7, sd=0.3),
    rnorm(n=16, mean=2.7, sd=0.3),
    rnorm(n=8, mean=3.1, sd=0.3)
  )
  ggplot(
    data.frame(x=x, y=y),
    aes(x=x, y=y)
  ) +
    geom_point(size=2, alpha=0.2, shape=16) +
    xlim(1,6) +
    ylim(0,4) +
    geom_point(
      data = data.frame(
        x = c(3,4,5,6),
        y = c(1.8,2.7,2.7,3.1)
      ),
      color = "turquoise",
      size = 3
    ) +
    labs(
      title = "Lagged effect stepped wedge model",
      y = "Treatment effect",
      x = "Time (steps)",
      color = "State"
    )

  # Plot 4
  # Exported as 600w X 400h
  ggplot(
    data.frame(
      x = c(3,4,5,6),
      y = c(1.8,2.7,2.7,3.1)
    ),
    aes(x=x, y=y)
  ) +
    stat_function(
      fun = function(x) { return(
        ifelse(x<=2, 0, 3*(1-exp(-(x-2)/1)))
      )},
      color="#333333"
    ) +
    xlim(1,6) +
    ylim(0,4) +
    geom_point(color="turquoise", size=3) +
    labs(
      title = "Lagged effect stepped wedge model",
      y = "Treatment effect",
      x = "Time (steps)",
      color = "State"
    )

  # Plot 5 new
  # Exported as 600w X 400h
  ggplot(
    data.frame(x=x, y=y),
    aes(x=x, y=y)
  ) +
    geom_point(size=2, alpha=0.2, shape=16) +
    stat_function(
      fun = function(x) { return(
        ifelse(x<=2, 0, ifelse(x>=5, 2.9,
                               2.1*(x-2) + (0.8/2 - 2.1)*pmax(0,(x-2)-1)
        ))
      )},
      color="#333333"
    ) +
    xlim(1,6) +
    ylim(0,4) +
    geom_point(aes(x=2, y=0), colour="black", size=3) +
    geom_point(aes(x=5, y=2.9), colour="green", size=3) +
    geom_point(aes(x=3, y=2.1), colour="purple", size=3) +
    labs(
      title = "Lagged effect stepped wedge model",
      y = "Treatment effect",
      x = "Time (steps)",
      color = "State"
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