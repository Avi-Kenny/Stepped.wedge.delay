# Title: "Stepped wedge lagged effect simulation"
# Author: Avi Kenny
# Date: 2020-10-25



#################.
##### Setup #####
#################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jim Hughes/Project - Stepped wedge lag/z.stepped.wedge/R")
} else {
  setwd("z.stepped.wedge/R")
}

# Load functions
{
  source("generate_dataset.R")
  source("one_simulation.R")
  source("plot_outcome.R")
  source("plot_sw_design.R")
  source("run_analysis.R")
  source("effect_curve.R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_process_results <- FALSE
  run_viz <- FALSE
  run_misc <- FALSE
  run_testing <- FALSE
}



########################################.
##### SETUP: Load packages locally #####
########################################.

if (FALSE) {

  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
  library(dplyr)
  library(stringr)
  library(lme4)
  library(rjags)
  # library(rstan)
  # library(INLA)
  library(sqldf)
  library(glmmTMB)
  library(mgcv)
  library(fastDummies)
  library(scales)
  library(car)
  library(ggplot2)
  library(parallel)
  library(viridis)

}



######################################################################.
##### TESTING: Generate dataset for testing new analysis methods #####
######################################################################.

if (FALSE) {

  # Generate dataset
  data <- generate_dataset(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 5, # 50
    theta = -0.5,
    tau = 1,
    alpha = -2,
    data_type = "normal",
    sigma = 0.5, # 2.1
    delay_model = list(type="spline", params=list(knots=c(0,1),slopes=1)),
    n_extra_time_points = 0
  )

  # Set number of time points
  J <- data$params$n_time_points

}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (run_main) {
  if (Sys.getenv("run") %in% c("first", "")) {

    # Simulation 1: compare all methods
    level_set_1 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = -0.5,
      tau = 1,
      sigma = 2.1,
      data_type = "normal",
      method = list(
        # "HH" = list(method="HH"),
        "ETI" = list(method="ETI"),
        # "MCMC-STEP (exp; N(1,10) prior)" = list(method="MCMC-STEP-MON",
        #                                     enforce="exp; N(1,10) prior"),
        "MCMC-STEP (exp; mix prior 0.1)" = list(method="MCMC-STEP-MON",
                                                enforce="exp; mix prior 0.1"),
        "MCMC-STEP (exp; mix prior 0.2)" = list(method="MCMC-STEP-MON",
                                                enforce="exp; mix prior 0.2"),
        "MCMC-STEP (exp; mix prior 0.3)" = list(method="MCMC-STEP-MON",
                                            enforce="exp; mix prior 0.3")
      ),
      delay_model = list(
        "Instantaneous" = list(
          type = "spline",
          params = list(knots=c(0,1), slopes=1)
        ),
        "Lagged" = list(
          type = "spline",
          params = list(knots=c(0,2,3), slopes=c(0,1))
        ),
        "Curved" = list(
          type = "exp",
          params = list(d=1.5)
        ),
        "Partially convex" = list(
          type = "spline",
          params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
        )
      ),
      n_extra_time_points = 0
    )

    # Simulation 2: compare all methods (power)
    level_set_2 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = seq(-0.5,0,0.1),
      tau = 1,
      sigma = 2.1,
      data_type = "normal",
      method = list(
        "HH" = list(method="HH"),
        "ETI" = list(method="ETI"),
        "MCMC-STEP (exp; mix prior 0.2)" = list(method="MCMC-STEP-MON",
                                                enforce="exp; mix prior 0.2")
      ),
      delay_model = list(
        "Instantaneous" = list(
          type = "spline",
          params = list(knots=c(0,1), slopes=1)
        ),
        "Lagged" = list(
          type = "spline",
          params = list(knots=c(0,2,3), slopes=c(0,1))
        ),
        "Curved" = list(
          type = "exp",
          params = list(d=1.5)
        ),
        "Partially convex" = list(
          type = "spline",
          params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
        )
      ),
      n_extra_time_points = 0
    )

    # Simulation 3: n_extra_time_points
    level_set_3 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = -0.5,
      tau = 1,
      sigma = 2.1,
      data_type = "normal",
      method = list("ETI" = list(method="ETI")),
      delay_model = list(
        "Instantaneous" = list(
          type = "spline",
          params = list(knots=c(0,1), slopes=1)
        ),
        "Lagged" = list(
          type = "spline",
          params = list(knots=c(0,2,3), slopes=c(0,1))
        ),
        "Curved" = list(
          type = "exp",
          params = list(d=1.5)
        ),
        "Partially convex" = list(
          type = "spline",
          params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
        )
      ),
      n_extra_time_points = c(0,1,2)
    )

    # Simulation 4: effect_reached
    level_set_4 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = -0.5,
      tau = 1,
      sigma = 2.1,
      data_type = "normal",
      method = list(
        "HH" = list(method="HH"),
        "ETI (effect_reached=0)" = list(method="ETI", effect_reached=0),
        "ETI (effect_reached=3)" = list(method="ETI", effect_reached=3),
        "ETI (effect_reached=4)" = list(method="ETI", effect_reached=4)
      ),
      delay_model = list(
        "Instantaneous" = list(
          type = "spline",
          params = list(knots=c(0,1), slopes=1)
        ),
        "Lagged" = list(
          type = "spline",
          params = list(knots=c(0,2,3), slopes=c(0,1))
        ),
        "Partially convex" = list(
          type = "spline",
          params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
        )
      ),
      n_extra_time_points = 0
    )

    # # Simulation 5: !!!!!
    # level_set_5 <- list(...)

  }
}



################################################.
##### MAIN: Choose which simulation to run #####
################################################.

if (run_main) {
  if (Sys.getenv("run") %in% c("first", "")) {

    # Set this manually
    level_set <- level_set_1

  }
}



#################################.
##### MAIN: Main simulation #####
#################################.

# Commands for job sumbission on Slurm:
# sbatch --export=run='first',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-16 --export=run='main',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

# Commands for job sumbission on SGE:
# qsub -v run='first',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992344 -t 1-3 -v run='main',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992345 -v run='last',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh

if (run_main) {

  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")

  run_on_cluster(

    first = {

      # Set up and configure simba object
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = 1000, # !!!!!
        parallel = "none",
        stop_at_error = FALSE, # !!!!!
        packages = c("dplyr", "magrittr", "stringr", "lme4", "rjags", # "geepack", "restriktor", "scam", "gamlss"
                     "sqldf", "glmmTMB", "mgcv", "fastDummies", "scales", "car")
        # packages = c("dplyr", "magrittr", "stringr", "lme4", "rjags", "rstan", # "geepack", "restriktor", "scam", "gamlss"
        #              "sqldf", "glmmTMB", "mgcv", "fastDummies", "scales", "car")
      )
      sim %<>% add_constants(alpha = -2)

      # Add functions to simba object
      sim %<>% add_creator(generate_dataset)
      sim %<>% set_script(one_simulation)
      sim %<>% add_method(run_analysis)
      sim %<>% add_method(effect_curve)

      # Set levels
      sim <- do.call(set_levels, c(list(sim), level_set))

    },

    main = {
      sim %<>% run()
    },

    last = {
      sim %>% summary() %>% print()
      sim$errors %>% print()
    },

    cluster_config = list(
      # js = "sge",
      # dir = "/home/students/avikenny/Desktop/z.stepped.wedge"
      js = "slurm",
      dir = "/home/akenny/z.stepped.wedge"
    )

  )

}



#####################################.
##### MAIN: Process sim results #####
#####################################.

if (run_process_results) {

  # Read in simulation object
  sim <- readRDS("../simba.out/sim_20210131_2.simba")

  # Generate true ATE values
  sim$results %<>% mutate(
    ate = case_when(
      delay_model=="Instantaneous" ~ theta,
      delay_model=="Lagged" ~ theta * mean(effect_curve(
        x=c(1:6), type="spline", params=list(knots=c(0,2,3),slopes=c(0,1))
      )),
      delay_model=="Curved" ~ theta * mean(effect_curve(
        x=c(1:6), type="exp", params=list(d=1.5)
      )),
      delay_model=="Partially convex" ~ theta * mean(effect_curve(
        x=c(1:6), type="spline", params=list(knots=c(0,2,4),slopes=c(0.1,0.4))
      ))
    )
  )

  # Summarize data
  summ <- summary(
    sim_obj = sim,
    mean = list(
      list(name="ate", x="ate"),
      list(name="mean_ate", x="ate_hat"),
      list(name="mean_lte", x="lte_hat")
    ),
    bias = list(
      list(name="bias_ate", estimate="ate_hat", truth="ate"),
      list(name="bias_lte", estimate="lte_hat", truth="theta")
    ),
    mse = list(
      list(name="mse_ate", estimate="ate_hat", truth="ate"),
      list(name="mse_lte", estimate="lte_hat", truth="theta")
    ),
    # quantile = list(
    #   list(name="q025_ate", x="ate_hat", prob=0.025, na.rm=TRUE),
    #   list(name="q975_ate", x="ate_hat", prob=0.975, na.rm=TRUE)
    # ),
    coverage = list(
      list(name="cov_ate", truth="ate", estimate="ate_hat", se="se_ate_hat"),
      list(name="cov_lte", truth="theta", estimate="lte_hat", se="se_lte_hat"),
      list(name="beta_ate", truth=0, estimate="ate_hat", se="se_ate_hat"),
      list(name="beta_lte", truth=0, estimate="lte_hat", se="se_lte_hat")
    )
  )

  # Drop some columns
  summ %<>% subset(select=-c(1:4,6:8))

  # Transform summary data
  summ %<>% mutate(
    # method = factor(method, levels=c("ETI","MCMC-STEP (exp; mix prior 0.1)","MCMC-STEP (exp; mix prior 0.2)","MCMC-STEP (exp; mix prior 0.3)")), # !!!!!
    method = factor(method, levels=c("HH","ETI (effect_reached=3)","ETI (effect_reached=4)","ETI (effect_reached=0)")), # !!!!!
    # method = factor(method, levels=c("HH","ETI","MCMC-STEP (exp; N(1,10) prior)","MCMC-STEP (exp; mix prior 0.2)")), # !!!!!
    delay_model = factor(delay_model, levels=c("Instantaneous","Lagged",
                                               "Curved","Partially convex")),
    power_ate = 1 - beta_ate,
    power_lte = 1 - beta_lte
  )
  p_data <- sqldf("
    SELECT method, delay_model, 'ATE' AS which, bias_ate AS bias,
    cov_ate AS coverage, power_ate AS power, mse_ate AS mse FROM summ
    UNION SELECT method, delay_model, 'LTE', bias_lte,
    cov_lte, power_lte, mse_lte FROM summ
  ")
  p_data <- sqldf("
    SELECT method, delay_model, which, 'bias' AS stat, bias AS value FROM p_data
    UNION SELECT method, delay_model, which, 'coverage', coverage FROM p_data
    UNION SELECT method, delay_model, which, 'power', power FROM p_data
    UNION SELECT method, delay_model, which, 'mse', mse FROM p_data
  ")

  # Plot results
  ggplot(
    filter(p_data, stat!="power"),
    aes(x=which, y=value, fill=method)
  ) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8, color="white") +
    facet_grid(cols=vars(delay_model), rows=vars(stat), scales="free") +
    theme(legend.position="bottom") +
    scale_fill_manual(values=viridis(5)) +
    labs(y=NULL, x=NULL)

}



############################################.
##### VIZ: Point estimates + quantiles #####
############################################.

if (run_viz) {

  # sim_main_1026.simba
  # Export: 800 x 400
  ggplot(
    data = summ,
    aes(
      x = method,
      y = mean_ate_hat,
      color = method,
      ymin = q025_ate,
      ymax = q975_ate
    )
  ) +
    geom_hline(
      mapping = aes(yintercept=mean_ate),
      linetype = 2
    ) +
    geom_point(size=2)+
    geom_errorbar(
      aes(
        ymin = q025_ate,
        ymax = q975_ate,
        color = method
      ),
      width = 0.2,
      cex = 1
    ) +
    labs(
      title = paste("Average point estimates and 2.5%/97.5% quantiles (1,000",
                    "sims per level): theta=log(0.5), tau=1"),
      x = "Analysis method",
      y = NULL
    ) +
    facet_wrap(~delay_model, ncol=3) +
    # facet_grid(rows=vars(date), cols=vars(dgm)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2", "green3", "darkorchid2"))

  # sim_main_1031.simba
  # Export: 800 x 400
  ggplot(
    data = summ,
    aes(
      x = method,
      y = mean_ate_hat,
      color = method,
      ymin = q025_ate,
      ymax = q975_ate
    )
  ) +
    geom_hline(
      mapping = aes(yintercept=mean_ate),
      linetype = 2
    ) +
    geom_point(size=2)+
    geom_errorbar(
      aes(
        ymin = q025_ate,
        ymax = q975_ate,
        color = method
      ),
      width = 0.2,
      cex = 1
    ) +
    labs(
      title = paste("Average point estimates and 2.5%/97.5% quantiles (1,000",
                    "sims per level)"),
      x = "Analysis method",
      y = NULL
    ) +
    coord_cartesian(ylim=c(-1,0)) +
    facet_grid(rows=vars(data_type), cols=vars(delay_model)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2"))

}



#########################.
##### VIZ: Coverage #####
#########################.

if (run_viz) {

  # sim_main_1026.simba
  # Export: 800 x 400
  ggplot(
    data = summ,
    aes(
      x = method,
      y = cov_ate,
      color = method
    )
  ) +
    geom_hline(
      mapping = aes(yintercept=0.95),
      linetype = 2
    ) +
    geom_point(size=2)+
    labs(
      title = "Coverage (1,000 sims per level)",
      x = "Analysis method",
      y = NULL
    ) +
    facet_wrap(~delay_model, ncol=3) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2", "green3", "darkorchid2"))

  # sim_main_1031.simba
  # Export: 800 x 400
  ggplot(
    data = summ,
    aes(
      x = method,
      y = cov_ate,
      color = method
    )
  ) +
    geom_hline(
      mapping = aes(yintercept=0.95),
      linetype = 2
    ) +
    geom_point(size=2)+
    labs(
      title = "Coverage (1,000 sims per level)",
      x = "Analysis method",
      y = NULL
    ) +
    facet_grid(rows=vars(data_type), cols=vars(delay_model)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2"))

}



######################.
##### VIZ: Power #####
######################.

if (run_viz) {

  # sim_main_1101.simba
  # Export: 800 x 400
  ggplot(
    data = summ,
    aes(
      x = theta_log,
      y = power,
      color = method
    )
  ) +
    geom_line(aes(group=method)) +
    geom_point(size=0.9) +
    labs(
      title = "Power of CI-based hypothesis test (500 sims per level)",
      x = "Theta",
      color = "Analysis method",
      y = NULL
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2")) +
    facet_grid(rows=vars(data_type), cols=vars(delay_model)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4)
    )

}



############################################.
##### MAIN: New visualizations (10/26) #####
############################################.

if (run_viz) {

  plot_data <- data.frame(
    "dgm" = character(),
    "date" = character(),
    "method" = character(),
    "estimate" = double(),
    "se" = double(),
    "power" = double(),
    stringsAsFactors = FALSE
  )

  sim <- readRDS("../simba.out/sim_main_1026.simba")
  summ <- summary(
    sim_obj = sim,
    mean = list(all=TRUE, na.rm=TRUE),
    coverage = list(
      name = "beta",
      truth = 0,
      estimate = "ate_hat",
      se = "se_ate_hat"
    )
  )

  for (i in 1:nrow(summ)) {
    plot_data[nrow(plot_data)+1,] <- list(
      "dgm" = summ[i,"delay_model"],
      "date" = "10/26",
      "method" = summ[i,"method"],
      "estimate" = summ[i,"mean_ate_hat"],
      "se" = summ[i,"mean_se_ate_hat"],
      "power" = 1 - summ[i,"beta"]
    )
  }

  plot_data %<>% mutate(
    method = factor(method, levels=c("HH","SS","MSS"))
  )

  # Plot: Estimates and CIs
  # Export: 800 x 400
  ggplot(
    data = plot_data,
    aes(
      x = method,
      y = estimate,
      color = method,
      ymin = estimate-(1.96*se),
      ymax = estimate-(1.96*se)
    )
  ) +
    geom_hline(
      mapping = aes(yintercept=y),
      data = data.frame(
        dgm = c("EXP (d=0)", "EXP (d=1.4)", "SPL (k=2,4 s=0.1,0.4)"),
        y = c(log(0.5)*1, log(0.5)*0.77, log(0.5)*0.57)
      ),
      linetype = 2
    ) +
    geom_point(aes(color=method), size=2)+
    geom_errorbar(
      aes(
        ymin = estimate-(1.96*se),
        ymax = estimate+(1.96*se),
        color = method
      ),
      width = 0.2,
      cex = 1
    ) +
    labs(
      title = paste0("Average point estimates and CIs (1,000 sims per level): ",
                     "theta=log(0.5), tau=1"),
      x = "Analysis method",
      y = NULL
    ) +
    coord_cartesian(ylim=c(-1,0)) +
    # facet_wrap(~dgm, ncol=3) +
    facet_grid(rows=vars(date), cols=vars(dgm)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    )
    # scale_color_manual(values = c("#F8766D", "#7CAE00", "#C77CFF", "green3", "dodgerblue2"))
    # scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2", "green3", "darkorchid2"))

  # !!!!! Combined plot
  # use plot_data2

}



####################################################.
##### MAIN: Visualize sim results (LTE, power) #####
####################################################.

if (run_viz) {

  # Set this manually
  # which <- "Power"
  which <- "Estimates and CIs"

  if (which=="Estimates and CIs") {
    sim_files <- c(
      "../simba.out/sim_main_0526.simba",
      "../simba.out/sim_main_0526_new.simba",
      "../simba.out/sim_main_0602.simba",
      "../simba.out/sim_main_0628_part1.simba",
      "../simba.out/sim_main_0630.simba",
      "../simba.out/sim_main_0702.simba",
      "../simba.out/sim_main_0708.simba"
    )
  }

  if (which=="Power") {
    # sim_files <- "../simba.out/sim_power_0715.simba"
    sim_files <- "../simba.out/sim_power_0728.simba"
  }

  plot_data <- data.frame(
    "dgm" = character(),
    "tau" = integer(),
    "theta" = double(),
    "method" = character(),
    "estimate" = double(),
    "se" = double(),
    "power" = double(),
    stringsAsFactors = FALSE
  )

  for (sims in sim_files) {

    sim <- readRDS(sims)
    summ <- summary(
      sim_obj = sim,
      mean = list(all=TRUE, na.rm=TRUE),
      coverage = list(
        name = "beta",
        truth = 0,
        estimate = "theta_hat",
        se = "se_theta_hat"
      )
    )

    for (i in 1:nrow(summ)) {
      plot_data[nrow(plot_data)+1,] <- list(
        "dgm" = summ[i,"delay_model"],
        "tau" = summ[i,"tau"],
        "theta" = summ[i,"theta"],
        "method" = summ[i,"method"],
        "estimate" = summ[i,"mean_theta_hat"],
        "se" = summ[i,"mean_se_theta_hat"],
        "power" = 1 - summ[i,"beta"]
      )
    }

  }

  # Transform/clean data
  plot_data %<>% filter(
    !( dgm %in% c("SPL (k=2,4 s=0.4,0.1)") ) &
      !(dgm=="SPL (k=2,4 s=0.1,0.4)" & method=="2S LMM") &
      !(method=="Last") &
      !(dgm=="EXP (d=0.5)") &
      !(is.na(method))
  )
  plot_data$method %<>% car::Recode(paste0(
    "'SPL (1-6) MONO'='SPL Mono';'WASH 1'='Wash 1';",
    "'WASH 2'='Wash 2'"
  ))

  plot_data$theta_log <- rep(NA, nrow(plot_data))
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==-0.7,
                            "log(0.5)", plot_data$theta_log)
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==-0.5,
                            "log(0.6)", plot_data$theta_log)
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==-0.4,
                            "log(0.7)", plot_data$theta_log)
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==-0.2,
                            "log(0.8)", plot_data$theta_log)
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==-0.1,
                            "log(0.9)", plot_data$theta_log)
  plot_data$theta_log <- ifelse(round(as.numeric(plot_data$theta),1)==0,
                            "log(1.0)", plot_data$theta_log)
  plot_data$theta_log <- as.factor(plot_data$theta_log)
  plot_data %<>% mutate(
    tau = ifelse(tau==1,"Tau=1","Tau=0")
  )
  if (which=="Estimates and CIs") {
    plot_data %<>% mutate(
      method = factor(method, levels=c(
        "HH", "Wash 1", "Wash 2", "SPL (1-6)", "Smooth 1",
        "Smooth 2", "SPL Mono", "SPL (1,6)", "2S LMM"
      ))
    )
  }
  if (which=="Power") {
    # plot_data %<>% filter(method %in% c(
    #   "HH", "Wash 2", "SPL (1-6)", "Smooth 2", "SPL (1,6)", "ETI"
    # ))
    # plot_data %<>% mutate(
    #   method = factor(method, levels=c(
    #     "HH", "Wash 2", "SPL (1-6)", "Smooth 2", "SPL (1,6)", "ETI"
    #   ))
    # )
    # plot_data %<>% mutate(
    #   method = factor(method, levels=c(
    #     "HH", "Wash 1", "Wash 2", "SPL (1-6)", "Smooth 1",
    #     "Smooth 2", "SPL Mono", "SPL (1,6)", "2S LMM", "ETI"
    #   ))
    # )
  }

  # Plot: Estimates and CIs
  # Export: 800 x 500
  if (which=="Estimates and CIs") {

    # # Plot for LTE
    # ggplot(
    #   data = plot_data,
    #   aes(
    #     x = method,
    #     y = estimate,
    #     color = method,
    #     ymin = estimate-(1.96*se),
    #     ymax = estimate-(1.96*se)
    #   )
    # ) +
    #   geom_point(aes(color=method), size=2)+
    #   geom_errorbar(
    #     aes(
    #       ymin = estimate-(1.96*se),
    #       ymax = estimate+(1.96*se),
    #       color = method
    #     ),
    #     width = 0.2,
    #     cex = 1
    #   ) +
    #   geom_point(
    #     data = data.frame(
    #       dgm = "SPL (k=2,4 s=0.1,0.4)",
    #       tau = c("Tau=0","Tau=1"),
    #       method = "2S LMM",
    #       estimate = log(0.5),
    #       se = 0
    #     ),
    #     pch = 7,
    #     color = "darkred"
    #   ) +
    #   geom_hline(
    #     yintercept = log(0.5),
    #     linetype = 2
    #   ) +
    #   labs(
    #     title = "Average point estimates and CIs (~1,000 sims per level)",
    #     x = "Analysis method",
    #     y = NULL
    #   ) +
    #   facet_grid(rows=vars(tau), cols=vars(dgm)) +
    #   theme(
    #     axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
    #     legend.position = "none"
    #   )

    # Plot for ATE

    for (th in c(0.5)) {
    # for (th in c(0.5, 0.9)) {

      plot <- ggplot(
        data = plot_data %>% filter(theta==log(th)),
        aes(
          x = method,
          y = estimate,
          color = method,
          ymin = estimate-(1.96*se),
          ymax = estimate-(1.96*se)
        )
      ) +
        geom_hline(
          mapping = aes(yintercept=y),
          data = data.frame(
            dgm = c("EXP (d=0)", "EXP (d=1.4)", "SPL (k=1,6 s=0.8,0.04)",
                    "SPL (k=2,4 s=0.1,0.4)"),
            y = c(log(th)*1, log(th)*0.77, log(th)*0.82, log(th)*0.57)
          ),
          linetype = 2
        ) +
        geom_point(aes(color=method), size=2)+
        geom_errorbar(
          aes(
            ymin = estimate-(1.96*se),
            ymax = estimate+(1.96*se),
            color = method
          ),
          width = 0.2,
          cex = 1
        ) +
        labs(
          title = paste0("Average point estimates and CIs (200 sims per level): ",
                         "theta=log(",th,")"),
          x = "Analysis method",
          y = NULL
        ) +
        facet_grid(rows=vars(tau), cols=vars(dgm)) +
        theme(
          axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
          legend.position = "none"
        )

      print(plot)

    }

  }

  # Plot: Power
  # Export: 800 x 500
  if (which=="Power") {

    ggplot(
      data = plot_data,
      aes(
        x = theta_log,
        y = power,
        color = method
      )
    ) +
      geom_line(aes(group=method)) +
      geom_point(size=0.9) +
      labs(
        title = "Power of CI-based hypothesis test (200 sims per level)",
        x = "Theta",
        color = "Analysis method",
        y = NULL
      ) +
      scale_color_manual(values=c("grey4","red","blue","green","purple","lightblue")) +
      facet_grid(rows=vars(tau), cols=vars(dgm)) +
      theme(
        axis.text.x = element_text(angle=90, hjust=0, vjust=0.4)
      )

  }

}



#########################################################.
##### MISC: Graphs of delay models (for manuscript) #####
#########################################################.

if (run_misc) {

  # Generate data
  d1 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,1), slopes=1))
  d2 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,2,3), slopes=c(0,1)))
  d3 <- effect_curve(seq(0,6,0.1), type="exp", params=list(d=1.4))
  d4 <- effect_curve(seq(0,6,0.1), type="non-monotonic", params=NULL)
  # d5 <- sapply(seq(0,6,0.1), function(x) { sin(((pi*x)/12)-pi/2)+1 })
  d5 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,2,4), slopes=c(0.1,0.4)))

  curve_labels <- c("(a) Instantaneous","(b) Lagged","(c) Curved",
                    "(d) Non-monotonic","(e) Partially convex")

  # Plot functions
  # Export: PDF 8"x3"
  ggplot(
    data.frame(
      x = rep(seq(0,6,0.1),5),
      y = c(d1,d2,d3,d4,d5),
      fn = factor(rep(curve_labels, each=61), levels=curve_labels)
    ),
    aes(x=x, y=y)
  ) +
    geom_line() +
    facet_wrap(~fn, ncol=5) +
    scale_y_continuous(labels=percent) +
    labs(x="Time (steps)", y="Percent of maximum effect achieved")

}

###########################################################.
##### MISC: Graphs to illustrate "modified X" approach #####
###########################################################.

if (run_misc) {

  # Generate data
  x <- rep(c(1,2,3,4,5,6), each=10)
  x_mod <- c( rep(c(1,2,3), each=10), rep(4, 30) )
  y <- pmin(x,4)^1.5+1 + rnorm(length(x))

  # Scatterplot with original X-values
  # Export: 500 x 400
  spl1 <- gam(y~s(x,k=5,bs="cr"))
  ggplot(data.frame(x=x,y=y), aes(x=x, y=y)) +
    geom_point(alpha=0.3) +
    xlim(c(1,6)) +
    labs(title="Original X coordinates") +
    geom_line(
      aes(x=x2, y=y2),
      data = data.frame(x2=x, y2=fitted(spl1)),
      color = "forestgreen"
    )

  # Scatterplot with modified X-values
  # Export: 500 x 400
  spl2 <- gam(y~s(x_mod,k=4,bs="cr"))
  ggplot(data.frame(x=x_mod, y=y), aes(x=x, y=y)) +
    geom_point(alpha=0.3) +
    xlim(c(1,6)) +
    labs(title="Modified X coordinates") +
    geom_line(
      aes(x=x2, y=y2),
      data = data.frame(x2=x, y2=fitted(spl2)),
      color = "forestgreen"
    )


}



#############################################.
##### MISC: SW design diagram for paper #####
#############################################.

if (run_misc) {

  library(ggpubr)

  data <- generate_dataset(
    alpha = log(0.1),
    tau = 1,
    theta = log(0.5),
    n_clusters = 8,
    n_time_points = 5,
    n_ind_per_cluster = 10,
    data_type = "normal",
    sigma = 0.01,
    delay_model = list(type="exp", params=list(d=1.4))
  )

  # Export: PDF 6"x3"
  plot_sw_design(
    data,
    title = c("Stepped wedge design", "Parallel design"),
    compare_to_parallel = TRUE
  )

}



################################.
##### MISC: Graphs for PPT #####
################################.

if (run_misc) {

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
      fun = function(x) {
        ifelse(x<=2, 0, ifelse(x>=3, 3, 3*(x-2)^2))
      },
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
      fun = function(x) {
        ifelse(x<=2, 0, ifelse(x>=5, 3, 1.507543*(1+sin((x-2-(pi/2))))))
      },
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
      fun = function(x) {
        ifelse(x<=2, 0, 3*(1-exp(-(x-2)/1)))
      },
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
      fun = function(x) {
        ifelse(x<=2, 0, ifelse(x>=5, 2.9,
          2.1*(x-2) + (0.8/2 - 2.1)*pmax(0,(x-2)-1)
        ))
      },
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



###################################################.
##### TESTING: SS vs ETI for pathological DGM #####
###################################################.
if (run_testing) {

  # !!!!!
  {source("generate_dataset.R")
  source("one_simulation.R")
  source("plot_outcome.R")
  source("plot_sw_design.R")
  source("run_analysis.R")
  source("effect_curve.R")}

  # Set up and configure simba object
  sim <- new_sim()
  sim %<>% set_config(
    num_sim = 10,
    stop_at_error = TRUE,
    packages = c("dplyr", "magrittr", "stringr", "lme4", "rjags", # "geepack", "restriktor", "scam", "gamlss"
                 "glmmTMB", "mgcv", "fastDummies", "scales", "car")
  )
  sim %<>% add_constants(
    alpha = log(0.1)
  )

  # Add functions to simba object
  sim %<>% add_creator(generate_dataset)
  sim %<>% set_script(one_simulation)
  sim %<>% add_method(run_analysis)
  sim %<>% add_method(effect_curve)

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    theta = -0.5,
    tau = 1,
    sigma = 2.1,
    data_type = "normal",
    method = c("ETI", "SS"),
    delay_model = list(
      "spl" = list(
        type = "spline",
        params = list(knots=c(0:6),slopes=c(1,-1,0.8,-0.8,1,-0.8))
      )
    )
  )

  sim %<>% run()
  s <- sim %>% summary()
  spl_true <- -0.5 * effect_curve(
    x = c(0:6),
    type = "spline",
    params = list(knots=c(0:6), slopes=c(1,-1,0.8,-0.8,1,-0.8))
  )
  spl_eti <- c(0,s[1,12],s[1,13],s[1,14],s[1,15],s[1,16],s[1,17])
  spl_ss <- c(0,s[2,12],s[2,13],s[2,14],s[2,15],s[2,16],s[2,17])

  ggplot(
    data.frame(
      x = rep(c(0:6),3),
      y = c(spl_true,spl_eti,spl_ss),
      which = rep(c("True","ETI","SS"), each=7)
    ),
    aes(x=x, y=y, color=which)) +
    geom_line()

}
