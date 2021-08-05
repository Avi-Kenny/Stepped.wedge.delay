# Title: "Stepped wedge lagged effect simulation"
# Author: Avi Kenny

##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  level_set_which = "level_set_test",
  run_or_update = "run",
  num_sim = 1000,
  pkgs = c("dplyr", "stringr", "lme4", "Iso", "sqldf", "mgcv", "MASS",
           "fastDummies", "car", "splines", "glmmTMB", "rstan"),
  pkgs_nocluster = c("ggplot2", "viridis", "scales", "facetscales", "glmmTMB", # devtools::install_github("zeehio/facetscales")
                     "readxl"),
  parallel = "none",
  stop_at_error = FALSE,
  # mcmc = list(n.adapt=300, n.iter=300, n.burn=300, n.chains=1, thin=1)
  mcmc = list(n.adapt=2000, n.iter=3000, n.burn=1000, n.chains=3, thin=1)
  # mcmc = list(n.adapt=2000, n.iter=5000, n.burn=3000, n.chains=3, thin=1)
)

# Set cluster config
cluster_config <- list(
  js = "slurm",
  dir = "/home/akenny/z.stepped.wedge"
  # js = "sge",
  # dir = "/home/users/avikenny/Desktop/z.stepped.wedge"
)
if (cluster_config$js=="sge") {
  cfg$pkgs <- cfg$pkgs[cfg$pkgs!="rstan"]
}



#################.
##### Setup #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Jim Hughes/Project - Stepped wedge lag/z.stepped.wedge/R"))
  load_pkgs_local <- TRUE
} else {
  # Cluster
  setwd("z.stepped.wedge/R")
  load_pkgs_local <- FALSE
}

# Load packages (if running locally)
if (load_pkgs_local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    do.call("library", list(pkg))
  }
}

# Load simba + functions
{
  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
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
  run_realdata <- FALSE
  run_misc <- FALSE
  run_testing <- FALSE
}



######################################################################.
##### TESTING: Generate dataset for testing new analysis methods #####
######################################################################.

if (FALSE) {

  # Generate dataset
  data <- generate_dataset(
    mu = 1,
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 20,
    theta = 0.5,
    tau = 0.25, # 0.5
    sigma = 1, # 2
    data_type = "normal",
    # delay_model = list(type="spline", params=list(knots=c(0,1),slopes=1)),
    # delay_model = list(type="exp", params=list(d=1.5)),
    delay_model = list(
      type = "spline",
      params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
    ),
    # n_extra_time_points = 2,
    n_extra_time_points = 0,
    rte = NA,
    # rte = list(type="height", rho=-0.2, nu=0.4)
    # rte = list(type="height+shape", rho1=-0.1, rho2=0.6, nu=0.4)
    time_trend = "none" # incr
  )

  # # Set variables needed in run_analysis
  J <- data$params$n_time_points
  data_type <- "normal"

}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (run_main) {
  if (Sys.getenv("simba_run") %in% c("first", "")) {

    delay_models <- list(
      "Instantaneous" = list(
        type = "spline",
        params = list(knots=c(0,0.1), slopes=10)
      ),
      "Lagged" = list(
        type = "spline",
        params = list(knots=c(0,2,2.1), slopes=c(0,10))
      ),
      "Curved" = list(
        type = "exp",
        params = list(d=1.5)
      ),
      "Partially convex" = list(
        type = "spline",
        params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
      )
    )

    # !!!!! Testing
    level_set_test <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = 0.5,
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI"),
        "MEC: Flat Dirichlet" = list(method = "MCMC-MON-Stan", mcmc = cfg$mcmc,
                                       enforce="Flat Dirichlet")
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("whole_curve"=list(whole_curve=TRUE))
    )

    # Simulation 1: pitfalls of the immediate treatment (IT) model
    # Simulation 2: estimation of the TATE and LTE
    # Simulation 3: estimation of the entire effect curve
    # 16 level combos
    level_set_123 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = 0.5,
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list(
        "IT" = list(method="IT"),
        "ETI" = list(method="ETI"),
        "NCS (4df)" = list(method="NCS-4df"),
        "MEC" = list(
          method = "MCMC-MON-Stan", enforce="simplex",
          mcmc = cfg$mcmc)
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("whole_curve"=list(whole_curve=TRUE))
    )

    # Simulation 4: power of Wald-type hypothesis tests
    # 96 level combos
    level_set_4 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = seq(0,0.5,0.1),
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list(
        "IT" = list(method="IT"),
        "ETI" = list(method="ETI"),
        "NCS (4df)" = list(method="NCS-4df"),
        "MEC" = list(
          method = "MCMC-MON-Stan", enforce="simplex",
          mcmc = cfg$mcmc)
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("none"=list())
    )

    # Simulation 5: performance of RETI models
    # 12 level combos
    level_set_5 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = 0.5,
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI", effect_reached=0),
        "RETI (3 steps)" = list(method="ETI", effect_reached=3),
        "RETI (4 steps)" = list(method="ETI", effect_reached=4)
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("none"=list())
    )

    # Simulation 6: performance of RTE models; data are generated with no RTEs
    # Simulation 7: performance of RTE models; data are generated with RTEs
    level_set_67 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = 0.5,
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI"),
        "ETI (RTE; height)" = list(method="ETI", re="height")
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = list(
        "none" = list(),
        "height" = list(type="height", nu=1, rho1=-0.2)
      ),
      return_extra = list("rte"=list(rte=TRUE))
    )

    # Simulation 8: Simulation results: effect of adding extra time points
    # 12 level combos
    level_set_8 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # 50
      theta = 0.5,
      tau = 0.25, # 0.5
      sigma = 1, # 2
      data_type = "normal",
      method = list("ETI" = list(method="ETI")),
      delay_model = delay_models,
      n_extra_time_points = c(0,1,2),
      rte = NA,
      return_extra = list("none"=list())
    )

    level_set <- eval(as.name(cfg$level_set_which))

  }
}



#################################.
##### MAIN: Main simulation #####
#################################.

# Commands for job sumbission on Slurm:
# sbatch --export=simba_run='first',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-16 --export=simba_run='main',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=simba_run='last',cluster='bionic',type='R',project='z.stepped.wedge' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

# Commands for job sumbission on SGE:
# qsub -v simba_run='first',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992344 -t 1-3 -v simba_run='main',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992345 -v simba_run='last',cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh

if (run_main) {

  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")

  if (cfg$run_or_update=="run") {

    run_on_cluster(

      first = {

        # Set up and configure simba object
        sim <- new_sim()
        sim %<>% set_config(
          num_sim = cfg$num_sim,
          parallel = cfg$parallel,
          stop_at_error = cfg$stop_at_error,
          seed = as.integer(1000*runif(1)),
          packages = cfg$pkgs
        )
        sim %<>% add_constants(mu = 1)

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
        # sim %>% summarize() %>% print()
        sim$results %>% head() %>% print()
        sim$errors %>% print()
      },

      cluster_config = cluster_config

    )

  }

  if (cfg$run_or_update=="update") {

    update_sim_on_cluster(

      first = {
        sim <- readRDS(paste0(cluster_config$dir,'/sim.simba'))
        sim <- do.call(set_levels, c(list(sim), level_set))
      },

      main = {
        sim %<>% update_sim()
      },

      last = {},

      cluster_config = cluster_config

    )

  }

}



############################################.
##### MAIN: Process simulation results #####
############################################.

if (run_process_results) {

  # Set simulation
  whichsim <- 2

  # Read in simulation object
  file <- case_when(
    whichsim %in% c(1,2,3) ~ "sim123_20210726.simba",
    whichsim==4 ~ "sim4_20210727.simba",
    whichsim==5 ~ "sim5_20210727.simba",
    whichsim %in% c(6,7) ~ "sim67_20210801.simba",
    whichsim==8 ~ "sim8_20210729.simba"
  )
  sim <- readRDS(paste0("../simba.out/",file))

  # Generate true TATE values
  # Note: now using step function approximations
  sim$results %<>% mutate(
    ate = case_when(
      delay_model=="Instantaneous" ~ theta,
      delay_model=="Lagged" ~ theta * mean(
        effect_curve(x = c(1:6),
                     type = "spline",
                     params = list(knots=c(0,2,2.1),slopes=c(0,10)))),
      delay_model=="Curved" ~ theta * mean(
        effect_curve(x = c(1:6),
                     type = "exp",
                     params = list(d=1.5))),
      delay_model=="Partially convex" ~ theta * mean(
        effect_curve(x = c(1:6),
                     type = "spline",
                     params = list(knots=c(0,2,4),slopes=c(0.1,0.4))))
    ),
    method = case_when(
      method=="NCS (4df)" ~ "NCS-4",
      method=="ETI (RTE; height)" ~ "RTE",
      TRUE ~ method,
    )
  )

  # Generate true theta values
  if (whichsim %in% c(1,3)) {
    theta <- sim$results[1,"theta"]
    thetas <- list(
      "I" = rep(theta,6),
      "L" = theta * effect_curve(x = c(1:6),
                                 type = "spline",
                                 params = list(knots=c(0,2,2.1),slopes=c(0,10))),
      "C" = theta * effect_curve(x = c(1:6),
                                 type = "exp",
                                 params = list(d=1.5)),
      "P" = theta * effect_curve(x = c(1:6),
                                 type = "spline",
                                 params = list(knots=c(0,2,4),slopes=c(0.1,0.4)))
    )
    thetas_df <- data.frame(
      delay_model = c("Instantaneous", "Lagged", "Curved", "Partially convex"),
      theta_1 = c(thetas$I[1],thetas$L[1],thetas$C[1],thetas$P[1]),
      theta_2 = c(thetas$I[2],thetas$L[2],thetas$C[2],thetas$P[2]),
      theta_3 = c(thetas$I[3],thetas$L[3],thetas$C[3],thetas$P[3]),
      theta_4 = c(thetas$I[4],thetas$L[4],thetas$C[4],thetas$P[4]),
      theta_5 = c(thetas$I[5],thetas$L[5],thetas$C[5],thetas$P[5]),
      theta_6 = c(thetas$I[6],thetas$L[6],thetas$C[6],thetas$P[6])
    )
    sim$results %<>% inner_join(thetas_df, by="delay_model")
    sim$results %<>% mutate(
      mse_wc = (1/6) * (
        (theta_1_hat-theta_1)^2+(theta_2_hat-theta_2)^2+(theta_3_hat-theta_3)^2+
          (theta_4_hat-theta_4)^2+(theta_5_hat-theta_5)^2+(theta_6_hat-theta_6)^2
      )
    )

    summ_mean <- list(
      list(name="ate", x="ate"),
      list(name="mean_ate", x="ate_hat"),
      list(name="mse_wc", x="mse_wc"),
      list(name="mean_lte", x="lte_hat"),
      list(name="theta_1", x="theta_1"),
      list(name="theta_2", x="theta_2"),
      list(name="theta_3", x="theta_3"),
      list(name="theta_4", x="theta_4"),
      list(name="theta_5", x="theta_5"),
      list(name="theta_6", x="theta_6")
    )

  } else {
    summ_mean <- list(
      list(name="ate", x="ate"),
      list(name="mean_ate", x="ate_hat"),
      list(name="mean_lte", x="lte_hat")
    )
  }

  # Summarize data
  summ <- summarize(
    sim_obj = sim,
    mean = summ_mean,
    bias_pct = list(
      list(name="bias_ate", estimate="ate_hat", truth="ate"),
      list(name="bias_lte", estimate="lte_hat", truth="theta")
    ),
    mse = list(
      list(name="mse_ate", estimate="ate_hat", truth="ate"),
      list(name="mse_lte", estimate="lte_hat", truth="theta")
    ),
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
  s_d_models <- c("(a) Instantaneous","(b) Lagged","(c) Curved",
                  "(d) Partially convex")
  summ %<>% mutate(
    delay_model = case_when(
      delay_model=="Instantaneous" ~ s_d_models[1],
      delay_model=="Lagged" ~ s_d_models[2],
      delay_model=="Curved" ~ s_d_models[3],
      delay_model=="Partially convex" ~ s_d_models[4]
    )
  )
  s_methods <- unique(sim$results$method)
  if (whichsim==5) s_methods <- c(s_methods[1],s_methods[3],s_methods[2])
  if (whichsim==1) summ %<>% filter(method %in% c("IT"))
  if (whichsim==2) summ %<>% filter(method!="IT")
  if (whichsim==6) summ %<>% filter(rte=="none")
  if (whichsim==7) summ %<>% filter(rte=="height")
  summ %<>% mutate(
    method = factor(method, levels=s_methods),
    delay_model = factor(delay_model, levels=s_d_models),
    power_ate = 1 - beta_ate,
    power_lte = 1 - beta_lte
  )
  summ %<>% rename("Method"=method, "n_extra"=n_extra_time_points)
  p_data <- sqldf("
    SELECT Method, n_extra, delay_model, 'TATE' AS which, bias_ate AS bias,
      theta, cov_ate AS Coverage, power_ate AS Power, mse_ate AS MSE FROM summ
    UNION SELECT Method, n_extra, delay_model, 'LTE', bias_lte,
      theta, cov_lte, power_lte, mse_lte FROM summ
  ")
  p_data <- sqldf("
    SELECT Method, n_extra, delay_model, which, 'Bias' AS stat, Bias AS value,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'Coverage', Coverage,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'Power', Power,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'MSE', MSE,
      theta FROM p_data
  ")
  p_data %<>% mutate(which = factor(which, levels=c("TATE", "LTE")))
  p_data %<>% mutate(n_extra = as.character(n_extra))
  p_data %<>% rename("Extra time points"=n_extra)

  if (whichsim %in% c(1,3)) {

    summ_it <- filter(summ, Method=="IT")

    p_data2 <- sqldf("
      SELECT delay_model, 'True effect curve' AS which, theta_1 AS value,
        1 AS time FROM summ_it
      UNION SELECT delay_model, 'True effect curve', theta_2, 2 FROM summ_it
      UNION SELECT delay_model, 'True effect curve', theta_3, 3 FROM summ_it
      UNION SELECT delay_model, 'True effect curve', theta_4, 4 FROM summ_it
      UNION SELECT delay_model, 'True effect curve', theta_5, 5 FROM summ_it
      UNION SELECT delay_model, 'True effect curve', theta_6, 6 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 1 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 2 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 3 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 4 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 5 FROM summ_it
      UNION SELECT delay_model, 'IT estimate', mean_ate, 6 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 1 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 2 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 3 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 4 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 5 FROM summ_it
      UNION SELECT delay_model, 'True TATE', ate, 6 FROM summ_it

    ")

  }

  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    IT = cb_colors[2],
    ETI = cb_colors[4],
    `NCS-4` = cb_colors[3],
    `MEC` = cb_colors[6],
    `RETI (3 steps)` = cb_colors[6],
    `RETI (4 steps)` = cb_colors[7],
    `RTE` = cb_colors[6],
    `0` = cb_colors[4],
    `1` = cb_colors[7],
    `2` = cb_colors[6]
  )

}



########################################.
##### VIZ: Figure for simulation 1 #####
########################################.

if (run_viz) {

  # Export: 8: x 3"
  ggplot(
    filter(p_data2, which!="True TATE"),
    aes(x=time, y=value, color=which, shape=which, group=which)
  ) +
    geom_hline(
      aes(yintercept=value, linetype="True TATE"),
      filter(p_data2, which=="True TATE" & time==1),
      color = cb_colors[6]
    ) +
    geom_line() +
    geom_point() +
    labs(
      x = "Study time",
      color = NULL,
      shape = NULL,
      linetype = NULL,
      y = "Treatment effect"
    ) +
    theme(legend.position="bottom") +
    scale_color_manual(values=cb_colors[c(2,3,4)]) +
    scale_linetype_manual(
      values = "longdash",
      guide = guide_legend(override.aes=list(color=cb_colors[6]))
    ) +
    scale_x_continuous(breaks=c(1:6)) +
    facet_grid(cols=dplyr::vars(delay_model)) +
    scale_y_continuous(labels=percent_format())

}



#############################################.
##### VIZ: Figure for simulations 2,5-7 #####
#############################################.

if (run_viz) {

  # # Table for sim1
  # if (FALSE) {
  #   p_data %>% filter(Method=="IT" & stat!="Power" & which=="LTE") %>%
  #     mutate(value=round(value,3)) %>%
  #     arrange(stat, delay_model, which)
  # }

  # Sim-specific ggplot config
  if (whichsim==2) {
    mc <- c(2,3,4)
    y_lims <- c(0.7,1)
  }
  if (whichsim==5) {
    mc <- c(2,5,6)
    y_lims <- c(0.5,1)
  }
  if (whichsim==6) {
    mc <- c(2,7)
    y_lims <- c(0.7,1)
  }
  if (whichsim==7) {
    mc <- c(2,7)
    y_lims <- c(0.4,1)
  }

  # Export: 8: x 4"
  ggplot(
    filter(p_data, stat!="Power"),
    aes(x=which, y=value, fill=Method)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8,
             color="#555555", size=0.35, alpha=0.8) +
    facet_grid_sc(
      cols = dplyr::vars(delay_model),
      rows = dplyr::vars(stat),
      scales=list(y=list(
        Bias = scale_y_continuous(labels=percent_format()),
        # Bias = scale_y_continuous(labels=percent_format(),
        #                           limits=c(-0.1,0.1)),
        Coverage = scale_y_continuous(labels=percent_format(),
                                      limits=y_lims),
        MSE = scale_y_continuous()
      ))
    ) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors[mc]) +
    labs(y=NULL, x=NULL, fill="Analysis model")

}



########################################.
##### VIZ: Figure for simulation 3 #####
########################################.

if (run_viz) {

  # Export: 8: x 3"
  summ2 <- filter(summ, Method!="IT")
  ggplot(
    summ2,
    aes(x=Method, y=mse_wc, fill=Method) # x=which,
  ) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8,
             color="#555555", size=0.35, alpha=0.8) +
    facet_grid(cols=dplyr::vars(delay_model)) +
    theme(
      legend.position="bottom",
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    ) +
    scale_fill_manual(values=m_colors[c(2,3,4)]) +
    labs(y="Average pointwise MSE", x=NULL, fill="Analysis model")

}



########################################.
##### VIZ: Figure for simulation 4 #####
########################################.

if (run_viz) {

  # Export: 8: x 4"
  ggplot(
    filter(p_data, stat=="Power"),
    aes(x=theta, y=value, color=Method, shape=Method, group=Method)
  ) +
    geom_line() +
    geom_point() +
    labs(
      x = unname(latex2exp::TeX("$\\delta$ (maximum effect size)")),
      color = "Analysis model",
      shape = "Analysis model",
      y = "Power"
    ) +
    theme(legend.position="bottom") +
    scale_color_manual(values=m_colors[c(1:4)]) +
    # scale_shape_manual(values=c(15,16,17)) +
    facet_grid(rows=dplyr::vars(which), cols=dplyr::vars(delay_model)) +
    scale_y_continuous(labels=percent_format()) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4)
    )

}



########################################.
##### VIZ: Figure for simulation 8 #####
########################################.

if (run_viz) {

  # Export: 8: x 4"
  ggplot(
    filter(p_data, stat!="Power"),
    aes(x=which, y=value, fill=`Extra time points`)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8,
             color="#555555", size=0.35, alpha=0.8) +
    facet_grid_sc(
      cols = dplyr::vars(delay_model),
      rows = dplyr::vars(stat),
        scales=list(y=list(
        Bias = scale_y_continuous(labels = percent_format(),
                                  limits=c(-0.1,0.1)),
        Coverage = scale_y_continuous(labels=percent_format(),
                                      limits=c(0.7,1)),
        MSE = scale_y_continuous()
      ))
    ) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors[c(8,9,10)]) +
    labs(y=NULL, x=NULL, fill="Extra time points")

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
    # facet_grid(rows=dplyr::vars(date), cols=dplyr::vars(dgm)) +
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
    facet_grid(rows=dplyr::vars(data_type), cols=dplyr::vars(delay_model)) +
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
    facet_grid(rows=dplyr::vars(data_type), cols=dplyr::vars(delay_model)) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("turquoise", "salmon", "dodgerblue2"))

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
  summ <- summarize(
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
    method = factor(method, levels=c("IT","SS","MSS"))
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
    facet_grid(rows=dplyr::vars(date), cols=dplyr::vars(dgm)) +
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
    summ <- summarize(
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
        "IT", "Wash 1", "Wash 2", "SPL (1-6)", "Smooth 1",
        "Smooth 2", "SPL Mono", "SPL (1,6)", "2S LMM"
      ))
    )
  }
  if (which=="Power") {
    # plot_data %<>% filter(method %in% c(
    #   "IT", "Wash 2", "SPL (1-6)", "Smooth 2", "SPL (1,6)", "ETI"
    # ))
    # plot_data %<>% mutate(
    #   method = factor(method, levels=c(
    #     "IT", "Wash 2", "SPL (1-6)", "Smooth 2", "SPL (1,6)", "ETI"
    #   ))
    # )
    # plot_data %<>% mutate(
    #   method = factor(method, levels=c(
    #     "IT", "Wash 1", "Wash 2", "SPL (1-6)", "Smooth 1",
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
    #   facet_grid(rows=dplyr::vars(tau), cols=dplyr::vars(dgm)) +
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
        facet_grid(rows=dplyr::vars(tau), cols=dplyr::vars(dgm)) +
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
      facet_grid(rows=dplyr::vars(tau), cols=dplyr::vars(dgm)) +
      theme(
        axis.text.x = element_text(angle=90, hjust=0, vjust=0.4)
      )

  }

}



######################################################.
##### MAIN: Real data analysis: helper functions #####
######################################################.

if (run_realdata) {

  #' Extract estimates and SEs
  #'
  #' @param model One of c("")
  #' @param type One of c("IT", "ETI", "NCS", "RETI")
  #' @param J Number of unique study time points
  #' @param len Number of unique exposure time points (excluding zero)
  #' @param R The "effect reached" value for the RETI model
  est <- function(model, type, J, len, R=NA) {

    # Extract treatment effect estimates
    {
      coeff_names <- names(summary(model)$coefficients$cond[,1])

      if (type=="IT") {

        delta <- as.numeric(summary(model)$coefficients$cond[,1]["x_ij"])
        se <- as.numeric(sqrt(diag(summary(model)$vcov[[1]])["x_ij"]))
        delta <- rep(delta,len)
        se <- rep(se,len)

      } else if (type=="ETI") {

        ind_tx <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
        delta <- as.numeric(summary(model)$coefficients$cond[,1][ind_tx])
        vcov <- summary(model)$vcov[[1]][ind_tx,ind_tx]
        se <- as.numeric(sqrt(diag(vcov)))

      } else if (type=="RETI") {

        ind_tx <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
        delta <- as.numeric(summary(model)$coefficients$cond[,1][ind_tx])
        vcov <- summary(model)$vcov[[1]][ind_tx,ind_tx]
        se <- as.numeric(sqrt(diag(vcov)))

      } else if (type=="NCS") {

        ns_basis <- ns(c(0:len), knots=c(len/4,(2*len)/4,(3*len)/4))
        ind_spl <- c(1:length(coeff_names))[coeff_names %in% c("c1","c2","c3","c4")]
        b_hat <- as.numeric(summary(model)$coefficients$cond[,1][ind_spl])
        sigma_b_hat <- summary(model)$vcov[[1]][ind_spl,ind_spl]
        B <- matrix(NA, nrow=len, ncol=4)
        for (i in 1:len) { for (j in 1:4) { B[i,j] <- ns_basis[i+1,j] }}
        delta <- as.numeric(B %*% b_hat)
        vcov <- B %*% sigma_b_hat %*% t(B)
        se <- sqrt(diag(vcov))

      } else if (type=="MEC") {

        # Extract beta_s means
        beta_s_hat <- c()
        for (i in 1:len) {
          beta_s_hat[i] <- mean(rstan::extract(model)[[paste0("beta_s_",i)]])
        }

        # Construct covariance matrix of s terms
        sigma_s_hat <- matrix(NA, nrow=len, ncol=len)
        for (i in 1:len) {
          for (j in 1:len) {
            sigma_s_hat[i,j] <- cov(
              rstan::extract(model)[[paste0("beta_s_",i)]],
              rstan::extract(model)[[paste0("beta_s_",j)]],
              use = "complete.obs"
            )
          }
        }

        # Calculate theta_l_hat vector and sigma_l_hat matrix
        B <- matrix(0, nrow=len, ncol=len)
        for (i in 1:len) {
          for (j in 1:len) {
            if (i>=j) B[i,j] <- 1
          }
        }
        delta <- B %*% beta_s_hat
        vcov <- B %*% sigma_s_hat %*% t(B)
        se <- as.numeric(sqrt(diag(vcov)))

      }
    }

    # Calculate time trend estimates and CIs
    if (type=="MEC") {
      # Extract beta_j means and sds
      time <- c()
      se2 <- c()
      for (i in 2:J) {
        time[i-1] <- mean(rstan::extract(model)[[paste0("beta_j_",i)]])
        se2[i-1] <- sd(rstan::extract(model)[[paste0("beta_j_",i)]])
      }
      time <- c(0,time)
      se2 <- c(0,se2)
    } else {

      ind_time <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(j)"]
      time <- c(0,as.numeric(summary(model)$coefficients$cond[,1][ind_time]))
      se2 <- c(0,as.numeric(sqrt(diag(summary(model)$vcov[[1]])[ind_time])))

    }

    # Calculate TATE and LTE
    if (type=="IT") {
      ate_hat <- delta[1]
      se_ate_hat <- se[1]
      lte_hat <- delta[1]
      se_lte_hat <- se[1]
    } else if (type=="RETI") {
      A <- (1/len) * matrix(c(rep(1,R-1),len+1-R), nrow=1)
      ate_hat <- (A %*% delta)[1,1]
      se_ate_hat <- sqrt(A %*% vcov %*% t(A))[1,1]
      lte_hat <- delta[length(delta)]
      se_lte_hat <- se[length(delta)]
    } else if (type %in% c("ETI","NCS","MEC")) {
      A <- (1/len) * matrix(rep(1,len), nrow=1)
      ate_hat <- (A %*% delta)[1,1]
      se_ate_hat <- sqrt(A %*% vcov %*% t(A))[1,1]
      lte_hat <- delta[len]
      se_lte_hat <- se[len]
    }

    # RETI modification
    if (type=="RETI") {
      len2 <- length(delta)
      delta <- c(delta[1:(len2-1)],rep(delta[len2],len+1-R))
      se <- c(se[1:(len2-1)],rep(se[len2],len+1-R))
    }

    return(list(
      delta = c(0,delta),
      ci_l = c(0,delta)-1.96*c(0,se),
      ci_u = c(0,delta)+1.96*c(0,se),
      time = time,
      ci2_l = time-1.96*se2,
      ci2_u = time+1.96*se2,
      ate_hat = ate_hat,
      ate_l = ate_hat-1.96*se_ate_hat,
      ate_u = ate_hat+1.96*se_ate_hat,
      lte_hat = lte_hat,
      lte_l = lte_hat-1.96*se_lte_hat,
      lte_u = lte_hat+1.96*se_lte_hat
    ))

  }

  #' Print model results (TATE and LTE)
  #'
  #' @param models A named list of estimate objects (each returned by `est()`)
  #' @param which Which data analysis; one of c("disinvestment","wa_state")
  print_results <- function(models, which) {

    # Code specific to analyses
    if (which=="disinvestment") {
      trans <- function(x) {x}
    } else if (which=="wa_state") {
      trans <- function(x) {exp(x)}
    }

    # Print results
    for (i in 1:length(models)) {

      print(paste("Model:",names(models)[i]))
      print(paste0("TATE: ", round(trans(models[[i]]$ate_hat),4),
                   " (", round(trans(models[[i]]$ate_l),4), " -- ",
                   round(trans(models[[i]]$ate_u),4), ")"))
      print(paste0("LTE: ", round(trans(models[[i]]$lte_hat),4),
                   " (", round(trans(models[[i]]$lte_l),4), " -- ",
                   round(trans(models[[i]]$lte_u),4), ")"))

    }

  }

  #' Print model results (TATE and LTE)
  #'
  #' @param models A named list of estimate objects (each returned by `est()`)
  #' @param which Which data analysis; one of c("disinvestment","wa_state")
  #' @param ncol Number of columns for facet_wrap()
  print_graphs <- function(models, which, ncol) {

    # Helper function to extract attributes; models accessed globally
    v <- function(attr) {
      as.numeric(sapply(models, function(m) { m[[attr]] }))
    }

    # Extract info
    labels <- names(models)
    n_models <- length(models)
    len_tx <- length(models[[1]]$delta)
    len_time <- length(models[[1]]$time)

    # Code specific to analyses
    if (which=="disinvestment") {
      trans <- function(x) {x}
      lab_y_tx <- "Treatment effect"
      lab_y_time <- "Time effect"
    } else if (which=="wa_state") {
      trans <- function(x) {exp(x)}
      lab_y_tx <- "Treatment effect (odds ratio)"
      lab_y_time <- "Time effect: (odds ratio)"
    }

    # Treatment effect estimates
    # Export: 7" x 4"
    plot_tx <- ggplot(
      data.frame(
        x = rep(c(0:(len_tx-1)),n_models),
        y = trans(v("delta")),
        ci_l = trans(v("ci_l")),
        ci_u = trans(v("ci_u")),
        model = rep(factor(labels, levels=labels), each=len_tx)
      ),
      aes(x=x, y=y, color=model)
    ) +
      facet_wrap(~model, ncol=ncol) +
      geom_line() +
      geom_ribbon(
        aes(ymin=ci_l, ymax=ci_u, fill=as.factor(model)),
        alpha = 0.2,
        linetype = "dotted"
      ) +
      theme(legend.position="none") + # bottom
      labs(y=lab_y_tx, x="Exposure time", color="Model", fill="Model")

    # Time trend estimates
    # Export: 7" x 4"
    plot_time <- ggplot(
      data.frame(
        x = rep(c(1:len_time),n_models),
        y = trans(v("time")),
        ci_l = trans(v("ci2_l")),
        ci_u = trans(v("ci2_u")),
        model = rep(factor(labels, levels=labels), each=len_time)
      ),
      aes(x=x, y=y, color=model)
    ) +
      facet_wrap(~model, ncol=ncol) +
      geom_line() +
      geom_ribbon(
        aes(ymin=ci_l, ymax=ci_u, fill=as.factor(model)),
        alpha = 0.2,
        linetype = "dotted"
      ) +
      theme(legend.position="none") + # bottom
      labs(y=lab_y_time, x="Study time", color="Model", fill="Model")

    # Display graphs
    print(plot_tx)
    print(plot_time)

  }

}



#######################################################.
##### MAIN: Real data analysis #1 (Disinvestment) #####
#######################################################.

if (run_realdata) {

  # Read/process data
  df <- read_excel("../realdata/disinvestment/S2 Data.xlsx")
  df %<>% filter(study1==1)
  df %<>% mutate(
    y = log(acute_los),
    i = factor(index_ward),
    j = factor(sw_step)
  )
  df %<>% rename(
    "x_ij" = no_we_exposure
  )
  df %<>% subset(select=c(y,i,j,x_ij))
  df2 <- df %>% group_by(i,j) %>%
    dplyr::summarize(x_ij=x_ij[1]) %>%
    arrange(i,j)
  df2$l <- c(0,rep(NA,nrow(df2)-1))
  for (row in 2:nrow(df2)) {
    if (df2[row,"i"]==df2[row-1,"i"]) {
      df2[row,"l"] <- df2[row-1,"l"] + df2[row,"x_ij"]
    } else {
      df2[row,"l"] <- 0
    }
  }
  df2$x_ij <- NULL
  df %<>% inner_join(df2, by=c("i","j"))

  # IT model
  {
    model_it <- glmmTMB(
      y ~ factor(j) + x_ij + (1|i),
      data = df
    )
    est_it <- est(model_it, "IT", J=7, len=length(unique(df$l))-1)
  }

  # ETI model
  {
    model_eti <- glmmTMB(
      y ~ factor(j) + factor(l) + (1|i),
      data = df
    )
    est_eti <- est(model_eti, "ETI", J=7, len=length(unique(df$l))-1)
  }

  # ETI (RTE; height) model
  {
    model_rte <- glmmTMB(
      y ~ factor(j) + factor(l) + (x_ij|i),
      data = df
    )
    est_rte <- est(model_rte, "ETI", J=7, len=length(unique(df$l))-1)
  }

  # RETI model (R=3)
  {
    r <- 3
    df_reti <- df %>% mutate(l = ifelse(l>=r,r,l))
    model_reti <- glmmTMB(
      y ~ factor(j) + factor(l) + (1|i),
      data = df_reti
    )
    est_reti <- est(model_reti, "RETI", J=7, len=length(unique(df$l))-1, R=r)
  }

  # Natural cubic spline (4df)
  {
    J <- 7
    len <- length(unique(df$l))-1
    ns_basis <- ns(c(0:len), knots=c(len/4,(2*len)/4,(3*len)/4))
    df$c1 <- ns_basis[df$l+1,1]
    df$c2 <- ns_basis[df$l+1,2]
    df$c3 <- ns_basis[df$l+1,3]
    df$c4 <- ns_basis[df$l+1,4]
    model_ncs <- glmmTMB(
      y ~ factor(j) + c1+c2+c3+c4 + (1|i),
      data = df
    )
    est_ncs <- est(model_ncs, "NCS", J=7, len=len)
  }

  # Monotone Effect Curve Bayesian model
  {
    # Setup
    mcmc <- cfg$mcmc
    df %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    df %<>% mutate(
      s_1 = as.integer(l>=1),
      s_2 = as.integer(l>=2),
      s_3 = as.integer(l>=3),
      s_4 = as.integer(l>=4),
      s_5 = as.integer(l>=5),
      s_6 = as.integer(l>=6)
    )
    options(mc.cores = parallel::detectCores()-1)
    # Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
    rstan_options(auto_write=TRUE)

    # Put data into Stan format
    stan_data <- list(
      I = length(unique(df$i)),
      N = nrow(df),
      y = df$y,
      i = as.integer(df$i),
      j_2 = df$j_2,
      j_3 = df$j_3,
      j_4 = df$j_4,
      j_5 = df$j_5,
      j_6 = df$j_6,
      j_7 = df$j_7,
      s_1 = df$s_1,
      s_2 = df$s_2,
      s_3 = df$s_3,
      s_4 = df$s_4,
      s_5 = df$s_5,
      s_6 = df$s_6
    )

    # Stan model
    stan_code <- quote("
      data {
        int I;
        int N;
        real y[N];
        int i[N];
        real j_2[N];
        real j_3[N];
        real j_4[N];
        real j_5[N];
        real j_6[N];
        real j_7[N];
        real s_1[N];
        real s_2[N];
        real s_3[N];
        real s_4[N];
        real s_5[N];
        real s_6[N];
      }
      parameters {
        real beta0;
        real beta_j_2;
        real beta_j_3;
        real beta_j_4;
        real beta_j_5;
        real beta_j_6;
        real beta_j_7;
        real delta;
        simplex[6] smp;
        real alpha[I];
        real<lower=0> sigma;
        real<lower=0> tau;
      }
      transformed parameters {
        real beta_s_1;
        real beta_s_2;
        real beta_s_3;
        real beta_s_4;
        real beta_s_5;
        real beta_s_6;
        vector[N] y_mean;
        beta_s_1 = delta * smp[1];
        beta_s_2 = delta * smp[2];
        beta_s_3 = delta * smp[3];
        beta_s_4 = delta * smp[4];
        beta_s_5 = delta * smp[5];
        beta_s_6 = delta * smp[6];
        for (n in 1:N) {
          y_mean[n] = beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] +
          beta_j_7*j_7[n] + delta*(
            smp[1]*s_1[n] + smp[2]*s_2[n] + smp[3]*s_3[n] +
            smp[4]*s_4[n] + smp[5]*s_5[n] + smp[6]*s_6[n]
          ) + alpha[i[n]];
        }
      }
      model {
        delta ~ normal(0,100);
        alpha ~ normal(0,tau);
        smp ~ dirichlet([1,1,1,1,1,1]');
        y ~ normal(y_mean,sigma);
    }")

    # Fit model
    model_mec <- stan(
      model_code = stan_code,
      data = stan_data,
      chains = mcmc$n.chains,
      iter = mcmc$n.burn+mcmc$n.iter,
      warmup = mcmc$n.burn,
      thin = mcmc$thin
    )
    est_mec <- est(model_mec, "MEC", J=7, len=6)

    # !!!!! MCMC diagnostics
    if (FALSE) {
      # vars <- c("beta_j_2", "beta_j_3", "beta_j_4",
      #           "beta_j_5", "beta_j_6", "beta_j_7")
      vars <- c("beta_s_1", "beta_s_2", "beta_s_3",
                "beta_s_4", "beta_s_5", "beta_s_6")
      # vars <- c("beta0", "delta", "omega")
      # vars <- "smp"
      # vars <- "alpha"
      # vars <- c("sigma", "tau")
      rstan::traceplot(fit, vars, inc_warmup=TRUE)
    }

  }

  # Display results
  models <- list("IT"=est_it, "ETI"=est_eti, "RTE"=est_rte,
                 "RETI-3"=est_reti, "NCS-4"=est_ncs, "MEC"=est_mec)
  print_results(models=models, which="disinvestment")
  print_graphs(models=models, which="disinvestment", ncol=3)

}



##################################################.
##### MAIN: Real data analysis #2 (WA state) #####
##################################################.

if (run_realdata) {

  # Read/process data
  df <- read.csv("../realdata/wa_state/wa_state.csv")
  df %<>% rename(
    "y" = ct,
    "i" = hdist,
    "j" = timex,
    "l" = timeontrtx,
    "x_ij" = treatx
  )
  df %<>% filter(!is.na(y))
  df %<>% group_by(i,j,l,x_ij,site)
  df %<>% dplyr::summarize(
    n = n(),
    y = sum(y)
  )
  df %<>% mutate(
    j = j+1
  )

  # IT model
  {
    model_it <- glmmTMB(
      cbind(y,n-y) ~ factor(j) + x_ij + (1|i) + (1|i:site),
      data = df,
      family = "binomial" # binomial(link="log")
    )
    est_it <- est(model_it, "IT", J=15, len=length(unique(df$l))-1)
  }

  # ETI model
  {
    model_eti <- glmmTMB(
      cbind(y,n-y) ~ factor(j) + factor(l) + (1|i) + (1|i:site),
      data = df,
      family = "binomial"
    )
    est_eti <- est(model_eti, "ETI", J=15, len=length(unique(df$l))-1)
  }

  # ETI (RTE; height) model
  {
    model_rte <- glmmTMB(
      cbind(y,n-y) ~ factor(j) + factor(l) + (x_ij|i),
      data = df,
      family = "binomial"
    )
    est_rte <- est(model_rte, "ETI", J=15, len=length(unique(df$l))-1)
  }

  # RETI model (R=8)
  {
    r <- 8
    df_reti <- df %>% mutate(l = ifelse(l>=r,r,l))
    model_reti <- glmmTMB(
      cbind(y,n-y) ~ factor(j) + factor(l) + (1|i) + (1|i:site),
      data = df_reti,
      family = "binomial"
    )
    est_reti <- est(model_reti, "RETI", J=15, len=length(unique(df$l))-1, R=r)
  }

  # Natural cubic spline (4df)
  {
    J <- 15
    len <- length(unique(df$l))-1
    ns_basis <- ns(c(0:len), knots=c(len/4,(2*len)/4,(3*len)/4))
    df$c1 <- ns_basis[df$l+1,1]
    df$c2 <- ns_basis[df$l+1,2]
    df$c3 <- ns_basis[df$l+1,3]
    df$c4 <- ns_basis[df$l+1,4]
    model_ncs <- glmmTMB(
      cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4 + (1|i) + (1|i:site),
      data = df,
      family = "binomial"
    )
    est_ncs <- est(model_ncs, "NCS", J=15, len=len)
  }

  # Monotone Effect Curve Bayesian model
  {
    # Setup
    mcmc <- cfg$mcmc
    df %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
    df %<>% mutate(
      s_1=as.integer(l>=1), s_2=as.integer(l>=2), s_3=as.integer(l>=3),
      s_4=as.integer(l>=4), s_5=as.integer(l>=5), s_6=as.integer(l>=6),
      s_7=as.integer(l>=7), s_8=as.integer(l>=8), s_9=as.integer(l>=9),
      s_10=as.integer(l>=10), s_11=as.integer(l>=11), s_12=as.integer(l>=12),
      s_13=as.integer(l>=13), s_14=as.integer(l>=14)
    )
    options(mc.cores = parallel::detectCores()-1)
    # Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
    rstan_options(auto_write=TRUE)

    # Put data into Stan format
    stan_data <- list(
      I = length(unique(df$i)),
      N = nrow(df),
      y = df$y,
      bin_n = df$n,
      i = as.integer(df$i),
      j_2 = df$j_2, j_3 = df$j_3, j_4 = df$j_4, j_5 = df$j_5, j_6 = df$j_6,
      j_7 = df$j_7, j_8 = df$j_8, j_9 = df$j_9, j_10 = df$j_10, j_11 = df$j_11,
      j_12 = df$j_12, j_13 = df$j_13, j_14 = df$j_14, j_15 = df$j_15,
      s_1 = df$s_1, s_2 = df$s_2, s_3 = df$s_3, s_4 = df$s_4, s_5 = df$s_5,
      s_6 = df$s_6, s_7 = df$s_7, s_8 = df$s_8, s_9 = df$s_9, s_10 = df$s_10,
      s_11 = df$s_11, s_12 = df$s_12, s_13 = df$s_13, s_14 = df$s_14
    )

    # Stan model
    stan_code <- quote("
      data {
        int I;
        int N;
        int y[N];
        int bin_n[N];
        int i[N];
        real j_2[N]; real j_3[N]; real j_4[N]; real j_5[N]; real j_6[N];
        real j_7[N]; real j_8[N]; real j_9[N]; real j_10[N]; real j_11[N];
        real j_12[N]; real j_13[N]; real j_14[N]; real j_15[N];
        real s_1[N]; real s_2[N]; real s_3[N]; real s_4[N]; real s_5[N];
        real s_6[N]; real s_7[N]; real s_8[N]; real s_9[N]; real s_10[N];
        real s_11[N]; real s_12[N]; real s_13[N]; real s_14[N];
      }
      parameters {
        real beta0;
        real beta_j_2; real beta_j_3; real beta_j_4; real beta_j_5;
        real beta_j_6; real beta_j_7; real beta_j_8; real beta_j_9;
        real beta_j_10; real beta_j_11; real beta_j_12; real beta_j_13;
        real beta_j_14; real beta_j_15;
        real delta;
        simplex[14] smp;
        real alpha[I];
        real<lower=0> sigma;
        real<lower=0> tau;
      }
      transformed parameters {
        real beta_s_1; real beta_s_2; real beta_s_3; real beta_s_4;
        real beta_s_5; real beta_s_6; real beta_s_7; real beta_s_8;
        real beta_s_9; real beta_s_10; real beta_s_11; real beta_s_12;
        real beta_s_13; real beta_s_14;
        beta_s_1 = delta * smp[1]; beta_s_2 = delta * smp[2];
        beta_s_3 = delta * smp[3]; beta_s_4 = delta * smp[4];
        beta_s_5 = delta * smp[5]; beta_s_6 = delta * smp[6];
        beta_s_7 = delta * smp[7]; beta_s_8 = delta * smp[8];
        beta_s_9 = delta * smp[9]; beta_s_10 = delta * smp[10];
        beta_s_11 = delta * smp[11]; beta_s_12 = delta * smp[12];
        beta_s_13 = delta * smp[13]; beta_s_14 = delta * smp[14];
      }
      model {
        delta ~ normal(0,100);
        alpha ~ normal(0,tau);
        smp ~ dirichlet([1,1,1,1,1,1,1,1,1,1,1,1,1,1]');
        for (n in 1:N) {
          y[n] ~ binomial_logit(bin_n[n],
            (beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
            beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] +
            beta_j_7*j_7[n] + beta_j_8*j_8[n] + beta_j_9*j_9[n] +
            beta_j_10*j_10[n] + beta_j_11*j_11[n] + beta_j_12*j_12[n] +
            beta_j_13*j_13[n] + beta_j_14*j_14[n] + beta_j_15*j_15[n] + delta*(
              smp[1]*s_1[n] + smp[2]*s_2[n] + smp[3]*s_3[n] + smp[4]*s_4[n] +
              smp[5]*s_5[n] + smp[6]*s_6[n] + smp[7]*s_7[n] + smp[8]*s_8[n] +
              smp[9]*s_9[n] + smp[10]*s_10[n] + smp[11]*s_11[n] +
              smp[12]*s_12[n] + smp[13]*s_13[n] + smp[14]*s_14[n]
            ) +
            alpha[i[n]])
          );
        }
    }")

    # Fit model
    model_mec <- stan(
      model_code = stan_code,
      data = stan_data,
      chains = mcmc$n.chains,
      iter = mcmc$n.burn+mcmc$n.iter,
      warmup = mcmc$n.burn,
      thin = mcmc$thin
    )
    est_mec <- est(model_mec, "MEC", J=15, len=14)

    # !!!!! MCMC diagnostics
    if (FALSE) {
      vars <- c("beta_j_2", "beta_j_3", "beta_j_4", "beta_j_5", "beta_j_6",
                "beta_j_7", "beta_j_8", "beta_j_9", "beta_j_10", "beta_j_11",
                "beta_j_12", "beta_j_13", "beta_j_14", "beta_j_15")
      # vars <- c("beta_s_1", "beta_s_2", "beta_s_3",
      #           "beta_s_4", "beta_s_5", "beta_s_6")
      # vars <- c("beta0", "delta", "omega")
      # vars <- "smp"
      # vars <- "alpha"
      # vars <- c("sigma", "tau")
      rstan::traceplot(model_mec, vars, inc_warmup=TRUE)
    }

  }

  # Display results
  models <- list("IT"=est_it, "ETI"=est_eti, "RTE"=est_rte,
                 "RETI-8"=est_reti, "NCS-4"=est_ncs, "MEC"=est_mec)
  print_results(models=models, which="wa_state")
  print_graphs(models=models, which="wa_state", ncol=3)

}



#########################################################.
##### MISC: Graphs of delay models (for manuscript) #####
#########################################################.

if (run_misc) {

  # Generate data
  d1 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,0.1), slopes=10))
  d2 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,2,2.1), slopes=c(0,10)))
  d3 <- effect_curve(seq(0,6,0.1), type="exp", params=list(d=1.5))
  d4 <- effect_curve(seq(0,6,0.1), type="spline",
                     params=list(knots=c(0,2,4), slopes=c(0.1,0.4)))
  d5 <- effect_curve(seq(0,6,0.1), type="non-monotonic", params=NULL)

  curve_labels <- c("(a) Instantaneous","(b) Lagged","(c) Curved",
                    "(d) Partially convex","(e) Non-monotonic")

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



##########################################################.
##### MISC: Graphs of WLS estimator (for manuscript) #####
##########################################################.

if (run_misc) {

  grid <- seq(0.01,0.99,0.01)
  weight <- function(Q,s,phi) {
    ( 6*(s-Q-1)*(s+2*phi*Q*s-Q*(1+phi+phi*Q)) ) /
      ( Q*(Q+1)*(phi*Q^2+2*Q-phi*Q-2) )
  }

  # Cases corresponding to a fixed value of S; weights as a function of rho
  # Export: 5" x 3"
  for (Q in c(3,5,8)) {
    assign(paste0("wts_",Q), c())
    assign(paste0("labels_",Q), c())
    for (s in c(1:Q)) {
      assign(paste0("wts_",Q), c(
        eval(as.name(paste0("wts_",Q))),
        sapply(grid, function(phi) { weight(Q=Q, s=s, phi=phi) })
      ))
      assign(paste0("labels_",Q), c(
        eval(as.name(paste0("labels_",Q))),
        paste0("s=", s)
      ))
    }
  }
  len <- length(grid)
  cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999")

  # Export: PDF 8"x3"
  ggplot(
    data.frame(
      x = c(rep(grid,3), rep(grid,5), rep(grid,8)),
      y = c(wts_3, wts_5, wts_8),
      which = c(rep(labels_3, each=len),
                rep(labels_5, each=len),
                rep(labels_8, each=len)),
      Q = c(rep("Q=3",3*len), rep("Q=5",5*len), rep("Q=8",8*len))
    ),
    aes(x=x, y=y, color=which)) +
    geom_line() +
    scale_color_manual(values=cb_colors) +
    facet_wrap(~Q, ncol=3) +
    labs(x=unname(latex2exp::TeX("$\\phi$")), y="Weight", color="")

  # TEMP GRAPH #1: Effect curve
  {

    which <- c()
    Qs <- c()
    x <- c()
    y <- c()
    for (Q in c(3,5,8)) {

      # Add "delay 1" points
      which <- c(which, rep("Delay 1", length(grid)))
      Qs <- c(Qs, rep(paste0("Q=",Q), length(grid)))
      x <- c(x, grid)
      y <- c(y, sapply(grid, function(phi) {
        wts <- sapply(c(1:Q), function(s) { weight(Q=Q, s=s, phi=phi) })
        return(sum(wts*c(0,rep(1,Q-1))))
      }))

      # Add "delay 2" points
      which <- c(which, rep("Delay 2", length(grid)))
      Qs <- c(Qs, rep(paste0("Q=",Q), length(grid)))
      x <- c(x, grid)
      y <- c(y, sapply(grid, function(phi) {
        wts <- sapply(c(1:Q), function(s) { weight(Q=Q, s=s, phi=phi) })
        return(sum(wts*c(0,0,rep(1,Q-2))))
      }))

      # Add "linear" points
      which <- c(which, rep("Linear", length(grid)))
      Qs <- c(Qs, rep(paste0("Q=",Q), length(grid)))
      x <- c(x, grid)
      y <- c(y, sapply(grid, function(phi) {
        wts <- sapply(c(1:Q), function(s) { weight(Q=Q, s=s, phi=phi) })
        return(sum(wts*(seq(0,1,length.out=(Q+1))[2:(Q+1)])))
      }))

    }
    cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")

    # Export: PDF 8"x3"
    ggplot(
      data.frame(which=which, Q=Qs, x=x, y=y),
      aes(x=x, y=y, color=which)) +
      geom_line() +
      scale_color_manual(values=cb_colors) +
      labs(x=unname(latex2exp::TeX("$\\phi$")), y="delta_hat", color="") +
      geom_hline(
        aes(yintercept=y, color=which),
        data = data.frame(
          y = c(c(2/3-0.015,4/5,7/8), c(1/3,3/5-0.015,6/8), c(2/3,3/5,9/16)),
          which = rep(c("Delay 1", "Delay 2", "Linear"), each=3),
          Q = rep(c("Q=3","Q=5","Q=8"), times=3)
        ),
        linetype = "longdash"
      ) +
      facet_wrap(~Q, ncol=3)

  }

  # TEMP GRAPH #2: Random time effect
  # (values from Mathematica)
  {
    rho_w <- c(0.2,0.5,0.8,0.2,0.5,0.8,0.2,0.5,0.8)
    rho_b <- c(0.2,0.2,0.2,0.5,0.5,0.5,0.8,0.8,0.8)
    w1 <- c(1,0.964,0.938,1.167,1.125,1.091,1.25,1.212,1.179)
    w2 <- c(0.167,0.179,0.188,0.111,0.125,0.136,0.083,0.096,0.107)
    w3 <- c(-0.167,-0.143,-0.125,-0.278,-0.25,-0.227,-0.333,-0.308,-0.286)
    p_data2 <- data.frame(
      which = rep(c("s=1","s=2","s=3"),each=9),
      rho_w = paste0("rho_w=",rep(rho_w, 3)),
      rho_b = rep(rho_b, 3),
      weight = c(w1,w2,w3)
    )
    ggplot(
      p_data2,
      aes(x=rho_b, y=weight, color=which)) +
      geom_line() +
      scale_color_manual(values=cb_colors) +
      facet_wrap(~rho_w, ncol=3) +
      labs(x="rho_b", y="Weight", color="")
  }

  # TEMP GRAPH #2: Random Tx effect
  # (values from Mathematica)
  {
    r0 <- rep(c(0.2,0.5,0.8,0.2,0.5,0.8,0.2,0.5,0.8),3)
    r1 <- rep(c(0.2,0.2,0.2,0.5,0.5,0.5,0.8,0.8,0.8),3)
    r2 <- c(rep(0.2,9),rep(0.5,9),rep(0.8,9))
    w1 <- c(
      c(0.923,0.872,0.804,1.088,1.037,0.956,1.267,1.224,1.141),
      c(0.985,0.864,0.736,1.191,1.071,0.934,1.400,1.291,1.154),
      c(1.025,0.547,0.529,1.547,1.000,0.717,1.581,1.541,1.159)
    )
    w2 <- c(
      c(0.192,0.221,0.265,0.125,0.153,0.203,0.041,0.062,0.112),
      c(0.170,0.224,0.284,0.087,0.143,0.212,-.011,0.039,0.108),
      c(0.150,0.377,0.534,-.004,0.500,0.279,0.026,-.061,0.114)
    )
    w3 <- c(
      c(-.115,-.093,-.069,-.213,-.190,-.159,-.308,-.286,-.253),
      c(-.155,-.087,-.020,-.278,-.214,-.147,-.389,-.330,-.262),
      c(-.175,0.075,-.063,-.543,-.500,0.004,-.607,-.480,-.273)
    )
    p_data3 <- data.frame(
      which = rep(c("s=1","s=2","s=3"),each=27),
      r0 = paste0("r0=",rep(r0, 3)),
      r1 = rep(r1, 3),
      r2 = paste0("r2=",rep(r2, 3)),
      weight = c(w1,w2,w3)
    )
    ggplot(
      p_data3,
      aes(x=r1, y=weight, color=which)) +
      geom_line() +
      scale_color_manual(values=cb_colors) +
      facet_grid(rows=vars(r2), cols=vars(r0))
      labs(x="r1", y="Weight", color="")
  }

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
    mu = log(0.1),
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



#################################.
##### MISC: Dirichlet plots #####
#################################.

if (run_misc) {

  sample1 <- rdirichlet(10,c(5,5,5,1,1,1))
  df1 <- data.frame(
    sample = factor(rep(c(1:10), each=6)),
    element = factor(rep(c(1:6), 10)),
    y = c(sample1[1,],sample1[2,],sample1[3,],sample1[4,],sample1[5,],
          sample1[6,],sample1[7,],sample1[8,],sample1[9,],sample1[10,])
  )
  sample2 <- rdirichlet(10,c(500,500,500,100,100,100))
  df2 <- data.frame(
    sample = factor(rep(c(1:10), each=6)),
    element = factor(rep(c(1:6), 10)),
    y = c(sample2[1,],sample2[2,],sample2[3,],sample2[4,],sample2[5,],
          sample2[6,],sample2[7,],sample2[8,],sample2[9,],sample2[10,])
  )
  sample3 <- rdirichlet(10,c(0.05,0.05,0.05,0.01,0.01,0.01))
  df3 <- data.frame(
    sample = factor(rep(c(1:10), each=6)),
    element = factor(rep(c(1:6), 10)),
    y = c(sample3[1,],sample3[2,],sample3[3,],sample3[4,],sample3[5,],
          sample3[6,],sample3[7,],sample3[8,],sample3[9,],sample3[10,])
  )

  ggplot(df1, aes(fill=element, y=y, x=sample)) +
    geom_bar(position="stack", stat="identity") +
    labs(title="Dirichlet(5,5,5,1,1,1)")
  ggplot(df2, aes(fill=element, y=y, x=sample)) +
    geom_bar(position="stack", stat="identity") +
    labs(title="Dirichlet(500,500,500,100,100,100)")
  ggplot(df3, aes(fill=element, y=y, x=sample)) +
    geom_bar(position="stack", stat="identity") +
    labs(title="Dirichlet(0.05,0.05,0.05,0.01,0.01,0.01)")

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
    packages = c("dplyr", "magrittr", "stringr", "lme4", "rjags", # "geepack", "restriktor", "scam", "gamlss", "glmmTMB"
                 "mgcv", "MASS", "fastDummies", "scales", "car")
  )
  sim %<>% add_constants(
    mu = log(0.1)
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
    theta = 0.5,
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
  s <- sim %>% summarize()
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
