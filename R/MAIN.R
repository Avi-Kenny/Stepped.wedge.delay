# Title: "Stepped wedge lagged effect simulation"
# Author: Avi Kenny



# install.packages(
#   pkgs = "fastDummies",
#   lib = "/home/students/avikenny/Desktop/R_lib", # UW - Bayes
#   repos = "http://cran.us.r-project.org",
#   dependencies = TRUE
# )



##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  level_set_which = "level_set_1",
  run_or_update = "run",
  num_sim = 1000,
  pkgs = c("dplyr", "stringr", "lme4", "rjags", "Iso", "sqldf", "mgcv", "MASS",
           "fastDummies", "car"),
  pkgs_nocluster = c("ggplot2", "viridis", "scales", "facetscales"), # devtools::install_github("zeehio/facetscales")
  parallel = "none",
  stop_at_error = FALSE
)

# Set cluster config
cluster_config <- list(
  # js = "sge",
  # dir = "/home/students/avikenny/Desktop/z.stepped.wedge"
  js = "slurm",
  dir = "/home/akenny/z.stepped.wedge"
)



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
library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
source("generate_dataset.R")
source("one_simulation.R")
source("plot_outcome.R")
source("plot_sw_design.R")
source("run_analysis.R")
source("effect_curve.R")

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
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 30,
    theta = 0.5,
    tau = 1,
    mu = 1,
    data_type = "binomial",
    sigma = 2,
    delay_model = list(type="spline", params=list(knots=c(0,1),slopes=1)),
    # delay_model = list(type="exp", params=list(d=1.5)),
    n_extra_time_points = 0,
    rte = NA
    # rte = list(type="height", rho=-0.2, nu=0.4)
    # rte = list(type="height+shape", rho1=-0.1, rho2=0.6, nu=0.4)
  )

  # # Set variables needed in run_analysis
  # J <- data$params$n_time_points
  # data_type <- "normal"

  # !!!!!

  df <- data$data # !!!!!
  df %<>% group_by(i,j,l,x_ij)
  df %<>% summarize(n=n(), y=sum(y))

  # ETI
  model_tmb <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + factor(l) + (1|i),
    data = df,
    family = "binomial"
  )
  coeffs_eti <- as.numeric(summary(model_tmb)$coefficients$cond[,1][8:13])

  # Manual cubic polynomial
  df %<>% mutate(
    c1 = l,
    c2 = l^2,
    c3 = l^3,
    c4 = pmax(0,(l-3)^3),
    # c5 = pmax(0,(l-4)^3)
  )
  model_tmb <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4 + (1|i),
    # cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4+c5 + (1|i),
    data = df,
    family = "binomial"
  )
  mb <- as.numeric(summary(model_tmb)$coefficients$cond[,1][8:12])
  coeffs_cub <- sapply(c(1:6), function(l) {
    mb[1]*l + mb[2]*l^2 + mb[3]*l^3 + mb[4]*(max(0,(l-3)^3))
    # + mb[4]*(max(0,(l-2)^3)) + mb[5]*(max(0,(l-4)^3))
  })

  # Smoothing spline
  model_ss <- gamm(
    cbind(y,n-y) ~ factor(j) + s(l, k=7, fx=FALSE, bs="cr", m=c(3,2), pc=0),
    random = list(i=~1),
    data = df,
    family = "binomial"
  )
  coeffs_ss <- sapply(c(1:(n_knots-1)), function(l) {
    predict(model_ss$gam, newdata=list(j=1, l=l), type = "terms")[2]
  })

  # Regression spline
  model_rs <- gamm(
    cbind(y,n-y) ~ factor(j) + s(l, k=4, fx=TRUE, bs="cr", pc=0),
    random = list(i=~1),
    data = df,
    family = "binomial"
  )
  coeffs_rs <- sapply(c(1:(n_knots-1)), function(l) {
    predict(model_rs$gam, newdata=list(j=1, l=l), type = "terms")[2]
  })

  # Plot
  ggplot(
    data.frame(
      x = rep(c(0:6),4),
      y = c(c(0,coeffs_eti), c(0,coeffs_ss), c(0,coeffs_rs), c(0,coeffs_cub)),
      which = rep(c("ETI","SS","RS","CUBE"), each=7)
    ),
    aes(x=x, y=y, color=which)
  ) +
    geom_line() +
    labs(y="log(OR)", x="Exposure time", color="Model")
  #



}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (run_main) {
  if (Sys.getenv("run") %in% c("first", "")) {

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

    # !!!!! Method archive
    # "MCMC (exp; N(1,10) mix 0.2)" = list(method="MCMC-STEP-MON", enforce="exp; N(1,10) mix (0.2)",mcmc=list(n.adapt=1000, n.burn=2000, n.iter=1000, n.chains=2)),
    # "PAVA (wts: equal)" = list(method="MCMC-STEP-PAVA", wts="equal",mcmc=list(n.adapt=1000, n.burn=1000, n.iter=2000, n.chains=2)),
    # "PAVA (wts: samp_size)" = list(method="MCMC-STEP-PAVA", wts="samp_size",mcmc=list(n.adapt=1000, n.burn=1000, n.iter=2000, n.chains=2)),
    # "PAVA (wts: sqrt_samp_size)" = list(method="MCMC-STEP-PAVA", wts="sqrt_samp_size",mcmc=list(n.adapt=1000, n.burn=1000, n.iter=2000, n.chains=2))
    # "MCMC-STEP (exp; mix prior 0.2)" = list(method="MCMC-STEP-MON",enforce="exp; mix prior 0.2")

    # Simulation 1: dangers of "immediate treatment" model
    # 12 level combos
    level_set_1 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = 0.5,
      tau = 1,
      sigma = 2,
      data_type = "normal",
      method = list(
        # "HH" = list(method="HH"),
        "ETI" = list(method="ETI"),
        "CUBIC-3df" = list(method="CUBIC-3df"), # !!!!!
        "CUBIC-4df" = list(method="CUBIC-4df"), # !!!!!
        "CUBIC-5df" = list(method="CUBIC-5df"), # !!!!!
        "SS" = list(method="SS")
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("none"=list())
    )

    # Simulation 2: power of CI-based hypothesis tests
    # 72 level combos
    level_set_2 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = seq(0,0.5,0.1),
      tau = 1,
      sigma = 2,
      data_type = "normal",
      method = list(
        "HH" = list(method="HH"),
        "ETI" = list(method="ETI"),
        "SS" = list(method="SS")
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("none"=list())
    )

    # Simulation 3: effect_reached
    # 12 level combos
    level_set_3 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = 0.5,
      tau = 1,
      sigma = 2,
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

    # Simulation 4: n_extra_time_points
    # 12 level combos
    level_set_4 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 50,
      theta = 0.5,
      tau = 1,
      sigma = 2,
      data_type = "normal",
      method = list("ETI" = list(method="ETI")),
      delay_model = delay_models,
      n_extra_time_points = c(0,1,2),
      rte = NA,
      return_extra = list("none"=list())
    )

    # Simulation 5: random treatment effects
    level_set_5 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20,
      theta = 0.5,
      tau = 0.5,
      sigma = 0.2,
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI"),
        # "ETI (RTE; height)" = list(method="ETI", re="height"),
        "ETI (RTE MCMC; height)"=list(
          method = "MCMC-RTE-height",
          mcmc = list(n.adapt=2000, n.iter=2000, n.burn=2000, n.chains=3)),
        "ETI (RTE MCMC; height+shape)" = list(
          method = "MCMC-RTE-height+shape",
          mcmc = list(n.adapt=2000, n.burn=2000))
      ),
      delay_model = delay_models,
      # delay_model = list("Curved"=list(type="exp",params=list(d=1.5))),
      n_extra_time_points = 0,
      rte = list(
        # "none" = NA,
        # "height" = list(type="height", nu=0.4, rho=-0.2),
        "height+shape" = list(type="height+shape", nu=0.4, rho1=-0.2, rho2=0.5)
      ),
      return_extra = list("rte"=list(rte=TRUE))
    )

    # Simulation 6: Monotone Effect Curve (MEC) model
    # 12 level combos
    level_set_6 <- list(
      n_clusters = 24,
      n_time_points = 7,
      n_ind_per_cluster = 20, # !!!!! 50
      theta = 0.5,
      tau = 1,
      sigma = 0.2, # !!!!! 2
      data_type = "normal",
      method = list(
        "ETI" = list(method="ETI"),
        "SS" = list(method="SS"),
        "MEC (0.1 mix)" = list(
          method = "MCMC-STEP-MON", enforce="exp; mix prior 0.1",
          mcmc = list(n.adapt=2000, n.iter=2000, n.burn=2000, n.chains=3)),
        "MEC (0.2 mix)" = list(
          method = "MCMC-STEP-MON", enforce="exp; mix prior 0.2",
          mcmc = list(n.adapt=2000, n.iter=2000, n.burn=2000, n.chains=3))
      ),
      delay_model = delay_models,
      n_extra_time_points = 0,
      rte = NA,
      return_extra = list("whole_curve"=list(whole_curve=TRUE))
    )

    # # Simulation 7: !!!!!
    # level_set_7 <- list(...)

    level_set <- eval(as.name(cfg$level_set_which))

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

  if (cfg$run_or_update=="run") {

    run_on_cluster(

      first = {

        # Set up and configure simba object
        sim <- new_sim()
        sim %<>% set_config(
          num_sim = cfg$num_sim,
          parallel = cfg$parallel,
          stop_at_error = cfg$stop_at_error,
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
        # print(paste("Check 1:",Sys.time())) # !!!!!
        sim %<>% run()
        # print(paste("Check 2:",Sys.time())) # !!!!!
      },

      last = {
        sim %>% summary() %>% print()
        # sim$results %>% print()
        sim$errors %>% print()
      },

      cluster_config = cluster_config

    )

  }

  if (cfg$run_or_update=="update") {

    update_on_cluster(

      first = {
        sim <- readRDS(paste0(cluster_config$dir,'/sim.simba'))
        sim <- do.call(set_levels, c(list(sim), level_set))
      },

      main = {
        sim %<>% update()
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
  whichsim <- 1

  # Read in simulation object
  sim <- readRDS("../simba.out/sim_newSmooths_20210426.simba")

  # Generate true TATE values
  sim$results %<>% mutate(
    ate = case_when(
      delay_model=="Instantaneous" ~ theta,
      delay_model=="Lagged" ~ theta * mean(
        effect_curve(x = seq(0.1,6,0.1),
                     type = "spline",
                     params = list(knots=c(0,2,2.1),slopes=c(0,10)))),
      delay_model=="Curved" ~ theta * mean(
        effect_curve(x = seq(0.1,6,0.1),
                     type = "exp",
                     params = list(d=1.5))),
      delay_model=="Partially convex" ~ theta * mean(
        effect_curve(x = seq(0.1,6,0.1),
                     type = "spline",
                     params = list(knots=c(0,2,4),slopes=c(0.1,0.4))))
    )
  )

  # Generate true theta values
  if (whichsim==6) {
    theta <- sim$results[1,"theta"]
    thetas <- list(
      "I" = rep(theta,6),
      "L" = theta * effect_curve(x = seq(1,6,1),
                                 type = "spline",
                                 params = list(knots=c(0,2,2.1),slopes=c(0,10))),
      "C" = theta * effect_curve(x = seq(1,6,1),
                                 type = "exp",
                                 params = list(d=1.5)),
      "P" = theta * effect_curve(x = seq(1,6,1),
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
      # pbias_wc = (1/6) * ((theta_1_hat-theta_1)/theta_1+(theta_2_hat-theta_2)/theta_2+(theta_3_hat-theta_3)/theta_3+(theta_4_hat-theta_4)/theta_4+(theta_5_hat-theta_5)/theta_5+(theta_6_hat-theta_6)/theta_6)
    )

    summ_mean <- list(
      list(name="ate", x="ate"),
      list(name="mean_ate", x="ate_hat"),
      list(name="mse_wc", x="mse_wc"),
      list(name="mean_lte", x="lte_hat")
    )

  } else {
    summ_mean <- list(
      list(name="ate", x="ate"),
      list(name="mean_ate", x="ate_hat"),
      list(name="mean_lte", x="lte_hat")
    )
  }

  # Summarize data
  summ <- summary(
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
  if (whichsim %in% c(1,2)) {
    summ %<>% mutate(
      method = ifelse(method=="HH", "IT", method)
    )
  }
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
  if (whichsim==6) {
    summ %<>% filter(method!="MEC (0.0 mix)")
    # summ %<>% filter(method!="SS") # !!!!!
  }
  if (whichsim %in% c(1,2)) {
    s_methods <- c("IT", "ETI", "SS")
  } else if (whichsim==3) {
    s_methods <- c("ETI", "RETI (3 steps)", "RETI (4 steps)")
  } else if (whichsim==6) {
    s_methods <- c("ETI", "SS", "MEC (0.1 mix)", "MEC (0.2 mix)")
  }
  summ %<>% mutate(
    method = factor(method, levels=s_methods),
    delay_model = factor(delay_model, levels=s_d_models),
    power_ate = 1 - beta_ate,
    power_lte = 1 - beta_lte
  )
  summ %<>% rename("Method"=method)
  p_data <- sqldf("
    SELECT Method, delay_model, 'TATE' AS which, bias_ate AS bias, theta,
    cov_ate AS Coverage, power_ate AS Power, mse_ate AS MSE FROM summ
    UNION SELECT Method, delay_model, 'LTE', bias_lte, theta,
    cov_lte, power_lte, mse_lte FROM summ
  ")
  p_data <- sqldf("
    SELECT Method, delay_model, which, 'Bias' AS stat, Bias AS value, theta FROM p_data
    UNION SELECT Method, delay_model, which, 'Coverage', Coverage, theta FROM p_data
    UNION SELECT Method, delay_model, which, 'Power', Power, theta FROM p_data
    UNION SELECT Method, delay_model, which, 'MSE', MSE, theta FROM p_data
  ")
  p_data %<>% mutate(which = factor(which, levels=c("TATE", "LTE")))

  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    IT = cb_colors[2],
    ETI = cb_colors[4],
    SS = cb_colors[3],
    `CUBIC-3df` = cb_colors[6],
    `CUBIC-4df` = cb_colors[7],
    `CUBIC-5df` = cb_colors[8],
    `MEC (0.1 mix)` = cb_colors[6],
    `MEC (0.2 mix)` = cb_colors[7],
    `RETI (3 steps)` = cb_colors[6],
    `RETI (4 steps)` = cb_colors[7]
  )
  # viridis(5)

}



###########################################.
##### VIZ: Figure for simulations 1+6 #####
###########################################.

if (run_viz) {

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
    geom_bar(stat="identity", position=position_dodge(),
             width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=vars(delay_model), rows=vars(stat), scales=list(y=list(
      Bias = scale_y_continuous(labels = percent_format()),
      Coverage = scale_y_continuous(labels = percent_format()),
      MSE = scale_y_continuous()
    ))) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors) +
    labs(y=NULL, x=NULL, fill="Analysis model")

}



########################################.
##### VIZ: Figure for simulation 2 #####
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
    scale_color_manual(values=m_colors) +
    scale_shape_manual(values=c(15,16,17)) +
    facet_grid(rows=vars(which), cols=vars(delay_model)) +
    scale_y_continuous(labels=percent_format()) +
    theme(
      axis.text.x = element_text(angle=90, hjust=0, vjust=0.4)
    )

}



##############################################.
##### VIZ: Figure (2nd) for simulation 6 #####
##############################################.

if (run_viz) {

  # Export: 8: x 4"
  ggplot(
    summ,
    aes(x=Method, y=mse_wc, fill=Method) # x=which,
  ) +
    geom_bar(stat="identity", position=position_dodge(),
             width=0.8, color="white", size=0.35) +
    facet_grid(cols=vars(delay_model)) +
    theme(
      legend.position="bottom",
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    ) +
    scale_fill_manual(values=m_colors) +
    labs(y="MSE", x=NULL, fill="Analysis model")

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



###############################################.
##### MAIN: Real data analysis (WA state) #####
###############################################.

if (run_realdata) {

  # !!!!! Also plot time trend estimates (with CI)

  # Read/process data
  df <- read.csv("../realdata/wa_state.csv")
  df %<>% rename(
    "y" = ct,
    "i" = hdist,
    "j" = timex,
    "l" = timeontrtx,
    "x_ij" = treatx
  )
  df %<>% filter(!is.na(y))
  df %<>% group_by(i,j,l,x_ij)
  df %<>% summarize(
    n = n(),
    y = sum(y)
  )

  # HH model
  model_hh <- glmer(
    cbind(y,n-y) ~ factor(j) + x_ij + (1|i),
    data = df,
    family = "binomial"
  )
  coeff_hh <- summary(model_hh)$coefficients["x_ij","Estimate"]

  # ETI model
  model_eti <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + factor(l) + (1|i),
    data = df,
    family = "binomial"
  )
  coeffs_eti <- as.numeric(summary(model_eti)$coefficients$cond[,1][16:29])

  # Cubic spline (3df)
  df %<>% mutate(
    c1 = l,
    c2 = l^2,
    c3 = l^3
  )
  model_cube3 <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + c1+c2+c3 + (1|i),
    data = df,
    family = "binomial"
  )
  mb <- as.numeric(summary(model_cube3)$coefficients$cond[,1][16:18])
  coeffs_cub3 <- sapply(c(1:14), function(l) {
    mb[1]*l + mb[2]*l^2 + mb[3]*l^3
  })

  # Cubic spline (7df)
  df %<>% mutate(
    c4 = pmax(0,(l-1*(14/5))^3),
    c5 = pmax(0,(l-2*(14/5))^3),
    c6 = pmax(0,(l-3*(14/5))^3),
    c7 = pmax(0,(l-4*(14/5))^3)
  )
  model_cube7 <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4+c5+c6+c7 + (1|i),
    data = df,
    family = "binomial"
  )
  mb <- as.numeric(summary(model_cube7)$coefficients$cond[,1][16:22])
  coeffs_cub7 <- sapply(c(1:14), function(l) {
    mb[1]*l + mb[2]*l^2 + mb[3]*l^3 + mb[4]*pmax(0,(l-1*(14/5))^3) +
      mb[5]*pmax(0,(l-2*(14/5))^3) + mb[6]*pmax(0,(l-3*(14/5))^3) +
      mb[7]*pmax(0,(l-4*(14/5))^3)
  })

  # Cubic spline (8df)
  df %<>% mutate(
    c4 = pmax(0,(l-1*(14/6))^3),
    c5 = pmax(0,(l-2*(14/6))^3),
    c6 = pmax(0,(l-3*(14/6))^3),
    c7 = pmax(0,(l-4*(14/6))^3),
    c8 = pmax(0,(l-5*(14/6))^3)
  )
  model_cube8 <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4+c5+c6+c7+c8 + (1|i),
    data = df,
    family = "binomial"
  )
  mb <- as.numeric(summary(model_cube8)$coefficients$cond[,1][16:23])
  coeffs_cub8 <- sapply(c(1:14), function(l) {
    mb[1]*l + mb[2]*l^2 + mb[3]*l^3 + mb[4]*pmax(0,(l-1*(14/6))^3) +
      mb[5]*pmax(0,(l-2*(14/6))^3) + mb[6]*pmax(0,(l-3*(14/6))^3) +
      mb[7]*pmax(0,(l-4*(14/6))^3) + mb[8]*pmax(0,(l-5*(14/6))^3)
  })

  # Cubic spline (9df)
  df %<>% mutate(
    c4 = pmax(0,(l-1*(14/7))^3),
    c5 = pmax(0,(l-2*(14/7))^3),
    c6 = pmax(0,(l-3*(14/7))^3),
    c7 = pmax(0,(l-4*(14/7))^3),
    c8 = pmax(0,(l-5*(14/7))^3),
    c9 = pmax(0,(l-6*(14/7))^3)
  )
  model_cube9 <- glmmTMB(
    cbind(y,n-y) ~ factor(j) + c1+c2+c3+c4+c5+c6+c7+c8+c9 + (1|i),
    data = df,
    family = "binomial"
  )
  mb <- as.numeric(summary(model_cube9)$coefficients$cond[,1][16:24])
  coeffs_cub9 <- sapply(c(1:14), function(l) {
    mb[1]*l + mb[2]*l^2 + mb[3]*l^3 + mb[4]*pmax(0,(l-1*(14/7))^3) +
      mb[5]*pmax(0,(l-2*(14/7))^3) + mb[6]*pmax(0,(l-3*(14/7))^3) +
      mb[7]*pmax(0,(l-4*(14/7))^3) + mb[8]*pmax(0,(l-5*(14/7))^3) +
      mb[9]*pmax(0,(l-6*(14/7))^3)
  })

  # Smoothing spline
  n_knots <- length(unique(df$l))
  model_ss <- gamm(
    cbind(y,n-y) ~ factor(j) + s(l, k=n_knots, fx=FALSE, bs="cr", m=c(3,2), pc=0),
    random = list(i=~1),
    data = df,
    family = "binomial"
  )
  coeffs_ss <- sapply(c(0:(n_knots-1)), function(l) {
    predict(model_ss$gam, newdata=list(j=1, l=l), type="terms")[2]
  })

  # Plot estimates
  ggplot(
    data.frame(
      x = rep(c(0:14),5),
      y = c(c(0,coeffs_eti), c(0,coeffs_cub3), c(0,coeffs_cub7),
            c(0,coeffs_cub9), coeffs_ss),
      which = rep(c("ETI","CUB3","CUB7","CUB9","SS"), each=15)
    ),
    aes(x=x, y=y, color=which)
  ) +
    geom_line() +
    geom_hline(yintercept=coeff_hh, linetype="longdash") +
    labs(y="log(OR)", x="Exposure time", color="Model")

  # !!!!! Testing smoothing spline
  n <- 20
  x <- seq(0,8, length.out=n)
  y <- sin(x) + rnorm(n, sd=1)
  model_ss <- gam(y ~ s(x, k=n, fx=FALSE, bs="cr", m=c(3,2)))
  ss <- as.numeric(predict(model_ss, newdata=list(x=x), type="terms"))
  ss <- ss + as.numeric(summary(model_ss)$p.coeff)
  model_ss2 <- smooth.spline(x=x, y=y, cv=F, all.knots=T, keep.data=F)
  ss2 <- list(x=predict(model_ss2)$x, y=predict(model_ss2)$y)
  ggplot(
    data.frame(
      x = c(x, x, ss2$x),
      y = c(y, ss, ss2$y),
      which = rep(c("Data","SS","SS2"), each=n)
    ),
    aes(x=x, y=y, color=which)
  ) +
    geom_line()
  #

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
