
# Setup
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") { setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jim Hughes/Project - Stepped wedge lag/z.stepped.wedge/R") } else { setwd("z.stepped.wedge/R") }
library(dplyr); library(magrittr); library(ggplot2); library(lme4); library(geepack); library(car); library(stringr); library(simba); library(parallel); library(glmmTMB); library(restriktor); library(mgcv); library(scam)
source("generate_dataset.R"); source("log_lik_spline.R"); source("one_simulation.R"); source("plot_outcome.R"); source("plot_sw_design.R"); source("run_analysis.R"); source("effect_curve.R")

# Generate dataset for testing
data <- generate_dataset(
  alpha = log(0.1),
  tau = 0.01, # 1
  theta = log(0.5),
  n_clusters = 12,
  n_time_points = 7,
  n_ind_per_cluster = 20,
  data_type = "normal",
  sigma = 0.01,
  delay_model = list(type="parabola", params=list(a=(-1/16), b=(1/2), c=0))
  # delay_model = list(type="exp", params=list(d=1.4))
)
J <- data$params$n_time_points

# Run model that is not shape-constrained (no random effect)
model_1 <- gam(
  y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
  data = data$data
)
plot(model_1)

model_2 <- scam(
  y ~ factor(j) + s(l, k=J, m=2, bs="mpd"),
  data = data$data
)
plot(model_2)

theta_hats <- sapply(c(1:(J-1)), function(l) {
  predict(model$gam, newdata=list(j=1, l=l), type = "terms")[2]
})

