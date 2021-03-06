#' Generate one stepped wedge dataset
#'
#' @param alpha Log baseline prevalence
#' @param tau Cluster random effect SD
#' @param theta Treatment effect
#' @param n_clusters Total number of clusters
#' @param n_time_points Total number of time points ("J")
#' @param n_ind_per_cluster Total number of individuals per cluster (assumes
#'     equal cluster sizes)
#' @param data_type Type of outcome; options include "binomial" or "normal"
#' @param sigma Standard deviation of the outcome (if data_type="normal"; omit
#'     otherwise)
#' @param delay_model A list containing `type` and `params`, that will be passed
#'     to effect_curve()
#' @param n_extra_time_points Number of extra time points at the end of the
#'     study (all clusters are in the treatment state)
#' @return A list containing the following: \cr
#'     * `params`: a list of the parameters supplied in the function call \cr
#'     * `data`: the resulting data frame

generate_dataset <- function(alpha, tau, theta, n_clusters, n_time_points,
                             n_ind_per_cluster, data_type, sigma=NA,
                             delay_model, n_extra_time_points) {

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
  if (n_clust_per_time %% 1 != 0) {
    stop("n_clusters must be divisible by n_time_points-1")
  }
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

    for (j in 1:(n_time_points+n_extra_time_points)) {

      x_ij <- ifelse(j<crossover_times[i], 0, 1)
      l <- ifelse(j<crossover_times[i], 0, (j-crossover_times[i])+1)

      if (l<=length(theta_ls)) {
        theta_l <- ifelse(l>0, theta_ls[l], 0)
      } else {
        # This will only be reached for j>n_time_points, i.e. when
        #     n_extra_time_points>0
        theta_l <- theta_ls[length(theta_ls)]
      }

      if (data_type=="normal") {
        y_ij <- alpha + beta_js[j] + theta_l*x_ij + v_i
      } else if (data_type=="binomial") {
        expit <- function(x) {1/(exp(-x)+1)}
        y_ij <- expit(alpha + beta_js[j] + theta_l*x_ij + v_i)
      } else {
        stop ("`data_type` must be either 'normal' or 'binomial'")
      }

      k <- n_ind_per_cluster
      if (data_type=="normal") {
        y <- y_ij + rnorm(k, mean=0, sd=sigma)
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=c(1:k), l=rep(l,k), v_i=rep(v_i,k),
          y_ij=rep(y_ij,k), x_ij=rep(x_ij,k), c_i=c_i, y=y
        )))
      } else if (data_type=="binomial") {
        # if (y_ij>1) {
        #   warning(paste0("Probability y_ij=",y_ij,", so y_ij was set to 1"))
        #   y_ij <- 1
        # }
        y <- rbinom(n=k, size=1, prob=y_ij)
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=c(1:k), l=rep(l,k), v_i=rep(v_i,k),
          y_ij=rep(y_ij,k), x_ij=rep(x_ij,k), c_i=c_i, y=y
        )))
      }

    }

  }

  params <- as.list(match.call())
  params$crossover_times <- crossover_times

  return (list(
    "params" = params,
    "beta_js" = beta_js,
    "theta_ls" = theta_ls,
    "data" = data
  ))

}
