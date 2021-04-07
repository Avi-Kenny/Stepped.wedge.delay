#' Generate one stepped wedge dataset
#'
#' @param mu Log baseline prevalence
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
#' @param rte Specification of random treatment effects. Either NA or a list
#'     containing `type`, `rho`, and `nu`. `type` can be either "height"
#'     (allowing the height of the effect curve to vary) or "height+shape"
#'     (allowing both the height and the shape of the effect curve to vary)).
#'     `rho` and `nu` correspond to the covariance and RTE sd; see the
#'     manuscript for details.
#' @return A list containing the following: \cr
#'     * `params`: a list of the parameters supplied in the function call \cr
#'     * `data`: the resulting data frame

generate_dataset <- function(mu, tau, theta, n_clusters, n_time_points,
                             n_ind_per_cluster, data_type, sigma=NA,
                             delay_model, n_extra_time_points, rte=NA) {

  # Generate data frame
  data <- data.frame(
    "i" = integer(), # cluster
    "j" = integer(), # time; 1=baseline, J=endline
    "k" = integer(), # individual
    "l" = integer(), # time since intervention
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

    if((is.na(rte))[[1]]) {
      alpha_i <- rnorm(1, mean=0, sd=tau)
      eta_i <- 0
    } else {
      if(rte$type=="height") {
        Sigma <- rbind(
          c(tau^2,rte$rho),
          c(rte$rho,rte$nu^2)
        )
        re <- mvrnorm(n=1, mu=c(0,0), Sigma=Sigma)
        alpha_i <- re[1]
        eta_i <- re[2]
      }
      if(rte$type=="height+shape") {
        J <- n_time_points
        Sigma <- rbind(
          c(tau^2,rep(rte$rho,J-1)),
          cbind(
            rep(rte$rho,J-1),
            diag(rep(rte$nu^2,J-1))
          )
        )
        re <- mvrnorm(n=1, mu=rep(0,J), Sigma=Sigma)
        alpha_i <- re[1]
        eta_it <- re[2:(J-1)]
      }
    }

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

      if ((is.na(rte))[[1]] || rte$type=="height") {
        mu_ij <- mu + beta_js[j] + (theta_l+eta_i)*x_ij + alpha_i
      } else if (rte$type=="height+shape") {
        if (l==0) {
          mu_ij <- mu + beta_js[j] + theta_l*x_ij + alpha_i
        } else {
          mu_ij <- mu + beta_js[j] + (theta_l+eta_it[l])*x_ij + alpha_i
        }
      }

      if (data_type=="binomial") {
        expit <- function(x) {1/(exp(-x)+1)}
        mu_ij <- expit(mu_ij)
      }

      k <- n_ind_per_cluster
      if (data_type=="normal") {
        y <- mu_ij + rnorm(k, mean=0, sd=sigma)
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=c(1:k), l=rep(l,k), # alpha_i=rep(alpha_i,k)
          x_ij=rep(x_ij,k), c_i=rep(c_i,k), y=y # mu_ij=rep(mu_ij,k)
        )))
      } else if (data_type=="binomial") {
        y <- rbinom(n=k, size=1, prob=mu_ij)
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=c(1:k), l=rep(l,k), # alpha_i=rep(alpha_i,k)
          x_ij=rep(x_ij,k), c_i=rep(c_i,k), y=y # mu_ij=rep(mu_ij,k)
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
