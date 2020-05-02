#' Generate one stepped wedge dataset
#'
#' @param alpha Log baseline prevalence
#' @param tau Cluster random effect SD
#' @param theta Treatment effect
#' @param d Treatment lag rate parameter
#' @param n_clusters Total number of clusters
#' @param n_time_points Total number of time points ("J")
#' @param n_ind_per_cluster Total number of individuals per cluster (assumes
#' @param type Type of outcome; options include "binomial" or "normal"
#' @param sigma Standard deviation of the outcome (if type="normal"; omit
#'     otherwise)
#'     equal cluster sizes)
#' @return A list containing the following: \cr
#'     * `params`: a list of the parameters supplied in the function call \cr
#'     * `data`: the resulting data frame
#' @export
generate_dataset <- function(alpha, tau, theta, d, n_clusters, n_time_points,
                             n_ind_per_cluster, type, sigma=NA) {

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
  if (type == "binomial") {
    data %<>% rename('p_ij' = `y_ij`)
  }

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

  # Create theta_ls (intervention effects) based on continuous function
  theta_ls <- sapply(1:(n_time_points-1), function(l){
    theta*(1-exp(-l/d))
  })

  # Loop through clusters, time, and individuals
  for (i in 1:n_clusters) {

    v_i <- rnorm(1, mean=0, sd=tau)
    c_i <- crossover_times[i]-1

    for (j in 1:n_time_points) {

      x_ij <- ifelse(j<crossover_times[i], 0, 1)
      l <- ifelse(j<crossover_times[i], 0, (j-crossover_times[i])+1)

      # !!!!! this currently excludes the MVN error term
      # !!!!! Check theta_l against estimates
      x_il <- ifelse(l>0, 1, 0)
      theta_l <- ifelse(l>0, theta_ls[l], 0)

      # !!!!! This is exp rather than expit; make sure it doesn't throw error
      if (type=="normal") {
        y_ij <- alpha + beta_js[j] + theta_l*x_il + v_i
      } else if (type=="binomial") {
        p_ij <- exp(alpha + beta_js[j] + theta_l*x_il + v_i)
      } else {
        stop ("`type` must be either 'normal' or 'binomial'")
      }

      k <- n_ind_per_cluster
      if (type=="normal") {
        y <- y_ij + rnorm(k, mean=0, sd=sigma)
        # !!!!! This table is too bulky
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=rep(k,k), l=rep(l,k), v_i=rep(v_i,k),
          y_ij=rep(y_ij,k), x_ij=rep(x_ij,k), c_i=c_i, y=y
        )))
      } else if (type=="binomial") {
        y <- rbinom(n=k, size=1, prob=p_ij)
        # !!!!! This table is too bulky
        data <- rbind(data, data.frame(cbind(
          i=rep(i,k), j=rep(j,k), k=rep(k,k), l=rep(l,k), v_i=rep(v_i,k),
          p_ij=rep(p_ij,k), x_ij=rep(x_ij,k), c_i=c_i, y=y
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
