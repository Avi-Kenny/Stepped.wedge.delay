#' Run a single simulation
#'
#' @return The list created by `run_analysis`

one_simulation <- function() {

  # Generate dataset
  data <- generate_dataset(
    mu = C$mu,
    tau = L$tau,
    theta = L$theta,
    n_clusters = L$n_clusters,
    n_time_points = L$n_time_points,
    n_ind_per_cluster = L$n_ind_per_cluster,
    data_type = L$data_type,
    sigma = L$sigma,
    delay_model = L$delay_model,
    n_extra_time_points = L$n_extra_time_points,
    rte = L$rte
  )

  results <- run_analysis(
    data = data,
    data_type = L$data_type,
    method = L$method
  )

  return (results)

}
