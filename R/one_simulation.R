#' Run a single simulation
#'
#' @param L Passed via simba; list of simulation levels
#' @param C Passed via simba; list of simulation constants
#' @return The list created by `run_analysis`
#' @export
one_simulation <- function(L,C) {

  # Generate dataset
  data <- generate_dataset(
    alpha = C$alpha,
    tau = L$tau,
    theta = L$theta,
    d = L$d,
    n_clusters = L$n_clusters,
    n_time_points = L$n_time_points,
    n_ind_per_cluster = L$n_ind_per_cluster,
    data_type = L$data_type,
    sigma = L$sigma
    # alpha = log(0.1),
    # tau = 0, # 0.25
    # theta = log(0.5),
    # d = 1.4,
    # n_clusters = 12,
    # n_time_points = 7,
    # n_ind_per_cluster = 20,
    # data_type = "binomial",
    # sigma = 3
  )

  results <- run_analysis(
    data = data,
    analysis_type = L$analysis_type,
    data_type = L$data_type,
    L = L,
    C = C
  )

  return (results)

}
