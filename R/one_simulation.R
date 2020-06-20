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
    n_clusters = L$n_clusters,
    n_time_points = L$n_time_points,
    n_ind_per_cluster = L$n_ind_per_cluster,
    data_type = L$data_type,
    sigma = L$sigma,
    delay_model = L$delay_model
    # alpha = log(0.1),
    # tau = 0,
    # theta = log(0.5),
    # n_clusters = 12,
    # n_time_points = 7,
    # n_ind_per_cluster = 20,
    # data_type = "normal",
    # sigma = 0.3,
    # delay_model = list(type="exp", params=list(d=1))
  )

  results <- run_analysis(
    data = data,
    analysis = L$analysis,
    data_type = L$data_type,
    L = L,
    C = C
  )

  return (results)

}
