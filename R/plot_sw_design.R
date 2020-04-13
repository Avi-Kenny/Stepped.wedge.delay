#' Plot the stepped wedge design
#'
#' @param data A dataset returned by generate_dataset()
#' @return A ggplot2 plot object representing the design
#' @export
plot_sw_design <- function(data) {

  n_clusters <- data$params$n_clusters
  crossover_times <- eval(data$params$crossover_times)

  design <- data.frame(
    cluster = rep(1:n_clusters,each=2),
    state_times = c(rbind(
      crossover_times,
      (max(crossover_times)+1) - crossover_times
    )),
    state = rep(c("Control","Treatment"),n_clusters)
  )

  return(
    ggplot(design, aes(x=state_times, y=as.factor(cluster))) +
      geom_col(
        aes(fill = as.factor(state)),
        position = position_stack(reverse=TRUE)
      ) +
      labs(title="SW design diagram",
           fill="State", x="Time", y="Cluster")
  )

}
