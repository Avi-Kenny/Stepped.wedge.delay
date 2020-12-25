#' Plot the stepped wedge design
#'
#' @param data A dataset returned by generate_dataset()
#' @param title Either a character string, or waiver() for no title
#' @return A ggplot2 plot object representing the design

plot_sw_design <- function(data, title="SW design diagram",
                           compare_to_parallel=FALSE) {

  n_clusters <- data$params$n_clusters
  crossover_times <- data$params$crossover_times
  J <- data$params$n_time_points

  design <- data.frame(
    cluster = rep(1:n_clusters,each=2),
    state_times = c(rbind(
      crossover_times-1,
      J - crossover_times + 1
    )),
    state = rep(c("Control","Treatment"),n_clusters)
  )

  sw_plot <- ggplot(design, aes(x=state_times, y=as.factor(cluster))) +
    geom_col(
      aes(fill = as.factor(state)),
      position = position_stack(reverse=TRUE)
    ) +
    labs(title=title[1], fill="State", x="Time", y="Cluster")

  if (compare_to_parallel) {

    design_par <- data.frame(
      cluster = rep(1:n_clusters,each=2),
      state_times = c(rep(c(0,J),n_clusters/2), rep(c(J,0),n_clusters/2)),
      state = rep(c("Control","Treatment"),n_clusters)
    )

    par_plot <- ggplot(design_par, aes(x=state_times, y=as.factor(cluster))) +
      geom_col(
        aes(fill = as.factor(state)),
        position = position_stack(reverse=TRUE)
      ) +
      labs(title=title[2], fill="State", x="Time", y="Cluster")

    return(

      ggarrange(
        par_plot,
        sw_plot,
        ncol = 2,
        nrow = 1,
        common.legend = TRUE,
        legend = "bottom"
      )

    )

  } else {

    return(sw_plot)

  }

}
