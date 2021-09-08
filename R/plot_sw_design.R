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
  design %<>% mutate(
    state_times = ifelse(state=="Control", state_times+0.05, state_times),
    cluster = factor(cluster, labels=sample(c(1:n_clusters)))
  )

  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  sw_plot <- ggplot(design, aes(x=state_times, y=as.factor(cluster))) +
    geom_col(
      aes(fill = as.factor(state)),
      position = position_stack(reverse=TRUE)
    ) +
    labs(title=title[1], fill="State", x="Time", y="Cluster") +
    # geom_vline(xintercept = c(0:4), linetype="longdash", color="lightgrey") +
    scale_fill_manual(values=cb_colors[c(2,3)]) +
    coord_cartesian(xlim=c(1.1,5))


  if (compare_to_parallel) {

    design_par <- data.frame(
      cluster = rep(1:n_clusters),
      state_times = rep(J+0.05,n_clusters),
      state = rep(c("Control","Treatment"), each=n_clusters/2)
    )
    design_par %<>% mutate(
      cluster = factor(cluster, labels=sample(c(1:n_clusters)))
    )

    par_plot <- ggplot(design_par, aes(x=state_times, y=as.factor(cluster))) +
      geom_col(
        aes(fill = as.factor(state)),
        position = position_stack(reverse=TRUE)
      ) +
      labs(title=title[2], fill="State", x="Time", y="Cluster") +
      scale_fill_manual(values=cb_colors[c(2,3)]) +
      coord_cartesian(xlim=c(1.1,5))

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
