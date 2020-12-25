#' Plot the outcome variable by cluster over time
#'
#' @param data A dataset returned by generate_dataset()
#' @param type Either "no error" or "realized"; "no error" gives the actual
#'     value of p_ij or y_ij, whereas "realized" gives the average of the y
#'     values
#' @return A ggplot2 plot object: spaghetti plot of mean cluster outcomes over
#'     time

plot_outcome <- function(data, type) {

  if (data$params$type == "binomial") {
    data$data %<>% rename("y_ij" = `p_ij`)
  }

  if (type=="no error") {
    d <- data$data %>% group_by(i,j) %>% summarize(y_grp=mean(y_ij))
    title <- "Outcome over time, by cluster (no error) (jittered)"
  }

  if (type=="realized") {
    d <- data$data %>% group_by(i,j) %>%
      summarize(y_grp=mean(y))
    title <- "Outcome over time, by cluster (realized) (jittered)"
  }

  return(
    ggplot(d, aes(x=j, y=y_grp, group=i, color=as.factor(i))) +
      geom_line(position=position_jitter(w=0.1)) +
      labs(title=title,
           x="Time", y="Outcome", color="Cluster") +
      theme(legend.position = "none")
  )

}
