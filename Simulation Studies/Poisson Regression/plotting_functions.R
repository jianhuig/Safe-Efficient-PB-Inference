plot_results <- function(df, x_title = " ", trans_type = "log"){
  
  plot_df <- df %>%
    pivot_longer(
      cols = c(mean_ci_width, cp),
      names_to = "metric",
      values_to = "value"
    )
  
  plot_df <- df %>%
    pivot_longer(
      cols = c(mean_ci_width, cp),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(!Method %in% c("naive", "true")) %>%
    mutate(
      Method = factor(
        Method,
        levels = c("classical", "ppi", "pdc", "chen-chen")
      )
    )
  
  facet_labels <- c(
    mean_ci_width = "95% Confidence Interval Width",
    cp = "Coverage Probability"
  )
  
  plot <- ggplot(plot_df, aes(x = param, y = value, color = Method)) +
    geom_line() +
    geom_point() +
    geom_hline(
      data = subset(plot_df, metric == "cp"),
      aes(yintercept = 0.95),
      linetype = "dashed",
      linewidth = 0.8,
      color = "black"
    ) +
    facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = facet_labels)) +
    facetted_pos_scales(
      y = list(
        metric == "cp" ~ scale_y_continuous(limits = c(0.9, 1)),
        metric == "mean_ci_width" ~ scale_y_continuous()
      )) + 
    labs(
      x = x_title,
      y = NULL,
      color = "Method"
    ) +
    theme_bw() +
    theme(
      # Axis labels
      axis.title.x = element_text(size = 16),   # x-axis label
      axis.title.y = element_text(size = 16),   # y-axis label
      
      # Axis tick labels
      axis.text.x = element_text(size = 14),    # x-axis tick numbers
      axis.text.y = element_text(size = 14),    # y-axis tick numbers
      
      # Facet strip labels
      strip.text = element_text(size = 16)  # facet labels
    ) +
    scale_color_discrete(
      limits = c("classical", "ppi", "pdc", "chen-chen"),
      labels = c(
        "classical" = "Classical",
        "ppi" = "PPI",
        "pdc" = "PDC",
        "chen-chen" = "CC"
      )
    )
  
  if(trans_type == "log"){
    
    plot <- plot + scale_x_log10(
      breaks = c(0.1, 0.5, 1, 10),           # only these will get labels
      minor_breaks = plot_df$param  # ~5 nicely spaced ticks
    ) 
    
  }
  
  return(plot)
  
}
