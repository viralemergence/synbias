#' This function creates a set of partial dependence plots
#' @param variables Vector of the variables to be plotted (if log transformed,
#' they need to have "log_" added to the front). They should be of length 16.
#' @param shapley_data These are the Shapley values produced from the
#' bootstrapping process. They do not need to be processed, but it should be a
#'  data.frame with columns for: variable (variable name), value (Shapley
#' value), rfvalue (training data value), stdfvalue (standardized training
#'  data value from 0 to 1), mean_value (mean Shapley value), ID
#' (bootstrap ID), row_index (row of training data associated with this
#' value)
#' @param train_data The data.frame used to train these models (used to grab
#' the distribution of values for each variable). The column names should
#' match what is seen in the variables vector as that is what is used to
#' grab the values.
#' @param colors_vec A vector of color values that should be equal in length
#'  to the number of variables given
#' @param data_labels A character vector of labels for the subplots (if NULL,
#' the column names will be used). The labels must be in the same order as
#' `variables`.
#' @param importance An optional argument that allows for the subplots to be
#' ordered from most to least important. If `importance` is provided in the
#' form of a data frame with a `character`/`factor` column with the `variables`
#' and a `double` column with importance scores, the subplots will be displayed
#' from most to least important.
#' @return Returns a partial dependence plot object that is actually 33
#' panels (16 panels of partial dependence plots, 16 panels of the density
#'  curves showing the distribution of values in the training data, panel of
#'  enlarged y-axis title from {cowplot})
#' @example partial.dependence.plot(CHC, SHAP, DATA, p_cols)
partial.dependence.plot <- function(
  variables,
  shapley_data,
  train_data,
  colors_vec,
  data_labels = NULL,
  importance = NULL
) {
  # Argument checking
  num_variables <- length(variables)
  if (num_variables != length(colors_vec)) {
    stop("The number of variables must match the number of colors provided.")
  }
  if (num_variables < 1) {
    stop("There must be at least one variable to plot.")
  }
  if (!is.data.frame(shapley_data)) {
    stop("shapley_data must be a data.frame.")
  }
  if (
    !all(
      c("variable", "value", "rfvalue", "ID") %in% colnames(shapley_data)
    )
  ) {
    stop(
      "shapley_data must contain the columns: variable, value, rfvalue, and ID."
    )
  }
  if (!all(variables %in% colnames(train_data))) {
    stop("All variables must be present in the train_data.")
  }
  if (!is.null(importance)) {
    if (!is.data.frame(importance)) {
      stop("`importance` must be a data frame.")
    }
    if (!all(complete.cases(importance))) {
      stop("`importance` must not contain NA values.")
    }
    if (ncol(importance) > 2) {
      stop(
        "`importance` must have 2 columns: a column with the variable names
       and a column with the importance scores."
      )
    }
    importance_class <- map_chr(importance, class)
    if (!any(c("character", "factor") %in% importance_class)) {
      stop(
        "`importance` must have a character or factor column defining the
        variable names"
      )
    }
    if (!("numeric" %in% importance_class)) {
      stop("`importance` must have a numeric column with importance scores")
    }
    score_index <- names(importance)[which(importance_class == "numeric")]
    name_index <- which(
      importance_class %in% c("character", "factor")
    )
    if (!all(variables %in% pull(importance, name_index))) {
      stop("All `variables` must be in the `importance` data frame")
    }
    importance <- importance[match(variables, pull(importance, name_index)), ]
    if (!is.null(data_labels)) {
      importance$labels <- data_labels
    }
    setorderv(importance, score_index, -1)
    variables <- pull(importance, name_index)
    data_labels <- pull(importance, labels)
  }
  if (is.null(data_labels)) {
    data_labels <- variables
  }
  if (length(data_labels) != num_variables) {
    stop("The number of labels must match the number of variables.")
  }
  if (!is.character(data_labels)) {
    stop("Labels must be a character vector.")
  }

  # a labeller for later
  dropLeadingZero <- function(l) {
    str_replace(l, '0(?=\\.)', '')
  }

  # Data manipulation
  MNVL <- shapley_data %>%
    filter(!is.na(rfvalue)) %>%
    group_by(variable, ID, rfvalue) %>%
    summarise(
      value = mean(value, na.rm = T),
      rfvalue = mean(rfvalue, na.rm = T)
    ) %>%
    filter(variable %in% variables) %>%
    ungroup() %>%
    data.frame()
  PDPR <- shapley_data %>%
    filter(variable %in% variables)

  # Plotting
  PDP <- lapply(1:num_variables, function(i) {
    PDPR %>%
      filter(variable == variables[i]) %>%
      ggplot(aes(x = rfvalue, y = value, group = ID)) +
      stat_smooth(geom = "line", se = F, alpha = 0.5, color = "grey70") +
      geom_smooth(
        data = MNVL %>% filter(variable == variables[i]),
        inherit.aes = F,
        method = "loess",
        aes(x = rfvalue, y = value),
        color = colors_vec[i],
        se = F,
        linewidth = 1.5
      ) +
      labs(y = NULL, x = data_labels[i]) +
      theme(
        aspect.ratio = 0.45,
        panel.background = element_rect(fill = "white", color = "grey50"),
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        plot.margin = margin(0, 0, 0, 0, "pt")
      )
  })
  p_lab <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Shapley Score",
      angle = 90,
      #fontface = "bold",
      size = 14
    )
  RUGV <- train_data[, variables] %>%
    reshape2::melt() %>%
    rename(rfvalue = value) %>%
    mutate(value = 0)
  TPLT <- map(1:num_variables, function(i) {
    RUGV %>%
      filter(variable == variables[i]) %>%
      ggplot() +
      geom_density(
        aes(x = round(rfvalue, digits = 2)),
        fill = colors_vec[i],
        alpha = 0.4,
        color = colors_vec[i]
      ) +
      scale_x_continuous(position = "bottom", labels = dropLeadingZero) +
      scale_y_reverse(expand = expansion(mult = c(0.1, 0))) +
      labs(x = data_labels[i], y = NULL) +
      theme(
        panel.background = element_rect(fill = "white", color = "grey50"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = "transparent"),
        aspect.ratio = 0.15,
        plot.margin = margin(0, 0, 0, 0, "pt")
      )
  })
  PTPL <- lapply(1:num_variables, function(i) PDP[[i]] / TPLT[[i]])
  cowplot::plot_grid(
    p_lab,
    patchwork::wrap_plots(PTPL, ncol = round(sqrt(num_variables))),
    nrow = 1,
    rel_widths = c(0.05, 1),
    rel_heights = c(1, 1)
  ) +
    theme(plot.background = element_rect(fill = "white", color = "transparent"))
}
