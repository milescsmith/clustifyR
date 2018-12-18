#' Plot a tSNE colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#' @param legend_name legend name to display, defaults to no name
#' @param c_cols character vector of colors to built color gradient
#' for continuous values. defaults to [`clustifyR::pretty_palette`]
#' @param d_cols character vector of colors for discrete values.
#' defaults to RColorBrewer paired palette
#' @param pt_size point size
#' @param scale_limits defaults to min = 0, max = max(data$x),
#' otherwise a two-element numeric vector indicating min and max to plot
#'
#' @importFrom glue glue
#' @importFrom dplyr arrange
#' @importFrom rlang !!
#' @importFrom cowplot theme_cowplot
#' @import ggplot2
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2",
                      feature,
                      legend_name = "",
                      c_cols = pretty_palette,
                      d_cols = NULL,
                      pt_size = 0.25,
                      scale_limits = NULL) {

  # sort data to avoid plotting null values over colors
  data <- arrange(data, !!sym(feature))

  p <- ggplot(data, aes_string(x, y)) +
    geom_point(aes_string(color = glue("`{feature}`")), # backticks protect special character gene names
      size = pt_size
    )

  if (typeof(data[[feature]]) %in% c("character","logical") | is.factor(data[[feature]])) {
    if (!is.null(d_cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = d_cols,
        name = legend_name
      )
    } else {
      p <- p + scale_color_brewer(
        palette = "Paired",
        name = legend_name
      )
    }
  } else {
    # continuous values
    if (is.null(scale_limits)) {
      scale_limits <- c(
        ifelse(min(data[[feature]]) < 0,
          min(data[[feature]]),
          0
        ),
        max(data[[feature]])
      )
    }
    p <- p + scale_color_gradientn(
      colors = c_cols,
      name = legend_name,
      limits = scale_limits
    )
  }

  p + theme_cowplot()
}
