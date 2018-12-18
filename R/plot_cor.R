#' Plot similarity measures on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with per cell tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param cluster_col colname of clustering data in metadata, defaults to rownames of the
#' metadata if not supplied.
#' @param dim1_col metadata column name with 1st axis dimension.
#' defaults to "tSNE_1".
#' @param dim2_col metadata column name with 2nd axis dimension.
#' defaults to "tSNE_2".
#' @param scale_legends if TRUE scale all legends to maximum values in entire
#' correlation matrix. if FALSE scale legends to maximum for each plot. A
#' two-element numeric vector can also be passed to supply custom values i.e. c(0, 1)
#' @param ... passed to plot_tsne
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join matches filter
#' @importFrom tidyr gather
#'
#' @export
plot_cor <- function(correlation_matrix,
                     metadata,
                     bulk_data_to_plot = colnames(correlation_matrix),
                     cluster_col = NULL,
                     dim1_col = "tSNE_1",
                     dim2_col = "tSNE_2",
                     scale_legends = FALSE,
                     ...) {
  if (!any(bulk_data_to_plot %in% colnames(correlation_matrix))) {
    stop("cluster ids not shared between metadata and correlation matrix")
  }

  if (is.null(cluster_col)) {
    cluster_col <- "rownames"
    metadata <- rownames_to_column(metadata, cluster_col)
  }

  cor_df <- as.data.frame(correlation_matrix)
  cor_df <- rownames_to_column(cor_df, cluster_col)
  cor_df_long <- gather(cor_df,
                        bulk_cluster,
                        expr,
                        -matches(cluster_col)
  )

  # checks matrix rownames, 2 branches for cluster number (avg) or cell bar code (each cell)
  if (cor_df[[cluster_col]][1] %in% metadata[[cluster_col]]) {
    plt_data <- left_join(cor_df_long,
      metadata,
      by = cluster_col
    )
  } else {
    plt_data <- left_join(cor_df_long,
                          metadata,
                          by = structure(names = cluster_col, "rn")
    )
  }

  # determine scaling method, either same for all plots, or per plot (default)
  if (typeof(scale_legends) == "logical" && scale_legends) {
    scale_limits <- c(
      ifelse(min(plt_data$expr) < 0,
        min(plt_data$expr),
        0
      ),
      max(max(plt_data$expr))
    )
  } else if (typeof(scale_legends) == "logical" && !scale_legends) {
    scale_limits <- NULL
  } else {
    scale_limits <- scale_legends
  }

  plts <- vector("list", length(bulk_data_to_plot))
  for (i in seq_along(bulk_data_to_plot)) {
    tmp_data <- filter(
      plt_data,
      bulk_cluster == bulk_data_to_plot[i]
    )
    plts[[i]] <- plot_tsne(tmp_data,
      x = dim1_col,
      y = dim2_col,
      feature = "expr",
      legend_name = bulk_data_to_plot[i],
      scale_limits = scale_limits,
      ...
    )
  }
  plts
}
