#' Plot called clusters on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#'
#' @export
plot_call <- function(correlation_matrix,
                      metadata,
                      bulk_data_to_plot = colnames(correlation_matrix)) {
  df_temp <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp == 0] <- "1"
  df_temp[df_temp != "1"] <- "0"
  plot_cor(df_temp, metadata, bulk_data_to_plot)
}
