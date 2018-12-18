#' Binarize scRNA seq data
#'
#' @param expr_mat single-cell expression matrix
#' @param n number of top expressing genes to keep
#'
#' @export
binarize_expr <- function(expr_mat,
                          n = 1000) {
  expr_df <- as.data.frame(as.matrix(expr_mat))
  df_temp <- apply(expr_df, 2, function(x) x - x[order(x, decreasing = TRUE)[n + 1]])
  df_temp[df_temp > 0] <- 1
  df_temp[df_temp < 0] <- 0
  df_temp
}
