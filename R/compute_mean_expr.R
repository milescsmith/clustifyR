#' compute mean of clusters
#' @noRd
#' @importFrom Matrix rowMeans
compute_mean_expr <- function(expr_mat, sc_assign, sc_clust) {
  return(sapply(sc_clust, function(x) rowMeans(expr_mat[, sc_assign == x])))
}
