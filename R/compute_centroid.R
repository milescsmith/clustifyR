#' @noRd
#' @importFrom Matrix rowMeans
compute_centroid <- function(expr_mat, sc_cluster) {
  num_cluster <- max(sc_cluster)
  num_genes <- nrow(expr_mat)
  mu_info <- matrix(NA, nrow = num_genes, ncol = num_cluster)
  sigma_info <- matrix(NA, nrow = num_genes, ncol = num_cluster)
  for (i in 1:num_cluster) {
    curr_cluster_expr <- expr_mat[, sc_cluster == i]
    mu_info[, i] <- rowMeans(curr_cluster_expr)
    sigma_info[, i] <- sqrt(rowMeans(curr_cluster_expr^2) - mu_info[, i]^2)
  }
  return(list(mu = mu_info, sigma = sigma_info))
}
