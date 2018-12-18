#' Compute the p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#' @param expr_mat single-cell expression matrix
#' @param bulk_mat bulk expression matrix
#' @param cluster_ids vector of cluster ids for each cell
#' @param compute_method method(s) for computing similarity scores
#' @param per_cell run per cell?
#' @param ... additional parameters not used yet
#' @noRd
get_similarity <- function(expr_mat,
                           bulk_mat,
                           cluster_ids,
                           compute_method,
                           per_cell = FALSE,
                           ...) {
  bulk_clust <- colnames(bulk_mat)

  if (!per_cell) {
    sc_clust <- sort(unique(cluster_ids))
    clust_avg <- compute_mean_expr(
      expr_mat,
      cluster_ids,
      sc_clust
    )
  } else {
    sc_clust <- cluster_ids
    clust_avg <- expr_mat
  }

  assigned_score <- calc_similarity(
    clust_avg,
    bulk_mat,
    compute_method,
    ...
  )

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- bulk_clust

  return(assigned_score)
}
