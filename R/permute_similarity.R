#' Compute a p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#'
#' @param expr_mat single-cell expression matrix
#' @param cluster_ids clustering info of single-cell data assume that
#'  genes have ALREADY BEEN filtered
#' @param bulk_mat bulk expression matrix
#' @param num_perm number of permutations
#' @param per_cell run per cell?
#' @param compute_method method(s) for computing similarity scores
#' @param ... additional parameters
#' @noRd
permute_similarity <- function(expr_mat,
                               bulk_mat,
                               cluster_ids,
                               num_perm,
                               per_cell = F,
                               compute_method, ...) {
  bulk_clust <- colnames(bulk_mat)

  if (!per_cell) {
    sc_clust <- sort(unique(cluster_ids))
    clust_avg <- compute_mean_expr(
      expr_mat,
      cluster_ids,
      sc_clust
    )
  } else {
    sc_clust <- colnames(expr_mat)
    clust_avg <- expr_mat
  }

  assigned_score <- calc_similarity(
    clust_avg,
    bulk_mat,
    compute_method,
    ...
  )

  # perform permutation
  sig_counts <- matrix(0L, nrow = length(sc_clust), ncol = length(bulk_clust))

  for (i in 1:num_perm) {
    resampled <- sample(cluster_ids,
      length(cluster_ids),
      replace = FALSE
    )

    if (!per_cell) {
      permuted_avg <- compute_mean_expr(
        expr_mat,
        resampled,
        sc_clust
      )
    } else {
      permuted_avg <- expr_mat[, resampled]
    }

    # permutate assignment
    new_score <- calc_similarity(
      permuted_avg,
      bulk_mat,
      compute_method,
      ...
    )
    sig_counts <- sig_counts + as.numeric(new_score > assigned_score)
  }

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- bulk_clust
  rownames(sig_counts) <- sc_clust
  colnames(sig_counts) <- bulk_clust

  return(list(
    score = assigned_score,
    p_val = sig_counts / num_perm
  ))
}
