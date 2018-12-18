#' @noRd
compute_bulk_score <- function(expr_mat, bulk_mat, num_cluster, default_similiarity, compute_method, ...) {
  # for clusters corresponding to bulk data: compute similarity score
  num_cells <- ncol(expr_mat)
  num_bulk <- ncol(bulk_mat)
  bulk_score <- matrix(NA, nrow = num_cells, ncol = num_bulk)
  for (i in 1:num_cells) {
    curr_cell_expr <- expr_mat[, i]
    bulk_score[i, ] <- sapply(1:num_bulk, function(x) calc_similarity(curr_cell_expr, bulk_mat[, x], compute_method, ...))
  }
  # for clusters without bulk data, use a default value
  if (num_cluster > num_bulk) {
    for (j in (num_bulk + 1):num_cluster) {
      bulk_score <- cbind(bulk_score, rep(default_similiarity, num_cells))
    }
  }
  return(bulk_score)
}
