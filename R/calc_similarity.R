#' compute similarity
#'
#' @importFrom magrittr %<>%
#' @noRd
calc_similarity <- function(sc_avg,
                            bulk_mat,
                            compute_method, ...) {

  # use stats::cor matrix method if possible
  if (any(compute_method %in% c("pearson", "spearman"))) {
    if (!is.matrix(sc_avg)) {
      sc_avg %<>% as.matrix()
    }
    if (!is.matrix(bulk_mat)) {
      bulk_mat %<>% as.matrix()
    }

    similarity_score <- cor(sc_avg,
      bulk_mat,
      method = compute_method
    )
    return(similarity_score)
  }

  sc_clust <- colnames(sc_avg)
  bulk_clust <- colnames(bulk_mat)
  similarity_score <- matrix(NA,
    nrow = length(sc_clust),
    ncol = length(bulk_clust)
  )
  for (i in seq_along(sc_clust)) {
    for (j in seq_along(bulk_clust)) {
      similarity_score[i, j] <- vector_similarity(
        sc_avg[, sc_clust[i]],
        bulk_mat[, bulk_clust[j]],
        compute_method, ...
      )
    }
  }
  return(similarity_score)
}
