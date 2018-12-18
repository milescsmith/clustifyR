#' Compute similarity between two vectors
#'
#' @description Compute the similarity score between two vectors using a
#' customized scoring function
#' Two vectors may be from either scRNA-seq or bulk RNA-seq data.
#' The lengths of vec1 and vec2 must match, and must be arranged in the
#' same order of genes.
#' Both vectors should be provided to this function after pre-processing,
#' feature selection and dimension reduction.
#'
#' @param vec1 test vector
#' @param vec2 reference vector
#' @param compute_method method to run i.e. corr_coef
#' @param ... arguments to pass to compute_method function
#' @export
vector_similarity <- function(vec1, vec2, compute_method, ...) {
  # examine whether two vectors are of the same size
  if (!is.numeric(vec1) || !is.numeric(vec2) || length(vec1) != length(vec2)) {
    stop("compute_similarity: two input vectors are not numeric or of different sizes.")
  }

  if (!(compute_method %in% c("cosine", "kl_divergence"))) {
    stop(paste(compute_method, "not implemented"))
  }

  if (compute_method == "kl_divergence") {
    res <- kl_divergence(vec1, vec2, ...)
  } else if (compute_method == "cosine") {
    res <- cosine(vec1, vec2, ...)
  }
  # return the similarity score, must be
  return(res)
}
