#' Cosine distance
#' @param vec1 test vector
#' @param vec2 reference vector
#' @export
cosine <- function(vec1, vec2) {
  return(sum(vec1 * vec2) / sqrt(sum(vec1^2) * sum(vec2^2)))
}
