#' Find entries shared in all vectors
#' @description return entries found in all supplied vectors. If the vector supplied
#' is NULL or NA, then it will be excluded from the comparision.
#' @param ... vectors
#' @examples
#' a <- rep(1:5)
#' b <- rep(4:10)
#' c <- rep(4:6)
#' get_common_elements(a, b, c)
#' @export
get_common_elements <- function(...) {
  vecs <- list(...)
  # drop NULL elements of list
  vecs <- vecs[!sapply(vecs, is.null)]
  # drop NA elements of list (NA values OK in a vector)
  vecs <- vecs[!is.na(vecs)]

  return(Reduce(intersect, vecs))
}
