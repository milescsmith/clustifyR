#' generate variable gene list from marker matrix
#'
#' @description Variable gene list is required for run_cor main function. This
#' function parses variables genes from a matrix input.
#'
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#'
#' @export
get_vargenes <- function(marker_m) {
  return(unique(unlist(marker_m, use.names = FALSE)))
}
