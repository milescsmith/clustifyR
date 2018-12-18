#' Average expression values per cluster
#'
#' @param mat expression matrix
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#'
#' @importFrom Matrix rowMeans
#' @export
average_clusters <- function(mat, cluster_info,
                             log_scale = T,
                             cluster_col = "cluster") {
  if (is.vector(cluster_info)) {
    cluster_ids <- split(colnames(mat), cluster_info)
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
  } else {
    stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
  }

  out <- lapply(
    cluster_ids,
    function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix")
      }
      if (log_scale) {
        mat_data <- expm1(mat[, cell_ids])
      } else {
        mat_data <- mat[, cell_ids]
      }
      res <- rowMeans(mat_data)
      if (log_scale) {
        res <- log1p(res)
      }
      res
    }
  )
  return(do.call(cbind, out))
}
