#' calculate adjusted p-values for hypergeometric test of gene lists
#' or jaccard index
#'
#' @param bin_mat binarized single-cell expression matrix
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#' @param n number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#'
#' @export
compare_lists <- function(bin_mat,
                          marker_m,
                          n = 30000,
                          metric = "hyper",
                          output_high = TRUE) {
  # check if matrix is binarized
  if (length(unique(bin_mat[, 1])) > 2) {
    metric <- "spearman"
  }

  # "expressed" genes per single cell data cluster
  if (metric == "hyper") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
            bin_temp <- bin_mat[, x][bin_mat[, x] == 1]
            list_top <- names(bin_temp)

            t <- length(intersect(list_top, marker_list))
            a <- max(length(list_top), length(marker_list))
            b <- min(length(list_top), length(marker_list))
            sum(dhyper(t:b, a, n - a, b))
          }
        )
        do.call(cbind, as.list(p.adjust(per_col)))
      }
    )
  }

  if (metric == "jaccard") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
            bin_temp <- bin_mat[, x][bin_mat[, x] == 1]
            list_top <- names(bin_temp)

            I <- length(intersect(list_top, marker_list))
            I / (length(list_top) + length(marker_list) - I)
          }
        )
        do.call(cbind, per_col)
      }
    )
  }

  if (metric == "spearman") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
            v1 <- marker_list
            bin_temp <- as.matrix(bin_mat)[, x]
            bin_temp <- bin_temp[order(bin_temp, decreasing = TRUE)]
            list_top <- names(bin_temp)
            v2 <- list_top[list_top %in% v1]
            sum(sapply(seq_along(v1), function(i) abs(i - (which(v2 == v1[i])))))
          }
        )
        do.call(cbind, per_col)
      }
    )
  }

  res <- do.call(rbind, out)
  rownames(res) <- colnames(bin_mat)
  colnames(res) <- colnames(marker_m)

  if (output_high == TRUE) {
    if (metric == "hyper") {
      res <- -log10(res)
    } else if (metric == "spearman") {
      res <- -res
    }
  }

  return(res)
}
