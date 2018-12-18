#' Plot gene expression on to tSNE
#'
#'
#' @param expr_mat input single cell matrix
#' @param metadata data.frame with tSNE coordinates
#' @param genes gene(s) to color tSNE
#' @param cell_col column name in metadata containing cell ids, defaults
#' to rownames if not supplied
#' @param ... additional arguments passed to `[clustifyR::plot_tsne()]`
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @export
plot_gene <- function(expr_mat,
                      metadata,
                      genes,
                      cell_col = NULL,
                      ...) {
  genes_to_plot <- genes[genes %in% rownames(expr_mat)]
  genes_missing <- setdiff(genes_to_plot, genes)

  if (length(genes_missing) != 0) {
    warning(paste0(
      "the following genes were not present in the input matrix ",
      paste(genes_missing, collapse = ",")
    ))
  }

  if (length(genes_to_plot) == 0) {
    stop("no genes present to plot")
  }
  expr_dat <- t(as.matrix(expr_mat[genes_to_plot, ]))
  expr_dat <- rownames_to_column(as.data.frame(expr_dat), "cell")

  if (is.null(cell_col)) {
    mdata <- rownames_to_column(metadata, "cell")
    cell_col <- "cell"
  } else {
    mdata <- metadata
  }

  if (!cell_col %in% colnames(mdata)) {
    stop("please supply a cell_col that is present in metadata")
  }

  plt_dat <- left_join(expr_dat, metadata,
    by = c("cell" = cell_col)
  )

  lapply(
    genes,
    function(gene) {
      plot_tsne(plt_dat,
        feature = gene,
        legend_name = gene,
        ...
      )
    }
  )
}
