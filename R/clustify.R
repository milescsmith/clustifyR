#' @title clustify
#' @description Main function to compare scRNA-seq data to bulk RNA-seq data.
#'
#' @param expr_mat Single-cell expression matrix or Seurat object
#' @param query_metadata Query cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set. Not required if running correlation per cell.
#' @param ref_mat Reference bulk expression matrix
#' @param ref_metadata Reference metadata.  If provided, will be used to rename study observations with cell names.
#' @param ref_columns A list of columns in the reference metadata to use for translating from study observaton to cell name.  Default: the first and second columns.
#' @param query_metadata
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col Column in the query metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param per_cell If true run per cell, otherwise per cluster.  Default: FALSE.
#' @param num_perm number of permutations.  Default: 0
#' @param compute_method Correlation method(s) to use in computing similarity scores.  Default: 'spearman'
#' @param use_var_genes If providing a seurat object, should only the variable genes be used?  Default: FALSE.
#' @param ... additional arguments to pass to compute_method function
#'
#' @rdname clustify
#' @export clustify
setGeneric(name = "clustify",
           def = function(expr_mat, ...) {
             standardGeneric("clustify")
             }
)
#'@export
clustifyr_methods <- c(
  "pearson",
  "spearman",
  "cosine",
  "kl_divergence"
)

#' @rdname clustify
#'
#' @importFrom plyr mapvalues
#'
#' @return \code{NULL}
#' @rdname clustify
#' @method clustify default
#'
setMethod('clustify',
          signature(expr_mat = "matrix"),
          function(expr_mat,
                   ref_mat,
                   ref_metadata = NULL,
                   ref_columns = NULL,
                   query_metadata = NULL,
                   query_genes = NULL,
                   cluster_col = NULL,
                   per_cell = FALSE,
                   num_perm = 0,
                   compute_method = "pearson",
                   ...) {
            if (!compute_method %in% clustifyr_methods) {
              stop(paste(compute_method, "correlation method not implemented"))
            }

            if (!is.null(ref_metadata)) {
              if (!is.null(ref_columns)){
                translation_table <- data.frame(ref_metadata[,ref_columns])
              } else {
                translation_table <- data.frame(ref_metadata[,c(1:2)])
              }
              colnames(ref_mat) <- mapvalues(x = colnames(ref_mat),
                                             from = translation_table[,1],
                                             to = translation_table[,2])
            }
            # select gene subsets
            gene_constraints <- get_common_elements(
              rownames(expr_mat),
              rownames(ref_mat),
              query_genes
            )

            expr_mat <- expr_mat[gene_constraints, , drop = FALSE]
            ref_mat <- ref_mat[gene_constraints, , drop = FALSE]


            if (!per_cell) {
              if (is.vector(query_metadata)) {
                cluster_ids <- query_metadata
              } else if (is.data.frame(query_metadata) & !is.null(cluster_col)) {
                cluster_ids <- query_metadata[[cluster_col]]
              } else {
                stop("query_metadata not formatted correctly,
           supply either a character vector or a dataframe")
              }
            }

            if (per_cell) {
              cluster_ids <- colnames(expr_mat)
            }

            if (num_perm == 0) {
              res <- get_similarity(
                expr_mat,
                ref_mat,
                cluster_ids = cluster_ids,
                per_cell = per_cell,
                compute_method = compute_method, ...
              )
            } else {
              # run permutation
              res <- permute_similarity(
                expr_mat,
                ref_mat,
                cluster_ids = cluster_ids,
                num_perm = num_perm,
                per_cell = per_cell,
                compute_method = compute_method,
                ...
              )
            }
            return(res)
          }
)

#' @rdname clustify
#'
#' @importFrom Seurat VariableFeatures GetAssayData
#' @importFrom magrittr "%>%"
#'
#' @return \code{NULL}
#' @rdname clustify
#' @method clustify Seurat
setMethod('clustify',
          signature(expr_mat = "Seurat"),
          function(expr_mat,
                   use_var_genes = FALSE,
                   ...) {
            expr_data <- GetAssayData(expr_mat) %>% as.matrix()
            query_metadata <- expr_mat@meta.data

            if (use_var_genes) {
              query_genes <- VariableFeatures(object = expr_mat)
            }

            res <- clustify(expr_mat = expr_data,
                            query_genes = query_genes,
                            query_metadata = query_metadata,
                            ...)
            return(res)
          }
)
