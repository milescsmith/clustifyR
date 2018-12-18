#' convert candidate genes list into matrix
#'
#' @param marker_df dataframe of candidate genes, must contain "gene" and "cluster" columns, or a matrix of gene names to convert to ranked
#' @param ranked unranked gene list feeds into hyperp, the ranked
#' gene list feeds into regular corr_coef
#' @param n number of genes to use
#' @param step_weight ranked genes are tranformed into pseudo expression by descending weight
#' @param background_weight ranked genes are tranformed into pseudo expression with
#' added weight
#' @param unique whether to use only unique markers to 1 cluster
#' @param labels vector or dataframe of cluster names
#'
#' @importFrom tidyr gather spread
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
matrixize_markers <- function(marker_df,
                              ranked = FALSE,
                              n = NULL,
                              step_weight = 1,
                              background_weight = 0,
                              unique = FALSE,
                              labels = NULL) {
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  marker_df <- as_tibble(marker_df)

  # if "gene" not present in column names, assume df is a matrix to be converted to ranked
  if (!("gene" %in% colnames(marker_df))) {
    marker_df <- data.frame(lapply(marker_df,
                                   as.character),
                            stringsAsFactors = FALSE) %>%
      gather(factor_key = TRUE,
             key = "cluster",
             value = "gene")
  }

  if (unique == TRUE) {
    nonunique <- marker_df %>%
      group_by(gene) %>%
      summarise(n = n()) %>%
      filter(n > 1)
    marker_df <- anti_join(marker_df,
                           nonunique,
                           by = "gene")
  }

  cut_num <- min((marker_df %>%
                    group_by(cluster) %>%
                    summarise(n = n()))$n)

  if (!is.null(n)) {
    if (n < cut_num) {
      cut_num <- n
    }
  }

  marker_temp <- marker_df %>%
    select(gene, cluster) %>%
    group_by(cluster) %>%
    slice(1:cut_num)
  if (ranked == TRUE) {
    marker_temp <- marker_temp %>%
      mutate(n = seq(step_weight * cut_num,
                     by = -step_weight,
                     length.out = cut_num) +
               background_weight)
    marker_temp2 <- as.data.frame(spread(marker_temp, key = "cluster", value = n) %>%
                                    replace(is.na(.), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- marker_temp2 %>% select(-gene)
  } else {
    marker_temp <- marker_temp %>% mutate(n = 1:cut_num)
    marker_temp2 <- as.data.frame(spread(marker_temp, key = "cluster", value = "gene") %>% select(-n))
  }

  # if labels is vector, adopt names in vector; if labels is a metadata dataframe, pulls names from "classified" column
  if (!is.null(labels)) {
    if (typeof(labels) != "character") {
      label_df <- labels
      labels <- left_join(data_frame(cluster = colnames(marker_temp2)),
        unique(data_frame(
          cluster = labels$cluster,
          classified = labels$classified)),
        by = "cluster") %>%
        pull(classified)
    }
    colnames(marker_temp2) <- labels
  }

  marker_temp2
}
