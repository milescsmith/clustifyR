#' dl_recount
#'
#' @param sra_id ID of SRA study to download
#'
#' @importFrom recount download_study scale_counts
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom dplyr mutate select mutate_at left_join group_by summarise filter
#' @importFrom stringr str_match str_c
#' @importFrom purrr map_chr
#' @importFrom tidyr separate spread gather
#' @importFrom tibble data_frame rownames_to_column column_to_rownames
#'
#' @return
#' @export dl_recount
#'
#' @examples
dl_recount <- function(sra_id){
  download_study(sra_id)
  load(file.path(sra_id, "rse_gene.Rdata"))
  # no longer need to downloaded data
  unlink(sra_id, recursive = TRUE)

  rse <- scale_counts(rse_gene)
  read_counts <- assay(rse, "counts")
  gene_ids <- rownames(read_counts)
  # get gene symbols, which are stored in rowData
  id2symbol <- data_frame(ids = rowData(rse_gene)$gene_id,
             symbols = rowData(rse_gene)$symbol@listData) %>%
    mutate(symbols = map_chr(symbols, ~.x[1]))

  # clean up metadata into a dataframe
  mdata <- colData(rse)
  mdata_cols <- lapply(mdata$characteristics,
      function(x){str_match(x, "^([^:]+):")[, 2]}) %>%
        unique() %>%
        unlist()

  mdata <- data_frame(run =  mdata$run,
                      all_data = as.list(mdata$characteristics)) %>%
    mutate(out = purrr::map_chr(all_data,
                                ~str_c(.x, collapse = "::"))) %>%
  separate(out,
           sep = "::",
           into = mdata_cols) %>%
    select(-all_data) %>%
    mutate_at(.vars = vars(-matches("run")),
              .funs = function(x) str_match(x, ": (.+)")[, 2])

  # convert ids to symbols
  row_ids_to_symbols <- left_join(data_frame(ids = gene_ids),
            id2symbol, by = "ids")

  if(length(gene_ids) != nrow(row_ids_to_symbols)) {
    warning("gene id mapping to symbols produce more or less ids")
  }

  row_ids_to_symbols <- filter(row_ids_to_symbols, !is.na(symbols))

  out_df <- read_counts %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    left_join(., row_ids_to_symbols,
              by = c("gene_id" = "ids")) %>%
    select(-gene_id) %>%
    select(symbols, everything()) %>%
    filter(!is.na(symbols))

  out_matrix <- gather(out_df, library, expr, -symbols) %>%
    group_by(symbols, library) %>%
    summarise(expr = sum(expr)) %>%
    spread(library, expr) %>%
    as.data.frame() %>%
    column_to_rownames("symbols") %>%
    as.matrix()

  list(readcounts = out_matrix,
       metadata = mdata)
}


