#' KL divergence
#'
#' @description Use package entropy to compute Kullback-Leibler divergence.
#' The function first converts each vector's reads to pseudo-number of
#' transcripts by normalizing the total reads to total_reads.
#' The normalized read for each gene is then rounded to serve as the
#' pseudo-number of transcripts.
#' Function [entropy::KL.shrink()] is called to compute the KL-divergence between
#' the two vectors, and the maximal allowed divergence is set to max_KL.
#' Finally, a linear transform is performed to convert the KL divergence,
#' which is between 0 and max_KL, to a similarity score between -1 and 1.
#'
#' @param vec1 Test vector
#' @param vec2 Reference vector
#' @param if_logcounts Whether the vectors are log-transformed. If so, the
#' raw count should be computed before computing KL-divergence.
#' @param total_reads Pseudo-library size
#' @param max_KL Maximal allowed value of KL-divergence.
#'
#' @importFrom entropy KL.shrink
#' @export
kl_divergence <- function(vec1, vec2, if_logcounts = FALSE,
                          total_reads = 1000, max_KL = 1) {
  if (if_logcounts) {
    vec1 <- expm1(vec1)
    vec2 <- expm1(vec2)
  }
  count1 <- round(vec1 * total_reads / sum(vec1))
  count2 <- round(vec2 * total_reads / sum(vec2))
  est_KL <- KL.shrink(count1, count2,
    unit = "log",
    verbose = FALSE
  )
  return((max_KL - est_KL) / max_KL * 2 - 1)
}
