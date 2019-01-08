#' mRNA expression data downsampled from GSE28490 (Roche)
#'
#' Three cell subtypes (B cells, CD4+ T cells, CD8+ T cells) were isolated from
#' 5 pools of 5 healthy donors each. RNA obtained from these 15 purified
#' populations were subsequently used for mRNA expression profiling by
#' HG-U133Plus2.0 microarrays. We downsample the original data to 5000
#' probes/probesets.
#' Subtype labels for purified populations are also included.
#' (Data generation script can be found in ./data_row folder.)
#' @docType data
#' @name RocheBT
#' @usage data(RocheBT)
#' @format A list with one expression mixture (y) and a categorical vector
#' giving subtypes (group).
#' @references Allantaz et al. PLoS One 2012;7(1):e29979. PMID: 22276136
NULL
