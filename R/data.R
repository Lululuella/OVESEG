#' mRNA expression data downsampled from GSE28490 (Roche)
#'
#' Three cell subtypes (B cells, CD4+ T cells, CD8+ T cells) were isolated from
#' 5 pools of 5 healthy donors each. RNA obtained from these 15 purified
#' populations were subsequently used for mRNA expression profiling by
#' HG-U133Plus2.0 microarrays. We downsample the original data to 5000
#' probes/probesets.
#' Subtype labels for purified populations are also included.
#' (Data generation script can be found in ./data_raw folder.)
#' @docType data
#' @name RocheBT
#' @usage data(RocheBT)
#' @format A list with one expression matrix (y) and a categorical vector
#' giving subtypes (group).
#' @references Allantaz et al. PLoS One 2012;7(1):e29979. PMID: 22276136
NULL



#' RNAseq count data downsampled from GSE60424
#'
#' Three cell subtypes (B cells, CD4+ T cells, CD8+ T cells) were purified from
#' 20 fresh blood samples. RNA was extracted from each of these cell subsets
#' and processed into RNA sequencing libraries (Illumina TruSeq). Sequencing
#' libraries were analyzed on an Illumina HiScan, with a target read depth
#' of ~20M reads. Reads were demultiplexed, mapped to human gene
#' models (ENSEMBL), and tabulated using HTSeq.
#' We downsample the original data to 10000 genes.
#' Subtype labels for purified populations are also included.
#' (Data generation script can be found in ./data_raw folder.)
#' @docType data
#' @name countBT
#' @usage data(countBT)
#' @format A list with one count mixture (count) and a categorical vector
#' giving subtypes (group).
#' @references Linsley et al. PLoS One 2014;9(10):e109760. PMID: 25314013
NULL
