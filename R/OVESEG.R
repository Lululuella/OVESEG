#' OVESEG: A package for marker gene test.
#'
#' Function \code{\link{OVESEGtest}} performs OVESEG-test for expression
#' profiles from multiple groups to detect subtype-specific marker genes.
#' While it may take a long time to execute permutations for p-value estimation,
#' users can apply \code{\link{OVESEGtstat}} to obtain OVESEG-test statistics
#' to rank genes and apply \code{\link{postProbNull}} to obtain each gene's
#' posterior probabilities of component null hypotheses.
#' \code{\link{nullDistri}} estimates probabilities of any one group being
#' upregulated under null hypotheses.
#' \code{\link{patternDistri}} estimates probabilities of all kinds of
#' upregulation patterns among groups.
#'
#' @references Chen, L., Herrington, D., Clarke, R., Yu, G., and Wang, Y.
#' (2019). “Data-Driven Robust Detection of Tissue/Cell-Specific Markers.”
#' bioRxiv. https://doi.org/10.1101/517961.
#'
#' @docType package
#' @name OVESEG-package
#' @aliases OVESEG
#'
#' @import BiocParallel
#' @importFrom methods is
#' @importFrom stats pf
#' @importFrom stats pchisq
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom fdrtool fdrtool
#' @importFrom Rcpp sourceCpp
#' @useDynLib OVESEG
NULL
