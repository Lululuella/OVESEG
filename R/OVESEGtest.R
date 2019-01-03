#' OVESEG test
#'
#' This function reduces data dimension by loading matrix and then project
#' dimension-reduced data to the hyperplane orthogonal to c(1,0,...,0),
#' i.e., the first axis in the new coordinate system..
#' @param y A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample. Missing values
#'     are not supported.
#'     All-zero rows will be removed internally.
#' @param group The matrix whose rows are loading vectors;
#' @details This function can project gene expression vectors to simplex
#' filtering. This function helps observe all genes in simplex plot.
#' @return  The data after perspective projection.
#' @export
#' @examples
#' #obtain data
#' data(RocheBT)
#' data <- RocheBT$y
#'
#' #preprocess data
OVESEGtest <- function(y, group, weights = NULL, alpha = 'moderated',
                       NumPerm = 999, seed = 111, detail.return=TRUE, cplusplus = TRUE,
                       BPPARAM=SerialParam()){
    K <- length(unique(group))
    if (K < 2) {
        stop("At least two groups are needed.")
    }
    #if (table(group))


    group <- factor(group)

    ppnull <- posteriorProbNull(y, group, weights = weights, alpha = alpha,
                                  detail.return = detail.return, cplusplus=cplusplus)

    ## permutation among top M groups
    tstat.perm <- c()
    topidx.perm <- c()
    for(M in seq.int(2,K)) {

        ovet.perm <- OVEtstatPermTopM(y, group, ppnull$groupOrder, M=M, weights = weights, alpha = alpha,
                         NumPerm = NumPerm, seed = seed, cplusplus=cplusplus, BPPARAM = BPPARAM)

        tstat.perm <- cbind(tstat.perm, ovet.perm$tstat.perm)
        topidx.perm <- cbind(topidx.perm, ovet.perm$topidx.perm)
    }

    wPerm <- ppnull$W[,rep(1:(K-1),each=NumPerm+1)]

    testing <- tstat.perm[,1]
    pv.overall <- pvalueWeightedEst(testing, tstat.perm, wPerm)
    # pv.oneside <- rep(0, nrow(y))
    # pv.oneside.max <- rep(0,K)
    # for(M in seq_len(K)) {
    #     idx <-
    #     tstat.perm.oneside <- tstat.perm
    #     tstat.perm.oneside[topidx.perm!=M] <- -1
    #     pvtmp <- pvalueWeightedEst(c(0,testing[idx]), tstat.perm.oneside, wPerm)
    #     pv.oneside[idx] <- pvtmp[-1]
    #     pv.oneside.max[M] <- pvtmp[1]
    # }

    gidx <- lapply(seq_len(K), function(M) which(ppnull$groupOrder[,1]==M))

    pvforOne <- function(M) {
        tstat.perm.oneside <- tstat.perm
        tstat.perm.oneside[topidx.perm!=M] <- -1
        pv <- pvalueWeightedEst(c(0,testing[gidx[[M]]]), tstat.perm.oneside, wPerm)
        pv
    }
    if (is.null(BPPARAM)) {
        pvK <- lapply(seq_len(K), pvforOne)
    } else {
        pvK <- bplapply(seq_len(K), pvforOne,
                                BPPARAM = BPPARAM)
    }
    pv.oneside <- pv.overall
    pv.oneside.max <- rep(0,K)
    for(M in seq_len(K)) {
        pv.oneside[gidx[[M]]] <- pvK[[M]][-1]
        pv.oneside.max[M] <- pvK[[M]][1]
    }

    minpv <- min(rowSums(ppnull$W)) / (NumPerm + 1) / sum(ppnull$W)
    # pv.multiside <- unlist(lapply(pv.oneside*K,
    #                             function(x) max(minpv, min(1, x))))

    pv.multiside <- pv.oneside*K
    pv.multiside[pv.multiside<0] <- 0
    pv.multiside[pv.multiside>1] <- 1



    return(c(list(pv.oneside=pv.oneside, pv.oneside.max=pv.oneside.max,
                pv.multiside=pv.multiside,
                pv.overall=pv.overall),ppnull))

}
