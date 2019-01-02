#' OVESEG test perm
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
OVEtstatPermTopM <- function(y, group, groupOrder, M, weights=NULL,
                        alpha='moderated', NumPerm=999, seed=111,
                        BPPARAM=SerialParam()) {

    K <- length(unique(group))
    if (K < M)
        stop("At least M groups are needed.")

    if (length(seed) == 1) {
        seeds <- unlist(lapply(seq_len(NumPerm), function(x) seed+x-1))
    } else if (length(seed) == NumPerm){
        seeds <- seed
    } else {
        stop('the length of seed should be one or the same as NumPerm!')
    }

    group <- factor(group)

    combM <- t(combn(K, M))
    ncombM <- nrow(combM)
    geneSubset <- integer(nrow(y))
    for(j in 1:ncombM) {
        geneidx <- apply(groupOrder[,1:M], 1,
                        function(x) setequal(x, combM[j,]))
        geneSubset[geneidx] <- j
    }



    permfunc <- function(p) {
        y.p <- y
        weights.p <- weights
        if(p > 1){
            set.seed(seeds[p-1])
            for(j in 1:ncombM) {
                gidx <- group %in% levels(group)[combM[j,]]
                sidx <- seq_along(group)
                sidx[gidx] <- sample(sidx[gidx])
                geneidx <- which(geneSubset==j)
                y.p[geneidx,] <- y.p[geneidx, sidx]
                if(!is.null(weights)) {
                    weights.p[geneidx,] <- weights.p[geneidx, sidx]
                }
            }
        }

        ovesegt <- OVESEGtstat(y.p, group, weights = weights.p, alpha = alpha,
                               order.return = 'top1', lmfit.return = FALSE)

        return(ovesegt)
    }



    NPerm <- NumPerm + 1

    if (is.null(BPPARAM)) {
        oveseg.perm <- lapply(seq_len(NPerm), permfunc)
    } else {
        oveseg.perm <- bplapply(seq_len(NPerm), permfunc,
                             BPPARAM = BPPARAM)
    }

    tstat.perm <- vapply(oveseg.perm, function(x) x$tstat,
                                FUN.VALUE= numeric(nrow(y)))

    topidx.perm <- vapply(oveseg.perm, function(x) x$groupOrder,
                                FUN.VALUE= integer(nrow(y)))


    return(list(tstat.perm=tstat.perm, topidx.perm=topidx.perm))
}
