#' Probabilities of all upregulation patterns
#'
#' This function estimates probabilities of all kinds of upregulation patterns
#' among subtypes.
#' @param ppnull a list returned by \code{\link{postProbNull}} or
#'     \code{\link{OVESEGtest}}.
#' @details The probability of each upregulation pattern is calculated by
#' accumulating and normalizing genewise posterior probability of null
#' hypotheses and of alternative hypotheses.
#' @return a data frame object containing all possible upregulation patterns
#' and corresponding probabilities.
#' @export
#' @examples
#' data(RocheBT)
#' ppnull <- postProbNull(RocheBT$y, RocheBT$group, alpha='moderated')
#' pd<- patternDistri(ppnull)
patternDistri <- function(ppnull)
{
    ngene <- nrow(ppnull$W)
    K <- ncol(ppnull$W) + 1

    W0 <- 1 - rowSums(ppnull$W)
    compWAll <- c()
    pattern <- c()
    for (M in seq_len(K)) {
        combM <- combn(K,M)
        ncombM <- ncol(combM)
        geneSubset <- integer(ngene)
        for (j in seq_len(ncombM)) {
            geneidx <- apply(ppnull$groupOrder[,seq_len(M), drop=FALSE], 1,
                            function(x) setequal(x, combM[,j]))
            geneSubset[geneidx] <- j
        }
        if (M == 1){
            compW <- unlist(lapply(seq_len(ncombM),
                            function(x) sum(W0[geneSubset==x])))
        } else {
            compW <- unlist(lapply(seq_len(ncombM),
                            function(x) sum(ppnull$W[geneSubset==x, M-1])))
        }
        compWAll <- c(compWAll, compW)

        for (j in seq_len(ncombM)) {
            cellc <- rep(0, K)
            cellc[c(combM[,j])] <- 1
            pattern <- rbind(pattern, cellc)
        }
    }
    Wpattern <- compWAll/sum(compWAll)
    colnames(pattern) <- ppnull$label
    rownames(pattern) <- NULL

    return(data.frame(pattern, Wpattern=Wpattern))
}
