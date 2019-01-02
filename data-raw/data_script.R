library(GEOquery)
gsm <- getGEO("GSE28490")

y <- exprs(gsm[[1]])
group <- pData(phenoData(gsm[[1]]))$characteristics_ch1


idx <- which(group %in% c('cell type: B cells', 'cell type: CD4+ T cells',
                          'cell type: CD8+ T cells'))

set.seed((111))
y <- y[sample(seq_len(nrow(y)), 5000),idx]

group <- group[idx]
group <- as.character(group)
group <- unlist(lapply(strsplit(group, ': '), function(x) x[2]))
group <- factor(group)

RocheBT <- list(y=y, group=group)
devtools::use_data(RocheBT)
