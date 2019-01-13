library(GEOquery)

####RocheBT
gsm <- getGEO("GSE28490")

y <- exprs(gsm[[1]])
group <- pData(phenoData(gsm[[1]]))$characteristics_ch1


idx <- which(group %in% c('cell type: B cells', 'cell type: CD4+ T cells',
                          'cell type: CD8+ T cells'))

set.seed((111))
y <- y[sample(seq_len(nrow(y)), 5000),idx]

group <- group[idx]
group <- as.character(group)
group[group=='cell type: B cells'] <- 'B'
group[group=='cell type: CD4+ T cells'] <- 'CD4'
group[group=='cell type: CD8+ T cells'] <- 'CD8'
group <- factor(group)

RocheBT <- list(y=y, group=group)
devtools::use_data(RocheBT)



####countBT
gsm <- getGEO("GSE60424")

mat <- as.matrix(read.table(
    gzfile("../GSE60424_GEOSubmit_FC1to11_normalized_counts.txt.gz"),
    header=T,row.names=1))

group <- pData(phenoData(gsm[[1]]))$characteristics_ch1.2
idx <- which(group %in% c('celltype: B-cells', 'celltype: CD4',
                          'celltype: CD8'))

set.seed((111))
mat <- mat[sample(seq_len(nrow(mat)), 10000),idx]

group <- group[idx]
group <- as.character(group)
group[group=='celltype: B-cells'] <- 'B'
group[group=='celltype: CD4'] <- 'CD4'
group[group=='celltype: CD8'] <- 'CD8'
group <- factor(group)

countBT <- list(count=mat, group=group)
devtools::use_data(countBT)
