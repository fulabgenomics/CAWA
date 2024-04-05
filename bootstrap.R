library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
set.seed(123)
input_cel <- readRDS('../data/CAWA.rds')
tissues<-c('Neuron','Muscle','Intestine','Hypodermis','Germline','Vulva_uterus','Pharynx')
levels(as.factor(input_cel[[]]$orig.ident))->batches
for (tissueinput in tissues) {
  tissue_sub = subset(input_cel, subset = tissue == tissueinput)
  for (batchinput in batches) {
    batch_sub = subset(tissue_sub, subset = orig.ident == batchinput)
    counts <- batch_sub@assays$RNA@counts
    metadata <- batch_sub@meta.data
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = metadata)
    groups <- colData(sce)[, c('orig.ident')]
    output_data<-c()
    for (i in 1:100) {
      group1 <- rep('neg', length(groups))
      group1[sample(c(1:length(groups)), size = 15, replace = FALSE)] <-'pos'
      pb2 <- Matrix.utils:::aggregate.Matrix(t(counts(sce)),
                                             groupings = group1, fun = 'sum')
      cbind(output_data, pb2['pos', ]) -> output_data
    }
    saveRDS(output_data,
            paste0('bootstrapcells/', tissueinput, '_', batchinput, '.RDS'))
  }
}
