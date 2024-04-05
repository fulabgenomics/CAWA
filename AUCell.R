counts <- N2all@assays$RNA@counts 
metadata <- N2all@meta.data
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
exprMatrix <- assay(sce)
exprMatrix <- as(exprMatrix, 'dgCMatrix')
half <- round(ncol(exprMatrix)/2)
ExprMat_1 <- exprMatrix[,1:half]
ExprMat_2 <- exprMatrix[,(half+1):ncol(exprMatrix)]
cells_rankings_1 <- AUCell_buildRankings(ExprMat_1)
cells_rankings_2 <- AUCell_buildRankings(ExprMat_2)
cells_rankings <- AUCell::cbind(cells_rankings_1, cells_rankings_2)
saveRDS(cells_rankings,'cells_rankings.N2.RDS')
my_txt <- readLines('TF2DNA_2018.txt')
data.frame(tissue=N2all$tissue)->out1
for (i in 1:length(my_txt)){
  unlist(strsplit(my_txt[i], '	'))->tmp
  tmp[1]->setname
  output_string <- strsplit(setname,'_')[[1]][1]
  print(output_string)
  genes<-tmp[3:length(tmp)]
  geneSets <- list(setname=genes)
  names(geneSets)<-c(setname)
  
  AuCell.Signature.all <- AUCell_calcAUC(geneSets = geneSets,
                                         rankings =cells_rankings,
                                         nCores = 4,
                                         aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
  getAUC(AuCell.Signature.all)->AucellScore
  data.frame(t(AucellScore))->auout
  merge(auout, out1, by.x = 0, by.y = 0)->out1
  row.names(out1)<-out1$Row.names
  out1<-subset(out1, select = -c(Row.names))
}
