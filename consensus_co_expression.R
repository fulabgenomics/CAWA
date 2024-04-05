for (tissueinput in tissues) {
  tissu_seurat<-subset(input_cel,subset = tissue == tissueinput)
  hd_in <- subset(tissu_seurat,subset = orig.ident %in% identsin)
  vector2=hd_in[[]]$orig.ident
  vector=genotypes[vector2]
  hd_in<-AddMetaData(object = hd_in, metadata = vector, col.name = 'genotype_wg')
  hd_in <- SetupForWGCNA(hd_in,
  gene_select = 'fraction,
  fraction = 0.05,
  wgcna_name = 'tissue')

  hd_in <- MetacellsByGroups(
    seurat_obj = hd_in,
    group.by = c('genotype_wg'), 
    k = 25, 
    max_shared = 10,
    ident.group = 'genotype_wg'
  )

  hd_in <- NormalizeMetacells(hd_in)
  
  hd_in <- SetMultiExpr(hd_in,
                            multi.group.by='genotype_wg', 
                            assay = 'RNA', 
                            slot = 'data')
  
  hd_in <- TestSoftPowersConsensus(hd_in, networkType = 'signed')
  
  
  hd_in <- ConstructNetwork(
    hd_in, 
    consensus=TRUE,
    tom_name = paste0(tissueinput,'_Genotype_Consensus'),
    overwrite_tom = TRUE
  )

  hd_in <- ModuleEigengenes(hd_in)
  hd_in <- ModuleConnectivity(hd_in)
  

  #module trait correlation
  vector2=hd_in[[]]$orig.ident
  vector=ages[vector2]
  hd_in<-AddMetaData(object = hd_in, metadata = vector, col.name = 'ages')
  hd_in$ages <- as.numeric(hd_in$ages)
  hd_in <- ModuleTraitCorrelation(
    hd_in,
    traits = 'ages',
    group.by='genotype_wg')
  mt_cor <- GetModuleTraitCorrelation(hd_in)

  saveRDS(hd_in,paste0('network/', tissueinput, '_',  'consensus.RDS'))
