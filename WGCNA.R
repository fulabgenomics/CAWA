for (genoinput in unique(genotypes) ){

  genotype_tissue <- subset(hd_in,subset = genotype_wg %in% genoinput)
  genotype_tissue <- SetupForWGCNA(
    genotype_tissue,
    gene_select = 'fraction',
    fraction = 0.05,
    wgcna_name = 'tissue')
  
  genotype_tissue <- MetacellsByGroups(
    seurat_obj = genotype_tissue,
    k = 25, 
    max_shared = 10)
  
  genotype_tissue <- NormalizeMetacells(genotype_tissue)
  genotype_tissue <- SetDatExpr(genotype_tissue,group_name =genoinput ,assay = 'RNA', slot = 'data')
  genotype_tissue <- TestSoftPowers(genotype_tissue)
  
  genotype_tissue <- ConstructNetwork(genotype_tissue, tom_name = paste0(tissueinput,'_',genoinput,'_seperate'),overwrite_tom = TRUE)
  genotype_tissue <- ModuleEigengenes(genotype_tissue)
  genotype_tissue <- ModuleConnectivity(genotype_tissue)
  vector2=genotype_tissue[[]]$orig.ident
  vector=ages[vector2]
  genotype_tissue<-AddMetaData(object = genotype_tissue, metadata = vector, col.name = 'ages')
  genotype_tissue$ages <- as.numeric(genotype_tissue$ages)
  genotype_tissue <- ModuleTraitCorrelation(genotype_tissue,traits = 'ages')
  mt_cor <- GetModuleTraitCorrelation(genotype_tissue)
  mt_cor$fdr$all_cells<0.01->positive_module
  data.frame(mt_cor$cor$all_cells)->correlationout
  correlationout[positive_module,]->correlationout_po
  names(correlationout_po)<-row.names(correlationout)[positive_module]
  GetModules(genotype_tissue)->gene_out
  saveRDS(gene_out,paste0('seperate/',tissueinput,'_',genoinput,'_raw.RDS'))
  gene_out[gene_out$module %in% names(correlationout_po),]->filtered
  filtered$up<-(correlationout_po>0)[as.character(filtered$module)]
  saveRDS(filtered,paste0('seperate/',tissueinput,'_',genoinput,'.RDS'))
}}
