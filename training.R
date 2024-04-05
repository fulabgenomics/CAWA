for (tissueinput in tissues) {
  training_x <- c()
  training_y <- c()
    for (batch_training in batches) {
      x <-readRDS(paste0('bootstrapcells/',tissueinput,'_',batch_training,'.RDS'))
      log1p(t(x) * 100000 / colSums(x)) -> training_x_in
      training_x <- rbind(training_x, training_x_in)
      y <- yforbatches[batch_training]
      training_y_in <- t(t(rep(y, 100)))
      training_y <- rbind(training_y, training_y_in)
    }
    model <- glmnet::cv.glmnet(
      x = training_x,
      y = training_y,
      type.measure = 'mae',
      standardize = F,
      relax = F,
      nfolds = 5
    )
    coef(model,s = 'lambda.min')-> coefout
    data.frame(coefout[,1])-> selectfeature
    data.frame(gene=rownames(x),value=selectfeature[2:nrow(selectfeature),])-> out
    out[out$value!=0,]->geneout
    saveRDS(model,paste0('training_predict/', tissueinput, '_model.RDS'))
    write.csv(geneout, paste0('training_predict/',tissueinput,'_coef.csv'),row.names = FALSE,quote=FALSE)
}
