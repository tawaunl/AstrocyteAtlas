require(ggplot2)
require(Seurat)


plot_reduc_library <- function(sobj,meta,sample,reduction){
  scvi_coords <- Embeddings(sobj[[reduction]])[,1:2]
  colnames(scvi_coords) <- c("UMAP1","UMAP2")
  scvi_coords <- as.data.frame(cbind(scvi_coords,sobj@meta.data))
  
  p <- ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2)) +
    geom_point(color="grey", size=1,alpha = 1) +
    geom_point(data = scvi_coords %>% dplyr::filter((!!sym(meta))==sample), color="red", size=0.5, alpha=0.5) +
    theme_classic() + ggtitle(paste(meta,"=",sample)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 24, face = "bold")) + xlab(paste0(toupper(reduction),"_1")) +
    ylab(paste0(toupper(reduction),"_2"))
  
  return(p)
  
}

plot_umap_library <- function(sobj,meta,sample,color="red"){
  scvi_coords <- Embeddings(sobj[["umap"]])
  colnames(scvi_coords) <- c("UMAP1","UMAP2")
  scvi_coords <- as.data.frame(cbind(scvi_coords,sobj@meta.data))
  
  p <- ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2)) +
    ggrastr::geom_point_rast(color="grey", size=1,alpha = 1) +
    ggrastr::geom_point_rast(data = scvi_coords %>% dplyr::filter((!!sym(meta))==sample), color=color, size=0.5, alpha=0.5) +
    theme_classic() + ggtitle(paste(sample)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 18))
  
  return(p)
  
}

plot_pca_library <- function(sobj,meta,sample){
  scvi_coords <- Embeddings(sobj[["umap"]])
  colnames(scvi_coords) <- c("UMAP1","UMAP2")
  scvi_coords <- as.data.frame(cbind(scvi_coords,sobj@meta.data))
  
  p <- ggplot(data=scvi_coords, aes(x=UMAP1, y=UMAP2)) +
    geom_point(color="grey", size=1,alpha = 1) +
    geom_point(data = scvi_coords %>% dplyr::filter((!!sym(meta))==sample), color="red", size=0.5, alpha=0.5) +
    theme_classic() + ggtitle(paste(meta,"=",sample)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 24, face = "bold"))
  
  return(p)
  
}

plot_harmony_library <- function(sobj,meta,sample){
  scvi_coords <- Embeddings(sobj[["harmony"]])
  scvi_coords <- as.data.frame(cbind(scvi_coords,sobj@meta.data))
  
  p <- ggplot(data=scvi_coords, aes(x=harmony_1, y=harmony_2)) +
    geom_point(color="grey", size=1,alpha = 1) +
    geom_point(data = scvi_coords %>% dplyr::filter((!!sym(meta))==sample), color="red", size=0.5, alpha=0.5) +
    theme_classic() + ggtitle(paste(meta,"=",sample)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 24, face = "bold"))
  
  return(p)
  
}

get_scvi_coords <- function(se,meta) {
  scvi_coords <- Embeddings(se[["umap"]])
  colnames(scvi_coords) <- c("UMAP1","UMAP2")
  scvi_coords <- as.data.frame(cbind(scvi_coords,se@meta.data))
  
  results <- fastDummies::dummy_cols(meta)
  colnames(results) <- unlist(lapply(colnames(results),function(x) strsplit(x,".data_")[[1]][2]))
  results <- results[,-1]
  scvi_coords <- cbind(scvi_coords, results)
  
  return(scvi_coords)
  
}

get_dimred_coords <- function(se,meta,reduc="umap") {
  scvi_coords <- Embeddings(se[[reduc]])
  colnames(scvi_coords) <- c(paste0(toupper(reduc),"1"),paste0(toupper(reduc),"2"))
  scvi_coords <- as.data.frame(cbind(scvi_coords,se@meta.data))
  
  results <- fastDummies::dummy_cols(meta)
  colnames(results) <- unlist(lapply(colnames(results),function(x) strsplit(x,".data_")[[1]][2]))
  results <- results[,-1]
  scvi_coords <- cbind(scvi_coords, results)
  
  return(scvi_coords)
  
}
runSeurat <- function(obj,res=0.5){
  obj <- NormalizeData(obj,verbose=FALSE)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(object = obj, verbose=FALSE)
  obj <- RunPCA(obj, npcs = 30,verbose=FALSE)
  obj <- FindNeighbors(obj, reduction = "pca",dims=1:30,verbose=FALSE)
  obj <- FindClusters(obj, resolution = res,verbose=FALSE)
  obj <- RunUMAP(obj,reduction= "pca",dims=1:30,verbose=FALSE)
  print("Seurat Pipeline Finished!")
  return(obj)
}


edgeRVolcano <- function(res=res, FCcutoff=1,
                         PvalCutoff = 0.05,plot.title="VolcanoPlot",
                         labels=NULL,legend="none") {
  
  require(ggrepel)
  res <- res %>% mutate(sig = ifelse((PValue <= PvalCutoff | PValue <= PvalCutoff) & abs(logFC) >= 1, "yes", "no"))
  
  
  res$diffexpressed <- "unchanged"
  
  res$diffexpressed[res$logFC >= FCcutoff&
                      res$PValue<=PvalCutoff &
                      res$sig=="yes" ] <- "up"
  
  res$diffexpressed[res$logFC <= -1 * FCcutoff &
                      res$PValue <=0.05 &
                      res$sig=="yes" ] <- "down"
  
  mycolors <- c("red", "dodgerblue", "black")
  names(mycolors) <- c("down", "up", "unchanged")
  
  p <-  ggplot(data=res, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(res))) +
    geom_point(size=2) +
    theme_classic()+
    scale_colour_manual(values = mycolors) + ggtitle(plot.title)+
    theme(legend.position=legend,legend.text = element_text(size=16,face = "bold"),
          legend.title = element_text(size=18,face='bold'),
          legend.key.size = unit(1 ,'cm'),
          axis.title = element_text(size=20,face='bold'),
          plot.title = element_text(size=24,face='bold',hjust = 0.5)) +
    xlab("Avg. LogFoldChange") + ylab(bquote(bold(-log[10](PValue))))
  
  if(is.null(labels)){
    return(p)
  }else{
    p <- p + geom_label_repel(data= subset(res, rownames(res) %in% labels),
                              aes(x=logFC, y=-log10(PValue),
                                  label=rownames(subset(res, rownames(res) %in% labels)),
                                  colour="lightblue")) +
      guides(fill = guide_legend(override.aes = aes(label = "")))
    return(p)
  }
}

VolcanoPlot <- function(res=res, FCcutoff=1,lfc = "logFC",Pval="pval",
                        PvalCutoff = 0.05,plot.title="VolcanoPlot",
                        labels=NULL,legend="none",NScolor="grey",
                        upcolor="dodgerblue",downcolor="red",labsize=6) {
  if(rownames(res)[1]=="1"){
    stop("Looks like genenames are not on rows, be sure to set rownames to gene names")
  }
  colnames(res)[which(colnames(res)==Pval)] <- "PValue"
  colnames(res)[which(colnames(res)==lfc)] <- "logFC"
  res$logFC <- as.numeric(res$logFC)
  res$PValue <- as.numeric(res$PValue)
  require(ggrepel)
  res <- res %>% mutate(sig = ifelse((PValue <= PvalCutoff) & abs(logFC) >= FCcutoff, "yes", "no"))
  
  
  res$diffexpressed <- "unchanged"
  
  res$diffexpressed[res$logFC >= FCcutoff&
                      res$PValue<=PvalCutoff &
                      res$sig=="yes" ] <- "up"
  
  res$diffexpressed[res$logFC <= -1 * FCcutoff &
                      res$PValue <=0.05 &
                      res$sig=="yes" ] <- "down"
  
  mycolors <- c(downcolor, upcolor, NScolor)
  names(mycolors) <- c("down", "up", "unchanged")
  
  p <-  ggplot(data=res, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(res))) +
    ggrastr::geom_point_rast(size=2) +
    theme_classic()+
    scale_colour_manual(values = mycolors) + ggtitle(plot.title) +
    theme(legend.position=legend,legend.text = element_text(size=16,face = "bold"),
          legend.title = element_text(size=18,face='bold'),
          legend.key.size = unit(1 ,'cm'),
          axis.title = element_text(size=20,face='bold'),
          plot.title = element_text(size=24,face='bold',hjust = 0.5)) +
    xlab("Avg. LogFoldChange") + ylab(bquote(bold(-log[10](PValue))))
  
  if(is.null(labels)){
    return(p)
  }else{
    res$lab <- rownames(res)
    if (!is.null(labels)) {
      names.new <- rep(NA, length(res$lab))
      indices <- which(res$lab %in% labels)
      names.new[indices] <- res$lab[indices]
      res$lab <- names.new
    }
    p <- p + geom_text_repel(data= subset(res, rownames(res) %in% labels),
                             aes(x=logFC, y=-log10(PValue),
                                 label=rownames(subset(res, rownames(res) %in% labels))),
                             size=labsize,max.overlaps = 25,min.segment.length = 0,na.rm = T) +
      guides(fill = guide_legend(override.aes = aes(label = "")))
    return(p)
  }
}


PCA_project <- function(ref.obj,query.obj,transferField="finalClusters"){
  require(stringr)
  transfer.integrated  <- FindTransferAnchors(ref.obj,query = query.obj,
                                              features = VariableFeatures(ref.obj))
  
  anchors <- FindTransferAnchors(reference = ref.obj, query = query.obj,
                                 dims = 1:30,
                                 reference.reduction = "pca")
  field <- t(ref.obj[[transferField]])
  nn<- colnames(field)
  field <- factor(field)
  names(field) <- nn
  predictions <- TransferData(anchorset = anchors, refdata =field,
                              dims = 1:30)
  
  query.obj <- AddMetaData(query.obj, metadata = predictions)
  
  #project into integrated reference UMAP
  query.obj <- MapQuery(anchorset = anchors, reference = ref.obj, query = query.obj,
                        refdata = field,
                        reference.reduction = "umap", reduction.model = "umap")
}

generatePlotTable <- function(res_dsl){
  require(dplyr)
  plotdata_up <- res_dsl %>% filter(LogFC > 0)
  res_dsl_plot <- data.frame(ID=plotdata_up$ID,
                             symbol=plotdata_up$symbol,
                             dl_mu=plotdata_up$dl_mu,
                             LogFC = plotdata_up$LogFC,
                             PValue=plotdata_up$metap_up,
                             FDR=plotdata_up$adj_metap_up,
                             AveExpr=plotdata_up$AveExpr,
                             n_tested=plotdata_up$n_tested,
                             n_up=plotdata_up$n_up,
                             n_down=plotdata_up$n_down)
  plotdata_down <- res_dsl %>% filter(LogFC < 0)
  res_dsl_plot <- rbind(res_dsl_plot,data.frame(ID=plotdata_down$ID,
                                                symbol=plotdata_down$symbol,
                                                dl_mu=plotdata_down$dl_mu,
                                                LogFC = plotdata_down$LogFC,
                                                PValue=plotdata_down$metap_down,
                                                FDR=plotdata_down$adj_metap_down,
                                                AveExpr=plotdata_down$AveExpr,
                                                n_tested=plotdata_down$n_tested,
                                                n_up=plotdata_down$n_up,
                                                n_down=plotdata_down$n_down))
  res_dsl_plot <- res_dsl_plot %>% mutate(sig = ifelse((FDR<=0.05 | FDR <= 0.05) & n_up | n_down >= n_tested/2 & abs(dl_mu) >= 1, "yes", "no"))
  return(res_dsl_plot)
}


df4way <- function(df1,df2,...){
  
  
}


get_harmony_coords <- function(se, harmony_umap_dim) {
  message(paste(names(reducedDims(se))[harmony_umap_dim]))
  harmony_coords_final <- reducedDims(se)[[harmony_umap_dim]]
  colnames(harmony_coords_final) <- c("UMAP1","UMAP2")
  harmony_coords_final <- as.data.frame(cbind(harmony_coords_final,colData(se)))
  
  results <- fastDummies::dummy_cols(se$studybatch)
  colnames(results) <- unlist(lapply(colnames(results),function(x) strsplit(x,".data_")[[1]][2]))
  results <- results[,-1]
  harmony_coords_final <- cbind(harmony_coords_final, results)
  
  return(harmony_coords_final)
  
}


plot_4way <- function(res1, res2, disease1, disease2,lfc="dl_mu",pval="Pvalue",FCcutoff=0.5) {
  require(tidyr)
  require(ggrepel)
  res1 <- data.frame(res1)
  res2 <- data.frame(res2)
  if(!"symbol" %in% colnames(res1)){
    res1$symbol <- rownames(res1)
  }
  if(!"symbol" %in% colnames(res2)){
    res2$symbol <- rownames(res2)
  }
  
  colnames(res1)[which(colnames(res1)==pval)] <- "PValue"
  colnames(res1)[which(colnames(res1)==lfc)] <- "dl_mu"
  colnames(res2)[which(colnames(res2)==pval)] <- "PValue"
  colnames(res2)[which(colnames(res2)==lfc)] <- "dl_mu"
  res1$dl_mu <- as.numeric(res1$dl_mu)
  res1$PValue <- as.numeric(res1$PValue)
  res2$dl_mu <- as.numeric(res2$dl_mu)
  res2$PValue <- as.numeric(res2$PValue)
  res1 <- as.data.frame(res1) %>% mutate(sig = ifelse((PValue<=0.05 | PValue <= 0.05) & abs(dl_mu) >= FCcutoff, "yes", "no"))
  res2 <- as.data.frame(res2) %>% mutate(sig = ifelse((PValue<=0.05 | PValue <= 0.05) & abs(dl_mu) >= FCcutoff, "yes", "no"))
  
  res1 <- as.data.frame(res1) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=FCcutoff & PValue<=0.05 , "up",
                                                                     ifelse(dl_mu<=-FCcutoff & PValue<=0.05 , "down", "unchanged")))
  
  res2 <- as.data.frame(res2) %>% dplyr::mutate(diffexpressed=ifelse(dl_mu>=FCcutoff & PValue<=0.05 , "up",
                                                                     ifelse(dl_mu<=-FCcutoff & PValue<=0.05 , "down", "unchanged")))
  res_merged <- merge(res1,res2,by="symbol")
  
  res_merged <- res_merged %>% dplyr::mutate(sig=
                                               ifelse(diffexpressed.x %in% c("up", "down") & diffexpressed.y %in% c("up", "down"), "Significant in both",
                                                      ifelse(diffexpressed.x %in% c("up", "down") & diffexpressed.y == "unchanged", paste("Changed only in ", disease1, sep=""),
                                                             ifelse(diffexpressed.y %in% c("up", "down") & diffexpressed.x == "unchanged", paste("Changed only in ", disease2, sep=""), "Unchanged in both"))))
  
  res_merged <- res_merged  %>% dplyr::select(dl_mu.x,dl_mu.y,diffexpressed.x, diffexpressed.y,sig,symbol) %>% unique %>% drop_na()
  
  p <- ggplot(res_merged, aes(x=dl_mu.x, y = dl_mu.y, color=sig)) +
    ggrastr::geom_point_rast(size=3,alpha=0.6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 24),
          legend.text = element_text(size=16),
          plot.title = element_text(hjust=0.5,size=28),
          legend.position = "right") +
    scale_color_manual(values=c("goldenrod1",  "red", "mediumblue", "gray67"), name = "") +
    geom_hline(yintercept =0) +
    geom_vline(xintercept=0) +
    geom_text_repel(size=8,
                    data = res_merged[res_merged$sig == "Significant in both",],
                    aes(x=dl_mu.x,y=dl_mu.y,label=symbol),
                    show.legend = FALSE,
                    box.padding = 0.8, max.overlaps = 20, fontface = "italic") +
    xlab(paste("LogFC in ", disease1, sep="")) + ylab(paste("LogFC in ", disease2, sep="")) + guides(fill=guide_legend(title=NULL))
  
  list_out <- list()
  list_out[[1]] <- res_merged
  list_out[[2]] <- p
  return(list_out)
}

GetAbundanceTable <- function(abundances,coldata,samplecol,subset.cols,calc.col,return_clr=F){
  require(edgeR)
  require(RColorBrewer)
  require(viridis)
  abundances <- unclass(abundances)
  
  extra.info <- coldata[match(colnames(abundances),coldata[[samplecol]]),]
  d <- DGEList(abundances, samples=extra.info)
  
  d = calcNormFactors(d)
  d= estimateCommonDisp(d, verbose=TRUE)
  
  norm_counts <- as.data.frame(t(d$counts)) 
  norm_counts <- norm_counts %>% mutate( !! samplecol := rownames(.))
  coldata <- as.data.frame(coldata)
  coldata_short <- coldata %>% dplyr::select(all_of(subset.cols)) %>% unique
  
  df_long_final <- left_join(norm_counts,coldata_short, by=samplecol)
  
  exclude <- (dim(df_long_final)[2]-dim(coldata_short)[2]+1):dim(df_long_final)[2]
  
  df_long_final[[samplecol]] <- factor(df_long_final[[samplecol]])
  percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
  sums <- data.frame(clustersums=rowSums(df_long_final[,-exclude]),ident=df_long_final[[samplecol]])
  
  for (clust in 1:length(levels(factor(coldata[[calc.col]])))) {
    for (sample in 1:length(df_long_final[,clust])) {
      percent <- (df_long_final[sample,clust]/sums[sample,1])*100
      percentages[sample,clust]<- percent
    }
  }
  
  
  percentages[,exclude] <- df_long_final[,exclude]
  
  colnames(percentages) <-  colnames(df_long_final)
  
  
  df_long <- percentages %>% 
    pivot_longer(unique(coldata[[calc.col]]),values_to = "Percent",names_to = "Cluster")
  
  
  if(return_clr){
    percentages <- percentages[,-exclude]
    value_to_add_to_offset_zeros <- min(percentages[percentages > 0],na.rm = T)
    percentages <- percentages + value_to_add_to_offset_zeros
    #Clr
    
    clr <- percentages
    for (sample in 1:dim(percentages)[1]) {
      tran <- log((percentages[sample,-exclude]/exp(mean(as.numeric(log( percentages[sample,-exclude]))))))
      #tran <- clr(percentages[sample,1:4])
      clr[sample,1:length(percentages[sample,-exclude])]<- tran
    }
    clr[,exclude] <- df_long_final[,exclude]
    colnames(clr) <-  colnames(df_long_final)
    df_transformed <- clr %>% 
      pivot_longer(unique(coldata[[calc.col]]),values_to = "Percent",names_to = "Cluster")
    
    return(df_transformed)
  }else{
    return(df_long)
  }
  
}


doMetaDE <- function(res,feat){
  require(metap)
  require(metafor)
  require(Matrix)
  require(dplyr)
  require(tidyr)
  require(tidyvers)
  
  
  #create empty matrix to hold values
  fc <- matrix(0, ncol=length(res), nrow=dim(feat)[1])
  rownames(fc) <- feat$symbol #fold changes
  
  pval <- matrix(NA, ncol=length(res), nrow=dim(feat)[1])
  rownames(pval) <- feat$symbol #pvalues
  
  mexp <- matrix(NA, ncol=length(res), nrow=dim(feat)[1])
  rownames(mexp) <- feat$symbol # mean expression
  
  sefc <- matrix(NA, ncol=length(res), nrow=dim(feat)[1])
  rownames(sefc) <- feat$symbol #standard error of fold change
  
  # for each study, at each gene,get FC, Avg expression, P-value and calculate standard error
  for(i in 1:length(res)){
    m   <- match(rownames(fc), rownames(res[[i]]))
    f.a <- !is.na(m)
    f.b <- m[f.a]
    fc[f.a,i] <- res[[i]][f.b,"LogFC"]
    pval[f.a,i] <- res[[i]][f.b,"PValue"]
    mexp[f.a,i] <- res[[i]][f.b,"AveExpr"]
    sefc[f.a,i] <- res[[i]][f.b,"LogFC"]/res[[i]][f.b,"t"] #standard error of effect size
  }
  
  # name columns of matricies
  x<-strsplit(names(res),".", fixed=T)
  x2 <- sapply(x, function(x) purrr::pluck(x,1))
  colnames(pval) = x2
  colnames(fc) = x2
  colnames(mexp) = x2
  colnames(sefc) = x2
  
  #combine matricies into list for calculations
  deg_master <- list(pval=pval, fc=fc, se_fc=sefc, mean_exp=mexp)
  
  f <- colnames(deg_master$pval)
  pval <- deg_master$pval[, f]
  fc <- deg_master$fc[, f]
  mexp <- deg_master$mean_exp[, f]
  sefc <- cbind(fc,deg_master$se_fc[,f])
  
  pval_up <- pval * 0
  pval_down <- pval * 0
  # for each gene calculate the get pvalue based on FC direction 
  for (i in 1:dim(fc)[2]) {
    f1 <- which(fc[, i] > 0)
    f2 <- which(fc[, i] < 0)
    pval_up[f1, i] <- pval[f1, i] * 0.5
    pval_up[f2, i] <- 1 - (pval[f2, i] * 0.5)
    pval_down[f1, i] <- 1 - (pval[f1, i] * 0.5)
    pval_down[f2, i] <- pval[f2, i] * 0.5
    
  }
  
  
  # Syntax
  x <- pval_down[rowSums(is.na(pval_down)) != ncol(pval_down), ]
  # Combine p-values by the sum of logs method, also known as Fisher's method 
  meta_up <- suppressWarnings(apply(pval_up, 1, function(y) sumlog(y[!is.na(y)])$p))
  meta_down <- suppressWarnings(apply(pval_down, 1, function(y) sumlog(y[!is.na(y)])$p))
  #average FC and expression 
  mean_fc <- rowMeans(fc, na.rm = T)
  mean_exp <- rowMeans(mexp, na.rm = T)
  
  
  #remove duplicates from directionality
  meta_down <- meta_down[-which(duplicated(names(meta_down))==TRUE)]
  meta_up <- meta_up[-which(duplicated(names(meta_up))==TRUE)]
  mean_fc <- mean_fc[-which(duplicated(names(mean_fc))==TRUE)]
  mean_exp <- mean_exp[-which(duplicated(names(mean_exp))==TRUE)]
  
  # find genes that are shared across studies
  f_se_fc <- apply(sefc[,1:length(res)],1,function(x) sum(!is.na(x))>=1)
  f_se_se <- apply(sefc[,(length(res)+1):ncol(sefc)],1,function(x) sum(!is.na(x))>=1)
  f_se <- f_se_fc & f_se_se
  
  #Function to fit meta-analytic random-effects model using DeSerminonian and Laird
  dsl_res <- apply(sefc[f_se,], 1,
                   function(row) rma(yi=as.vector(row[1:length(res)]),
                                     vi=as.vector(row[(length(res)+1):ncol(sefc)]),method="DL"))
  
  # Get columns  
  betas <- sapply(dsl_res,function(x) purrr::pluck(x,"beta"))
  ses <- sapply(dsl_res,function(x) purrr::pluck(x,"se"))
  dl_ps <- sapply(dsl_res,function(x) purrr::pluck(x,"pval"))
  dl_ps <- sapply(dsl_res,function(x) purrr::pluck(x,"pval"))
  
  # Create dataframe of results
  dl_sub <- data.frame(ID=rownames(sefc)[f_se],beta=betas,se=ses,pval=dl_ps)
  
  dl_all <- data.frame(ID=rownames(sefc))
  dl_all <- left_join(dl_all,dl_sub,by=c("ID"="ID"),relationship = "many-to-many")
  
  dl_all <- dl_all[-which(duplicated(dl_all$ID)==TRUE),]
  n_tested = apply(pval_up, 1, function(y)
    sum(!is.na(y))) # get number of studies gene shows up in
  n_tested <- n_tested[-which(duplicated(names(n_tested))==TRUE)]
  n_up = apply(fc, 1, function(y)
    sum(y > 0, na.rm = T)) # get number of studies gene FC is up-regulated
  n_up <- n_up[-which(duplicated(names(n_up))==TRUE)]
  n_down = apply(fc, 1, function(y)
    sum(y < 0, na.rm = T))# get number of studies gene FC is down-regulated
  n_down <- n_down[-which(duplicated(names(n_down))==TRUE)]
  
  #Create final Dataframe of results
  temp <-
    data.frame(
      ID = dl_all$ID,
      AveExpr = as.numeric(mean_exp),
      PValue = dl_all$pval,
      FDR = p.adjust(dl_all$pval,method="BH"),
      LogFC = as.numeric(mean_fc),
      dl_mu = dl_all$beta,
      dl_s = dl_all$se,
      metap_up = meta_up,
      metap_down = meta_down,
      adj_metap_up = p.adjust(meta_up, method = "BH"),
      adj_metap_down = p.adjust(meta_down, method = "BH"),
      n_tested = n_tested,
      n_up = n_up,
      n_down = n_down
    )
  
  # Clean dataframe of results
  x <- match(temp$ID,feat$symbol)
  features<- feat[x,]
  temp <- left_join(temp, features,by=c("ID"="symbol"))
  temp<- temp[which(is.nan(temp$AveExpr)==FALSE),]
  ord <- order(temp$metap_up)
  
  temp <- temp[ord, ]
  x <-
    c(
      "ID.y",
      "ID",
      "AveExpr",
      "PValue",
      "FDR",
      "LogFC",
      "dl_mu",
      "dl_s",
      "metap_up",
      "metap_down",
      "adj_metap_up",
      "adj_metap_down",
      "n_tested",
      "n_up",
      "n_down","desc"
    )
  
  res_dsl <- temp %>% dplyr::select(all_of(x))
  colnames(res_dsl)[1:2] <- c("ID","symbol")
  return(res_dsl)
}

DotPlotCompare <- function(
    gsea_list,
    n = 5,
    size_col = "NES",
    color_col = "p.adjust",
    size_cutoff = NULL,
    color_cutoff = NULL,
    direction = NULL,
    ySize=12
) {
  #' Create a Dot Plot of GSEA Results
  #'
  #' @param gsea_list A list where each element is a data frame of GSEA results for a specific cluster.
  #' Each data frame should have columns like `pathway`, size-related column (e.g., `NES`), p-value (`p.adjust`), and `setSize`.
  #' @param n Integer. Number of "top pathways" to select per cluster.
  #' @param size_col The column to use for scaling the point size in the plot.
  #' @param color_col The column to use for coloring the points.
  #' @param size_cutoff Numeric. Filter pathways based on the specified threshold for `size_col`.
  #'  - If positive (e.g., 2), includes values `>= size_cutoff`.
  #'  - If negative (e.g., -2), includes values `<= size_cutoff`.
  #'  - Default is `NULL` (no filtering by size_col).
  #' @param color_cutoff Numeric. Maximum `color_col` value to keep (e.g., a significance threshold). Default is `NULL` (no cutoff).
  #' @param direction Character. Direction of `size_col` values to filter:
  #'  - `"positive"`: Select top pathways with positive size_col (e.g., top NES).
  #'  - `"negative"`: Select top pathways with negative size_col (e.g., most negative NES).
  #'  - `NULL` (default): Select the top pathways with the smallest values for `color_col` (e.g., smallest p.adjust).
  #' @return A `ggplot2` dot plot showing pathways across clusters.
  
  # Load libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Validate the "direction" parameter
  if (!is.null(direction) && !direction %in% c("positive", "negative")) {
    stop("`direction` must be one of 'positive', 'negative', or NULL.")
  }
  
  # Step 1: Process the top n pathways for each cluster
  top_pathways <- lapply(names(gsea_list), function(cluster) {
    # Retrieve data for the specific cluster
    cluster_data <- gsea_list[[cluster]]
    
    # # Optional size cutoff filtering
    # if (!is.null(size_cutoff)) {
    #   if (size_cutoff > 0) {
    #     # Keep size_col >= size_cutoff
    #     cluster_data <- cluster_data %>%
    #       dplyr::filter(.data[[size_col]] >= size_cutoff)
    #   } else {
    #     # Keep size_col <= size_cutoff (negative cutoff for more negative values)
    #     cluster_data <- cluster_data %>%
    #       dplyr::filter(.data[[size_col]] <= size_cutoff)
    #   }
    # }
    
    # Optional color cutoff filtering
    if (!is.null(color_cutoff)) {
      cluster_data <- cluster_data %>%
        dplyr::filter(.data[[color_col]] <= color_cutoff)
    }
    
    # Logic for selecting pathways based on the "direction"
    if (!is.null(direction)) {
      if (direction == "positive") {
        # Top n most positive size_col values
        cluster_data <- cluster_data %>%
          dplyr::filter(.data[[size_col]] > 0) %>%
          dplyr::arrange(desc(.data[[size_col]])) %>%  # Sort by descending size_col
          dplyr::slice(1:n)  # Top n
      } else if (direction == "negative") {
        # Top n most negative size_col values
        cluster_data <- cluster_data %>%
          dplyr::filter(.data[[size_col]] < 0) %>%
          dplyr::arrange(.data[[size_col]]) %>%  # Sort by ascending size_col
          dplyr::slice(1:n)  # Top n
      } 
    } else {
      # Default: Top n pathways with the smallest color_col values
      cluster_data <- cluster_data %>%
        dplyr::arrange(.data[[color_col]]) %>%  # Sort by ascending color_col
        dplyr::slice(1:n)  # Top n
    }
    
    # Add cluster column
    cluster_data %>% dplyr::mutate(cluster = cluster)
  }) %>%
    dplyr::bind_rows()  # Combine data from all clusters

  
  # Step 2: Merge results for all clusters and prepare for plotting
  enriched_pathways <- lapply(names(gsea_list), function(cluster) {
    gsea_list[[cluster]] %>%
      dplyr::mutate(
        cluster = cluster
      )
  }) %>%
    dplyr::bind_rows()
  
  # Filter to keep only pathways that were in the top n in any cluster
  filtered_pathways <- top_pathways %>%
    dplyr::pull(Description) %>%
    unique()
  
  # Shape plot data
  plot_data <- enriched_pathways %>%
    dplyr::filter(Description %in% filtered_pathways) %>%  # Retain only relevant pathways
    dplyr::select(ID,Description, cluster, .data[[size_col]], .data[[color_col]]) %>%
    dplyr::rename(size_value = !!size_col, color_value = !!color_col) %>%
    tidyr::complete(Description, cluster, fill = list(size_value = NA, color_value = NA)) # Ensure all combinations exist
  
  # Arrange for clean grouping and ordered plotting
  plot_data <- plot_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(cluster, desc(size_value), .by_group = TRUE) %>%
    dplyr::ungroup()
  plot_data <- plot_data[complete.cases(plot_data),]
  
  plot_data$Description <- stringr::str_wrap(plot_data$Description,width = 40)
  
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(stringr::str_wrap(
                                    unique(top_pathways$Description),
                                    width = 40)
                                  ))
  
  plot_data$cluster <- factor(plot_data$cluster,levels = names(gsea_list))
  plot_data$size_value<- round(plot_data$size_value,2)
  # Step 3: Create the dot plot
  dot_plot <- ggplot(plot_data, aes(x = cluster, y = Description)) +
    geom_point(aes(size = abs(size_value), color = color_value), alpha = 0.95) +
    scale_colour_distiller(palette = "RdYlBu",direction = 1,
                           name = color_col) +
    scale_size_continuous(range = c(3, 10), name = size_col) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = ySize),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste0("Top Pathway Enrichment (Top ", n, " Pathways per Cluster)"),
    )+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(dot_plot)
}
