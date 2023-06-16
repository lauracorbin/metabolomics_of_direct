## contains functions to make figures

## load libraries
library('heatmap3')
library('RColorBrewer')
library('pcaMethods')
library('psych')
library('ggplot2')
library('patternplot')
library('stringr')


#################################################
################# volcano plots #################
#################################################

make.volcano = function(dtst, title, filename, plot_col, hline, legend_one,legend_col,plot_shape, legend_two, legend_shape){
  # pdf
  pdf(paste0(filename,".pdf"))
  with(dtst, plot(raw_fc_log2medianfoldchange_int_over_control, -log10(rnt_lm_treat_HolmAdjP), 
                              main=title,xlab = "log2 median fold change (intervention/control)",
                              ylab="-log10(Holm-corrected p)", col=plot_col, pch=plot_shape, xlim=(c(-2,2)) ))
  abline(h=hline,lty=2, v=0)
  legend(x="topright", inset=c(0,0),legend=legend_one, pch=18, bty="n", xjust=1, col=legend_col, cex=0.85)
  legend(x="topleft", inset=c(0,0), legend=str_to_title(legend_two), pch=legend_shape, bty="n", xjust=1, cex=0.85)
  text(x=-1.8,y=hline+0.3,"p=0.05")
  invisible(dev.off())
}

#######################################
############## heatmaps  ##############
#######################################

make.heatmap = function(dtst,feature_ids,filename,colcolours){
  # replace column names in data
  feature_id_list <- as.data.frame(names(dtst))
  feature_id_list$feature_label_id <- feature_ids$feature_label_id[match(feature_id_list$`names(dtst)`,feature_ids$feature_id)]
  names(dtst) <- feature_id_list$feature_label_id
  # order by name
  dtst <- dtst[,order(names(dtst))]
  
  # count no. of features per class - needs to be in alphabetical order
  AA_ncols <- length(grep( "AA" , names(dtst ) ) )
  CHO_ncols <- length(grep( "CHO" , names(dtst ) ) )
  Vit_ncols <- length(grep( "CoFac & Vit" , names(dtst) ) )
  Energy_ncols <- length(grep( "Energy" , names(dtst ) ) )
  Lipid_ncols <- length(grep( "Lipid" , names(dtst ) ) )
  NMR.derived_nolc <- length(grep( "NMR.derived" , names(dtst ) ) )
  Nucleotide_ncols <- length(grep( "Nucleotide" , names(dtst ) ) )
  Peptide_ncols <- length(grep( "Peptide" , names(dtst ) ) )
  unknown_ncols <- length(grep( "Unclassified" , names(dtst ) ) )
  Xeno_ncols <- length(grep( "Xeno" , names(dtst ) ) )
  
  #define row colours
  pcol <- brewer.pal(10,"Set3")
  RowSideColors <- cbind(Super.pathway=c(rep(pcol[1],AA_ncols),rep(pcol[2],CHO_ncols),rep(pcol[3],Vit_ncols),rep(pcol[4],Energy_ncols),
                                   rep(pcol[5],Lipid_ncols),rep(pcol[6],NMR.derived_nolc),rep(pcol[7],Nucleotide_ncols),
                                   rep(pcol[8],Peptide_ncols),rep(pcol[9],unknown_ncols),rep(pcol[10],Xeno_ncols) ))
  
  # transpose
  dtst_t <- t(dtst)
  dim(dtst_t)
  
  # set up empty vectors to use in place of row & column names
  empty.cols = unlist(lapply(colnames(dtst),function(x){
    a = " "
  }))
  empty.rows = unlist(lapply(row.names(dtst),function(x){
    a = " "
  }))
  
  # make dendograms
  features_dist <- dist(dtst_t)
  features_den <- as.dendrogram(hclust(features_dist, method = "ward.D2"))
  samples_dist <- dist(t(dtst_t))
  samples_den <- as.dendrogram(hclust(samples_dist, method = "ward.D2"))
  
  # make pdf plot
  pdf(paste0(filename,".pdf"))
  heatmap3(x=as.matrix(dtst_t),Rowv = features_den, Colv = samples_den,
           ColSideColors = colcolours, RowSideColors = RowSideColors,labCol = empty.cols,labRow = empty.rows, scale="none",
           legendfun=function() showLegend(legend = c("Amino Acid","Carbohydrate","CoFac & Vit","Energy","Lipid",
                                                      "NMR derived","Nucleotide","Peptide","Unclassified","Xenobiotics"),
                                           col=c(pcol[1],pcol[2],pcol[3],pcol[4],pcol[5],pcol[6],
                                                 pcol[7],pcol[8],pcol[9],pcol[10]), cex=0.8),
           margins=c(8,8))
  invisible(dev.off())
  # make postscript plot
  postscript(paste0(filename,".eps"))
  heatmap3(x=as.matrix(dtst_t),Rowv = features_den, Colv = samples_den,
           ColSideColors = colcolours, RowSideColors = RowSideColors,labCol = empty.cols,labRow = empty.rows, scale="none",
           legendfun=function() showLegend(legend = c("Amino Acid","Carbohydrate","CoFac & Vit","Energy","Lipid",
                                                      "NMR derived","Nucleotide","Peptide","Unclassified","Xenobiotics"),
                                           col=c(pcol[1],pcol[2],pcol[3],pcol[4],pcol[5],pcol[6],
                                                 pcol[7],pcol[8],pcol[9],pcol[10]), cex=0.8),
           margins=c(8,8))
  invisible(dev.off())
  
  # return modified dataset
  return(dtst)
}


########################################
############## pca plots  ##############
########################################

make.pca.plot = function(pca.data,colour_spec,shape_spec,scree.data,filename) {
  pdf(paste0(filename,".pdf"))
  par(mfrow=c(2,1))
  plot(scree.data, type = "b", pch = 21, bg = "blue", col = "grey40", lwd = 3, cex = 2,
       xlab = "PC", ylab = "% Variance explained | R2", main = "Scree Plot")
  # pca plots
  plot(pca.data@scores[,1],pca.data@scores[,2],
       pch = shape_spec$sex.grp, col=colour_spec$Treat.grp, cex = 0.8,
       xlab = paste0("PC1 VarExp = ", varexp[1], "%"),
       ylab = paste0("PC2 VarExp = ", varexp[2], "%"),
       main="Probabilistic PCA - by group", 
       ylim=c(min(pca.data@scores[,2])*1.7,max(pca.data@scores[,2])),xlim=c(min(pca.data@scores[,1])*1.3,max(pca.data@scores[,1])*1.3))
  legend(x="bottomleft",legend=as.character(unique(colour_spec$treat)),col=unique(colour_spec$Treat.grp),pch=21,bty="n", xjust=1,cex=1.0)
  legend(x="bottomright",ncol=1,
         legend=c(as.character(unique(shape_spec$sex))),col="black",pch = c(unique(shape_spec$sex.grp)),bty="n", xjust=1,cex=1.0)
  invisible(dev.off())
}

lm.pcs = function(pca.data,trial_data,id_list) {
  # make output file
  lm_results = as.data.frame(matrix(data = NA, nrow = ncol(trial_data)-1, ncol = 11))
  # set up data for model
  pc_data <- as.data.frame(pca.data@scores[,1])
  pc_data$pc2 <- pca.data@scores[,2]
  pc_data$pc3 <- pca.data@scores[,3]
  pc_data$pc4 <- pca.data@scores[,4]
  pc_data$pc5 <- pca.data@scores[,5]
  #names(pc_data)[1] <- "pc1"
  pc_data$id <- id_list
  data_for_model <- merge(pc_data,trial_data,all.x=T,by="id")
  data_for_model$site <- as.factor(data_for_model$site)
  # run model for each covariate
  for (i in 7:ncol(data_for_model)) {
    #print(i)
    n_levels <- 100
    if (class(data_for_model[,i]) == "factor"){
      n_levels <- length(unique(data_for_model[,i]))
    }
    #print(n_levels)
    if (n_levels > 1 & (length(which(is.na(data_for_model[,i]))) <= nrow(data_for_model)-1 ) ) {
      fit1 <- lm(data_for_model[,2]~data_for_model[,i])
      fit2 <- lm(data_for_model[,3]~data_for_model[,i])
      fit3 <- lm(data_for_model[,4]~data_for_model[,i])
      fit4 <- lm(data_for_model[,5]~data_for_model[,i])
      fit5 <- lm(data_for_model[,6]~data_for_model[,i])
      lm_results[i-4,1] <- names(data_for_model)[i] 
      lm_results[i-4,2] <- summary(fit1)$r.squared
      lm_results[i-4,3] <- lmp(fit1)
      lm_results[i-4,4] <- summary(fit2)$r.squared
      lm_results[i-4,5] <- lmp(fit2)
      lm_results[i-4,6] <- summary(fit3)$r.squared
      lm_results[i-4,7] <- lmp(fit3)
      lm_results[i-4,8] <- summary(fit4)$r.squared
      lm_results[i-4,9] <- lmp(fit4)
      lm_results[i-4,10] <- summary(fit5)$r.squared
      lm_results[i-4,11] <- lmp(fit5)
    }
    else {
      lm_results[i-4,1] <- names(data_for_model)[i]
      lm_results[i-4,2:11] <- NA
    }
  }
  names(lm_results)[1] <- "pheno"
  names(lm_results)[2] <- "pc1_R2"
  names(lm_results)[3] <- "pc1_p"
  names(lm_results)[4] <- "pc2_R2"
  names(lm_results)[5] <- "pc2_p"
  names(lm_results)[6] <- "pc3_R2"
  names(lm_results)[7] <- "pc3_p"
  names(lm_results)[8] <- "pc4_R2"
  names(lm_results)[9] <- "pc4_p"
  names(lm_results)[10] <- "pc5_R2"
  names(lm_results)[11] <- "pc5_p"
  return(lm_results)
}

lmp = function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
