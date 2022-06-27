# functions for analysis script 06.0

# load libraries
library('nFactors')
library('ggrepel')
library('ggplot2')
library('ggforce')
library('ggpubr')

# functions
# hypergeometric function for enrichment - from David Hughes
hgtest = function(allfeatures, selectedfeatures) {
  ## unique list of categories in the selected feature list
  cats = na.omit( unique(selectedfeatures) )
  ## iterate over those cat and perform test
  
  HG_F_test = sapply(cats, function(i){
    g = length(allfeatures)
    ## total number of tested features with some annotation
    N = length(allfeatures) - ( sum(is.na(allfeatures) | allfeatures == "") )
    ## number of features in cat i
    m = sum(allfeatures == i, na.rm = TRUE )
    ## number of features NOT in cat i
    n = N - m
    ## number of features selected by X criteria, that have some annotation
    k = length(selectedfeatures) - ( sum(is.na(selectedfeatures) | selectedfeatures == "") )
    ## number of selected features in cat i
    x = sum(selectedfeatures == i, na.rm = TRUE)
    
    
    ## estimate fold enrichment and p.value
    fold.enrichment <- (x / k ) / (m / N)
    p.value <- phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    
    ## FISHER EXACT TEST
    ## BUILD 2 X 2 contingency table
    dmat = matrix( c(x,k-x,m-x,n-(k-x)),2,2, byrow=TRUE, dimnames = list(c("Selected","NotSelected"), c("FocusedCat", "AllOtherCats") ))
    print(i)
    print(dmat)
    ftest <- fisher.test(x=dmat, alternative = "greater")
    
    
    ## data out
    out = c(fold.enrichment, p.value, ftest$estimate, ftest$p.value)
    names(out) = c("fold.enrichment", "HG_pval", "Ftest_OR", "Ftest_pval")
    return(out)
  })
  ## return to user
  return(t(HG_F_test))
  
  
}

# modified version of moose_biplot by David Hughes to deal with ppca data output
corbin_biplot_for_ppca = function(PCA, dataframe_of_phenotypes, pheno_list, npcs, class_info, annot_data,ordered_id,filename){
  
  ##############################
  ## build correlation matrix
  ## pf PCs against external variables
  ##############################
  #varexp = summary(PCA)[[6]][2,]
  varexp = summary(PCA)
  for(pc_x in 1:(npcs-1)){
    for(pc_y in 2:npcs){
      if ( (pc_x != pc_y) & (pc_x < pc_y) ) {
        print(paste0("Run for PC: ", pc_x, " and PC: ",pc_y))
        
        # plot original PCA loadings
        # extract loadings as a dataframe
        loading_for_plot <- as.data.frame(PCA@loadings)
        loading_for_plot$feature_id <- rownames(loading_for_plot)
        loading_for_plot$biochemical <- class_info$biochemical[match(loading_for_plot$feature_id,class_info$feature_id)]
        # look at main driving metabolites
        loading_for_plot[which(loading_for_plot$V1 <= -0.2),"biochemical"]
        loading_for_plot[which(loading_for_plot$V1 >= 0.2),"biochemical"]
        loading_for_plot[which(loading_for_plot$V2 <=- 0.2),"biochemical"]
        loading_for_plot[which(loading_for_plot$V2 >= 0.2),"biochemical"]
        # mark subset for labelling
        loading_for_plot$ToLabel <- 0
        loading_for_plot[loading_for_plot$V1 < -0.2,c("ToLabel")] <- 1
        loading_for_plot[loading_for_plot$V1 > 0.2,c("ToLabel")] <- 1
        loading_for_plot[loading_for_plot$V2 > 0.2,c("ToLabel")] <- 1
        loading_for_plot[loading_for_plot$V2 < -0.2,c("ToLabel")] <- 1
        # make plot
        graphfile <- paste0(results_dir,"Figures\\ForPaper\\",filename,"_biplot_loadings_",pc_x,"_",pc_y)
        pdf(paste0(graphfile,".pdf"))
        plotout = 
          ggplot(loading_for_plot, aes(x = V1, y = V2) ) +
          geom_point() +
          #geom_text(aes(label=biochemical)) + 
          coord_fixed(ratio = 1) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          geom_hline(yintercept = 0, linetype="dashed", colour="black") +
          geom_vline(xintercept = 0, linetype="dashed", colour="black") +
          geom_text_repel(aes(label=ifelse(ToLabel == 1,as.character(biochemical),'')),hjust=1,vjust=0.5,cex=2.5) +
          labs(x = "PC1",y="PC2")
        print(plotout)
        invisible(dev.off())
        
        #mydata = data.frame( PC1 = PCA@scores[,pc_x], PC2 = PCA@scores[,pc_y], class = groupfactor)
        mydata = data.frame( id = ordered_id,PC1 = PCA@scores[,pc_x], PC2 = PCA@scores[,pc_y])
        
        # order phenotype by id
        dataframe_of_phenotypes <- dataframe_of_phenotypes[order(dataframe_of_phenotypes$id),]
        ## empty Correlation Matrix
        cormat = matrix(NA, ncol(dataframe_of_phenotypes) , npcs*4)
        rownames(cormat) <- colnames(dataframe_of_phenotypes)
        est_col <- c(1,5)
        lci_col <- c(2,6)
        uci_col <- c(3,7)
        p_col <- c(4,8)
        col_names_vec <- vector()
        col_names_vec <- c(col_names_vec,1:(npcs*4))
        col_names_vec[est_col] <- paste0("PC", 1:npcs, ".r")
        col_names_vec[lci_col] <- paste0("PC", 1:npcs, ".lower.95CI")
        col_names_vec[uci_col] <- paste0("PC", 1:npcs, ".upper.95CI")
        col_names_vec[p_col] <- paste0("PC", 1:npcs, ".p")
        colnames(cormat) <- col_names_vec
        ## estimate pearson with pheno of interest
        for(i in 1:ncol(dataframe_of_phenotypes)){
          for(j in 1:npcs){
            est = cor.test(dataframe_of_phenotypes[,i], PCA@scores[,j], method = "pearson", use="pairwise.complete.obs")$estimate
            lci = cor.test(dataframe_of_phenotypes[,i], PCA@scores[,j], method = "pearson", use="pairwise.complete.obs")$conf.int[1]
            uci = cor.test(dataframe_of_phenotypes[,i], PCA@scores[,j], method = "pearson", use="pairwise.complete.obs")$conf.int[2]
            pval = cor.test(dataframe_of_phenotypes[,i], PCA@scores[,j], method = "pearson", use="pairwise.complete.obs")$p.value
            cormat[i,(j*4)-3] = est
            cormat[i,(j*4)-2] = lci
            cormat[i,(j*4)-1] = uci
            cormat[i,(j*4)] = pval
          }
        }

        #k = order(rhosum, decreasing = TRUE)[1:plot_top_N_phenotypes]
        cormat = cormat[pheno_list, ]
        print(cormat)

        # write out
        cormat_to_write <- as.data.frame(cormat)
        print(row.names(cormat_to_write))
        cormat_to_write$clinical.variable <- row.names(cormat_to_write)
        cormat_to_write <- cormat_to_write[,c(ncol(cormat_to_write),1:(ncol(cormat_to_write)-1))]
        print(cormat_to_write)
        write.table(cormat_to_write,file=paste0(results_dir,"Tables\\TableS5_",filename,"_pca_pheno_cor.txt"),
                    sep="\t", quote=F,row.names = FALSE)
        
        # keep only those that have p<0.05 for PC1 and or 2
        # NB currently works for PC1/PC2
        cormat <- as.data.frame(cormat)
        cormat <- cormat[(cormat[4] < 0.05 | cormat[8] < 0.05), ]
        print(cormat)
        
        # restrict to cor estimates (drop pval)
        # NB currently works for PC1/PC2
        cormat <- as.matrix(cormat[,c(1,5)])
        print(cormat)
        
        ## redefine the rotation | loadings in the PCA
        #PCA[[2]] = cormat
        PCA@loadings = cormat
        
        ## stealing  code from ggbiplot
        pcobj = PCA
        
        #nobs.factor <- sqrt(nrow(pcobj$x) - 1)
        nobs.factor <- sqrt(nrow(pcobj@scores) - 1)
        #d <- pcobj$sdev
        #d <- pcobj@sDev
        d <-apply(scores(pcobj),2,sd)
        #u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
        u <- sweep(pcobj@scores, 2, 1/(d * nobs.factor), FUN = "*")
        #v <- pcobj$rotation
        v <- pcobj@loadings
        
        choices = c(pc_x,pc_y)
        obs.scale = 0
        var.scale = 1
        #circle.prob = 0.69
        circle.prob = 0.99
        df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
        #v <- sweep(v, 2, d^var.scale, FUN = "*")
        df.v <- as.data.frame(v[, choices])
        names(df.u) <- c("xvar", "yvar")
        names(df.v) <- names(df.u)
        
        df.u <- df.u * nobs.factor
        
        r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
        v.scale <- rowSums(v^2)
        df.v <- r * df.v/sqrt(max(v.scale))
        df.v$labels = rownames(df.v)
        
        # merge annotations with PC data
        mydata <- merge(mydata,annot_data,by="id")
        
        # define shape var
        mydata$subgroup_v1 <- paste0(mydata$treat,"_",mydata$reversal)
        mydata$subgroup <- NA
        mydata[which(mydata$subgroup_v1 == "Control_No"),c("subgroup")] <- "Control, no remission"
        mydata[which(mydata$subgroup_v1 == "Control_Yes"),c("subgroup")] <- "Control, remission"
        mydata[which(mydata$subgroup_v1 == "Intervention_No"),c("subgroup")] <- "Intervention, no remission"
        mydata[which(mydata$subgroup_v1 == "Intervention_Yes"),c("subgroup")] <- "Intervention, remission"
        
        ### GGPLOT with labels
        graphfile <- paste0(results_dir,"Figures\\",filename,"_biplot_",pc_x,"_",pc_y)
        pdf(paste0(graphfile,".pdf"))
        plotout = mydata %>% ggplot( aes(x = PC1, y = PC2) ) +
          geom_point(aes(colour=weight.change, shape=subgroup)) +
          coord_fixed(ratio = 1) +
          geom_text(data = df.v, aes(label = labels,
                                     x = xvar, y = yvar), 
                    color = "black", size = 2) +
          geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar),
                       arrow = arrow(length = unit(1/2, "picas")),
                       size = 1, color = "grey50") +
          scale_shape_manual(values=c(1,16,2,17),name="Subgroup (allocation, remission status)") +
          scale_color_gradient(low="blue",high="red",name="Weight change (kg)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          geom_hline(yintercept = 0, linetype="dashed", colour="black") +
          geom_vline(xintercept = 0, linetype="dashed", colour="black") +
          labs(x = paste0("PC",pc_x,";  variance explained = ", signif(varexp[1,pc_x], d = 2) ),
               y = paste0("PC",pc_y,";  variance explained = ", signif(varexp[1,pc_y], d = 2) ) )
        print(plotout)
        dev.off()

        
        ### GGPLOT without labels
        graphfile <- paste0(results_dir,"Figures\\ForPaper\\",filename,"_biplot_",pc_x,"_",pc_y,"_nolabels")
        pdf(paste0(graphfile,".pdf"))
        plotout = mydata %>% ggplot( aes(x = PC1, y = PC2) ) +
          geom_point(aes(colour=weight.change, shape=subgroup)) +
          coord_fixed(ratio = 1) +
          #scale_x_continuous(limits=c(axis_min,axis_max)) +
          #scale_y_continuous(limits=c(axis_min,axis_max)) +
          geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar),
                       arrow = arrow(length = unit(1/2, "picas")),
                       size = 1, color = "grey50") +
          scale_shape_manual(values=c(1,16,2,17),name="Subgroup (allocation, remission status)") +
          scale_color_gradient(low="blue",high="red",name="Weight change (kg)") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          geom_hline(yintercept = 0, linetype="dashed", colour="black") +
          geom_vline(xintercept = 0, linetype="dashed", colour="black") +
          labs(x = paste0("PC",pc_x,";  variance explained = ", signif(varexp[1,pc_x], d = 2) ),
               y = paste0("PC",pc_y,";  variance explained = ", signif(varexp[1,pc_y], d = 2) ) )
        print(plotout)
        dev.off()
      }
    }
  }
}

