# This script runs analyses to extract a subset of representative metabolites based on correlative structure

# last run: 31st March 2022

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/05.1_variable_selection.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/05.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )

## call R scripts containing functions
source(paste0(scripts_dir,"02.0_qc_functions.R"))
source(paste0(scripts_dir,"05.0_functions.R"))
source(paste0(scripts_dir,"RPackages/growthPheno/R/PVA.R"))

# open pdf to print figures to
pdf(paste0(results_dir,"Figures\\variable_selection.pdf"))

# set cut height for tree
h = 0.80

# missing threshold for high missingness
miss_thresh = 0.40

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in  datasets
# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))
# clinic data
clinic_data_only <- readRDS(file=paste0(processeddata_dir,"clinic_data_input.rds"))
# summary stats
ds_results_list <- readRDS(file=paste0(results_dir,"metabolite_summaries.rds"))
# linear model residuals (w/o adjustment for treatment group)
lm_residuals_list <- list(metab = readRDS(file=paste0(processeddata_dir,"metabolon_lm_rnt_resid.rds")), nmr = readRDS(file=paste0(processeddata_dir,"nmr_lm_rnt_resid.rds"))) 


## metabolite feature metadata (post cleaning)
clean_feature_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_feature_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_feature_metadata.rds"))) 

## add column for feature id matching
clean_feature_list[[1]]$id_for_matching <- clean_feature_list[[1]]$compid_to_match
clean_feature_list[[2]]$id_for_matching <- clean_feature_list[[2]]$shortname

## sample metadata (post cleaning)
clean_sample_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_sample_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_sample_metadata.rds")))

####################################################################################
####################################################################################
## process data prior to data reduction
####################################################################################
####################################################################################

# define raw data list
raw_data_list <- list(metab=dataset_list$raw$metab, nmr=dataset_list$raw$nmr)

## identify derived vs raw measures in NMR
percent_measures <- names(raw_data_list$nmr[,grep("%",names(raw_data_list$nmr))])
print("No. of percent measures in NMR (baseline + endpoint):")
length(percent_measures)
ratio_measures <- names(raw_data_list$nmr[,grep("/",names(raw_data_list$nmr))])
print("No. of ratio measures in NMR (baseline + endpoint):")
length(ratio_measures)
derived_measures <- c(percent_measures, ratio_measures)
raw_measures <- names(raw_data_list$nmr[,! names(raw_data_list$nmr) %in% derived_measures]) # also captures id column
print("No. of raw measures in NMR (baseline + endpoint) - to be included in selection step:")
length(raw_measures) - 1

# get list of raw measures w/o extension
raw_measures_clean <- gsub("[.]e","",raw_measures)
raw_measures_clean <- gsub("[.]b","",raw_measures_clean)
raw_measures_clean <- unique(raw_measures_clean)

# get list of derived measures w/o extension
derived_measures_clean <- gsub("[.]e","",derived_measures)
derived_measures_clean <- gsub("[.]b","",derived_measures_clean)
derived_measures_clean <- unique(derived_measures_clean)

# write out lists 
write.table(derived_measures_clean, file=paste0(results_dir,"NMRderived.txt"),row.names = F, col.names = F, quote = F)
write.table(raw_measures_clean, file=paste0(results_dir,"NMRraw.txt"),row.names = F, col.names = F, quote = F)

# add IDs to raw data (residual data already has)
raw_data_list$nmr$id <- rownames(raw_data_list$nmr)
raw_data_list$metab$id <- rownames(raw_data_list$metab)

# merge nmr and metabolon
raw_data <- merge(raw_data_list$nmr,raw_data_list$metab,by="id",all=T)
resid_data <- merge(lm_residuals_list$nmr,lm_residuals_list$metab, by="id",all=T)
  
# drop clincal assay output
raw_data <- raw_data[, -grep("mmol",colnames(raw_data))]
raw_data <- raw_data[, -grep("insulin",colnames(raw_data))]
resid_data <- resid_data[, -grep("mmol",colnames(resid_data))]
resid_data <- resid_data[, -grep("insulin",colnames(resid_data))]

print("No. of features across two panels going into selection step:")
(ncol(raw_data)-1)/2

# drop endpoint measures
raw_data.b <- raw_data[, -grep("[.]e",colnames(raw_data))]

# drop baseline measures
raw_data.e <- raw_data[, -grep("[.]b",colnames(raw_data))]

####################################################################################
####################################################################################
## prep data for making tree
####################################################################################
####################################################################################

# remove ID and transpose dataset
raw_data_t <- as.data.frame(t(raw_data[,-c(1)]))
dim(raw_data_t)

raw_data.e_t <- as.data.frame(t(raw_data.e[,-c(1)]))
dim(raw_data.e_t)

raw_data.b_t <- as.data.frame(t(raw_data.b[,-c(1)]))
dim(raw_data.b_t)

resid_data_t <- as.data.frame(t(resid_data[,-c(1)]))
dim(resid_data_t)

# calculate missingness based on endpoint data
feature_missing <- feature.missing(dtst=raw_data.e_t)

# find those with less than threshold % missing at endpoint 
featuresIN.e <- names(which(feature_missing < miss_thresh))
print(paste0("No. of metabolite features with <", miss_thresh*100, "% missing at the 12m timepoint (to be used in cluster tree based selection): ", length(featuresIN.e)))

# get equivalent list labelled as baseline
featuresIN.b <- gsub("[.]e",".b",featuresIN.e)
print(paste0("Check no. of metabolite features with <", miss_thresh*100, "% missing transformed to the baseline timepoint (should be the same as above): ", length(featuresIN.b)))

# get equivalent with no extension
featuresIN <- gsub("[.]e","",featuresIN.e)
print(paste0("Check no. of metabolite features with <", miss_thresh*100, "% missing transformed to the baseline timepoint (should be the same as above): ", length(featuresIN)))

# write out list of features with less missing than threshold 
write.table(as.data.frame(c(featuresIN)),file=paste0(results_dir,"lm_feature_list.txt"),row.names=F,col.names=F,quote=F,sep="\t")


####################################################################################
####################################################################################
## make tree 
####################################################################################
####################################################################################

# generate spearmans correlation based tree using features (at baseline) with less than threshold % missing at 12mth
SpearTree <- make.tree(dtst=resid_data_t[featuresIN,] ,cor_method="spearman", hclust_method="complete")

# restrict based on cut off
k = cutree(SpearTree, h = h)
print("Chosen cut height is:")
h

# count clusters
print(paste0("There are ", length(unique(k)), " clusters at cut height ", h))

# save out groupings
k_out <- as.data.frame(k)
k_out$feature_id <- row.names(k_out)
k_out <- k_out[,c(2,1)]
names(k_out)[2] <- "k_group"

write.table(k_out,file=paste0(results_dir,"feature_group_allocation.txt"),row.names=F,quote=F,sep="\t")

# identify clusters/groups
k_group = table(k)

# plot group distributions
barplot(k_group, main="Distribution of features across clusters")
hist(k_group, main="No. of features per cluster")

# identify one feature classes)
w = which(k_group == 1); ind_k = names( k_group[w] ); w = which(k %in% ind_k)
ind = names(k)[w]
print("Total no. of one feature classes:")
length(ind)

# for remaining classes (with >1 metabolite), evaluate variance explained by each metabolite in the cluster and id lead metabolite (by PVA) 
w = which(k_group > 1); k_group_mult = names( k_group[w] )
k_group_summary <- as.data.frame(matrix(data=NA, nrow = length(k_group_mult), ncol = 7))
pva_output <- NULL
for (i in 1:length(k_group_mult)) {
  x <- k_group_mult[i]
  k_group_summary[i,1] = x
  w = which(k %in% x); n=names(k[w])
  pva1 <- PVA(responses=n,data=resid_data,plot=FALSE)
  # add cluster ref
  pva1$cluster <- x
  # add last var that's missing from pva list
  final_var <- which(!n %in% pva1$Selected)
  final_var_row <- data.frame("Variable" = length(n), "Selected" = n[final_var], 
                                "h.partial" = NA, "Added.Propn" = NA, "Cumulative.Propn" = NA, "cluster" = x)
  pva1$Selected <- as.character(pva1$Selected)
  pva1 <- rbind(pva1,final_var_row)
  # merge output
  if (i==1) {
    pva_output <- pva1
  } else {
    pva_output <- rbind(pva_output,pva1)
  }
  k_group_summary[i,2] = as.character(pva1[1,2])
  k_group_summary[i,3] = as.numeric(pva1[1,4])
  cor_matrix1 <- cor(t(resid_data_t[n,]),method="spearman", use="pairwise.complete.obs")
  # flatten cor matrix
  ut <- upper.tri(cor_matrix1)
  cor_output1 <- data.frame(
    row = rownames(cor_matrix1)[row(cor_matrix1)[ut]],
    column = rownames(cor_matrix1)[col(cor_matrix1)[ut]],
    cor = (cor_matrix1)[ut]
  )
  # extract mean, min and max
  k_group_summary[i,4] = mean(abs(cor_output1$cor))
  k_group_summary[i,5] = max(abs(cor_output1$cor))
  k_group_summary[i,6] = min(abs(cor_output1$cor))
  # no. in group
  k_group_summary[i,7] = length(n)
}
names(k_group_summary)[1] <- "k_group"
names(k_group_summary)[2] <- "lead_rep_feature"
names(k_group_summary)[3] <- "var_exp_by_rep"
names(k_group_summary)[4] <- "mean_within_grp_cor"
names(k_group_summary)[5] <- "max_within_grp_cor"
names(k_group_summary)[6] <- "min_within_grp_cor"
names(k_group_summary)[7] <- "n_features_in_group"

# add single var clusters to pva output
# id rows for single var clusters from cluster output file 
single_var_clusters <- which(k_out$feature_id %in% ind)
single_var_clusters_out <- NULL
for (i in 1:length(single_var_clusters)) {
  var_row <- data.frame("Variable" = 1, "Selected" = k_out[single_var_clusters[i],"feature_id"], "h.partial" = NA, 
                          "Added.Propn" = NA, "Cumulative.Propn" = NA, "cluster" = k_out[single_var_clusters[i],"k_group"])
  if (i==1) {
    single_var_clusters_out <- var_row
  } else {
    single_var_clusters_out <- rbind(single_var_clusters_out,var_row)
  }
}

# add on to pva out
pva_output_final <- rbind(pva_output,single_var_clusters_out)
# save out
write.table(pva_output_final,file=paste0(results_dir,"pva_output_by_feature.txt"),row.names=F,quote=F,sep="\t")


# summarise groups
print(paste0("Total no. of classes with 2 or more metabolites: ", nrow(k_group_summary)))
print(paste0("Summary of class size: "))
summary(k_group_summary$n_features_in_group)
print(paste0("Summary of variance explained by lead representative feature: "))
summary(k_group_summary$var_exp_by_rep)
print(paste0("Summary of mean correlation within cluster: "))
summary(k_group_summary$mean_within_grp_cor)
print(paste0("Summary of max. correlation within cluster: "))
summary(k_group_summary$max_within_grp_cor)
print(paste0("Summary of min. correlation within cluster: "))
summary(k_group_summary$min_within_grp_cor)

# join the two lists
rep_features = paste( c(ind, as.character(k_group_summary$rep_feature)))
#print(paste0("Total no. of metabolites selected to analyse in next step: ", length(rep_features)))

## plot tree
#pdf(paste0(results_dir,"Figures/hierarchical_tree.pdf"))
plot(SpearTree, hang = -1, cex = 0.5, main = paste0("Spearman Cluster Dendrogram with cut height of ",h), col = "royalblue", lwd = 2)
rect.hclust( SpearTree , h = h, border="red")
#invisible(dev.off())

# check corr
# 'representative' subset
cor_matrix_ind <- cor(t(resid_data_t[rep_features,]),method="spearman", use="pairwise.complete.obs")
# flatten cor matrix
ut <- upper.tri(cor_matrix_ind)
cor_output_ind <- data.frame(
  row = rownames(cor_matrix_ind)[row(cor_matrix_ind)[ut]],
  column = rownames(cor_matrix_ind)[col(cor_matrix_ind)[ut]],
  cor = (cor_matrix_ind)[ut]
)
hist(abs(cor_output_ind$cor),main="Histogram of correlations across lead features")
print("Summary of correlations across selected features")
summary(abs(cor_output_ind$cor))

# write out list of features with >set % missing to run through presence absence test
high_missing_features <- names(which(feature_missing >= miss_thresh))
# remove .e from metabolite name 
high_missing_features <- gsub("[.]e","",high_missing_features)
print(paste0("No. of features with >= ",miss_thresh*100,"% missing at the 36mth timepoint:"))
length(high_missing_features)

write.table(as.data.frame(high_missing_features),file=paste0(results_dir,"high_missingness_feature_list.txt"),row.names=F,col.names=F,quote=F,sep="\t")

# write out cluster summary
write.table(k_group_summary,file=paste0(results_dir,"clustering_summary.txt"),row.names=F,quote=F,sep="\t")


##################
## QUIT SESSION ##
##################

# close pdf
invisible(dev.off())

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()


