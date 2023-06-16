# This script makes results tables from analyses relating to the EFFECT OF THE INTERVENTION ON METABOLITES

# last run: 13th June 2023

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/08.0_make_tables.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/08.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )
rawdata_dir <- as.character( pfile[5,2] )

# load libraries
library("stargazer")
library("ggplot2")

####################################################################################
####################################################################################
## read in data files
####################################################################################
####################################################################################

# read in trial data
trial_data <- readRDS(paste0(rawdata_dir,"DiRECT_fu1_20171113.rds"))

# read in data used as model input
model_data_list <- readRDS(file=paste0(processeddata_dir,"model_input.rds"))

# separate metabolites from other pheno
raw_data_list <- list(metab=model_data_list$metab[,17:ncol(model_data_list$metab)],nmr=model_data_list$nmr[,17:ncol(model_data_list$nmr)])
model_data_only <- list(metab=model_data_list$metab[,1:16],nmr=model_data_list$nmr[,1:16])

# read in results files
lm_results_list <- readRDS(file=paste0(results_dir,"lm_results_list.rds"))
metab_lm_results <- lm_results_list$metab
metab_pa_results <- readRDS(file=paste0(results_dir,"metab_pa_results.rds"))
nmr_lm_results <- lm_results_list$nmr

# read in table containing labels
class_info <- readRDS(file=paste0(results_dir,"class_info_for_enrichment.rds"))

## read in clustering groups
clusters <- read.table(file=paste0(results_dir,"feature_group_allocation.txt"),sep="\t", stringsAsFactors = F,h=T)

## bring in feature lists from selection step
pa_features <- as.character((read.table(file=paste0(results_dir,"high_missingness_feature_list.txt"),sep="\t", stringsAsFactors = F))$V1)

## bring in info about clusters
k_group_summary <- read.table(file=paste0(results_dir,"clustering_summary.txt"),h=T,sep="\t")

# read in by cluster enrichment results
clusters_with_enrich <- readRDS(file=paste0(results_dir,"clusters_with_enrich.rds"))

# read in list of selected associated features 
selected_assoc_features <- read.table(file=paste0(results_dir,"selected_assoc_features.txt"),h=F, stringsAsFactors = F)
selected_assoc_features <- selected_assoc_features$V1

# read in lists of derived versus raw 
derived_measures <- read.table(file=paste0(results_dir,"NMRderived.txt"),h=F)
raw_measures <- read.table(file=paste0(results_dir,"NMRraw.txt"),h=F)

derived_measures$fullname <- class_info$biochemical[match(derived_measures$V1,class_info$feature_id)]
raw_measures$fullname <- class_info$biochemical[match(raw_measures$V1,class_info$feature_id)]

# read in feature metadata for metab features
metab_feature_metadata <- readRDS(paste0(processeddata_dir,"OrigScale_feature_metadata_with_QCstats.rds"))
nmr_feature_metadata <- readRDS(paste0(processeddata_dir,"nmr_feature_metadata_with_QCstats.rds"))

####################################################################################
## Table 1 - cohort characteristics
####################################################################################

# add data fields
model_data_only$metab$`Fasting glucose (mmol/L)` <- trial_data$glucose.mmol.l.b[match(model_data_only$metab$id,trial_data$id)]
model_data_only$metab$`Total cholesterol (mmol/L)` <- trial_data$chol.mmol.l.b[match(model_data_only$metab$id,trial_data$id)]
model_data_only$metab$`HDL cholesterol (mmol/L)` <- trial_data$hdl.mmol.l.b[match(model_data_only$metab$id,trial_data$id)]
model_data_only$metab$`Triglycerides (mmol/L)` <- trial_data$trig.mmol.l.b[match(model_data_only$metab$id,trial_data$id)]

# drop those not needed
model_data_only$metab <- model_data_only$metab[,c(1,4,6,7,8,11,17:20)]

# rename
names(model_data_only$metab)[3] <- "Age (years)"
names(model_data_only$metab)[5] <- "BMI (kg/m2)"
names(model_data_only$metab)[6] <- "Weight (kg)"

# count by class
table(model_data_only$metab$treat)
table(model_data_only$nmr$treat)

# restrict to intervention group
intervention_summary <- model_data_only$metab[model_data_only$metab$treat == "Intervention",]
filename <- paste0(results_dir,"Tables\\Table1_intervention.html")
stargazer(intervention_summary, out = filename, summary.stat = c("mean","sd"), digits=1)
# sex distribution
table(intervention_summary$sex)

# restrict to control group
control_summary <- model_data_only$metab[model_data_only$metab$treat == "Control",]
filename <- paste0(results_dir,"Tables\\Table1_control.html")
stargazer(control_summary, out = filename, summary.stat = c("mean","sd"), digits=1)
# sex distribution
table(control_summary$sex)


####################################################################################
## Table S5 - class specifications for enrichment analyis plus cluster membership
####################################################################################

# generate copy to edit
new_class_info <- class_info

# rename column headers
names(new_class_info)[1:5] <- c("Feature_ID","NMR_class","Metabolon_SuperPathway","Biochemical_name", "Source")

# add subpathway
new_class_info$Metabolon_SubPathway <- metab_feature_metadata$sub.pathway[match(new_class_info$Feature_ID,metab_feature_metadata$compid_to_match)]
  
# re-order columns 
new_class_info <- new_class_info[,c(1,4,5,2,3,6)]

# add cluster membership
new_class_info$cluster.group <- clusters$k_group[match(new_class_info$Feature_ID,clusters$feature_id)]

# write out 
write.table(new_class_info[,-c(1)],file=paste0(results_dir,"Tables\\TableS5_class_and_cluster_info.txt"),row.names=F,sep="\t")


####################################################################################
## Table S2 & S3 - NMR and Metabolon linear model full results
####################################################################################

## iterate over metabolon and nmr
lm_results_tables <- list(metab=NA, nmr=NA)

for (i in 1:length(lm_results_list)) {
  print(paste0("Processing data from: ", names(lm_results_list)[[i]]))
  # check model fails
  print("Linear model failures:")
  print(table(lm_results_list[[i]]$rnt_lm_model_fail_info,useNA="always"))
  # restrict results to columns of interest
  lm_results_tables[[i]] <- lm_results_list[[i]][,!grepl("*_warning",names(lm_results_list[[i]]))] # drop mixed model warnings
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("pathway.sortorder",names(lm_results_tables[[i]]))] # drop sort order var
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("comp.id",names(lm_results_tables[[i]]))] # drop id
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("chemical.id",names(lm_results_tables[[i]]))] # drop id
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("natlog_*",names(lm_results_tables[[i]]))] # drop nat log output
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("pa_*",names(lm_results_tables[[i]]))] # drop presence absence output
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("percent_*",names(lm_results_tables[[i]]))] # drop % success
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("*Qu.*",names(lm_results_tables[[i]]))] # drop quartiles
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("*var*",names(lm_results_tables[[i]]))] # drop variances
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("rnt_lm_n_intervention",names(lm_results_tables[[i]]))] #  drop lm extra info
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("rnt_lm_n_control",names(lm_results_tables[[i]]))] #  drop lm extra info
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("rnt_lm_n_imputed_baselines",names(lm_results_tables[[i]]))] #  drop lm extra info
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("*model_fail*",names(lm_results_tables[[i]]))] # drop model fails from lm
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("win_*",names(lm_results_tables[[i]]))] # drop winsorization fc results
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("raw_fc_mean_*",names(lm_results_tables[[i]]))] # drop mean based fc stats
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("*meanfold*",names(lm_results_tables[[i]]))] # drop mean based fc stats
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("raw_fc_geom_*",names(lm_results_tables[[i]]))] # drop geometric mean based fc stats
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("*geomfold*",names(lm_results_tables[[i]]))] # drop geometric mean based fc stats
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("high_missingness_flag",names(lm_results_tables[[i]]))] 
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("id_for_matching",names(lm_results_tables[[i]]))] 
  
  # add flag for if metabolite was selected as a rep feature
  lm_results_tables[[i]]$representative_feature <- "No"
  lm_results_tables[[i]][which(lm_results_tables[[i]]$feature_id %in% selected_assoc_features), "representative_feature"] <- "Yes"
  lm_results_tables[[i]][which(lm_results_tables[[i]]$rnt_lm_treat_HolmAdjP_flag == 0), "representative_feature"] <- NA
  
  # summarise selected
  print(paste0("No. of associated features that are selected for follow up based on clustering analysis:"))
  print(table(lm_results_tables[[i]]$rnt_lm_treat_HolmAdjP_flag,lm_results_tables[[i]]$representative_feature))
  
  # add in 95% CI for rnt betas
  lm_results_tables[[i]]$rnt_lm_treat_beta_lci <- lm_results_tables[[i]]$rnt_lm_treat_beta - (lm_results_tables[[i]]$rnt_lm_treat_SE * 1.96)
  lm_results_tables[[i]]$rnt_lm_treat_beta_uci <- lm_results_tables[[i]]$rnt_lm_treat_beta + (lm_results_tables[[i]]$rnt_lm_treat_SE * 1.96)
  lm_results_tables[[i]]$rnt_lm_treat_beta_lci_after_adj_for_weight.change <- lm_results_tables[[i]]$rnt_lm_treat_beta_after_adj_for_weight.change - (lm_results_tables[[i]]$rnt_lm_treat_SE_after_adj_for_weight.change * 1.96)
  lm_results_tables[[i]]$rnt_lm_treat_beta_uci_after_adj_for_weight.change <- lm_results_tables[[i]]$rnt_lm_treat_beta_after_adj_for_weight.change + (lm_results_tables[[i]]$rnt_lm_treat_SE_after_adj_for_weight.change * 1.96)
  
  # sort by p
  lm_results_tables[[i]] <- lm_results_tables[[i]][order(lm_results_tables[[i]]$rnt_lm_treat_p_rank),]
  
  # drop feature id
  lm_results_tables[[i]] <- lm_results_tables[[i]][,!grepl("feature_id",names(lm_results_tables[[i]]))] 
  
  # no. of associated features after weight adj
  print(paste0("No. of associated features from weight adjusted model after Holm correction: ", length(which(lm_results_tables[[i]]$rnt_lm_treat_HolmAdjp_after_adj_for_weight.change < 0.05 ))))
  
  # format (small) p values to scientific notation
  lm_results_tables[[i]][,c("rnt_lm_treat_p")] <- format(lm_results_tables[[i]][,c("rnt_lm_treat_p")],format="e",digits=3)
  lm_results_tables[[i]][,c("rnt_lm_treat_p_after_adj_for_weight.change")] <- format(lm_results_tables[[i]][,c("rnt_lm_treat_p_after_adj_for_weight.change")],format="e",digits=3)
  lm_results_tables[[i]][,c("rnt_lm_treat_HolmAdjP")] <- format(lm_results_tables[[i]][,c("rnt_lm_treat_HolmAdjP")],format="e",digits=3)
  #lm_results_tables[[i]][,c("rnt_lm_treat_HolmAdjp_after_adj_for_weight.change")] <- format(lm_results_tables[[i]][,c("rnt_lm_treat_HolmAdjp_after_adj_for_weight.change")],format="e",digits=3)
  
  if (i==1){
    # add panel
    lm_results_tables[[i]]$source <- new_class_info$Source[match(lm_results_tables[[i]]$biochemical,new_class_info$Biochemical_name)]
    # add cluster membership
    lm_results_tables[[i]]$cluster.group <- new_class_info$cluster.group[match(lm_results_tables[[i]]$biochemical,new_class_info$Biochemical_name)]
    # add superpathway
    lm_results_tables[[i]]$super.pathway <- new_class_info$Metabolon_SuperPathway[match(lm_results_tables[[i]]$biochemical,new_class_info$Biochemical_name)]
    # add subpathway
    lm_results_tables[[i]]$sub.pathway <- metab_feature_metadata$sub.pathway[match(lm_results_tables[[i]]$biochemical,metab_feature_metadata$biochemical)]
  } else {
    # add panel
    lm_results_tables[[i]]$source <- new_class_info$Source[match(lm_results_tables[[i]]$fullname,new_class_info$Biochemical_name)]
    # add cluster membership
    lm_results_tables[[i]]$cluster.group <- new_class_info$cluster.group[match(lm_results_tables[[i]]$fullname,new_class_info$Biochemical_name)]
    # add superpathway
    lm_results_tables[[i]]$super.pathway <- new_class_info$Metabolon_SuperPathway[match(lm_results_tables[[i]]$fullname,new_class_info$Biochemical_name)]
    # add subpathway
    lm_results_tables[[i]]$sub.pathway <- NA
  }
  
  # replace eta2 with VarExp
  names(lm_results_tables[[i]]) <- gsub(x=names(lm_results_tables[[i]]),pattern="eta2",replacement = "VarExp")
  
  # calculate change in variance explained after weight adjustment
  lm_results_tables[[i]]$Percent.change.in.VE.by.allocation.after.weight.adj <- 
    100*((lm_results_tables[[i]]$rnt_lm_VarExp_treat_after_adj_for_weight.change - lm_results_tables[[i]]$rnt_lm_VarExp_treat)/lm_results_tables[[i]]$rnt_lm_VarExp_treat)

  ## set weight adj results to NA if metabolite not associated in primary model
  # find weight adj columns
  weightadj.cols <- which(grepl("weight.",names(lm_results_tables[[i]])))
  # set to NA
  # set col to use as flag
  flag_col <- which(names(lm_results_tables[[i]]) == "rnt_lm_treat_HolmAdjP_flag")
  
  for (m in 1:nrow(lm_results_tables[[i]])){
    if (lm_results_tables[[i]][m,flag_col] == 0) {
      lm_results_tables[[i]][m,weightadj.cols] <- NA
    }
  }
  
  # post weight adjustment summary
  print("Number of associated metabolites for which allocation VE decreases (or is unchanged) after adjusting for weight change:")
  print(length(which(lm_results_tables[[i]]$rnt_lm_VarExp_treat >= lm_results_tables[[i]]$rnt_lm_VarExp_treat_after_adj_for_weight.change)))
  print("Number of associated metabolites for which allocation VE increases after adjusting for weight change:")
  print(length(which(lm_results_tables[[i]]$rnt_lm_VarExp_treat < lm_results_tables[[i]]$rnt_lm_VarExp_treat_after_adj_for_weight.change)))
  
  print("Number of associated metabolites for which allocation VE is greater than or equal to weight change VE:")
  print(length(which(lm_results_tables[[i]]$rnt_lm_VarExp_treat_after_adj_for_weight.change >= lm_results_tables[[i]]$rnt_lm_VarExp_weight.change)))
  print("Number of associated metabolites for which allocation VE is less than weight change VE:")
  print(length(which(lm_results_tables[[i]]$rnt_lm_VarExp_treat_after_adj_for_weight.change < lm_results_tables[[i]]$rnt_lm_VarExp_weight.change)))
  
  # limit numeric cols to 3dp
  num_cols <- which(unlist(lapply(lm_results_tables[[i]], is.numeric)))
  for (j in 1:length(num_cols)){
    lm_results_tables[[i]][,num_cols[j]] <- round(lm_results_tables[[i]][,num_cols[j]],3)
  }
  
  # check column names
  print(names(lm_results_tables[[i]]))
    
  # reorder
  if (i==1){
    lm_results_tables[[i]] <-  lm_results_tables[[i]][,c(1,55,56,2:8,57,58,22:24,51,52,25:35,53,54,36:39,59,40:50)]
  } else {
    lm_results_tables[[i]] <-  lm_results_tables[[i]][,c(1,49,2,50,51,52,16:18,45,46,19:29,47,48,30:33,53,34:44)]
  }
  
  ## extract some summary info for text
  assoc_subset <- lm_results_tables[[i]][lm_results_tables[[i]]$rnt_lm_treat_HolmAdjP_flag == 1,]
  
  # no. of features
  print(paste0("No. of features in LM results: ", nrow(lm_results_tables[[i]]) ))
  # no. of associated features
  print(paste0("No. of associated features after Holm correction: ", nrow(assoc_subset) ))
  # no. increased/decreased
  print(paste0("No. of associated features increased in intervention group: ", nrow(assoc_subset[assoc_subset$rnt_lm_treat_beta > 0,]) ))
  print(paste0("No. of associated features decreased in intervention group: ", nrow(assoc_subset[assoc_subset$rnt_lm_treat_beta < 0,]) ))
  # distribution by superpathway and direction of effect
  print("By class distribution for increasing metabolites:")
  print(table(assoc_subset[assoc_subset$rnt_lm_treat_beta > 0,"super.pathway"]))
  print("By class distribution for decreasing metabolites:")
  print(table(assoc_subset[assoc_subset$rnt_lm_treat_beta < 0,"super.pathway"]))
  # data regarding sample size
  print("Sample size summary:")
  print(summary(lm_results_tables[[i]]$rnt_lm_n_samples_in_model))
  # calculate 10% missing threshold
  thresh <- 0.9*max(lm_results_tables[[i]]$rnt_lm_n_samples_in_model)
  print(paste0("The 10% missingness threshold is ", thresh , " samples."))
  print("The number of metabolites with missingness less than 10% (ie more samples than this) is: ")
  print(nrow(lm_results_tables[[i]][lm_results_tables[[i]]$rnt_lm_n_samples_in_model > thresh,]))
  # for NMR only, add print statement re. no. of associated that are derived
  print(paste0("No. of associated features classed as derived (NMR): ", nrow(assoc_subset[assoc_subset$fullname %in% derived_measures$fullname,]) ))
  print(paste0("No. of associated features classed as raw (NMR): ", nrow(assoc_subset[assoc_subset$fullname %in% raw_measures$fullname,]) ))
}

# save out required var
# convert to reader friendly variable names
names(lm_results_tables$nmr)
# add shortname
lm_results_tables$nmr$shortname <- nmr_feature_metadata$shortname[match(lm_results_tables$nmr$fullname,nmr_feature_metadata$fullname)]
names(lm_results_tables$nmr) <- c("Biochemical.name","Source.platform","Units.of.measurement","Cluster.group","Super.pathway","Sub.pathway",
                                  "N.samples.in.model","Allocation.beta","Allocation.SE","Allocation.lower.95CI","Allocation.upper.95CI","Allocation.tstat","Allocation.p",
                                  "VE.age","VE.sex","VE.centre","VE.list.size","VE.metab.at.baseline","VE.allocation","VE.residuals",
                                  "Allocation.beta.with.weight.adj","Allocation.SE.with.weight.adj","Allocation.lower.95CI.with.weight.adj","Allocation.upper.95CI.with.weight.adj",
                                  "Allocation.t.with.weight.adj","Allocation.p.with.weight.adj","Allocation.VE.with.weight.adj","VE.weight.change","Percent.change.in.VE.by.allocation.after.weight.adj",
                                  "N.samples.in.FC.analysis","N.controls.in.FC.analysis","N.intervention.in.FC.analysis","Median.in.controls","Median.in.intervention",
                                  "FC.of.median.intervention.over.control","log2.FC.of.median","Allocation.p.Holm.adj","Allocation.assoc.flag","Allocation.p.rank",
                                  "Representative.feature.?","shortname")

filename <- paste0(results_dir,"Tables\\TableS2_nmr_lm_results.txt")
write.table(lm_results_tables$nmr[,c(1,ncol(lm_results_tables$nmr),2:(ncol(lm_results_tables$nmr)-1))],file=filename,sep="\t",row.names=F,quote=F)

# convert to reader friendly variable names
names(lm_results_tables$metab)
# add chem ids
lm_results_tables$metab$Metabolon.chem.id <- metab_feature_metadata$chemical.id[match(lm_results_tables$metab$biochemical,metab_feature_metadata$biochemical)]
names(lm_results_tables$metab) <- c("Biochemical.name","Source.platform","Cluster.group","Assay.run","Retention.index","Metabolite.mass",
                                    "CAS.number","Pubchem.id","Kegg.id","HMDB.id","Super.pathway","Sub.pathway",
                                  "N.samples.in.model","Allocation.beta","Allocation.SE","Allocation.lower.95CI","Allocation.upper.95CI","Allocation.tstat","Allocation.p",
                                  "VE.age","VE.sex","VE.centre","VE.list.size","VE.metab.at.baseline","VE.allocation","VE.residuals",
                                  "Allocation.beta.with.weight.adj","Allocation.SE.with.weight.adj","Allocation.lower.95CI.with.weight.adj","Allocation.upper.95CI.with.weight.adj",
                                  "Allocation.t.with.weight.adj","Allocation.p.with.weight.adj","Allocation.VE.with.weight.adj","VE.weight.change","Percent.change.in.VE.by.allocation.after.weight.adj",
                                  "N.samples.in.FC.analysis","N.controls.in.FC.analysis","N.intervention.in.FC.analysis","Median.in.controls","Median.in.intervention",
                                  "FC.of.median.intervention.over.control","log2.FC.of.median","Allocation.p.Holm.adj","Allocation.assoc.flag","Allocation.p.rank",
                                  "Representative.feature.?","Metabolon.chem.id")
filename <- paste0(results_dir,"Tables\\TableS3_metab_lm_results.txt")
write.table(lm_results_tables$metab[,c(1,ncol(lm_results_tables$metab),2:(ncol(lm_results_tables$metab)-1))],file=filename,sep="\t",row.names=F,quote=F)

####################################################################################
## Run some summaries of weight adjusted analysis results (in lm_results)
####################################################################################

# How many metabolites had reduced VE by allocation when weight change added?
n_reduced <- length(which(lm_results_tables$nmr$Percent.change.in.VE.by.allocation.after.weight.adj < 0 & lm_results_tables$nmr$Allocation.assoc.flag == 1))
n_total <- length(which( lm_results_tables$nmr$Allocation.assoc.flag == 1))
print("Percent of associated NMR metabolites that saw a reduction in VE by allocation after weight change fitted:")
100*(n_reduced/n_total)
n_reduced
n_total

n_reduced <- length(which(lm_results_tables$metab$Percent.change.in.VE.by.allocation.after.weight.adj < 0 & lm_results_tables$metab$Allocation.assoc.flag == 1))
n_total <- length(which( lm_results_tables$metab$Allocation.assoc.flag == 1))
print("Percent of associated MS metabolites that saw a reduction in VE by allocation after weight change fitted:")
100*(n_reduced/n_total)
n_reduced
n_total

# How many metabolites had greater VE by allocation than weight change?
lm_results_tables$nmr$VE_diff <- lm_results_tables$nmr$Allocation.VE.with.weight.adj - lm_results_tables$nmr$VE.weight.change
summary(lm_results_tables$nmr$VE_diff)
n_win <- length(which(lm_results_tables$nmr$VE_diff >= 0 & lm_results_tables$nmr$Allocation.assoc.flag == 1))
n_total <- length(which(lm_results_tables$nmr$Allocation.assoc.flag == 1))
print("Percent of associated NMR metabolites that had greater VE by allocation than weight change:")
100*(n_win/n_total)
n_win
n_total

lm_results_tables$metab$VE_diff <- lm_results_tables$metab$Allocation.VE.with.weight.adj - lm_results_tables$metab$VE.weight.change
summary(lm_results_tables$metab$VE_diff)
n_win <- length(which(lm_results_tables$metab$VE_diff >= 0 & lm_results_tables$metab$Allocation.assoc.flag == 1))
n_total <- length(which(lm_results_tables$metab$Allocation.assoc.flag == 1))
print("Percent of associated MS metabolites that had greater VE by allocation than weight change:")
100*(n_win/n_total)
n_win
n_total

####################################################################################
## Table S9 - Metabolon PA full results
####################################################################################

# check model fails
print("Model fails in logistic regression:")
table(metab_pa_results$pa_model_fail_info,useNA="always")

# restrict results to columns of interest
metab_pa_output <- metab_pa_results[,!grepl("*_warning",names(metab_pa_results))] # drop mixed model warnings
metab_pa_output <- metab_pa_output[,!grepl("pathway.sortorder",names(metab_pa_output))] # drop sort order var
metab_pa_output <- metab_pa_output[,!grepl("comp.id",names(metab_pa_output))] # drop id
metab_pa_output <- metab_pa_output[,!grepl("natlog_*",names(metab_pa_output))] # drop nat log output
metab_pa_output <- metab_pa_output[,!grepl("rnt_*",names(metab_pa_output))] # drop rnt output
metab_pa_output <- metab_pa_output[,!grepl("percent_*",names(metab_pa_output))] # drop % success
metab_pa_output <- metab_pa_output[,!grepl("*Qu.*",names(metab_pa_output))] # drop quartiles
metab_pa_output <- metab_pa_output[,!grepl("*var*",names(metab_pa_output))] # drop variances
metab_pa_output <- metab_pa_output[,!grepl("win_*",names(metab_pa_output))] # drop winsorization fc results
metab_pa_output <- metab_pa_output[,!grepl("raw_fc_mean_*",names(metab_pa_output))] # drop mean based fc stats
metab_pa_output <- metab_pa_output[,!grepl("*meanfold*",names(metab_pa_output))] # drop mean based fc stats
metab_pa_output <- metab_pa_output[,!grepl("raw_fc_geom_*",names(metab_pa_output))] # drop geometric mean based fc stats
metab_pa_output <- metab_pa_output[,!grepl("*geomfold*",names(metab_pa_output))] # drop geometric mean based fc stats
metab_pa_output <- metab_pa_output[,!grepl("high_missingness_flag",names(metab_pa_output))] 
metab_pa_output <- metab_pa_output[,!grepl("id_for_matching",names(metab_pa_output))] 

# add in 95% CI for rnt beta
metab_pa_output$pa_logit_treat_beta_lci <- metab_pa_output$pa_logit_treat_beta - (metab_pa_output$pa_logit_treat_SE * 1.96)
metab_pa_output$pa_logit_treat_beta_uci <- metab_pa_output$pa_logit_treat_beta + (metab_pa_output$pa_logit_treat_SE * 1.96)

# sort by p
metab_pa_output <- metab_pa_output[order(metab_pa_output$pa_logit_treat_p_rank),]

# format p values to scientific notation
metab_pa_output[,c("pa_logit_treat_p")] <- format(metab_pa_output[,c("pa_logit_treat_p")],format="e",digits=3)

# restrict to 2dp
num_cols <- which(unlist(lapply(metab_pa_output, is.numeric)))
for (j in 1:length(num_cols)){
  metab_pa_output[,num_cols[j]] <- round(metab_pa_output[,num_cols[j]],2)
}

# drop feature id
metab_pa_output <- metab_pa_output[,!grepl("feature_id",names(metab_pa_output))] 

# summary stats for manuscript text
assoc_subset <- metab_pa_output[metab_pa_output$pa_logit_treat_p_flag == 1,]
# no. of features
print(paste0("No. of features in MS logistic model results: ", nrow(metab_pa_output) ))
# no. of associated features
print(paste0("No. of associated features (p<0.05): ", nrow(assoc_subset) ))
# no. increased/decreased
print(paste0("No. of associated features enriched/depleted in intervention group: "))
print(table(assoc_subset$state_in_intervention))

# restrict to vars for output
names(metab_pa_output)
metab_pa_output_to_export <- metab_pa_output[,c(1:11,25:31,45,46,32:33,42:44)]

# convert to reader friendly variable names
names(metab_pa_output_to_export)
names(metab_pa_output_to_export) <- c("Biochemical.name","Super.pathway","Sub.pathway","Assay.run","Metabolon.chem.id","Retention.index","Metabolite.mass",
                                    "CAS.number","Pubchem.id","Kegg.id","HMDB.id",
                                    "N.samples.in.model","Count.at.baseline.in.control","Count.at.1year.in.control","Count.at.baseline.in.intervention","Count.at.1year.in.intervention",
                                    "Allocation.beta","Allocation.SE","Allocation.lower.95CI","Allocation.upper.95CI","Allocation.zstat","Allocation.p",
                                    "State.in.intervention.group","Allocation.assoc.flag","Allocation.p.rank")

# save out required var
filename <- paste0(results_dir,"Tables\\TableS9_metab_pa_results.txt")
write.table(metab_pa_output_to_export,file=filename,sep="\t",row.names=F,quote=F)

###################################################################################
## Table S6 - associated representative from LM
####################################################################################

## metabolon ##
# start with suppl table version
main_metab_assoc <- lm_results_tables$metab

# restrict to representative and associated
main_metab_assoc <- main_metab_assoc[main_metab_assoc$Allocation.assoc.flag == 1 & main_metab_assoc$`Representative.feature.?` == "Yes",]

# restrict to relevant results columns
main_metab_assoc <- main_metab_assoc[,c("Source.platform", "Biochemical.name", "Super.pathway", "Sub.pathway","Cluster.group", "Allocation.beta", "Allocation.lower.95CI", 
                                        "Allocation.upper.95CI","Allocation.p.Holm.adj","Percent.change.in.VE.by.allocation.after.weight.adj")]

# add info re. enrichment/clustering
main_metab_assoc$No.in.cluster <- k_group_summary$n_features_in_group[match(main_metab_assoc$Cluster.group,k_group_summary$k_group)]
# if only one metabolite in cluster, wont appear in k_group_summary so need to insert no. in cluster (which will be 1) 
main_metab_assoc$No.in.cluster <- replace(main_metab_assoc$No.in.cluster, is.na(main_metab_assoc$No.in.cluster), 1)


## NMR ##
# start with suppl table version
main_nmr_assoc <- lm_results_tables$nmr

# restrict to representative and associated
main_nmr_assoc <- main_nmr_assoc[main_nmr_assoc$Allocation.assoc.flag == 1 & main_nmr_assoc$`Representative.feature.?` == "Yes",]
main_nmr_assoc <- main_nmr_assoc[,c("Source.platform", "Biochemical.name", "Super.pathway", "Sub.pathway","Cluster.group", "Allocation.beta", "Allocation.lower.95CI", 
                                    "Allocation.upper.95CI","Allocation.p.Holm.adj","Percent.change.in.VE.by.allocation.after.weight.adj")]

# add info re. enrichment/clustering
main_nmr_assoc$No.in.cluster <- k_group_summary$n_features_in_group[match(main_nmr_assoc$Cluster.group,k_group_summary$k_group)]
# if only one metabolite in cluster, wont appear in k_group_summary so need to insert no. in cluster (which will be 1) 
main_nmr_assoc$No.in.cluster <- replace(main_nmr_assoc$No.in.cluster, is.na(main_nmr_assoc$No.in.cluster), 1)


# extract info on results for others in class
# metab
lm_results_list$metab$cluster.group <- new_class_info$cluster.group[match(lm_results_list$metab$biochemical,new_class_info$Biochemical_name)]
main_metab_assoc$Proportion.nom.assoc.in.cluster <- 0
for (i in 1:nrow(main_metab_assoc)) {
  class.no <- main_metab_assoc[i,c("Cluster.group")]
  restricted_set1 <- lm_results_list$metab[lm_results_list$metab$cluster.group == class.no,]
  restricted_set2 <- lm_results_list$nmr[lm_results_list$nmr$cluster.group == class.no,]
  total_length <- length(restricted_set1[as.numeric(restricted_set1$rnt_lm_treat_p) < 0.05,c("feature_id")]) + 
                    length(restricted_set2[as.numeric(restricted_set2$rnt_lm_treat_p) < 0.05,c("feature_id")])
  main_metab_assoc[i,ncol(main_metab_assoc)] <- round(total_length/main_metab_assoc[i,c("No.in.cluster")],2)
}

# nmr
lm_results_list$nmr$cluster.group <- new_class_info$cluster.group[match(lm_results_list$nmr$feature_id,new_class_info$Feature_ID)]
main_nmr_assoc$Proportion.nom.assoc.in.cluster <- 0
for (i in 1:nrow(main_nmr_assoc)) {
  class.no <- main_nmr_assoc[i,c("Cluster.group")]
  restricted_set1 <- lm_results_list$nmr[lm_results_list$nmr$cluster.group == class.no,]
  restricted_set2 <- lm_results_list$metab[lm_results_list$metab$cluster.group == class.no,]
  total_length <- length(restricted_set1[as.numeric(restricted_set1$rnt_lm_treat_p) < 0.05,c("feature_id")]) + 
                    length(restricted_set2[as.numeric(restricted_set2$rnt_lm_treat_p) < 0.05,c("feature_id")])
  main_nmr_assoc[i,ncol(main_nmr_assoc)] <- round(total_length/main_nmr_assoc[i,c("No.in.cluster")],2)
}

## merge the two sets of results
main_assoc <- rbind(main_metab_assoc,main_nmr_assoc)

## add in info. re. cluster enrichment
main_assoc$Cluster.super.pathway <- NA
main_assoc$Cluster.super.pathway.fold.enrich <- NA
main_assoc$Cluster.super.pathway.enrichP <- NA

## loop over assoc features and extract cluster enrichment details
for (i in 1:nrow(main_assoc)) {
  # extract cluster no.
  cluster_no <- main_assoc[i,"Cluster.group"]
  # extract n in cluster
  n_in_cluster <- main_assoc[i,"No.in.cluster"]
  # extract relevant cluster enrichment results
  cluster_enrich <- as.data.frame(clusters_with_enrich[[cluster_no]])
  # restrict to case where more than 5 features in cluster
  if (n_in_cluster > 5) {
    # restrict to the case where enrichment data exist (if super-pathways of all features in cluster are unknown then no enrichment results)
    if (ncol(cluster_enrich) > 0) {
      # sort by enrich p's
      cluster_enrich_sorted <- cluster_enrich[order(cluster_enrich$Ftest_pval),]
      # extract top class and enrich results
      main_assoc[i,13] <- row.names(cluster_enrich_sorted)[1]
      main_assoc[i,14] <- cluster_enrich_sorted[1,1]
      main_assoc[i,15] <- cluster_enrich_sorted[1,2]
    }
    else {
      main_assoc[i,13] <- NA
      main_assoc[i,14] <- NA
      main_assoc[i,15] <- NA
    }
  }
  else {
    main_assoc[i,13] <- NA
    main_assoc[i,14] <- NA
    main_assoc[i,15] <- NA   
  }
}

# sort by p
main_assoc <- main_assoc[order(as.numeric(main_assoc$Allocation.p.Holm.adj)),]

# re-order columns
main_assoc <- main_assoc[,c("Source.platform", "Biochemical.name", "Super.pathway","Sub.pathway","Cluster.group","No.in.cluster","Allocation.beta",
                            "Allocation.lower.95CI","Allocation.upper.95CI", "Allocation.p.Holm.adj", "Percent.change.in.VE.by.allocation.after.weight.adj",
                            "Proportion.nom.assoc.in.cluster","Cluster.super.pathway","Cluster.super.pathway.fold.enrich","Cluster.super.pathway.enrichP")]

dim(main_assoc)

# restrict to identified metabolites (ie exclude Metabolon X's)
print("Distribution of selected associated metabolites:")
table(main_assoc$Super.pathway,useNA="always")

# print statements re. no. of features from each source (nightingale/metabolon)
print("Distribution of associated representative by source:")
table(main_assoc$Source.platform)

# restrict to identified only
main_assoc_identified <- main_assoc[main_assoc$Super.pathway != "Unclassified",]

# save out required var
filename <- paste0(results_dir,"Tables\\TableS6_primary_results.txt")
write.table(main_assoc_identified,file=filename,sep="\t",row.names=F,quote=F)


##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
rm(list=ls())
