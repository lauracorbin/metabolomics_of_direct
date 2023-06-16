# This script processes results from analyses relating to the EFFECT OF THE INTERVENTION ON METABOLITES

# last run: 13th June 2023

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/06.0_results_process.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/06.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )

## call R script containing functions
source(paste0(scripts_dir,"04.0_functions.R"))
source(paste0(scripts_dir,"02.0_qc_functions.R"))

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in  datasets
# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))

# read in results files
results_list <- list(metab=readRDS(file=paste0(results_dir,"metabolon_results.rds")),nmr=readRDS(file=paste0(results_dir,"nmr_results.rds")))

# bring in list of features with less missing than threshold
lm_features <- as.character((read.table(file=paste0(results_dir,"lm_feature_list.txt"), h=F, sep="\t"))$V1)

## bring in feature lists from selection step
pa_features <- as.character((read.table(file=paste0(results_dir,"high_missingness_feature_list.txt"),sep="\t", stringsAsFactors = F))$V1)

## read in clustering groups
clusters <- read.table(file=paste0(results_dir,"feature_group_allocation.txt"),sep="\t", stringsAsFactors = F)

# read in residuals from mixed model
metab_lm_resid <- readRDS(file=paste0(processeddata_dir,"metabolon_lm_rnt_resid.rds"))
nmr_lm_resid <- readRDS(file=paste0(processeddata_dir,"nmr_lm_rnt_resid.rds"))

####################################################################################
####################################################################################
## split metab results by primary analysis
####################################################################################
####################################################################################

# Metabolon - linear model
clinical_measures <- results_list$metab$feature_id[which(is.na(results_list$metab$comp.id))]
# create results subsets
metab_lm_results <- results_list$metab[results_list$metab$feature_id %in% lm_features,]
dim(metab_lm_results)
clinical_lm_results <- results_list$metab[results_list$metab$feature_id %in% clinical_measures,]
dim(clinical_lm_results)

# Metabolon - presence absence
metab_pa_results <- results_list$metab[results_list$metab$feature_id %in% pa_features,]
dim(metab_pa_results)

# NMR - mixed model
# create results subset (currently all go through lm pipeline)
nmr_lm_results <- results_list$nmr[results_list$nmr$feature_id %in% lm_features,]
dim(nmr_lm_results)

####################################################################################
####################################################################################
## summarise pa results  
####################################################################################
####################################################################################

length(which(metab_pa_results$pa_logit_treat_p < 0.01))
# add var to describe direction of effect
metab_pa_results$state_in_intervention <- NA
metab_pa_results[which(metab_pa_results$pa_logit_treat_beta <= 0), "state_in_intervention"] <- "depleted"
metab_pa_results[which(metab_pa_results$pa_logit_treat_beta > 0), "state_in_intervention"] <- "enriched"
metab_pa_results$pa_logit_treat_p_flag <- 0
metab_pa_results[which(metab_pa_results$pa_logit_treat_p < 0.01), "pa_logit_treat_p_flag"] <- 1
table(metab_pa_results$pa_logit_treat_p_flag,useNA="always")

# evaluate model failures
print(paste0("Of those in the PA analysis, " , length(which(metab_pa_results$pa_model_fail_info == "model failed to run")), " did not run."))
print(metab_pa_results[which(metab_pa_results$pa_model_fail_info == "model failed to run"), "biochemical"])

# add p-value based rank 
# sort by p
metab_pa_results <- metab_pa_results[order(metab_pa_results$pa_logit_treat_p),]
# add rank
metab_pa_results$pa_logit_treat_p_rank <- rank(metab_pa_results$pa_logit_treat_p)

pa_sig <- metab_pa_results[metab_pa_results$pa_logit_treat_p_flag == 1,]
table(pa_sig$state_in_intervention,pa_sig$super.pathway, useNA="ifany")

####################################################################################
####################################################################################
## do p-value adjustment for linear model output
####################################################################################
####################################################################################

lm_results_list <- list(metab=metab_lm_results, nmr=nmr_lm_results)
  
for (i in 1:length(lm_results_list)){
  print(paste0("Processing: ",names(lm_results_list)[[i]]))
  # adjust p for multiple testing using holm
  lm_results_list[[i]]$rnt_lm_treat_HolmAdjP <- p.adjust(lm_results_list[[i]]$rnt_lm_treat_p, method=c("holm"))
  #lm_results_list[[i]]$rnt_lm_treat_HolmAdjp_after_adj_for_weight.change <- p.adjust(lm_results_list[[i]]$rnt_lm_treat_p_after_adj_for_weight.change, method=c("holm"))
  # replace values >1 with =1
  lm_results_list[[i]][which(lm_results_list[[i]]$rnt_lm_treat_HolmAdjP > 1), "rnt_lm_treat_HolmAdjP"] <- 1
  #lm_results_list[[i]][which(lm_results_list[[i]]$rnt_lm_treat_HolmAdjp_after_adj_for_weight.change > 1), "rnt_lm_treat_HolmAdjp_after_adj_for_weight.change"] <- 1
  
  # flag associated based on BF p<0.05
  lm_results_list[[i]]$rnt_lm_treat_HolmAdjP_flag <- 0
  lm_results_list[[i]][which(lm_results_list[[i]]$rnt_lm_treat_HolmAdjP < (0.05) ),c("rnt_lm_treat_HolmAdjP_flag")] <- 1
  print(table(lm_results_list[[i]]$rnt_lm_treat_HolmAdjP_flag))
  
  # add p-value based rank 
  # sort by p
  lm_results_list[[i]] <- lm_results_list[[i]][order(lm_results_list[[i]]$rnt_lm_treat_p),]
  # add rank
  lm_results_list[[i]]$rnt_lm_treat_p_rank <- rank(lm_results_list[[i]]$rnt_lm_treat_p)
  # print out details of metabolite with smallest p
  print(head(lm_results_list[[i]],1))
}

# evaluate model failures
print(paste0("Of those Metabolon metabolites in the linear analysis, " , length(which(lm_results_list$metab$rnt_lm_model_fail_info == "model failed to run")), " did not run."))
print(lm_results_list[[i]][which(lm_results_list$metab$rnt_lm_model_fail_info == "model failed to run"), "biochemical"])

print(paste0("Of those NMR metabolites in the linear analysis, " , length(which(lm_results_list$nmr$rnt_lm_model_fail_info == "model failed to run")), " did not run."))
print(lm_results_list[[i]][which(lm_results_list$nmr$rnt_lm_model_fail_info == "model failed to run"), "biochemical"])

####################################################################################
####################################################################################
## save out modified results files
####################################################################################
####################################################################################

write.table(lm_results_list$metab,file=paste0(results_dir,"metabolon_lm_results.txt"), sep="\t", row.names=F,quote=F)
write.table(lm_results_list$nmr,file=paste0(results_dir,"nmr_lm_results.txt"), sep="\t", row.names=F,quote=F)
saveRDS(lm_results_list,file=paste0(results_dir,"lm_results_list.rds"))

write.table(metab_pa_results,file=paste0(results_dir,"metab_pa_results.txt"), sep="\t", row.names=F,quote=F)
saveRDS(metab_pa_results,file=paste0(results_dir,"metab_pa_results.rds"))

write.table(clinical_lm_results,file=paste0(results_dir,"clinical_lm_results.txt"), sep="\t", row.names=F,quote=F)
saveRDS(clinical_lm_results,file=paste0(results_dir,"clinical_lm_results.rds"))

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

