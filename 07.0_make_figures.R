# This script makes figures from analyses relating to the EFFECT OF THE INTERVENTION ON METABOLITES

# last run: 7th September 2023

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/07.0_make_figures.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/07.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )
rawdata_dir <- as.character( pfile[5,2] )
ext_data_dir <- as.character( pfile[6,2] )

## call R script containing functions
source(paste0(scripts_dir,"02.0_qc_functions.R"))
source(paste0(scripts_dir,"04.0_functions.R"))
source(paste0(scripts_dir,"07.0_functions.R"))


####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in  datasets
# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))

# read in complete results files
results_list <- list(metab=readRDS(file=paste0(results_dir,"metabolon_results.rds")),nmr=readRDS(file=paste0(results_dir,"nmr_results.rds")))

# read in results files separated by primary analysis (based on missingness)
lm_results_list <- readRDS(file=paste0(results_dir,"lm_results_list.rds"))
metab_pa_results <- readRDS(file=paste0(results_dir,"metab_pa_results.rds"))
clinical_lm_results <- readRDS(file=paste0(results_dir,"clinical_lm_results.rds"))

# bring in list of features with less missing than threshold
lm_features <- as.character((read.table(file=paste0(results_dir,"lm_feature_list.txt"), h=F, sep="\t"))$V1)

## bring in feature lists from selection step
pa_features <- as.character((read.table(file=paste0(results_dir,"high_missingness_feature_list.txt"),sep="\t", stringsAsFactors = F))$V1)

## read in clustering groups
clusters <- read.table(file=paste0(results_dir,"feature_group_allocation.txt"),sep="\t", stringsAsFactors = F,h=T)

# read in labels for NMR and attach
load(file=paste0(ext_data_dir,"metabolite_labels.Rda"))

# read in feature metadata for metab features
metab_feature_metadata <- readRDS(paste0(processeddata_dir,"OrigScale_feature_metadata_with_QCstats.rds"))
nmr_feature_metadata <- readRDS(paste0(processeddata_dir,"nmr_feature_metadata_with_QCstats.rds"))

# read in residuals from linear model for selected associated features
selected_resid <- readRDS(file = paste0(processeddata_dir,"resid_selected_assoc.rds")) 

# read in trial data
trial_data <- readRDS(paste0(rawdata_dir,"DiRECT_fu1_20171113.rds"))

# read in data used as model input
model_data_list <- readRDS(file=paste0(processeddata_dir,"model_input.rds"))

## read in var classes
date_vars <- read.table(file=paste0(processeddata_dir,"trial_date_vars.txt"),h=F)
integer_vars <- read.table(file=paste0(processeddata_dir,"trial_integer_vars.txt"),h=F)
numeric_vars <- read.table(file=paste0(processeddata_dir,"trial_numeric_vars.txt"),h=F)
factor_vars <- read.table(file=paste0(processeddata_dir,"trial_factor_vars.txt"),h=F)

# read in labels 
class_info <- readRDS(file=paste0(results_dir,"class_info_for_enrichment.rds"))

# read in list of associated representative features
selected_assoc_features <- read.table(file=paste0(results_dir,"selected_assoc_features.txt"),h=F)

# read in PVA results
pva_output <- read.table(file=paste0(results_dir,"pva_output_by_feature.txt"),sep="\t", stringsAsFactors = F,h=T)

# clinic data
clinic_data_only <- readRDS(file=paste0(processeddata_dir,"clinic_data_input.rds"))

## metabolite feature metadata (post cleaning)
clean_feature_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_feature_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_feature_metadata.rds"))) 

## add column for feature id matching
clean_feature_list[[1]]$id_for_matching <- clean_feature_list[[1]]$compid_to_match
clean_feature_list[[2]]$id_for_matching <- clean_feature_list[[2]]$shortname

####################################################################################
####################################################################################
## process data
####################################################################################
####################################################################################

## merge mm results for NMR and Metabolon
# extract details of features to merge from mm analysis
# find col to keep and extract
pval <- which(names(lm_results_list$metab) == "rnt_lm_treat_HolmAdjP")
followup <- which(names(lm_results_list$metab) == "rnt_lm_treat_HolmAdjP_flag")
foldchange <- which(names(lm_results_list$metab) == "raw_fc_log2medianfoldchange_int_over_control")
ranks <- which(names(lm_results_list$metab) == "rnt_lm_treat_p_rank")
metab_results_to_merge <- lm_results_list$metab[,c(1:14,pval,foldchange,followup,ranks)]

pval <- which(names(lm_results_list$nmr) == "rnt_lm_treat_HolmAdjP")
followup <- which(names(lm_results_list$nmr) == "rnt_lm_treat_HolmAdjP_flag")
foldchange <- which(names(lm_results_list$nmr) == "raw_fc_log2medianfoldchange_int_over_control")
ranks <- which(names(lm_results_list$nmr) == "rnt_lm_treat_p_rank")
nmr_results_to_merge <- lm_results_list$nmr[,c(1:8,pval,foldchange,followup,ranks)]

# merge across two panels
merged_lm_results <- merge(metab_results_to_merge,nmr_results_to_merge,by="feature_id",all=T)

# join common columns
merged_lm_results$raw_fc_log2medianfoldchange_int_over_control <- merged_lm_results$raw_fc_log2medianfoldchange_int_over_control.x
to_insert <- which(is.na(merged_lm_results$raw_fc_log2medianfoldchange_int_over_control))
merged_lm_results$raw_fc_log2medianfoldchange_int_over_control[to_insert] <- merged_lm_results$raw_fc_log2medianfoldchange_int_over_control.y[to_insert]

merged_lm_results$rnt_lm_treat_HolmAdjP <- merged_lm_results$rnt_lm_treat_HolmAdjP.x
to_insert <- which(is.na(merged_lm_results$rnt_lm_treat_HolmAdjP))
merged_lm_results$rnt_lm_treat_HolmAdjP[to_insert] <- merged_lm_results$rnt_lm_treat_HolmAdjP.y[to_insert]

merged_lm_results$rnt_lm_treat_HolmAdjP_flag <- merged_lm_results$rnt_lm_treat_HolmAdjP_flag.x
to_insert <- which(is.na(merged_lm_results$rnt_lm_treat_HolmAdjP_flag))
merged_lm_results$rnt_lm_treat_HolmAdjP_flag[to_insert] <- merged_lm_results$rnt_lm_treat_HolmAdjP_flag.y[to_insert]

merged_lm_results$super.pathway <- merged_lm_results$super.pathway.x
to_insert <- which(is.na(merged_lm_results$super.pathway))
merged_lm_results$super.pathway[to_insert] <- merged_lm_results$super.pathway.y[to_insert]

# add panel and superpathway info from class info
merged_lm_results$panel <- class_info$panel[match(merged_lm_results$feature_id,class_info$feature_id)]
merged_lm_results$super.pathway <- class_info$super.pathway[match(merged_lm_results$feature_id,class_info$feature_id)]

# rename
metabprank <- which(names(merged_lm_results) == "rnt_lm_treat_p_rank.x")
names(merged_lm_results)[metabprank] <- "rnt_lm_treat_p_rank.metab"
nmrprank <- which(names(merged_lm_results) == "rnt_lm_treat_p_rank.y")
names(merged_lm_results)[nmrprank] <- "rnt_lm_treat_p_rank.nmr"

# save out merged results
saveRDS(merged_lm_results,file=paste0(results_dir,"merged_lm_results.rds"))

####################################################################################
####################################################################################
## volcano plot - Figure S2
####################################################################################
####################################################################################

# set up colours
safe_palette <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499",
                  "#44AA99","#999933","#882255","#661100","#6699CC","#888888")
merged_lm_results$super.pathway <- as.factor(merged_lm_results$super.pathway)
metab_col <- with(merged_lm_results, data.frame(pathway = levels(super.pathway),colour = safe_palette[1:nlevels(super.pathway)] ))

# set up shapes
merged_lm_results$panel <- as.factor(merged_lm_results$panel)
metab_shape <- with(merged_lm_results, data.frame(panel = levels(panel),shape = c(18,20) ))

# # make volcano
# filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS2_lm_results_all_volcano")
# make.volcano(dtst=merged_lm_results, title="", filename=filename, plot_col=metab_col$colour[match(merged_lm_results$super.pathway,metab_col$pathway)], 
#              hline=-log10(0.05), legend_one=as.character(metab_col$pathway),legend_col=metab_col$colour,
#              plot_shape=metab_shape$shape[match(merged_lm_results$panel,metab_shape$panel)],
#              legend_two=as.character(metab_shape$panel), legend_shape=metab_shape$shape)
# 
# ## histogram (aim to illustrate density of points up/down)
# hist(-log10(merged_lm_results$rnt_lm_treat_HolmAdjP),main="",ylab="",xlab="",yaxt = "n")
# 
####################################################################################
####################################################################################
## presentation of enrichment - Figure S3
####################################################################################
####################################################################################

# add flag to class info that indicates associated metabolites
class_info$associated <- merged_lm_results$rnt_lm_treat_HolmAdjP_flag[match(class_info$feature_id,merged_lm_results$feature_id)]

# restrict to features that went through lm model primary analysis
class_info_lm_only <- class_info[!is.na(class_info$associated),]
dim(class_info_lm_only)

# separate out associated
class_info_associated <- class_info_lm_only[class_info_lm_only$associated ==1,]

# extract numbers for stack
table_for_all <- as.data.frame(table(class_info_lm_only$super.pathway))
# add %
table_for_all$class_pct <- (table_for_all$Freq / sum(table_for_all$Freq))*100

# extract numbers for stack
table_for_assoc <- as.data.frame(table(class_info_associated$super.pathway))
# add %
table_for_assoc$class_pct <- (table_for_assoc$Freq / sum(table_for_assoc$Freq))*100

# add label for all versus assoc
table_for_assoc$outcome <- paste0("Associated","\n","metabolites")
table_for_all$outcome <- paste0("All","\n","metabolites")

# join to make graph input
stacked_bar_data <- rbind(table_for_all,table_for_assoc)
# update var labels
names(stacked_bar_data)[1] <- "super.pathway"
names(stacked_bar_data)[2] <- "freq_overall"
names(stacked_bar_data)[3] <- "Percentage.in.class"

# # plot stacked bar
# filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS3_class_distribution.pdf")
# pdf(filename)
# ggplot(stacked_bar_data, aes(fill=super.pathway,y=Percentage.in.class,x=outcome)) +
#   geom_bar(position="stack",stat="identity",colour="black",size=0.1) +
#   theme_bw() + 
#   theme(axis.title.x = element_blank()) +
#   labs(y="Percentage in class") +
#   scale_fill_manual(values=c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499",
#                              "#44AA99","#999933","#882255","#661100","#6699CC","#888888"))
# invisible(dev.off())


####################################################################################
####################################################################################
## presentation of N samples/beta SE - Figure S1
####################################################################################
####################################################################################
# 
# vline_val <- 0.4*max(results_list$metab$rnt_lm_n_samples_in_model)
# 
# filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS1_N_SE_plot.pdf")
# pdf(filename)
# plot(results_list$metab$rnt_lm_n_samples_in_model,results_list$metab$rnt_lm_treat_SE,
#      main="",xlab="Number of samples in linear model", ylab="Standard error of effect estimate (beta)")
# abline(v=vline_val,lty=2)
# invisible(dev.off())

####################################################################################
####################################################################################
## histogram of sample sizes 
####################################################################################
####################################################################################

# restrict to only those metabolites for which lm is main result
# 
# filename <- paste0(results_dir,"Figures\\NMR_N_hist.pdf")
# pdf(filename)
# hist(lm_results_list$nmr$rnt_lm_n_samples_in_model,main="",xlab="No. of samples in model")
# invisible(dev.off())
# 
# filename <- paste0(results_dir,"Figures\\Metabolon_N_hist.pdf")
# pdf(filename)
# hist(lm_results_list$metab$rnt_lm_n_samples_in_model,main="",xlab="No. of samples in model")
# invisible(dev.off())

#################################################################################### 
####################################################################################
## heatmap with additional factors - Figure 2
####################################################################################
####################################################################################

# extract trial data for annotations
annot_var <- c("weight.change","reversal","treat")
annot_data <- trial_data[trial_data$id %in% selected_resid$id,c("id",annot_var)]

# merge annotation data with residuals
heatmap_dataset <- merge(annot_data,selected_resid,by="id")

# order by treatment
heatmap_dataset <- heatmap_dataset[order(as.character((heatmap_dataset$treat))),]

# set up colour categories
colfunc_blue <- colorRampPalette(c("light blue", "dark blue"))
colfunc_green <- colorRampPalette(c("light green", "dark green"))

# start with relevant trial data
colour_spec <- heatmap_dataset[,1:4]
colour_spec$weight.change.quartiles <- ntile(colour_spec$weight.change, 4)
# define colour by intervention group
colour_spec$Treatment.group <- "blue"
colour_spec[colour_spec$treat == "Intervention", c("Treatment.group")] <- "goldenrod"
# define colour by % weight change
colour_spec$weight.change.grp <- (colfunc_blue(4))[colour_spec$weight.change.quartiles]
# define colour by diabetes reversal status
colour_spec$Remission.group <- "black"
colour_spec[colour_spec$reversal == "No", c("Remission.group")] <- "white"

# specify heatmap column colours
ColSideColors <- cbind(Treatment.group=colour_spec$Treatment.group,weight.grp=colour_spec$weight.change.grp,Remission.group=colour_spec$Remission.group)

# restrict to data only
heatmap_data_only <- heatmap_dataset[,c(5:ncol(heatmap_dataset) )]
dim(heatmap_data_only)
heatmap_data_only[1:5,1:5]

# generate id based on Super.pathway (instead of comp.id)
lm_feature_ids <- as.data.frame(names(heatmap_data_only))
names(lm_feature_ids)[1] <- "feature_id"
# create ordering var
lm_feature_ids$data_order <- seq(1,nrow(lm_feature_ids),1)
# bring in Super.pathway and mass info and numeric compid
lm_feature_ids$Super.pathway <- class_info_associated$super.pathway[match(lm_feature_ids$feature_id,class_info_associated$feature_id)]

# order by Super.pathway and id
lm_feature_ids <- lm_feature_ids[order(lm_feature_ids$Super.pathway,lm_feature_ids$feature_id),]
# make short name
lm_feature_ids$shortname <- NA
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Amino Acid"), "shortname" ] <- "AA"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Carbohydrate"), "shortname" ] <- "CHO"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Cofactors and Vitamins"), "shortname" ] <- "CoFac & Vit"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Energy"), "shortname" ] <- "Energy"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Lipid"), "shortname" ] <- "Lipid"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Nucleotide"), "shortname" ] <- "Nucleotide"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Peptide"), "shortname" ] <- "Peptide"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Xenobiotics"), "shortname" ] <- "Xeno"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "NMR ratio/percentage"), "shortname" ] <- "NMR.derived"
lm_feature_ids[which(lm_feature_ids$Super.pathway == "Unclassified"), "shortname" ] <- "Unclassified"

# create unique IDs for each metabolite in pathway
lm_feature_ids$class_count <- as.numeric(ave(lm_feature_ids$feature_id, lm_feature_ids$shortname, FUN = seq_along))
lm_feature_ids$feature_label_id <- "NA"
for (i in 1:nrow(lm_feature_ids)) {
  if (lm_feature_ids[i,5] < 10) {
    lm_feature_ids[i,6] <- paste0(lm_feature_ids[i,4], "_00", lm_feature_ids[i,5])
  }
  else if (lm_feature_ids[i,5] >= 10 & lm_feature_ids[i,5] < 100 ){
    lm_feature_ids[i,6] <- paste0(lm_feature_ids[i,4], "_0", lm_feature_ids[i,5])
  }
  else {
    lm_feature_ids[i,6] <- paste0(lm_feature_ids[i,4], "_", lm_feature_ids[i,5])
  }
}

#### plot with all follow up features ####
# make plot as pdf & eps
filename <- paste0(results_dir,"Figures\\ForPaper\\Figure2_heatmap_all_followup_features")
temp_data <- make.heatmap(dtst=heatmap_data_only,feature_ids=lm_feature_ids,filename=filename,colcolours=ColSideColors[,c(1,3)])

# make plot with default legend
filename <- paste0(results_dir,"Figures\\ForPaper\\Figure2_heatmap_all_followup_features_default_legend")
temp_data <- make.heatmap.with.legend(dtst=heatmap_data_only,feature_ids=lm_feature_ids,filename=filename,colcolours=ColSideColors[,c(1,3)])

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
#rm(list=ls())

