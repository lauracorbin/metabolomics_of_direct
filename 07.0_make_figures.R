# This script makes figures from analyses relating to the EFFECT OF THE INTERVENTION ON METABOLITES

# last run: 1st April 2022

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
## volcano plot - Figure 2
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

# make volcano
filename <- paste0(results_dir,"Figures\\ForPaper\\Figure2_lm_results_all_volcano")
make.volcano(dtst=merged_lm_results, title="", filename=filename, plot_col=metab_col$colour[match(merged_lm_results$super.pathway,metab_col$pathway)], 
             hline=-log10(0.05), legend_one=as.character(metab_col$pathway),legend_col=metab_col$colour,
             plot_shape=metab_shape$shape[match(merged_lm_results$panel,metab_shape$panel)],
             legend_two=as.character(metab_shape$panel), legend_shape=metab_shape$shape)

filename <- paste0(results_dir,"Figures\\Figure2_lm_results_all_volcano_dh")
make.volcano.dh(dtst=merged_lm_results, title="", filename=filename, plot_col=metab_col$colour[match(merged_lm_results$super.pathway,metab_col$pathway)], 
             hline=-log10(0.05), legend_one=as.character(metab_col$pathway),legend_col=metab_col$colour,
             plot_shape=metab_shape$shape[match(merged_lm_results$panel,metab_shape$panel)],
             legend_two=as.character(metab_shape$panel), legend_shape=metab_shape$shape)

## histogram (aim to illustrate density of points up/down)
hist(-log10(merged_lm_results$rnt_lm_treat_HolmAdjP),main="",ylab="",xlab="",yaxt = "n")

####################################################################################
####################################################################################
## presentation of enrichment - Figure 3
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

# plot stacked bar
filename <- paste0(results_dir,"Figures\\ForPaper\\Figure3_class_distribution.pdf")
pdf(filename)
ggplot(stacked_bar_data, aes(fill=super.pathway,y=Percentage.in.class,x=outcome)) +
  geom_bar(position="stack",stat="identity",colour="black",size=0.1) +
  theme_bw() + 
  theme(axis.title.x = element_blank()) +
  labs(y="Percentage in class") +
  scale_fill_manual(values=c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499",
                             "#44AA99","#999933","#882255","#661100","#6699CC","#888888"))
invisible(dev.off())


####################################################################################
####################################################################################
## presentation of N samples/beta SE - Figure S1
####################################################################################
####################################################################################

vline_val <- 0.4*max(results_list$metab$rnt_lm_n_samples_in_model)

filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS1_N_SE_plot.pdf")
pdf(filename)
plot(results_list$metab$rnt_lm_n_samples_in_model,results_list$metab$rnt_lm_treat_SE,
     main="",xlab="Number of samples in linear model", ylab="Standard error of effect estimate (beta)")
abline(v=vline_val,lty=2)
invisible(dev.off())

####################################################################################
####################################################################################
## histogram of sample sizes 
####################################################################################
####################################################################################

# restrict to only those metabolites for which lm is main result

filename <- paste0(results_dir,"Figures\\NMR_N_hist.pdf")
pdf(filename)
hist(lm_results_list$nmr$rnt_lm_n_samples_in_model,main="",xlab="No. of samples in model")
invisible(dev.off())

filename <- paste0(results_dir,"Figures\\Metabolon_N_hist.pdf")
pdf(filename)
hist(lm_results_list$metab$rnt_lm_n_samples_in_model,main="",xlab="No. of samples in model")
invisible(dev.off())

#################################################################################### 
####################################################################################
## heatmap with additional factors - Figure S3
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
colour_spec$treat.grp <- "blue"
colour_spec[colour_spec$treat == "Intervention", c("treat.grp")] <- "goldenrod"
# define colour by % weight change
colour_spec$weight.change.grp <- (colfunc_blue(4))[colour_spec$weight.change.quartiles]
# define colour by diabetes reversal status
colour_spec$reversal.grp <- "black"
colour_spec[colour_spec$reversal == "No", c("reversal.grp")] <- "white"

# specify heatmap column colours
ColSideColors <- cbind(treat.grp=colour_spec$treat.grp,weight.grp=colour_spec$weight.change.grp,reversal.grp=colour_spec$reversal.grp)

# restrict to data only
heatmap_data_only <- heatmap_dataset[,c(5:ncol(heatmap_dataset) )]
dim(heatmap_data_only)
heatmap_data_only[1:5,1:5]

# generate id based on super.pathway (instead of comp.id)
lm_feature_ids <- as.data.frame(names(heatmap_data_only))
names(lm_feature_ids)[1] <- "feature_id"
# create ordering var
lm_feature_ids$data_order <- seq(1,nrow(lm_feature_ids),1)
# bring in super.pathway and mass info and numeric compid
lm_feature_ids$super.pathway <- class_info_associated$super.pathway[match(lm_feature_ids$feature_id,class_info_associated$feature_id)]

# order by super.pathway and id
lm_feature_ids <- lm_feature_ids[order(lm_feature_ids$super.pathway,lm_feature_ids$feature_id),]
# make short name
lm_feature_ids$shortname <- NA
lm_feature_ids[which(lm_feature_ids$super.pathway == "Amino Acid"), "shortname" ] <- "AA"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Carbohydrate"), "shortname" ] <- "CHO"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Cofactors and Vitamins"), "shortname" ] <- "CoFac & Vit"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Energy"), "shortname" ] <- "Energy"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Lipid"), "shortname" ] <- "Lipid"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Nucleotide"), "shortname" ] <- "Nucleotide"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Peptide"), "shortname" ] <- "Peptide"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Xenobiotics"), "shortname" ] <- "Xeno"
lm_feature_ids[which(lm_feature_ids$super.pathway == "NMR ratio/percentage"), "shortname" ] <- "NMR.derived"
lm_feature_ids[which(lm_feature_ids$super.pathway == "Unclassified"), "shortname" ] <- "Unclassified"

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
# make plot
filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS4_heatmap_all_followup_features")
temp_data <- make.heatmap(dtst=heatmap_data_only,feature_ids=lm_feature_ids,filename=filename,colcolours=ColSideColors[,c(1,3)])


####################################################################################
####################################################################################
### make plots for graphical abstract
####################################################################################
####################################################################################

# need to calculate fold change at baseline

# set up list for results output
fc_interim_results_list <- list(metab=as.data.frame(matrix(data = NA, nrow = 1300, ncol = 18)),
                                nmr=as.data.frame(matrix(data = NA, nrow = 300, ncol = 18)))

# implement median fold change on raw data
#datasets_to_run <- c(1,5)
#for (j in 1:length(datasets_to_run)){
  # set dataset
dataset_to_use <- dataset_list$raw
for (k in 1:length(dataset_to_use)){
  print(paste0("Running fold change for: ",names(dataset_to_use)[k]))
  # set warnings to print after each model run
  op <- options("warn")
  on.exit(options(op))
  options(warn=1)
  # define dataset - merge trial data with metabolite data
  dtst <- merge(clinic_data_only[[k]],dataset_to_use[[k]],by="row.names")
  rownames(dtst) <- dtst$Row.names
  dtst <- dtst[,-1]
  # define var list (ie list of metabolites)
  depvar <- c(clean_feature_list[[k]]$id_for_matching,"glucose_mmol_l","insulin_uu_ml","hdl_mmol_l","trig_mmol_l","chol_mmol_l")
  ## set up reporting var frequency
  s = base::seq(from=1,to=length(depvar), by=100)
  for (i in 1:length(depvar)) {
    if(i %in% s){
      print(paste0( "Now processing metabolite ", i, " of ", length(depvar) ))
    }
    mtb_name <- depvar[i]
    mtb_name_e <- paste0(mtb_name,".e")
    mtb_name_b <- paste0(mtb_name,".b")
    fc_interim_results_list[[k]][i,1] <- mtb_name
    mtb_b_col_ref <- which( colnames(dtst) == mtb_name_b)
    # subset data to single feature
    mtb <- dtst[!is.na(dtst[mtb_b_col_ref]), c("id","treat",mtb_name_b)]
    # count non-missing 
    fc_interim_results_list[[k]][i,2] <- nrow(mtb)
    # count per group
    fc_interim_results_list[[k]][i,4] <- nrow(mtb[mtb$treat == "Intervention",])
    fc_interim_results_list[[k]][i,3] <- nrow(mtb[mtb$treat == "Control",])
    # calculate means by group
    mean_control <- mean(mtb[mtb$treat == "Control",3])
    var_control <- var(mtb[mtb$treat == "Control",3])
    med_control <- median(mtb[mtb$treat == "Control",3])
    geom_control <- exp(mean(log(mtb[mtb$treat == "Control",3])))
    fc_interim_results_list[[k]][i,5] <- mean_control
    fc_interim_results_list[[k]][i,6] <- var_control
    fc_interim_results_list[[k]][i,7] <- med_control
    fc_interim_results_list[[k]][i,8] <- geom_control
    mean_intervention <- mean(mtb[mtb$treat == "Intervention",3])
    var_intervention <- var(mtb[mtb$treat == "Intervention",3])
    med_intervention <- median(mtb[mtb$treat == "Intervention",3])
    geom_intervention <- exp(mean(log(mtb[mtb$treat == "Intervention",3])))
    fc_interim_results_list[[k]][i,9] <- mean_intervention
    fc_interim_results_list[[k]][i,10] <- var_intervention
    fc_interim_results_list[[k]][i,11] <- med_intervention
    fc_interim_results_list[[k]][i,12] <- geom_intervention
    # proceed to FC calculation only if metabolite is in the subset of selected metabolites from step 1
    fc <- mean_intervention/mean_control
    fc_interim_results_list[[k]][i,13] <- fc
    fc_interim_results_list[[k]][i,14] <- log2(fc)
    fc_med <- med_intervention/med_control
    fc_interim_results_list[[k]][i,15] <- fc_med
    fc_interim_results_list[[k]][i,16] <- log2(fc_med)
    fc_geom <- geom_intervention/geom_control
    fc_interim_results_list[[k]][i,17] <- fc_geom
    fc_interim_results_list[[k]][i,18] <- log2(fc_geom)
  }
  # name column headers
  names(fc_interim_results_list[[k]])[1] <- "feature_id"
  names(fc_interim_results_list[[k]])[2] <- "fc_n_samples"
  names(fc_interim_results_list[[k]])[3] <- "fc_n_control"
  names(fc_interim_results_list[[k]])[4] <- "fc_n_intervention"
  
  names(fc_interim_results_list[[k]])[5] <- "fc_mean_in_controls"
  names(fc_interim_results_list[[k]])[6] <- "fc_var_in_controls"
  names(fc_interim_results_list[[k]])[7] <- "fc_median_in_controls"
  names(fc_interim_results_list[[k]])[8] <- "fc_geom_in_controls"
  
  names(fc_interim_results_list[[k]])[9] <- "fc_mean_in_intervention"
  names(fc_interim_results_list[[k]])[10] <- "fc_var_in_intervention"
  names(fc_interim_results_list[[k]])[11] <- "fc_median_in_intervention"
  names(fc_interim_results_list[[k]])[12] <- "fc_geom_in_intervention"
  
  names(fc_interim_results_list[[k]])[13] <- "fc_meanfoldchange_int_over_control"
  names(fc_interim_results_list[[k]])[14] <- "fc_log2meanfoldchange_int_over_control"
  names(fc_interim_results_list[[k]])[15] <- "fc_medianfoldchange_int_over_control"
  names(fc_interim_results_list[[k]])[16] <- "fc_log2medianfoldchange_int_over_control"
  names(fc_interim_results_list[[k]])[17] <- "fc_geomfoldchange_int_over_control"
  names(fc_interim_results_list[[k]])[18] <- "fc_log2geomfoldchange_int_over_control"   
  
  names(fc_interim_results_list[[k]])[2:18] <- paste0("raw_", names(fc_interim_results_list[[k]])[2:18])
}

## put all baseline foldchange values in a matrix
# combine nmr and metab
fc_results_combined <- rbind(fc_interim_results_list$metab,fc_interim_results_list$nmr)
# combine with endpoint
fc_both <- merged_lm_results[,c("feature_id","raw_fc_log2medianfoldchange_int_over_control","rnt_lm_treat_HolmAdjP_flag")]
fc_both$raw_fc_log2medianfoldchange_int_over_control_baseline <- fc_results_combined$raw_fc_log2medianfoldchange_int_over_control[match(fc_both$feature_id,fc_results_combined$feature_id)]

# restrict to associated
fc_both <- fc_both[fc_both$rnt_lm_treat_HolmAdjP_flag ==1,]

# add blank row
#fc_both <- rbind(fc_both,data.frame(feature_id ="blank",raw_fc_log2medianfoldchange_int_over_control=0,raw_fc_log2medianfoldchange_int_over_control_baseline=0,rnt_lm_treat_HolmAdjP_flag=0))
# make into matrices
baseline_matrix <- matrix(fc_both$raw_fc_log2medianfoldchange_int_over_control_baseline,nrow=6,ncol=31)                                                                                                                                    
endpoint_matrix <- matrix(fc_both$raw_fc_log2medianfoldchange_int_over_control,nrow=6,ncol=31)                                                                                                                                    
overall <- rbind(endpoint_matrix,baseline_matrix)

# heatmaps
#paletteLength=50
#myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
#myBreaks <- c(seq(-0.90,0,length.out=ceiling(paletteLength/2) + 1),
            #  seq(1.4/paletteLength,1.8,length.out=floor(paletteLength/2)))

#heatmap3(baseline_matrix,cor=color.palette,breaks=palette.breaks)

#heatmap3(baseline_matrix,balanceColor = T, Rowv = NA, Colv = NA, scale="none",labRow = F,labCol = F)
#heatmap3(endpoint_matrix,balanceColor = T, Rowv = NA, Colv = NA, scale="none",labRow = F,labCol = F)

filename <- paste0(results_dir,"Figures\\ForPaper\\GA_heatmap.pdf")
pdf(filename)
heatmap3(overall,balanceColor = T, Rowv = NA, Colv = NA, scale="none",labRow = F,labCol = F)
dev.off()

#hist(endpoint_matrix)
#hist(baseline_matrix)

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
