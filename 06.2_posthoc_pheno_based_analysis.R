# This script runs posthoc analyses using trial phenos following the EFFECT OF THE INTERVENTION ON METABOLITES analysis

# last run: 9th Nov 2022

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/06.2_posthoc_pheno_based_analysis.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/06.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )
rawdata_dir <- as.character( pfile[5,2] )
ext_data_dir <- as.character( pfile[6,2] )

## call R script containing functions
source(paste0(scripts_dir,"04.0_functions.R"))
source(paste0(scripts_dir,"02.0_qc_functions.R"))
source(paste0(scripts_dir,"05.0_functions.R"))
source(paste0(scripts_dir,"06.0_functions.R"))
source(paste0(scripts_dir,"07.0_functions.R"))


####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in trial data
trial_data <- readRDS(paste0(rawdata_dir,"DiRECT_fu1_20171113.rds"))

# read in feature metadata for metab features
metab_feature_metadata <- readRDS(paste0(processeddata_dir,"OrigScale_feature_metadata_with_QCstats.rds"))
nmr_feature_metadata <- readRDS(paste0(processeddata_dir,"nmr_feature_metadata_with_QCstats.rds"))

# read in residuals from linear model for selected associated features
selected_resid <- readRDS(file = paste0(processeddata_dir,"resid_selected_assoc.rds")) 

# read in results files
lm_results_list <- readRDS(file=paste0(results_dir,"lm_results_list.rds"))
metab_pa_results <- readRDS(file=paste0(results_dir,"metab_pa_results.rds"))
clinical_lm_results <- readRDS(file=paste0(results_dir,"clinical_lm_results.rds"))

# read in list of selected associated features 
selected_assoc_features <- read.table(file=paste0(results_dir,"selected_assoc_features.txt"),h=F, stringsAsFactors = F)

# read in labels 
class_info <- readRDS(file=paste0(results_dir,"class_info_for_enrichment.rds"))

# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))

####################################################################################
####################################################################################
## summarise trial phenos
####################################################################################
####################################################################################

# work out pheno class
pheno_list <- as.data.frame(names(trial_data))
names(pheno_list)[1] <- "pheno_name"
pheno_list$pheno_name <- as.character(pheno_list$pheno_name)
pheno_list$class1 <- "unknown"
pheno_list$class2 <- NA
for (i in 1:ncol(trial_data)) {
  if (length(class(trial_data[,i])) == 1) {
    pheno_list[i,2] <- class(trial_data[,i])
  }
  else {
    pheno_list[i,2] <- class(trial_data[,i])[1]
    pheno_list[i,3] <- class(trial_data[,i])[2]
  }
}
table(pheno_list$class1)
table(pheno_list$class2)

date_vars <- pheno_list[pheno_list$class1 == "dates",c("pheno_name")]
print(paste0(length(date_vars), " Date vars:"))
date_vars
write.table(date_vars,file=paste0(processeddata_dir,"trial_date_vars.txt"),quote=F,row.names=F,col.names=F)

factor_vars <- pheno_list[pheno_list$class1 == "factor",c("pheno_name")]
print(paste0(length(factor_vars), " Factor vars:"))
factor_vars
write.table(factor_vars,file=paste0(processeddata_dir,"trial_factor_vars.txt"),quote=F,row.names=F,col.names=F)

integer_vars <- pheno_list[pheno_list$class1 == "integer",c("pheno_name")]
print(paste0(length(integer_vars), " Integer vars:"))
integer_vars
write.table(integer_vars,file=paste0(processeddata_dir,"trial_integer_vars.txt"),quote=F,row.names=F,col.names=F)

numeric_vars <- pheno_list[pheno_list$class1 == "numeric",c("pheno_name")]
print(paste0(length(numeric_vars), " Numeric vars:"))
numeric_vars
write.table(numeric_vars,file=paste0(processeddata_dir,"trial_numeric_vars.txt"),quote=F,row.names=F,col.names=F)

# join numeric and integer vars
number_vars <- c(integer_vars,numeric_vars)


## summarise according to var class
# generate tables for factors
for (i in 1:length(factor_vars)) {
  col_n <- which(names(trial_data) == factor_vars[i])
  print(factor_vars[i])
  print(table(trial_data[,col_n]))
}

# generate summaries for integers and numerics
filename <- paste0(results_dir,"Figures\\trial_pheno_summaries.pdf")
pdf(filename)
for (i in 1:length(number_vars)) {
  col_n <- which(names(trial_data) == number_vars[i])
  print(number_vars[i])
  print(summary(trial_data[,col_n]))
  hist(trial_data[,col_n],main=number_vars[i])
}
invisible(dev.off())


####################################################################################
####################################################################################
## build PCA from linear model residuals 
####################################################################################
####################################################################################

# order by id
selected_resid <- selected_resid[order(selected_resid$id),]

# make id list in order
ordered_id <- as.character(selected_resid$id)

# PCA
# pca model
pca_model <- ppca(as.matrix(selected_resid[,-c(1)]),seed=12345,nPcs=(ncol(selected_resid[,-c(1)])-1))

# scree
varexp= round( (summary(pca_model)[1,])*100, d=2)

## calculate no. of PCs to investigate based on acceleration factor, etc
pca_data <- selected_resid[,-c(1)]
ev <- eigen(stats::cor(pca_data, use="pairwise.complete.obs"))
ap <- nFactors::parallel(subject=nrow(pca_data),var=ncol(pca_data),rep=100,cent=0.05)
ns <- nFactors::nScree(x=ev$values, aparallel=ap$eigen$qevpea)
accelerationfactor = as.numeric( ns[[1]][2] )
if(accelerationfactor < 2) { accelerationfactor = 2 }
nsig_parallel = as.numeric(ns[[1]][3])
print(paste0("Number of PCs to retain is: ", nsig_parallel))
# plot scree
filename <- paste0(results_dir,"Figures\\pca_scree")
pdf(paste0(filename,".pdf"))
plotnScree(ns)
invisible(dev.off())

####################################################################################
####################################################################################
## restrict pheno to a subset of interest
####################################################################################
####################################################################################

# make lists of vars to include in PCA plots
annot_var <- c("weight.change","reversal","treat")
checklist_change_var <- c("sbp.change","dbp.change","chol.mmol.l.change","creat.umol.l.change","crp.mg.l.change",
                          "eq5d.change","eq5d.vas.change","weight.change","lab.hba1c.mmol.mol.change")
nafld_change_var <- c("ast.u.l.change","alt.u.l.change","acr.mg.mmol.change","ggt.u.l.change","liver.fat.change")

# extract data
annot_data <- trial_data[trial_data$id %in% selected_resid$id,c("id",annot_var)]
checklist_change_data <- trial_data[trial_data$id %in% selected_resid$id,c("id",checklist_change_var)]
nafld_change_data <- trial_data[trial_data$id %in% selected_resid$id,c("id",nafld_change_var)]

# check correlations in checklist change var
cor_matrix_checklist <- cor(checklist_change_data[,c(2:ncol(checklist_change_data))],method="spearman", use="pairwise.complete.obs")

# flatten cor matrix
ut <- upper.tri(cor_matrix_checklist)
cor_output_checklist <- data.frame(
  row = rownames(cor_matrix_checklist)[row(cor_matrix_checklist)[ut]],
  column = rownames(cor_matrix_checklist)[col(cor_matrix_checklist)[ut]],
  cor = (cor_matrix_checklist)[ut]
)
hist(abs(cor_output_checklist$cor),main="Histogram of correlations across checklist change vars")
print("Summary of correlations across checklist change vars")
summary(abs(cor_output_checklist$cor))


# check correlations in nafld change var
cor_matrix_nafld <- cor(nafld_change_data[,c(2:ncol(nafld_change_data))],method="spearman", use="pairwise.complete.obs")

# flatten cor matrix
ut <- upper.tri(cor_matrix_nafld)
cor_output_nafld <- data.frame(
  row = rownames(cor_matrix_nafld)[row(cor_matrix_nafld)[ut]],
  column = rownames(cor_matrix_nafld)[col(cor_matrix_nafld)[ut]],
  cor = (cor_matrix_nafld)[ut]
)
hist(abs(cor_output_nafld$cor),main="Histogram of correlations across nafld change vars")
print("Summary of correlations across nafld change vars")
summary(abs(cor_output_nafld$cor))


####################################################################################
####################################################################################
## explore relationship between PCs and selected trial vars for Figure S8 & 9
####################################################################################
####################################################################################

#NB This function currently only works for npcs=2 because of selection step that reduces pheno to those that are meaningfully correlated with the PCs
corbin_biplot_for_ppca(PCA=pca_model, dataframe_of_phenotypes=checklist_change_data, pheno_list=checklist_change_var, class_info = class_info,
                       npcs = 2, annot_data = annot_data, ordered_id = ordered_id, filename="FigureS8A_checklist_change")

corbin_biplot_for_ppca(PCA=pca_model, dataframe_of_phenotypes=nafld_change_data, pheno_list=nafld_change_var, class_info = class_info, 
                       npcs = 2, annot_data = annot_data, ordered_id = ordered_id, filename="FigureS8B_nafld_change")

# for interactive testing of biplot function only
#PCA=pca_model; dataframe_of_phenotypes=nafld_change_data; pheno_list=nafld_change_var; class_info = class_info
#npcs = 2; annot_data = annot_data; ordered_id = ordered_id; filename="Figure4_nafld_change"
#pc_x = 1; pc_y = 2
####################################################################################
####################################################################################
## plot eigenvalues for Figure S6
####################################################################################
####################################################################################

filename <- paste0(results_dir,"Figures\\ForPaper\\FigureS6_scree.pdf")
pdf(filename)
plot(ev$values, main="", xlab = "Principal component", ylab = "Eigenvalue", pch = 4, col="black", frame.plot=F)
invisible(dev.off())

####################################################################################
####################################################################################
## explore relationship between weight loss and GP practice (site)
####################################################################################
####################################################################################

# extract data
var_data <- trial_data[trial_data$id %in% selected_resid$id,c("id","site","treat","weight.change")]

# identify which sites are control and which intervention
site_def <- as.data.frame(table(var_data$site,var_data$treat))
# limit to control counts only
site_def <- site_def[site_def$Var2 == "Control",]
site_def$site_designation <- NA
site_def$Var2 <- as.character(site_def$Var2)
for (i in 1:nrow(site_def)) {
  if (site_def[i,3] != 0) {
    site_def[i,4] <- site_def[i,2]
  }
  else if (site_def[i,3] == 0 & site_def[i,2] == "Control") {
    site_def[i,4] <- "Intervention"
  }
  else {
    site_def[i,4] <- "Control"
  }
}

control_sites <- as.character(site_def[site_def$site_designation == "Control",c("Var1")])
intervention_sites <- as.character(site_def[site_def$site_designation == "Intervention",c("Var1")])

# look to see if weight change varies by site (within intervention group)
var_data_controls <- var_data[var_data$site %in% as.numeric(control_sites),]
table(var_data_controls$treat)
table(var_data_controls$site)
fitA <- lm(weight.change ~ as.factor(site), data=var_data_controls)
summary(fitA)

var_data_intervention <- var_data[var_data$site %in% as.numeric(intervention_sites),]
table(var_data_intervention$treat)
table(var_data_intervention$site)
fitB <- lm(weight.change ~ as.factor(site), data=var_data_intervention)
summary(fitB)

## look at relationship between NEFA and associated metabolites
trial_data_ids_to_keep <- trial_data[which(trial_data$id %in% selected_resid$id),]
nefa_pheno <- trial_data_ids_to_keep[,c("id","nefa.10.b","nefa.70.b","nefa.10","nefa.70","nefa.10.change","nefa.70.change")]

# check correlations 
cor_matrix_nefa <- cor(nefa_pheno[,c(2:ncol(nefa_pheno))],selected_resid[,c(2:ncol(selected_resid))], method="spearman", use="pairwise.complete.obs")

# transpose
cor_matrix_nefa <- as.data.frame(t(cor_matrix_nefa))

# add metabolite names
cor_matrix_nefa$feature_id <- rownames(cor_matrix_nefa)
cor_matrix_nefa$biochemical <- class_info$biochemical[match(cor_matrix_nefa$feature_id,class_info$feature_id)]

# write out
write.table(cor_matrix_nefa,file=paste0(results_dir,"Tables\\Table_of_nefa_correlations.txt"),row.names=F,sep="\t",quote=F)


######################################################################################################
######################################################################################################
## explore relationship between presence of metformin in metabolite data and clinical reporting
######################################################################################################
######################################################################################################

# extract clinical metformin data
var_data <- trial_data[trial_data$id %in% selected_resid$id,c("id","metformin","metformin.b")]

# extract metabolite metformin data (P/A)
# find compid for metformin
metformin_compid <- metab_feature_metadata[which(metab_feature_metadata$biochemical == "metformin"),c("compid_to_match")]
met_b <- paste0(metformin_compid,".b")
met_e <- paste0(metformin_compid,".e")

metformin_metab_data <- dataset_list$pa$metab[,c(met_b,met_e)]
metformin_metab_data$id <- row.names(metformin_metab_data)

var_data_plus_metab <- merge(var_data,metformin_metab_data,by="id",all.x=T)

## compare
print("Comparison of clinical metformin use (yes/no) and presence (1) or absence (0) of metformin in the sample, baseline")
table(var_data_plus_metab$metformin.b,var_data_plus_metab$compid38306.b)

print("Comparison of clinical metformin use (yes/no) and presence (1) or absence (0) of metformin in the sample, 12 months")
table(var_data_plus_metab$metformin,var_data_plus_metab$compid38306.e)

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
