# This script pulls together the trial phenotype data needed for the analysis
# And combines it with the metabolite data
# last run: 31st March 2022

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/03.0_pheno_data_prep.txt",split=TRUE)

# read in parameter file - contains file paths and file names
pfile = read.table("U:/Scripts/03.0_pheno_data_prep_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
rawdata_dir <- as.character( pfile[4,2] )
id_correction <- as.character( pfile[5,2] )
id_correction <- strsplit(id_correction,split=",")[[1]] 
exclusions <- as.character( pfile[6,2] )
exclusions <- strsplit(exclusions,split=",")[[1]]

####################################################################################
####################################################################################
# read in trial dataset
####################################################################################
####################################################################################

all_data <- readRDS(paste0(rawdata_dir,"DiRECT_fu1_20171113.rds"))
#names(all_data)
print(paste("Number of individuals in trial dataset:", nrow(all_data)))

####################################################################################
####################################################################################
# restrict trial data to required variables - separate pheno from clinical bloods
####################################################################################
####################################################################################

# retain variables for model
model_data <- all_data[,c("id","site","centre","treat","list.size","age","sex","bmi.b","bmi","bmi.change","weight.b","weight","weight.change")]
#dim(model_data)
#head(model_data)
#tail(model_data)

# rename weight and bmi endpoint var so .e added to signify endpoint
names(model_data)[9] <- "bmi.e"
names(model_data)[12] <- "weight.e"

# check class of factors for model
class(model_data$centre)
class(model_data$treat)
class(model_data$list.size)
class(model_data$site)
class(model_data$age)
class(model_data$sex)

# reformat as required
model_data$site <- as.factor(model_data$site)


# save out model data
saveRDS(model_data,file=paste0(processeddata_dir,"model_data.rds"))

# retain clinical bloods variables for analysis
bloods_data <- all_data[,c(1,126,127,129,130,139,140,141,152,153,154)]

print("Blood variables retained for analysis alongside metabolites:")
names(bloods_data)

# rename var - remove '.'
names(bloods_data)[2] <- "glucose_mmol_l.b"
names(bloods_data)[3] <- "glucose_mmol_l.e"
names(bloods_data)[4] <- "insulin_uu_ml.b"
names(bloods_data)[5] <- "insulin_uu_ml.e"
names(bloods_data)[6] <- "hdl_mmol_l.b"
names(bloods_data)[9] <- "hdl_mmol_l.e"
names(bloods_data)[7] <- "trig_mmol_l.b"
names(bloods_data)[10] <- "trig_mmol_l.e"
names(bloods_data)[8] <-"chol_mmol_l.b"
names(bloods_data)[11] <- "chol_mmol_l.e"

# save out variable data
saveRDS(bloods_data,file=paste0(processeddata_dir,"bloods_data.rds"))

# add flag for presence of individual in model_data
model_data$in_trial_data <- "yes"

####################################################################################
####################################################################################
## read in metabolite data files
####################################################################################
####################################################################################

## qc'd metabolon data
metab_data <- readRDS(paste0(processeddata_dir,"OrigScale_cleaned_data.rds"))
metab_feature <- readRDS(paste0(processeddata_dir,"OrigScale_cleaned_feature_metadata.rds"))
metab_sample <- readRDS(paste0(processeddata_dir,"OrigScale_sample_metadata_with_QCstats.rds"))

## replace client.identifier that seems to be incorrect
metab_sample[which(metab_sample$sample.name == id_correction[1] ),c("client.identifier")] <- id_correction[2]

## replace sample id with client identifier (cid) in metab_data
metab_sample_list <- as.data.frame(names(metab_data))
names(metab_sample_list)[1] <- "sample.name"
metab_sample_list$cid <- metab_sample$client.identifier[match(metab_sample_list$sample.name,metab_sample$sample.name)]
metab_sample_list$cid <- paste0("cid.",metab_sample_list$cid)
names(metab_data) <- metab_sample_list$cid

## qc'd nmr data
nmr_data <- readRDS(paste0(processeddata_dir,"nmr_cleaned_data.rds"))
nmr_feature <- readRDS(paste0(processeddata_dir,"nmr_cleaned_feature_metadata.rds"))
nmr_sample <- readRDS(paste0(processeddata_dir,"nmr_sample_metadata_with_QCstats.rds"))

## replace sample id with client identifier (cid) in nmr_data
nmr_sample_list <- as.data.frame(names(nmr_data))
names(nmr_sample_list)[1] <- "sampleid"
nmr_sample$client.identifier <- row.names(nmr_sample)
nmr_sample_list$cid <- nmr_sample$client.identifier[match(nmr_sample_list$sampleid,nmr_sample$sampleid)]
nmr_sample_list$cid <- paste0("cid.",nmr_sample_list$cid)
names(nmr_data) <- nmr_sample_list$cid

## add prefixes to cid
names(metab_sample)[3] <- "cid"
metab_sample$cid <- paste0("cid.",metab_sample$cid)
names(nmr_sample)[15] <- "cid"
nmr_sample$cid <- paste0("cid.",nmr_sample$cid)

## add prefixes to sample metadata files
names(metab_sample) <- paste0("metabolon_",colnames(metab_sample))
names(nmr_sample) <- paste0("nmr_",colnames(nmr_sample))

####################################################################################
####################################################################################
# add trial data for model to sample metadata
####################################################################################
####################################################################################

# merge
orig_sample <- merge(nmr_sample,metab_sample,by.x="nmr_cid",by.y="metabolon_cid",all=T)

# check qc exclusions
print("Check QC exclusions across two datasets:")
table(orig_sample$metabolon_qc_flag,orig_sample$nmr_qc_flag,useNA="always")

# add trial data to metabolite sample metadata file
orig_sample$treat <- model_data$treat[match(orig_sample$metabolon_id,model_data$id)]

# rename client id column
names(orig_sample)[1] <- "cid"

# rename trial id column
names(orig_sample)[20] <- "trial_id"

# rename visit column
names(orig_sample)[21] <- "visit"

# check sample numbers at time points and by treatment group
print("Check visit numbers by treatment status for n=574 samples:")
table(orig_sample$visit,orig_sample$treat,useNA="always")
# restrict to metab qc pass
metab_pass <- orig_sample[orig_sample$metabolon_qc_flag == "pass",]
print("Check visit numbers by treatment status after metabolon sample qc (duplicates still in):")
table(metab_pass$visit,metab_pass$treat,useNA="always")
# restrict to nmr qc pass
nmr_pass <- orig_sample[orig_sample$nmr_qc_flag == "pass",]
print("Check visit numbers by treatment status after nmr sample qc (duplicates still in):")
table(nmr_pass$visit,nmr_pass$treat,useNA="always")

# replace missing in 'visit' with alternative value to enable selection
orig_sample[which(is.na(orig_sample$visit)),c("visit")] <- "unknown"
print("Check visit numbers after NA replaced with unknown:")
table(orig_sample$visit, useNA = "always")

####################################################################################
####################################################################################
## make data files for running models on
####################################################################################
####################################################################################

## edit metabolite data so baseline/end added to metab name
# ID sample sets based on client identifier
baseline_samples <- orig_sample[orig_sample$visit == "V1",c("cid")]
end_samples <- orig_sample[orig_sample$visit == "V2",c("cid")]

## add flag for presence of baseline/endpoint samples for individuals to model data
ids_with_baseline <- orig_sample$trial_id[which(orig_sample$cid %in% baseline_samples)]
model_data$has_baseline_sample <- "no"
model_data[which(model_data$id %in% ids_with_baseline),c("has_baseline_sample")] <- "yes"
print("Sample numbers at baseline (before QC) when combined with trial data:")
table(model_data$treat,model_data$has_baseline_sample,useNA="always")

ids_with_end <- orig_sample$trial_id[which(orig_sample$cid %in% end_samples)]
model_data$has_end_sample <- "no"
model_data[which(model_data$id %in% ids_with_end),c("has_end_sample")] <- "yes"
print("Sample numbers at year 1 (before QC) when combined with trial data (includes two duplicates):")
table(model_data$treat,model_data$has_end_sample,useNA="always")


# extract data and transform
for (i in 1:2) {
  # set dataset
  if (i==1) {
    dataset <- metab_data
    print("Processing Metabolon data")}
  else {
    dataset <- nmr_data
    print("Processing NMR data")}
  
  # extract list of baseline samples
  baseline_data <- as.data.frame(t(dataset[,names(dataset) %in% baseline_samples]))
 
  # extract list of year 1 samples
  end_data <- as.data.frame(t(dataset[,names(dataset) %in% end_samples]))
  
  # add suffix
  names(baseline_data)[1:ncol(baseline_data)] <- paste0(names(baseline_data),".b")
  names(end_data)[1:ncol(end_data)] <- paste0(names(end_data),".e")
  
  # move sample ID to column in dataframe
  baseline_data$cid <- row.names(baseline_data)
  end_data$cid <- row.names(end_data)
  
  # bring in trial id
  baseline_data$trial_id <- orig_sample$trial_id[match(baseline_data$cid,orig_sample$cid)]
  end_data$trial_id <- orig_sample$trial_id[match(end_data$cid,orig_sample$cid)]

  # drop duplicate samples in end_data (occurs twice)
  # drop sample with higher nmr missingness
  exc1 <- which(end_data$cid == exclusions[1])
  exc2 <- which(end_data$cid ==  exclusions[2])
  end_data <- end_data[-c(exc1,exc2),]
  
  # drop cid
  baseline_data <- baseline_data[,-(ncol(baseline_data)-1)]
  end_data <- end_data[,-(ncol(end_data)-1)]

  # check numbers
  print("Baseline dataset dimensions:")
  print(dim(baseline_data))
  print("Year one dataset dimensions:")
  print(dim(end_data))
  
  # merge trial data and metabolite data
  phenofile_plus_b <- merge(model_data,baseline_data,by.x="id",by.y="trial_id",all = T)
  dim(phenofile_plus_b)
  
  phenofile_plus_e <- merge(phenofile_plus_b,end_data,by.x="id",by.y="trial_id",all = T)
  dim(phenofile_plus_e)
  
  # merge bloods data
  final_phenofile <- merge(phenofile_plus_e,bloods_data,by="id",all=T)
  
  # check numbers
  table(final_phenofile$treat,useNA = "always")
  
  # save out data
  if (i==1) {
    saveRDS(final_phenofile,file=paste0(processeddata_dir,"metab_combined_data.rds"))
    # write out ids that should be analysed
    write.table(baseline_data$trial_id,file=paste0(processeddata_dir,"metab_baseline_samples.txt"),row.names=F,col.names=F,quote=F)
    write.table(end_data$trial_id,file=paste0(processeddata_dir,"metab_end_samples.txt"),row.names=F,col.names=F,quote=F)
    }
  else {
    saveRDS(final_phenofile,file=paste0(processeddata_dir,"nmr_combined_data.rds"))
    # write out ids that should be analysed
    write.table(baseline_data$trial_id,file=paste0(processeddata_dir,"nmr_baseline_samples.txt"),row.names=F,col.names=F,quote=F)
    write.table(end_data$trial_id,file=paste0(processeddata_dir,"nmr_end_samples.txt"),row.names=F,col.names=F,quote=F)
    }
}


##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()


