# This script takes the raw metabolite data and processes it ahead of the analysis
# Removes individuals with missing group allocation and qc fails

# last run: 31st March 2022

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/04.0_metabolite_data_prep.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/04.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )

## call R scripts containing functions
source(paste0(scripts_dir,"02.0_qc_functions.R"))
source(paste0(scripts_dir,"04.0_functions.R"))

####################################################################################
####################################################################################
## read in and process data
####################################################################################
####################################################################################

# merged trial pheno and metabolite data
model_data_list <- list(metab=readRDS(file=paste0(processeddata_dir,"metab_combined_data.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_combined_data.rds")))

# read in lists of samples that passed QC and can be analysed
nmr_end_samples <- read.table(file=paste0(processeddata_dir,"nmr_end_samples.txt"),h=F,stringsAsFactors = F)
nmr_baseline_samples <- read.table(file=paste0(processeddata_dir,"nmr_baseline_samples.txt"),h=F,stringsAsFactors = F)
metab_end_samples <- read.table(file=paste0(processeddata_dir,"metab_end_samples.txt"),h=F,stringsAsFactors = F)
metab_baseline_samples <- read.table(file=paste0(processeddata_dir,"metab_baseline_samples.txt"),h=F,stringsAsFactors = F)

baseline_samples <- list(metab=metab_baseline_samples,nmr=nmr_baseline_samples)
end_samples <- list(metab=metab_end_samples,nmr=nmr_end_samples)


# drop those individuals with missing group allocation and qc fails and no baseline data
for (i in 1:length(model_data_list)) {
  print(paste0("Processing dataset: ",names(model_data_list)[[i]]))
  model_data_list[[i]] <- model_data_list[[i]][!is.na(model_data_list[[i]]$treat),]
  # restrict to only samples that passed qc
  index_list <- which(model_data_list[[i]]$id %in% as.character(end_samples[[i]]$V1))
  model_data_list[[i]] <- model_data_list[[i]][index_list,]
  # check numbers
  print("Those with year 1 data:")
  print(table(model_data_list[[i]]$treat,model_data_list[[i]]$has_end_sample))
  print("Check for baseline data in those with year 1 data")
  print(table(model_data_list[[i]]$has_end_sample,model_data_list[[i]]$has_baseline_sample))
  # remove those without a baseline sample that passed QC
  index_list <- which(model_data_list[[i]]$id %in% as.character(baseline_samples[[i]]$V1))
  model_data_list[[i]] <- model_data_list[[i]][index_list,]
  # check numbers
  print("Final numbers")  
  print(nrow(model_data_list[[i]]))
  print(table(model_data_list[[i]]$treat))
}

# add row names
row.names(model_data_list$metab) <- model_data_list$metab$id
row.names(model_data_list$nmr) <- model_data_list$nmr$id

# write out datasets
saveRDS(model_data_list,file=paste0(processeddata_dir,"model_input.rds"))

# separate metabolites from other pheno
raw_data_list <- list(metab=model_data_list$metab[,17:ncol(model_data_list$metab)],nmr=model_data_list$nmr[,17:ncol(model_data_list$nmr)])
clinic_data_only <- list(metab=model_data_list$metab[,1:16],nmr=model_data_list$nmr[,1:16])

# write out clinic dataset [raw metabolite data included in dataset_list below]
saveRDS(clinic_data_only,file=paste0(processeddata_dir,"clinic_data_input.rds"))

## metabolite feature metadata (post cleaning)
clean_feature_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_feature_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_feature_metadata.rds"))) 

## add column for feature id matching
clean_feature_list[[1]]$id_for_matching <- clean_feature_list[[1]]$compid_to_match
clean_feature_list[[2]]$id_for_matching <- clean_feature_list[[2]]$shortname

## sample metadata (post cleaning)
clean_sample_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_sample_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_sample_metadata.rds")))


####################################################################################
####################################################################################
### derive transformed datasets
####################################################################################
####################################################################################

# set up lists for new versions of data
rnt_data_list <- list(metab=NA,nmr=NA)
natlog_data_list <- list(metab=NA,nmr=NA)
pa_data_list <- list(metab=NA,nmr=NA)
win_data_list <- list(metab=NA,nmr=NA)

## create transformed datasets (metabolon then nmr)
## transformations done within timepoint since in raw data there is a column per metab at a timepoint
for (i in 1:length(raw_data_list)){
  print(paste0("Processing dataset: ",names(raw_data_list)[i]))
  # define dataset
  dtst <- raw_data_list[[i]]
  # rank normal transform
  print("To RNT ...")
  rnt_data_list[[i]] <- apply(dtst,2,function(x) qnorm((rank(x,na.last="keep",ties.method = "random")-0.5)/sum(!is.na(x))))
  # natural log transform
  print("To natural log ...")
  natlog_data_list[[i]] <- apply(dtst,2,function(x) log(x+1))
  # presence absence 
  print("To presence/absence ...")
  pa_data_list[[i]] <- pa.convert(dtst)
  # winsorise
  print("To winsorized ...")
  win_data_list[[i]] <- apply(dtst,2,function(x) winsorize_x(x))
}

dataset_list <- list(raw=raw_data_list,rnt=rnt_data_list,natlog=natlog_data_list,pa=pa_data_list,win=win_data_list)

# save out
saveRDS(dataset_list,file=paste0(processeddata_dir,"alternative_formats.rds"))


####################################################################################
####################################################################################
### run data summary on ALL features
####################################################################################
####################################################################################

# set up list for output
ds_results_list <- list(metab=NA,nmr=NA)

# implement summary
for (i in 1:length(raw_data_list)){
  print(paste0("Processing feature summary stats on raw data for: ",names(raw_data_list)[i]))
  # define dataset
  dtst <- raw_data_list[[i]]
  ds_results_list[[i]] <- as.data.frame(t(apply(dtst,2,function(x) { 
    o = summary(x)
    if(length(o) == 6){
      o = c(o, 0)
      names(o)[7] = "NA's"
    }
    return(o)
  }
  )))
  # add comp id to column
  ds_results_list[[i]]$full_feature_id <- row.names(ds_results_list[[i]])
  ds_results_list[[i]]$feature_id <- sapply(1:nrow(ds_results_list[[i]]), function(x) {
    strsplit(ds_results_list[[i]]$full_feature_id[x],split = "[.]")[[1]][1]
  })
  ds_results_list[[i]]$timepoint <- sapply(1:nrow(ds_results_list[[i]]), function(x) {
    strsplit(ds_results_list[[i]]$full_feature_id[x],split = "[.]")[[1]][2]
  })
  # split into baseline/end
  baseline_stats <- ds_results_list[[i]][ds_results_list[[i]]$timepoint == "b",]
  names(baseline_stats)[1:7] <- paste0(names(baseline_stats)[1:7],"_baseline.samples")
  baseline_stats <- baseline_stats[,c(1:7,9)]
  end_stats <- ds_results_list[[i]][ds_results_list[[i]]$timepoint == "e",]
  names(end_stats)[1:7] <- paste0(names(end_stats)[1:7],"_endpoint.samples")
  end_stats <- end_stats[,c(1:7,9)]
  # merge together
  all_stats <- merge(baseline_stats, end_stats, by="feature_id")
  ds_results_list[[i]] <- all_stats
}

# save out
saveRDS(ds_results_list,file=paste0(results_dir,"metabolite_summaries.rds"))

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()





