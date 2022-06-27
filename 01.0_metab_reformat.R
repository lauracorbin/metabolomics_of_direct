# script to reformat metabolite data
# original format is as extracted from the Metabolon excel file by staff in Glasgow
# derived format is that required to run existing metabolite data qc scripts  

# last run: 31st March 2022

# capture output
sink(file="U:/Logs/01.0_metab_reformat.txt",split=TRUE)

# read in parameter file - contains file paths and file names
pfile = read.table("U:/Scripts/01.0_metab_reformat_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
dataset <- as.character( pfile[2,2] )
processseddata_dir <- as.character( pfile[3,2] )
scripts_dir <- as.character( pfile[4,2] )
rawdata_dir <- as.character( pfile[5,2] )
rawdata_dir_metab <- as.character( pfile[6,2] )
externaldata_dir <- as.character( pfile[7,2] )

####################################################################################
####################################################################################
### read in data
####################################################################################
####################################################################################

# read in orig scale metabolite data plus sample info
alldata <- readRDS(paste0(rawdata_dir_metab,"metabolon_originalscale_20190807.rds"))
print("Dimensions of original scale data: ")
dim(alldata)
print(paste("No. of samples in original scale data is: ", nrow(alldata)))
print(paste("No. of features in original scale data is: ", ncol(alldata)-7))
alldata[1:5,1:6]
alldata[1:5,(ncol(alldata)-5):ncol(alldata)]

# read in ScaledImp metabolite data plus sample info
alldata_imp <- readRDS(paste0(rawdata_dir_metab,"metabolon_scaledimputed_20190807.rds"))
print("Dimensions of imputed data: ")
dim(alldata_imp)
print(paste("No. of samples in imputed data is: ", nrow(alldata_imp)))
print(paste("No. of features in imputed data is: ", ncol(alldata_imp)-7))
alldata_imp[1:5,1:6]
alldata_imp[1:5,(ncol(alldata_imp)-5):ncol(alldata_imp)]

# read in feature metadata
feature_data <- readRDS(paste0(rawdata_dir_metab,"metabolon_info_20190807.rds"))
dim(feature_data)
feature_data[1:5,1:13]

# read in corrected metabolite IDs - based on information received from Metabolon in Jan 2022
new_labels <- read.table(paste0(externaldata_dir,"metabolon_id_update.csv"),sep=",",h=T)

####################################################################################
####################################################################################
### update selected metabolite IDs (following library correction at Metabolon)
####################################################################################
####################################################################################

new_labels$old_compid <- NA

# replace details in feature data
for (i in 1:nrow(new_labels)) {
  # find target row in feature data
  to_replace <- which(feature_data$biochemical == new_labels[i,c("IncorrectID")])
  # continue only if found
  if (length(to_replace) > 0) {
    print(paste0("replacing ", new_labels[i,2], " with ", new_labels[i,3]))
    print(paste0("compid is: ", feature_data[to_replace,c("comp.id")]))
    # add old compid to new_labels so metab data can be updated
    new_labels[i,9] <- feature_data[to_replace,c("comp.id")]
    # update info in feature_data
    feature_data[to_replace,c("comp.id")] <- new_labels[i,5]
    feature_data[to_replace,c("cas")] <- new_labels[i,6]
    feature_data[to_replace,c("pubchem")] <- new_labels[i,7]
    feature_data[to_replace,c("hmdb")] <- new_labels[i,8]
    feature_data[to_replace,c("biochemical")] <- new_labels[i,3]
    feature_data[to_replace,c("kegg")] <- NA
  }
  else {
    print(paste0(new_labels[i,c("IncorrectID")], " not found"))
  }
}

# replace details in metabolite data files (column headers - compids)
# make compids in the right format
new_labels$old_compid <- paste0("compid", new_labels$old_compid)
new_labels$new_compid <- paste0("compid", new_labels$COMP_ID)

# replace in alldata
for (i in 1:nrow(new_labels)) {
  # find target column
  to_replace <- which(names(alldata) == new_labels[i,c("old_compid")])
  # continue if found
  if (length(to_replace) > 0) {
    print(paste0("replacing ", new_labels[i,9], " with ", new_labels[i,10]))
    names(alldata)[to_replace] <- new_labels[i,10]
  }
  else {
    print(paste0(new_labels[i,c("old_compid")], " not found"))
  }
}

# replace in alldata_imp
for (i in 1:nrow(new_labels)) {
  # find target column
  to_replace <- which(names(alldata_imp) == new_labels[i,c("old_compid")])
  # continue if found
  if (length(to_replace) > 0) {
    print(paste0("replacing ", new_labels[i,9], " with ", new_labels[i,10]))
    names(alldata_imp)[to_replace] <- new_labels[i,10]
  }
  else {
    print(paste0(new_labels[i,c("old_compid")], " not found"))
  }
}

####################################################################################
####################################################################################
### make raw data file (orig)
####################################################################################
####################################################################################

# restrict to metab data only - strip out columns at start and end that contain sample metadata
mydata <- alldata[,c(6:1281)] 
dim(mydata)
mydata[1:5,1:5]
mydata[1:5,(ncol(mydata)-5):ncol(mydata)]

# recode to numeric
mydata <- data.frame(lapply(mydata,as.numeric))

# make row names the sample names
row.names(mydata) <- alldata$sample.name

# transpose
mydata_t <- as.data.frame(t(mydata))
dim(mydata_t)
mydata_t[1:5,1:5]
mydata_t[1:5,(ncol(mydata_t)-5):ncol(mydata_t)]

# save out
saveRDS(mydata_t,file=paste0(processseddata_dir,"OrigScale_raw_data.rds"))

####################################################################################
####################################################################################
### make feature metadata file
####################################################################################
####################################################################################

myfeatures <- feature_data
dim(myfeatures)
myfeatures$compid_to_match <- paste0("compid",myfeatures$comp.id)
head(myfeatures)

# save out, assuming scaled metadata the same
saveRDS(myfeatures,file=paste0(processseddata_dir,"OrigScale_feature_metadata.rds"))
saveRDS(myfeatures,file=paste0(processseddata_dir,"ScaledImpData_feature_metadata.rds"))

####################################################################################
####################################################################################
### make sample metadata file
####################################################################################
####################################################################################

# restrict to columns at start and end that contain sample metadata
mysamples <- alldata[,c(1:5,(ncol(alldata)-1),ncol(alldata))]
dim(mysamples)
head(mysamples,20)
table(mysamples$visit)
row.names(mysamples) <- seq(1,nrow(mysamples),1)

# count no. of individuals represented
print(paste0("There are ", length(unique(mysamples$id)), " individuals represented in the raw data."))

# save out, assuming scaled metadata the same
saveRDS(mysamples,file=paste0(processseddata_dir,"OrigScale_sample_metadata.rds"))
saveRDS(mysamples,file=paste0(processseddata_dir,"ScaledImpData_sample_metadata.rds"))

####################################################################################
####################################################################################
### make raw data file (scaled imputed)
####################################################################################
####################################################################################

# restrict to metab data only
mydata_imp <- alldata_imp[,c(6:1281)] 
dim(mydata_imp)
mydata_imp[1:5,1:5]
mydata_imp[1:5,(ncol(mydata_imp)-5):ncol(mydata_imp)]

# recode to numeric
mydata_imp <- data.frame(lapply(mydata_imp,as.numeric))

# make row names the sample names
row.names(mydata_imp) <- alldata_imp$sample.name

# transpose
mydata_imp_t <- as.data.frame(t(mydata_imp))
dim(mydata_imp_t)
mydata_imp_t[1:5,1:5]
mydata_imp_t[1:5,(ncol(mydata_imp_t)-5):ncol(mydata_imp_t)]

# save out
saveRDS(mydata_imp_t,file=paste0(processseddata_dir,"ScaledImpData_raw_data.rds"))

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()