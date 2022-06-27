## preparation of metabolomics data from Nightingale to run existing metabolite data qc scripts ##
## input = raw data as received from lab in excel file ##
## output = three .rds files per data worksheet in excel file. Need to be read in with readRDS command. ##

# last run: 31st March 2022

# capture output
sink(file="U:/Logs/01.0_nmr_reformat.txt",split=TRUE)

## load required libraries
library(readxl)

#########################
## PREPARE ENVIRONMENT ##
#########################

# read in parameter file - contains file paths and file names
pfile = read.table("U:/Scripts/01.0_nmr_reformat_param_file.R", sep = "=", as.is=T)

# set paramaters
project <- as.character( pfile[1,2] )
dataset <- as.character( pfile[2,2] )
processseddata_dir <- as.character( pfile[3,2] )
scripts_dir <- as.character( pfile[4,2] )
rawdata_dir <- as.character( pfile[5,2] )
n_feature_meta_rows <- as.numeric( pfile[6,2] )
n_sample_meta_col <- as.numeric( pfile[7,2] )
feature_id_row <- as.numeric( pfile[8,2] ) # this needs to be a unique feature identifier - must be present for all metabolites
sample_id_col <- as.character( pfile[9,2] ) # this needs to be a unique sample identifier 
last_sample <- as.numeric( pfile[10,2] )
last_feature <- as.character( pfile[11,2] )
filename <- paste0(rawdata_dir,as.character( pfile[12,2] ))


#######################
## METABOLOMICS DATA ##
#######################

## find start row of data
start_row <- n_feature_meta_rows + 1

## find start column of data
start_col <- LETTERS[n_sample_meta_col + 1]

## define first cell of data
start_cell <- paste0(start_col,start_row)
finish_cell <- paste0(last_feature,last_sample)
cell_range <- paste(start_cell,finish_cell,sep=":")
print(paste0("Metabolite data contained in the following grid area in data worksheets:",cell_range, ". PLEASE CONFIRM THIS IS CORRECT BY EXAMINING THE EXCEL FILE!!" ))

## read in main body (i.e. metabolomics data only) of 'n' versions of data (appear on worksheets 2,3,4...)
sheet_names <- excel_sheets(filename)
n_sheets <- length(sheet_names)
n_datasets <- n_sheets
print(paste0("This data consists of ", n_datasets, " datasets. PLEASE CONFIRM THIS IS CORRECT BY EXAMINING THE EXCEL FILE!!"))

mydatalist <- lapply(sheet_names[1:n_sheets], function(f){
  out = as.data.frame(read_excel(filename, sheet = f,range = cell_range,col_names=F,col_types="text"))
  return(out)
})

## name each dataset according to sheet name
names(mydatalist) <- sheet_names[1:n_sheets]

## add row (feature) and column (sample) names
## grab relevant row/column from excel file as specified in arguments
start_feature_ids <- paste0(start_col,as.character(feature_id_row))
finish_feature_ids <- paste0(last_feature,as.character(feature_id_row))
cell_range_feature_ids <- paste(start_feature_ids,finish_feature_ids,sep=":")

start_sample_ids <- paste0(sample_id_col,start_row)
finish_sample_ids <- paste0(sample_id_col,last_sample)
cell_range_sample_ids <- paste(start_sample_ids,finish_sample_ids,sep=":")

mysampleids <- lapply(sheet_names[1:n_sheets], function(f){
  out = as.data.frame(read_excel(filename, sheet = f,range = cell_range_sample_ids,col_names=F))
  return(out)
})

myfeatureids <- lapply(sheet_names[1:n_sheets], function(f){
  out = as.data.frame(read_excel(filename, sheet = f,range = cell_range_feature_ids,col_names=F))
  return(out)
})

##############################
## PROCESS FEATURE METADATA ##
##############################  

## feature metadata
## get coordinates for grid of sample information
start_feature_meta <- paste0(LETTERS[n_sample_meta_col+1],"9")
finish_feature_meta <- paste0(last_feature, as.character(n_feature_meta_rows))
cell_feature_meta <- paste(start_feature_meta, finish_feature_meta, sep=":")

myfeaturemetadata <- lapply(sheet_names[1:n_sheets], function(f){
  out = as.data.frame(read_excel(filename, sheet = f,range = cell_feature_meta, col_names=F))
  transpose <- as.data.frame(t(out))
  rownames(transpose) <- transpose[,2]
  # convert all columns from factor to character
  transpose[] <- lapply(transpose, as.character)
  # convert success % to numeric
  transpose$V4 <- as.numeric(transpose$V4)
  # name columns
  names(transpose)[1] <- "fullname"
  names(transpose)[2] <- "shortname"
  names(transpose)[3] <- "units"
  names(transpose)[4] <- "percent_success"
  return(transpose)
})

## name each dataset according to sheet name
names(myfeaturemetadata) <- sheet_names[1:n_sheets]

## save out dataframes individually
for(i in 1:n_datasets){
  # grab label for dataset
  label = names(myfeaturemetadata)[i]
  filename_temp = paste0(processseddata_dir,"nmr_feature_metadata.rds")
  saveRDS(myfeaturemetadata[[i]],file = filename_temp)
}

## provide summary of sample metadata
print("The following column names exist in the feature metadata files: ")
names(myfeaturemetadata[[1]])
#print(paste0("The unique sample identifier being used is: ", names(myfeaturemetadata[[1]][feature_id_row])))
print(paste0("The data contains: ", nrow(myfeaturemetadata[[1]]), " features."))

#############################
## PROCESS SAMPLE METADATA ##
#############################

## sample metadata
## get coordinates for grid
start_sample_meta <- paste0("A",n_feature_meta_rows+1)
finish_sample_meta <- paste0(LETTERS[n_sample_meta_col],last_sample)
cell_sample_meta <- paste(start_sample_meta,finish_sample_meta,sep=":")

mysamplemetadata <- lapply(sheet_names[1:length(sheet_names)], function(f){
  out = as.data.frame(read_excel(filename, sheet = f,range = cell_sample_meta, col_names=F))
  rownames(out) <- out[,5]
  # name columns
  names(out)[1]<- "High pyruvate" 
  names(out)[2]<- "High lactate" 
  names(out)[3]<- "Low glutamine/high glutamate" 
  names(out)[4]<- "Plasma sample"
  names(out)[5]<- "sampleid"
  # name rows
  rownames(out) <- out[,5]
  # add sampleis to ids
  out$sampleid <- paste0("sampleid.",out$sampleid)
  return(out)
})

## name each dataset
names(mysamplemetadata) <- sheet_names[1:n_sheets]

## save out dataframes individually
for(i in 1:n_datasets){
  # grab label for dataset
  label = names(mysamplemetadata)[i]
  filename_temp = paste0(processseddata_dir,"nmr_sample_metadata.rds")
  saveRDS(mysamplemetadata[[i]],file = filename_temp)
}

## provide summary of feature metadata
print("The following column names exist in the sample metadata files: ")
names(mysamplemetadata[[1]])
print(paste0("The data contains: ", nrow(mysamplemetadata[[1]]), " samples."))

#########################
## FINALISE METAB DATA ##
#########################

## make feature IDs the row names and sample IDs the col names
for(i in 1:n_datasets){
  # generate row names
  row.names(mydatalist[[i]]) = paste0("sampleid.",row.names(mysamplemetadata[[i]]))
  # generate column names
  names(mydatalist[[i]]) = row.names(myfeaturemetadata[[i]])
}

## save out dataframes individually
for(i in 1:n_datasets){
  # grab label for dataset
  label = names(mydatalist)[i]
  data <- mydatalist[[i]]
  data <- as.data.frame(lapply(data, function(x) {
    gsub("NDEF", "-9", x) }),stringsAsFactors = F)
  data <- data.frame(lapply(data, function(x) {
    gsub("TAG", "-99", x) }),stringsAsFactors = F)
  data <- data.frame(lapply(data, function(x) {
    gsub("NA", "-999", x) }),stringsAsFactors = F)
  # change to numeric
  data <- data.frame(lapply(data,as.numeric))
  # transpose
  transposed <- as.data.frame(t(data))
  # replace sample ids
  names(transposed) <- paste0("sampleid.",row.names(mysamplemetadata[[i]]))
  # check dimensions
  print(paste0("The ",label," dataset contains: ", ncol(transposed), " samples and ", nrow(transposed), " features/metabolites."))
  filename_temp = paste0(processseddata_dir,"nmr_raw_data.rds")
  saveRDS(transposed,file = filename_temp)
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
