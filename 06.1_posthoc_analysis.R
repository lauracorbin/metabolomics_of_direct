# This script runs posthoc analyses following the EFFECT OF THE INTERVENTION ON METABOLITES analysis

# last run: 19th July 2023

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/06.1_posthoc_analysis.txt",split=TRUE)

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
source(paste0(scripts_dir,"02.0_qc_functions.R"))
source(paste0(scripts_dir,"04.0_functions.R"))
source(paste0(scripts_dir,"RPackages/growthPheno/R/PVA.R"))
source(paste0(scripts_dir,"05.0_functions.R"))
source(paste0(scripts_dir,"06.0_functions.R"))

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in  datasets
# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))

# read in results files
lm_results_list <- readRDS(file=paste0(results_dir,"lm_results_list.rds"))
metab_pa_results <- readRDS(file=paste0(results_dir,"metab_pa_results.rds"))
clinical_lm_results <- readRDS(file=paste0(results_dir,"clinical_lm_results.rds"))

# bring in list of features with less missing than threshold
lm_features <- as.character((read.table(file=paste0(results_dir,"lm_feature_list.txt"), h=F, sep="\t"))$V1)

## bring in feature lists from selection step
pa_features <- as.character((read.table(file=paste0(results_dir,"high_missingness_feature_list.txt"),sep="\t", stringsAsFactors = F))$V1)

## read in clustering groups
clusters <- read.table(file=paste0(results_dir,"feature_group_allocation.txt"),sep="\t", stringsAsFactors = F,h=T)

# read in PVA results
pva_output <- read.table(file=paste0(results_dir,"pva_output_by_feature.txt"),sep="\t", stringsAsFactors = F,h=T)

# read in residuals from mixed model
metab_lm_resid <- readRDS(file=paste0(processeddata_dir,"metabolon_lm_rnt_resid.rds"))
nmr_lm_resid <- readRDS(file=paste0(processeddata_dir,"nmr_lm_rnt_resid.rds"))

# read in feature metadata for metab features
metab_feature_metadata <- readRDS(paste0(processeddata_dir,"OrigScale_feature_metadata_with_QCstats.rds"))
nmr_feature_metadata <- readRDS(paste0(processeddata_dir,"nmr_feature_metadata_with_QCstats.rds"))

# read in labels for NMR 
load(file=paste0(ext_data_dir,"metabolite_labels.Rda"))

# read in trial data
trial_data <- readRDS(paste0(rawdata_dir,"DiRECT_fu1_20171113.rds"))

####################################################################################
####################################################################################
## Make data file with features to follow up (from linear model)
####################################################################################
####################################################################################

# extract list of lm features to follow up
# (p<0.05)
metab_lm_assoc_features <- lm_results_list$metab[lm_results_list$metab$rnt_lm_treat_HolmAdjP_flag == 1, "feature_id"]
# extract from pva_output
metab_pva_output_assoc <- pva_output[pva_output$Selected %in% metab_lm_assoc_features,] 

nmr_lm_assoc_features <- lm_results_list$nmr[lm_results_list$nmr$rnt_lm_treat_HolmAdjP_flag == 1, "feature_id"]
# extract from pva_output
nmr_pva_output_assoc <- pva_output[pva_output$Selected %in% nmr_lm_assoc_features,] 

#merge metab and nmr
# pva output for all associated
pva_output_assoc <- rbind(metab_pva_output_assoc, nmr_pva_output_assoc)
# list of all associated
all_assoc_features <- c(metab_lm_assoc_features,nmr_lm_assoc_features)

# select one var per cluster to take forward
# select the metabolite ranked highest by the PVA (i.e. the one that explains most variance in the cluster)
assoc_clusters <- unique(pva_output_assoc$cluster)
selected_assoc_features <- NULL
for (i in 1:length(assoc_clusters)) {
  pva_sub <- pva_output_assoc[pva_output_assoc$cluster == assoc_clusters[i],]
  pva_sub <- pva_sub[order(pva_sub$Variable),]
  var_to_keep <- pva_sub[1,"Selected"]
  if (i == 1) {
    selected_assoc_features <- var_to_keep
  } else {
    selected_assoc_features <- c(selected_assoc_features, var_to_keep)
  }
}


# extract residuals for selected associated features
cols_to_keep <- which(names(metab_lm_resid) %in% selected_assoc_features)
metab_resid_selected <- metab_lm_resid[,c(1,cols_to_keep)]

cols_to_keep <- which(names(nmr_lm_resid) %in% selected_assoc_features)
nmr_resid_selected <- nmr_lm_resid[,c(1,cols_to_keep)]

# merge nmr and metabolon
resid_selected_assoc <- merge(metab_resid_selected, nmr_resid_selected, by="id", all=T)

# save out
saveRDS(resid_selected_assoc, file = paste0(processeddata_dir,"resid_selected_assoc.rds")) 

# save out list of associated representative features 
write.table(as.data.frame(selected_assoc_features), file=paste0(results_dir,"selected_assoc_features.txt"),col.names=F, row.names=F,quote=F)


####################################################################################
####################################################################################
## Establish labels
####################################################################################
####################################################################################

## METABOLON ##
## sort out labelling of super-pathways (metabolon) and classes (NMR) and combine within a single nomenclature

## pull out superpathway info for features that passed QC
metab_classes <- metab_feature_metadata[metab_feature_metadata$qc_flag == "pass",c("compid_to_match", "super.pathway", "biochemical")]
# add class column to enable merge with NMR
metab_classes$class <- NA

## NMR ##
# attach NMR labels
#remove units
metabolite_labels$raw.label <- gsub("\\(.*","",metabolite_labels$raw.label)
# remove trailing white space
metabolite_labels$raw.label <- trimws(metabolite_labels$raw.label,which=c("right"))
# edit remnant-c label
metabolite_labels[which(metabolite_labels$raw.label == "Remnant cholesterol"), "raw.label"] <- "Remnant cholesterol (non-HDL, non-LDL -cholesterol)"
# rename col to match source data
names(metabolite_labels)[names(metabolite_labels)=="raw.label"] <- "fullname"
# merge with feature file
nmr_classes <- merge(metabolite_labels,nmr_feature_metadata,by="fullname",all.y=T)
# replace NA with other ratio designation
nmr_classes[is.na(nmr_classes$class),c("class")] <- "Other ratio/percentage"

## identify derived and raw measures
percent_measures <- nmr_classes[grep("%",nmr_classes$shortname),"shortname"]
ratio_measures <- nmr_classes[grep("/",nmr_classes$shortname),"shortname"]
derived_measures <- c(percent_measures, ratio_measures)
raw_measures <- nmr_classes[!nmr_classes$shortname %in% derived_measures,"shortname"]

# save out lists of derived versus raw
write.table(derived_measures, file=paste0(results_dir,"NMRderived.txt"),col.names=F, row.names=F,quote=F)
write.table(raw_measures, file=paste0(results_dir,"NMRraw.txt"),col.names=F, row.names=F,quote=F)

## restrict to QC passes and raw measures
nmr_classes <- nmr_classes[nmr_classes$qc_flag == "pass",]
#nmr_classes <- nmr_classes[nmr_classes$shortname %in% raw_measures,]
nmr_classes <- nmr_classes[,c("shortname","class")]

# summarise within panel
cat(paste("Distribution of metabolites across classes in NMR: "))
table(nmr_classes$class,useNA="always")
cat(paste("  ", "Distribution of metabolites across super.pathways in Metabolon: ", sep="\n"))
table(metab_classes$super.pathway,useNA="always")

## add labels to match Metabolon (based on Wahl et al 2015 designations)
nmr_classes$super.pathway <- NA
nmr_classes[which(nmr_classes$class == "Amino acids"),c("super.pathway")] <- "Amino Acid"
nmr_classes[which(nmr_classes$class == "Apolipoproteins"),c("super.pathway")] <- "Peptide"
nmr_classes[which(nmr_classes$class == "Cholesterol"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$class == "Fatty acids"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$class == "Fatty acids ratios"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$shortname == "Alb"),c("super.pathway")] <- "Peptide"
nmr_classes[which(nmr_classes$shortname == "Crea"),c("super.pathway")] <- "Amino Acid"
nmr_classes[which(nmr_classes$class == "Glycerides and phospholipids"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$shortname == "Cit"),c("super.pathway")] <- "Energy"
nmr_classes[which(nmr_classes$shortname == "Glc"),c("super.pathway")] <- "Carbohydrate"
nmr_classes[which(nmr_classes$shortname == "Glol"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$shortname == "Lac"),c("super.pathway")] <- "Carbohydrate"
nmr_classes[which(nmr_classes$shortname == "Pyr"),c("super.pathway")] <- "Carbohydrate"
nmr_classes[which(nmr_classes$shortname == "Gp"),c("super.pathway")] <- "Peptide"
nmr_classes[which(nmr_classes$shortname == "bOHBut"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$shortname == "Ace"),c("super.pathway")] <- "Energy"
nmr_classes[which(nmr_classes$class == "Lipoprotein particle size"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$class == "Lipoprotein subclasses"),c("super.pathway")] <- "Lipid"
nmr_classes[which(nmr_classes$class == "Other ratio/percentage"),c("super.pathway")] <- "NMR ratio/percentage"

cat(paste("  ", "Distribution of metabolites across super.pathways in NMR (after manual labelling procedure): ", sep="\n"))
table(nmr_classes$super.pathway,useNA="always")

# add biochem ids
nmr_classes$biochemical <- nmr_feature_metadata$fullname[match(nmr_classes$shortname,nmr_feature_metadata$shortname)]

# rename feature id columns to enable merging
names(metab_classes)[1] <- "feature_id"
names(nmr_classes)[1] <- "feature_id"

#add panel ID column
metab_classes$panel <- "metabolon"
nmr_classes$panel <- "nightingale"

# merge
class_info <- rbind(nmr_classes,metab_classes)

# summarise
cat(paste("  ", "Distribution of metabolites across super.pathways in NMR & Metabolon combined: ", sep="\n"))
table(class_info$super.pathway,useNA="always")

#replace NA with unclassified
class_info[which(is.na(class_info$super.pathway)),c("super.pathway")] <- "Unclassified"

# save out
saveRDS(class_info, file=paste0(results_dir,"class_info_for_enrichment.rds"))

####################################################################################
####################################################################################
## Run enrichment analyses
####################################################################################
####################################################################################

## 1. look at enrichment of associated features from LM to all features run in LM
print("Enrichment of associated in linear model to all linear model features")

# extract superpathways for list of associated
lm_all_assoc_spp <- class_info[class_info$feature_id %in% all_assoc_features, c("super.pathway")]
print("Distribution of classes in associated MM features:")
table(lm_all_assoc_spp, useNA="always")

# extract superpathways for all features that went through linear model analysis
all_lm_spp <- class_info[class_info$feature_id %in% lm_features, c("super.pathway")]
print("Distribution of classes in all mixed model features:")
table(all_lm_spp, useNA="always")

# apply enrichment function
lm_enrich <- hgtest(all_lm_spp, lm_all_assoc_spp)
print("Output summary for 'Enrichment of associated in linear model to all linear model features':")
lm_enrich


## look at enrichment of metabolite clusters (i.e. those metabolites in the same cluster) to full dataset
## aim is to characterise each cluster

# make empty list to store lists of features in each associated cluster
cluster_groups <- vector(mode = "list", length=length(unique(clusters$k_group)))

# populate list with feature lists
for (i in 1:length(unique(clusters$k_group)) ) {
  #metab_group <- clusters[which(clusters$feature_id == selected_assoc_features[i]),"k_group"]
  cluster_groups[[i]] <- as.character(clusters[clusters$k_group == i,c("feature_id")])
  names(cluster_groups)[[i]] <- paste0(i,"_grp")
}

# make empty list to add enrichment
clusters_with_enrich <- vector(mode = "list", length=length(cluster_groups))

# loop over clusters
for (i in 1:length(cluster_groups) ) {
  # extract superpathways for metabolites in each cluster
  cluster_spp <- class_info[class_info$feature_id %in% cluster_groups[[i]], c("super.pathway")]
  cat(paste("Starting next cluster ...", "  ", "  ", sep="\n"))
  print(paste0("Enrichment results for group of: ",names(cluster_groups)[[i]]))
  print(paste0("There are ", length(cluster_spp), " metabolites in this cluster."))
  # run enrichment analysis against all features that went through mixed model analysis (and were therefore used in the clustering exercise)
  cluster_enrich <- hgtest(all_lm_spp, cluster_spp)
  clusters_with_enrich[[i]] <- cluster_enrich
  names(clusters_with_enrich)[[i]] <- names(cluster_groups)[[i]]
  print(clusters_with_enrich[[i]])
}


# save out cluster enrichment results
saveRDS(clusters_with_enrich, file=paste0(results_dir,"clusters_with_enrich.rds"))


####################################################################################
####################################################################################
## Check for baseline differences and make boxplots
####################################################################################
####################################################################################

# vectors to collect p's
ranktest_output = NULL
metab_output = NULL
nmr_output = NULL

for (j in 1:length(dataset_list$raw)){
#for (j in 2:2){
  # pick associated list
  if (j==1) {
    assoc_features_list <- metab_lm_assoc_features
    class_file = metab_classes
  }
  else {
    assoc_features_list <- nmr_lm_assoc_features
    class_file = nmr_classes
  }
  # extract baseline
  baseline_data <- select(dataset_list$raw[[j]],contains(".b"))
  # remove .b
  names(baseline_data) <- gsub("[.]b","",names(baseline_data))
  # restrict to associated
  baseline_data_assoc <- baseline_data[,assoc_features_list] 
  #zscore
  baseline_data_assoc_z <- as.data.frame(apply(baseline_data_assoc,2,function(x) scale(x,center=F,scale=F)))
  # add timepoint label
  baseline_data_assoc_z$timepoint <- "baseline"                             
  # add IDs
  baseline_data_assoc_z$id <- rownames(baseline_data)
  
  # extract endpoint
  endpoint_data <- select(dataset_list$raw[[j]],contains(".e"))
  # remove .e
  names(endpoint_data) <- gsub("[.]e","",names(endpoint_data))
  # restrict to associated
  endpoint_data_assoc <- endpoint_data[,assoc_features_list]
  #zscore
  endpoint_data_assoc_z <- as.data.frame(apply(endpoint_data_assoc,2,function(x) scale(x,center=F,scale=F)))
  # add timepoint label
  endpoint_data_assoc_z$timepoint <- "endpoint (12 months)"
  # add IDs
  endpoint_data_assoc_z$id <- rownames(endpoint_data)
  
  # join by rbind
  data_assoc_z <- rbind(baseline_data_assoc_z,endpoint_data_assoc_z)
  
  # add allocation
  data_assoc_z$allocation <- trial_data$treat[match(data_assoc_z$id,trial_data$id)]
  
  # check numbers
  print(table(data_assoc_z$timepoint,data_assoc_z$allocation))

  metabolite_plots = lapply(assoc_features_list, function(feature){
    metab <- data_assoc_z[,c(feature, "timepoint", "allocation")]
    metab_id <- names(metab[1])
    metab_name <- class_file[which(class_file$feature_id == metab_id),"biochemical"]
    names(metab)[1] <- "metabolite"
    # ttest for difference at baseline and end
    baseline_data <- metab[metab$timepoint == "baseline",]
    ranktest_result_baseline <- wilcox.test(x=baseline_data[baseline_data$allocation == "Intervention",c("metabolite")],
                                   y=baseline_data[baseline_data$allocation == "Control",c("metabolite")])
    baseline_p <- format(ranktest_result_baseline$p.value,format="e",digits=3)
    if (j==1){
      metab_output <- c(metab_output,baseline_p)
      ylabel <- "Metabolite peak area abundance"
    }
    else {
      nmr_output <- c(nmr_output,baseline_p)
      unit_measure <- nmr_feature_metadata[which(nmr_feature_metadata$shortname == feature),c("units")]
      ylabel <- paste0("Metabolite concentration (", unit_measure, ")")
    }
    endpoint_data <- metab[metab$timepoint == "endpoint (12 months)",]
    ranktest_result_endpoint <- wilcox.test(x=endpoint_data[endpoint_data$allocation == "Intervention",c("metabolite")],
                                           y=endpoint_data[endpoint_data$allocation == "Control",c("metabolite")])
    endpoint_p <- format(ranktest_result_endpoint$p.value,format="e",digits=3)
    baseline_text <- data.frame(timepoint = c("baseline"),metabolite = c(max(metab$metabolite,na.rm=T)+(0.1*max(metab$metabolite,na.rm=T))),
                                allocation=as.factor(c("Control"))) 
    endpoint_text <- data.frame(timepoint = c("endpoint (12 months)"),metabolite = c(max(metab$metabolite,na.rm=T)+(0.1*max(metab$metabolite,na.rm=T))),
                                allocation=as.factor(c("Control"))) 
    # plot
    p1 <- ggplot(metab, aes(y=metabolite, x=timepoint)) +
      geom_boxplot(aes(col=allocation),outlier.shape=NA) +
      labs(x="Time point", y=ylabel, fill="Allocation group", title=metab_name) +
      theme(axis.text = element_text(size=8, colour = "black"),
            axis.title = element_text(size=8, colour = "black"),
            plot.title = element_text(size=8, colour = "black"),
            legend.title = element_text(size=8, colour = "black"),
            legend.text = element_text(size=8, colour = "black"),
            legend.position = "bottom") +
      geom_point(position=position_jitterdodge(),aes(col=allocation),size=0.3)
    p2 <- p1 + geom_text(data=baseline_text,label=paste0("p=",baseline_p),size=2.5) +
                geom_text(data=endpoint_text,label=paste0("p=",endpoint_p),size=2.5)
    return(p2)
  })
  if (j==1) {
     metab_plots <- metabolite_plots 
     }
  else {
     nmr_plots <- metabolite_plots
  }
  
  metabolite_ranktest <- lapply(assoc_features_list, function(feature){
    metab <- data_assoc_z[,c(feature, "timepoint", "allocation")]
    metab_id <- names(metab[1])
    metab_name <- class_file[which(class_file$feature_id == metab_id),"biochemical"]
    names(metab)[1] <- "metabolite"
    # ttest for difference at baseline and end
    baseline_data <- metab[metab$timepoint == "baseline",]
    ranktest_result_baseline <- wilcox.test(x=baseline_data[baseline_data$allocation == "Intervention",c("metabolite")],
                                            y=baseline_data[baseline_data$allocation == "Control",c("metabolite")])
    baseline_p <- format(ranktest_result_baseline$p.value,format="e",digits=3)
    ranktest_output <- c(ranktest_output,baseline_p)
    return(ranktest_output)
  })
    if (j==1){
      metab_output <- metabolite_ranktest
    }
    else {
      nmr_output <- metabolite_ranktest
    }

}

# calculate bonferroni correction N
n_adjust <- length(nmr_output) + length(metab_output)

print("Number of associated in NMR (linear model):")
length(nmr_output)
print("Number associated at baseline (p<0.05/186)")
length(which(nmr_output <= (0.05/n_adjust)))
print("Number associated at baseline (p<0.05)")
length(which(nmr_output <= (0.05)))

print("Number of associated in Metabolon (linear model):")
length(metab_output)
print("Number associated at baseline (p<0.05/186)")
length(which(metab_output <= (0.05/n_adjust)))
print("Number associated at baseline (p<0.05)")
length(which(metab_output <= (0.05)))

## write plots to PDF
filename = paste0(results_dir,"Figures\\ForPaper\\GitHub_Figures_metab_boxplots.pdf")
outputfig <- ggpubr::ggarrange(plotlist=metab_plots,
                  ncol=2, nrow=3,
                  common.legend = T,
                  legend="top")
ggpubr::ggexport(outputfig, filename=filename)

filename = paste0(results_dir,"Figures\\ForPaper\\GitHub_Figures_nmr_boxplots.pdf")
outputfig <- ggpubr::ggarrange(plotlist=nmr_plots,
                               ncol=2, nrow=3,
                               common.legend = T,
                               legend="top")
ggpubr::ggexport(outputfig, filename=filename)

#############################################################################################
#############################################################################################
## Calculate metab change 
#############################################################################################
#############################################################################################

# list of associated metabolites to include
compids_to_include <- c(metab_lm_assoc_features, nmr_lm_assoc_features)
  
# add timepoint extensions
compids_to_include_e <- paste0(compids_to_include, ".e")
compids_to_include_b <- paste0(compids_to_include, ".b")
compids_to_include_all <- c(compids_to_include_e,compids_to_include_b)

# restrict raw data to selected metabs
metab_features_for_boxplots <- as.data.frame(dataset_list$raw$metab[,names(as.data.frame(dataset_list$raw$metab)) %in% compids_to_include_all])
metab_features_for_boxplots$id <- row.names(metab_features_for_boxplots)
nmr_features_for_boxplots <- as.data.frame(dataset_list$raw$nmr[,names(as.data.frame(dataset_list$raw$nmr)) %in% compids_to_include_all])
nmr_features_for_boxplots$id <- row.names(nmr_features_for_boxplots)

features_for_boxplots <- merge(nmr_features_for_boxplots,metab_features_for_boxplots,by="id",all=T)

# calculate delta
delta_phenos <- data.frame(matrix(data=NA,nrow=nrow(features_for_boxplots),ncol=2*length(compids_to_include)))
delta_phenos$id <- features_for_boxplots$id

for (i in 1:length(compids_to_include)) {
  print(compids_to_include[i])
  # add column headers
  names(delta_phenos)[i] <- compids_to_include[i]
  names(delta_phenos)[i+length(compids_to_include)] <- paste0(compids_to_include[i],"_scaled")
  # calculate delta
  baseline_id <- compids_to_include_b[i]
  endpoint_id <- compids_to_include_e[i]
  baseline_col <- which(names(features_for_boxplots) == baseline_id)
  baseline_data <- features_for_boxplots[,c(1,baseline_col)]
  endpoint_col <- which(names(features_for_boxplots) == endpoint_id)
  endpoint_data <- features_for_boxplots[,c(1,endpoint_col)]
  baseline_data$time <- "b"
  endpoint_data$time <-"e"
  names(baseline_data)[2] <- "metabolite"
  names(endpoint_data)[2] <- "metabolite"
  all_data <- rbind(baseline_data,endpoint_data)
  all_data$rnt_metabolite <- qnorm((rank(all_data$metabolite,na.last="keep",ties.method = "random")-0.5)/sum(!is.na(all_data$metabolite)))
  baseline_data <- all_data[all_data$time == "b",]
  endpoint_data <- all_data[all_data$time == "e",]
  delta_phenos[,i] <- endpoint_data$rnt_metabolite - baseline_data$rnt_metabolite
  delta_phenos[,i+length(compids_to_include)] <- delta_phenos[,i]/sd(delta_phenos[,i],na.rm=T)
}


# add weight change
delta_phenos$weight.change <- trial_data$weight.change[match(delta_phenos$id,trial_data$id)]
# add weight change quantiles
weight_quantiles <- quantile(delta_phenos$weight.change,c(0:4/4))
weight_quantiles
delta_phenos$weight.change.quantiles <- with(delta_phenos, cut(weight.change,weight_quantiles,include.lowest=T,
                                                              labels=c("Q1:\n-31.6 to -9.0kg","Q2:\n-9.0 to -3.8kg","Q3:\n-3.8 to 0.1kg","Q4:\n0.1 to 13.7kg")))

# add remission status 
delta_phenos$Remission.status <- trial_data$reversal[match(delta_phenos$id,trial_data$id)]
delta_phenos$Remission.status.lm <- 0
delta_phenos <- delta_phenos %>% mutate(`Remission.status.lm` = replace(`Remission.status.lm`, `Remission.status` == "Yes",1 ))


# check frequencies
#table(delta_phenos$weight.change.quantiles,delta_phenos$Remission.status)

# add group/class
delta_phenos$Allocation <- trial_data$treat[match(delta_phenos$id,trial_data$id)]


#############################################################################################
#############################################################################################
## ttest of metabolite change within weight quantiles
#############################################################################################
#############################################################################################

# set up results dataframe
tt_results <- as.data.frame(matrix(data = NA, nrow = length(compids_to_include), ncol = 11))

for (i in 1:length(compids_to_include)){
  metab_id <- compids_to_include[i]
  metab <- delta_phenos[,c(metab_id,"weight.change.quantiles", "Remission.status","weight.change")]
  if (identical(grep("compid",compids_to_include[i]), integer(0)) == TRUE){
    metab_name <- nmr_classes[which(nmr_classes$feature_id == metab_id),"biochemical"]
  }
  else {
    metab_name <- metab_classes[which(metab_classes$feature_id == metab_id),"biochemical"]
  }
  names(metab)[1] <- "metabolite"
  # separate Q2 and Q2
  q1data <- metab[metab$weight.change.quantiles == "Q1:\n-31.6 to -9.0kg",]
  q2data <- metab[metab$weight.change.quantiles == "Q2:\n-9.0 to -3.8kg",]
  # run ttest
  t.test(q1data$metabolite ~ q1data$Remission.status,var.equal=F)
  t.test(q2data$metabolite ~ q2data$Remission.status,var.equal=F)
  # check correlation between weight change and metabolite change
  cor.test(q1data$metabolite, q1data$weight.change)
  cor.test(q2data$metabolite, q2data$weight.change)
  # put results in dataframe
  tt_results[i,1] <- metab_name
  tt_results[i,2] <- t.test(q1data$metabolite ~ q1data$Remission.status,var.equal=F)[3]
  tt_results[i,3] <- cor.test(q1data$metabolite, q1data$weight.change)[[4]]
  tt_results[i,4] <- cor.test(q1data$metabolite, q1data$weight.change)[[9]][1]
  tt_results[i,5] <- cor.test(q1data$metabolite, q1data$weight.change)[[9]][2]
  tt_results[i,6] <- cor.test(q1data$metabolite, q1data$weight.change)[[3]]
  tt_results[i,7] <- t.test(q2data$metabolite ~ q2data$Remission.status,var.equal=F)[3]
  tt_results[i,8] <- cor.test(q2data$metabolite, q2data$weight.change)[[4]]
  tt_results[i,9] <- cor.test(q2data$metabolite, q2data$weight.change)[[9]][1]
  tt_results[i,10] <- cor.test(q2data$metabolite, q2data$weight.change)[[9]][2]
  tt_results[i,11] <- cor.test(q2data$metabolite, q2data$weight.change)[[3]]
  if (metab_name == "1,5-anhydroglucitol (1,5-AG)" | metab_name == "fructose" | metab_name == "mannose" | metab_name == "glucose" ) {
    print(metab_name)
    print(summary(glm(q1data$Remission.status ~ q1data$weight.change + q1data$metabolite,family = "binomial")))
    }
}
names(tt_results)[1] <- "Biochemical.name"
names(tt_results)[2] <- "Q1.weight.change_ttest_p"
names(tt_results)[3] <- "Q1.weight.change_cortest_r"
names(tt_results)[4] <- "Q1.weight.change_cortest_r_lci"
names(tt_results)[5] <- "Q1.weight.change_cortest_r_uci"
names(tt_results)[6] <- "Q1.weight.change_cortest_p"
names(tt_results)[7] <- "Q2.weight.change_ttest_p"
names(tt_results)[8] <- "Q2.weight.change_cortest_r"
names(tt_results)[9] <- "Q2.weight.change_cortest_r_lci"
names(tt_results)[10] <- "Q2.weight.change_cortest_r_uci"
names(tt_results)[11] <- "Q2.weight.change_cortest_p"

# add source info
tt_results$Source.platform <- class_info$panel[match(lm_results$Biochemical.name,class_info$biochemical)]

# save out
filename <- paste0(results_dir,"Tables\\TableS4_weight_and_metabolite_on_remission.txt")
write.table(tt_results,file=filename,sep="\t",row.names=F,quote=F)

#############################################################################################
#############################################################################################
## plots of remission on metabolite change and weight change to show interactions
#############################################################################################
#############################################################################################

# list of metabolites to include
metabs_to_include <- c("glucose","1,5-anhydroglucitol (1,5-AG)","fructose","mannose")

# extract comp ids
compids_to_include <- class_info[class_info$biochemical %in% metabs_to_include,c("feature_id")]

# find columns to extract 
extract_cols <- c("id", "weight.change", "Remission.status", "Remission.status.lm", "Allocation", "weight.change.quantiles")
extract_cols <- c(extract_cols, compids_to_include)
delta_pheno_cols <- which(names(delta_phenos) %in% extract_cols)
delta_pheno_selected <- delta_phenos[,delta_pheno_cols]

# boxplots by remission and weight change
metabolite_plots_by_weight = lapply(compids_to_include, function(feature){
  metab <- delta_pheno_selected[,c(feature, "weight.change.quantiles", "Remission.status")]
  metab_id <- names(metab[1])
  metab_name <- metab_classes[which(metab_classes$feature_id == metab_id),"biochemical"]
  
  names(metab)[1] <- "metabolite"
  p <- ggboxplot(metab,x="weight.change.quantiles",y="metabolite",
               color="Remission.status",
               add="jitter",font.label=list(size=4),
               ylab = "Change in metabolite\n(normalised SD units)", xlab="Weight change quantiles", title="",
               palette = gray.colors(n=2, start=0,end=0.5)) + 
        theme(legend.position = "none")
  p1 <- p + stat_compare_means(aes(group=Remission.status),method="t.test",label="p.format",method.args = list(var.equal=F)) +
            scale_y_continuous(limits = c(-4,4), expand = expansion(mult = c(0.05,0.15)))
  return(p1)
})


# write out as pdf
filename1 = paste0(results_dir,"Figures\\ForPaper\\FigureS4a_metab_boxplots_by_weight.pdf")
pdf(filename1)
metabolite_plots_by_weight[1]
invisible(dev.off())
filename2 = paste0(results_dir,"Figures\\ForPaper\\FigureS4b_metab_boxplots_by_weight.pdf")
pdf(filename2)
metabolite_plots_by_weight[2]
invisible(dev.off())
filename3 = paste0(results_dir,"Figures\\ForPaper\\FigureS4c_metab_boxplots_by_weight.pdf")
pdf(filename3)
metabolite_plots_by_weight[3]
invisible(dev.off())
filename4 = paste0(results_dir,"Figures\\ForPaper\\FigureS4d_metab_boxplots_by_weight.pdf")
pdf(filename4)
metabolite_plots_by_weight[4]
invisible(dev.off())

# write out as postscript
filename1 = paste0(results_dir,"Figures\\ForPaper\\Figure2a_metab_boxplots_by_weight.eps")
postscript(filename1)
metabolite_plots_by_weight[1]
invisible(dev.off())
filename2 = paste0(results_dir,"Figures\\ForPaper\\Figure2b_metab_boxplots_by_weight.eps")
postscript(filename2)
metabolite_plots_by_weight[2]
invisible(dev.off())
filename3 = paste0(results_dir,"Figures\\ForPaper\\Figure2c_metab_boxplots_by_weight.eps")
postscript(filename3)
metabolite_plots_by_weight[3]
invisible(dev.off())
filename4 = paste0(results_dir,"Figures\\ForPaper\\Figure2d_metab_boxplots_by_weight.eps")
postscript(filename4)
metabolite_plots_by_weight[4]
invisible(dev.off())

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
