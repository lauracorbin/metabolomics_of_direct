# This script runs analyses relating to the EFFECT OF THE INTERVENTION ON METABOLITES
# both linear model and logistics model (on P/A) included

# last run: 9th Nov 2022

####################################################################################
####################################################################################
# set up
####################################################################################
####################################################################################

# capture output
sink(file="U:/Logs/05.0_intervention_effect.txt",split=TRUE)

# read in parameter file
pfile = read.table("U:/Scripts/05.0_param_file.R", sep = "=", as.is=T)

# set parameters
project <- as.character( pfile[1,2] )
processeddata_dir <- as.character( pfile[2,2] )
scripts_dir <- as.character( pfile[3,2] )
results_dir <- as.character( pfile[4,2] )
  
## call R script containing functions
source(paste0(scripts_dir,"05.0_functions.R"))

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in  datasets
# metabolite data in various forms 
dataset_list <- readRDS(file=paste0(processeddata_dir,"alternative_formats.rds"))
# clinic data
clinic_data_only <- readRDS(file=paste0(processeddata_dir,"clinic_data_input.rds"))
# summary stats
ds_results_list <- readRDS(file=paste0(results_dir,"metabolite_summaries.rds"))

## metabolite feature metadata (post cleaning)
clean_feature_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_feature_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_feature_metadata.rds"))) 

## add column for feature id matching
clean_feature_list[[1]]$id_for_matching <- clean_feature_list[[1]]$compid_to_match
clean_feature_list[[2]]$id_for_matching <- clean_feature_list[[2]]$shortname

## sample metadata (post cleaning)
clean_sample_list <- list(metab=readRDS(file=paste0(processeddata_dir,"OrigScale_cleaned_sample_metadata.rds")),nmr=readRDS(file=paste0(processeddata_dir,"nmr_cleaned_sample_metadata.rds")))

####################################################################################
####################################################################################
### run linear model 
####################################################################################
####################################################################################

# set up list for results output
lm_results_list <- list(rnt=NA,natlog=NA)
lm_interim_results_list <- list(metab=as.data.frame(matrix(data = NA, nrow = 1300, ncol = 17)),
                        nmr=as.data.frame(matrix(data = NA, nrow = 300, ncol = 17)))

# set up list for residuals output 
lm_residuals_list <- list(rnt=NA,natlog=NA)
lm_interim_residuals_list <- list(metab=as.data.frame(matrix(data=NA, nrow = nrow(dataset_list$raw$metab), ncol = nrow(clean_feature_list$metab)+1 )),
                                  nmr=as.data.frame(matrix(data=NA, nrow = nrow(dataset_list$raw$nmr), ncol = nrow(clean_feature_list$nmr)+1 )) )

# generate list of traditional blood measures to include in model
trad_bloods <- c("glucose_mmol_l","insulin_uu_ml","hdl_mmol_l","trig_mmol_l","chol_mmol_l")
trad_bloods.b <- paste0(trad_bloods,".b")

# implement linear model on rnt data and natlog data
datasets_to_run <- c(2,3)
for (j in 1:length(datasets_to_run)){
  # set dataset
  dataset_to_use <- dataset_list[[datasets_to_run[j]]]
  for (k in 1:length(dataset_to_use)){
    print(paste0("Running linear model on ", names(dataset_list)[datasets_to_run[j]] ," for: ",names(dataset_to_use)[k]))
    # set up pdf file for residual checks
    filename = paste0(results_dir, "Figures\\residual_plots_", names(dataset_list)[datasets_to_run[j]], " for ", names(dataset_to_use)[k], ".pdf")
    pdf(file = filename)
    par(mfrow=c(2,1))
    # define dataset - merge trial data with metabolite data
    dtst <- merge(clinic_data_only[[k]],dataset_to_use[[k]],by="row.names")
    rownames(dtst) <- dtst$Row.names
    dtst <- dtst[,-1]
    # define var list (ie list of metabolites)
    depvar <- c(clean_feature_list[[k]]$id_for_matching,trad_bloods)
    # add patient ids to residual file
    lm_interim_residuals_list[[k]][,1] <- dtst$id
    names(lm_interim_residuals_list[[k]])[1] <- "id"
    # set controls to be reference
    if (levels(dtst$treat)[1] == "Intervention") {
      dtst <- within(dtst,treat <- relevel(treat, ref = "Control"))
    }
    ## set up reporting var frequency
    s = base::seq(from=1,to=length(depvar), by=100)
    for (i in 1:length(depvar)) {
      if(i %in% s){
        print(paste0( "Now processing metabolite ", i, " of ", length(depvar) ))
      }
      mtb_name <- depvar[i]
      mtb_name_b <- paste0(mtb_name,".b")
      mtb_name_e <- paste0(mtb_name,".e")
      lm_interim_results_list[[k]][i,1] <- mtb_name
      mtb_e_col_ref <- which( colnames(dtst) == mtb_name_e)
      mtb_b_col_ref <- which( colnames(dtst) == mtb_name_b)
      # subset data to single feature (restrict to those with non-missing)
      mtb <- dtst[which(!is.na(dtst[mtb_e_col_ref]) ),
                  c("id","treat","site","centre","list.size","age","sex",mtb_name_b,mtb_name_e)]
      # subset data to those with non-missing basline
      mtb <- mtb[which(!is.na(mtb[8])),]
      # record no. of observations
      lm_interim_results_list[[k]][i,2] <- nrow(mtb)
      # count missing baseline
      lm_interim_results_list[[k]][i,5] <- length(which(is.na(mtb[,8])))
      # count per group
      lm_interim_results_list[[k]][i,4] <- nrow(mtb[mtb$treat == "Intervention",])
      lm_interim_results_list[[k]][i,3] <- nrow(mtb[mtb$treat == "Control",])
      # replace missing baseline values with mean
      #mtb[is.na(mtb[,8]),8] <- mean(mtb[,8],na.rm=T)
        # fit model
          fitA <- tryCatch(lm(mtb[,mtb_name_e] ~ age + sex + centre + list.size + mtb[,mtb_name_b] + treat , data=mtb), error = function(e) e )
          fitB <- tryCatch(lm(mtb[,mtb_name_e] ~ age + sex + centre + list.size + mtb[,mtb_name_b] , data=mtb), error = function(e) e )
          if (class(fitA)[1] == "lm") {
            # extract coefficients
            coef = summary(fitA)$coefficients
            # check for model fail because one of coef missing
            if ( (nrow(coef) == 7) ) {
              # run anova to calculate variance explained by each factor individually
              # alternative calc using Anova from car package (type 2) - need to restrict to cases where estimation possible
              if ( fitA$df.residual > 0 & (deviance(fitA) > sqrt(.Machine$double.eps)) ) {
                a = Anova(fitA)
                eta = a[,1]/sum(a[,1],na.rm=T)
              } else {
                eta = c(NA,NA,NA,NA,NA,NA,NA)
              }
              #a = anova(fitA)
              #eta = a[,2]/sum(a[,2],na.rm=T)
              # extract model results for treatment allocation
              treat_coef = coef[nrow(coef),]
              lm_interim_results_list[[k]][i,6:9] <- treat_coef
              # extract variance explained results
              lm_interim_results_list[[k]][i,10:16] <- eta
              # extract residuals (without intervention fitted)
              model_resid <- as.data.frame(resid(fitB))
              model_resid$id <- mtb$id
              lm_interim_residuals_list[[k]][,i+1] <- model_resid$`resid(fitB)`[match(lm_interim_residuals_list[[k]]$id,model_resid$id)]
              names(lm_interim_residuals_list[[k]])[i+1] <- mtb_name
              # check residuals for main model
              plot(fitted(fitA),residuals(fitA), main = mtb_name)
              qqnorm(residuals(fitA), main = mtb_name)
              qqline(residuals(fitA), col="red")
              abline(a=0,b=1)
            } else {
              lm_interim_results_list[[k]][i,17] <- "model failed to run - not all coefficients estimated"
              lm_interim_results_list[[k]][i,6:16] <- NA
            }
          } else {
            lm_interim_results_list[[k]][i,17] <- "model failed to run"
            lm_interim_results_list[[k]][i,6:16] <- NA
          }
    }
    # name column headers
    model_components <- c("age","sex","centre","list.size","mtb_at_baseline","treat","Residuals")
    names(lm_interim_results_list[[k]])[1] <- "feature_id"
    names(lm_interim_results_list[[k]])[2] <- "lm_n_samples_in_model"
    names(lm_interim_results_list[[k]])[3] <- "lm_n_control"
    names(lm_interim_results_list[[k]])[4] <- "lm_n_intervention"
    names(lm_interim_results_list[[k]])[5] <- "lm_n_imputed_baselines"
    names(lm_interim_results_list[[k]])[6] <- "lm_treat_beta"
    names(lm_interim_results_list[[k]])[7] <- "lm_treat_SE"
    names(lm_interim_results_list[[k]])[8] <- "lm_treat_t"
    names(lm_interim_results_list[[k]])[9] <- "lm_treat_p"
    names(lm_interim_results_list[[k]])[10:16] <- paste0("lm_eta2_", model_components)
    names(lm_interim_results_list[[k]])[17] <- "lm_model_fail_info"
    names(lm_interim_results_list[[k]])[2:17] <- paste0(names(dataset_list)[datasets_to_run[j]], "_", names(lm_interim_results_list[[k]])[2:17])
    # close pdf
    invisible(dev.off())
  }
  lm_results_list[[j]] <- lm_interim_results_list
  lm_residuals_list[[j]] <- lm_interim_residuals_list
}

# write out residuals
write.table(lm_residuals_list$rnt$metab,file=paste0(processeddata_dir,"metabolon_lm_rnt_resid.txt"),row.names=F,quote=F,sep="\t")
saveRDS(lm_residuals_list$rnt$metab,file=paste0(processeddata_dir,"metabolon_lm_rnt_resid.rds"))

write.table(lm_residuals_list$rnt$nmr,file=paste0(processeddata_dir,"nmr_lm_rnt_resid.txt"),row.names=F,quote=F,sep="\t")
saveRDS(lm_residuals_list$rnt$nmr,file=paste0(processeddata_dir,"nmr_lm_rnt_resid.rds"))

####################################################################################
####################################################################################
### run logistic model on presence/absence
####################################################################################
####################################################################################

# set up list for results output
logit_results_list <- list(pa=NA)
logit_interim_results_list <- list(metab=as.data.frame(matrix(data = NA, nrow = 1300, ncol = 11)),
                                nmr=as.data.frame(matrix(data = NA, nrow = 300, ncol = 11)))


# implement logistic model on rnt data 
datasets_to_run <- c(4)
for (j in 1:length(datasets_to_run)){
  # set dataset
  dataset_to_use <- dataset_list[[datasets_to_run[j]]]
  for (k in 1:length(dataset_to_use)){
    print(paste0("Running logistic model on ", names(dataset_list)[datasets_to_run[j]] ," for: ",names(dataset_to_use)[k]))
    # define dataset - merge trial data with metabolite data
    dtst <- merge(clinic_data_only[[k]],dataset_to_use[[k]],by="row.names")
    rownames(dtst) <- dtst$Row.names
    dtst <- dtst[,-1]
    # define var list (ie list of metabolites)
    depvar <- c(clean_feature_list[[k]]$id_for_matching,"glucose_mmol_l","insulin_uu_ml","hdl_mmol_l","trig_mmol_l","chol_mmol_l")
    ## set controls to be reference
    if (levels(dtst$treat)[1] == "Intervention") {
      dtst <- within(dtst,treat <- relevel(treat, ref = "Control"))
    }
    ## set up reporting var frequency
    s = base::seq(from=1,to=length(depvar), by=100)
    for (i in 1:length(depvar)) {
      if(i %in% s){
        print(paste0( "Now processing metabolite ", i, " of ", length(depvar) ))
      }
      mtb_name <- depvar[i]
      mtb_name_b <- paste0(mtb_name,".b")
      mtb_name_e <- paste0(mtb_name,".e")
      logit_interim_results_list[[k]][i,1] <- mtb_name
      mtb_e_col_ref <- which( colnames(dtst) == mtb_name_e)
      mtb_b_col_ref <- which( colnames(dtst) == mtb_name_b)
      # subset data to single feature
      mtb <- dtst[!is.na(dtst[mtb_e_col_ref]),
                  c("id","age","sex","treat",mtb_name_b,mtb_name_e)]
      # record no. of observations
      logit_interim_results_list[[k]][i,2] <- nrow(mtb)
      # count present by group and timepoint (i.e. equal to 1)
      logit_interim_results_list[[k]][i,3] <- sum(mtb[mtb$treat == "Control",5],na.rm=T) # baseline control
      logit_interim_results_list[[k]][i,4] <- sum(mtb[mtb$treat == "Control",6],na.rm=T) # endpoint control
      logit_interim_results_list[[k]][i,5] <- sum(mtb[mtb$treat == "Intervention",5],na.rm=T) # baseline intervention
      logit_interim_results_list[[k]][i,6] <- sum(mtb[mtb$treat == "Intervention",6],na.rm=T) # endpoint intervention
        # run model
        fitA <- try(glm(mtb[,mtb_name_e] ~ mtb[,mtb_name_b] + treat, data=mtb,family="binomial"))
        # check for model fail
        test <- !is.null(attr(fitA,"try-error"))
        if (test | !fitA$converged) {
          logit_interim_results_list[[k]][i,11] <- "model failed to run"
          logit_interim_results_list[[k]][i,7:10] <- NA
        } else {
          # extract results
          coef = summary(fitA)$coefficients
          treat_coef = coef[nrow(coef) ,]
          logit_interim_results_list[[k]][i,7:10] <- treat_coef
        }
    }
    # name column headers
    names(logit_interim_results_list[[k]])[1] <- "feature_id"
    names(logit_interim_results_list[[k]])[2] <- "logit_n_samples_in_model"
    names(logit_interim_results_list[[k]])[3] <- "logit_n_control_baseline"
    names(logit_interim_results_list[[k]])[4] <- "logit_n_control_end"
    names(logit_interim_results_list[[k]])[5] <- "logit_n_intervention_baseline"    
    names(logit_interim_results_list[[k]])[6] <- "logit_n_intervention_end"    
    names(logit_interim_results_list[[k]])[7] <- "logit_treat_beta"
    names(logit_interim_results_list[[k]])[8] <- "logit_treat_SE"
    names(logit_interim_results_list[[k]])[9] <- "logit_treat_z"
    names(logit_interim_results_list[[k]])[10] <- "logit_treat_p"
    names(logit_interim_results_list[[k]])[11] <- "model_fail_info"
    names(logit_interim_results_list[[k]])[2:11] <- paste0(names(dataset_list)[datasets_to_run[j]], "_", names(logit_interim_results_list[[k]])[2:11])
  }
  logit_results_list[[j]] <- logit_interim_results_list
}

####################################################################################
####################################################################################
### run fold change on raw, rnt and winsorized data 
####################################################################################
####################################################################################

# set up list for results output
fc_results_list <- list(raw=NA,win=NA)
fc_interim_results_list <- list(metab=as.data.frame(matrix(data = NA, nrow = 1300, ncol = 18)),
                                nmr=as.data.frame(matrix(data = NA, nrow = 300, ncol = 18)))

# implement fold change on raw, rnt and winsorized data
datasets_to_run <- c(1,5)
for (j in 1:length(datasets_to_run)){
  # set dataset
  dataset_to_use <- dataset_list[[datasets_to_run[j]]]
  for (k in 1:length(dataset_to_use)){
    print(paste0("Running fold change on ", names(dataset_list)[datasets_to_run[j]] ," for: ",names(dataset_to_use)[k]))
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
      mtb_e_col_ref <- which( colnames(dtst) == mtb_name_e)
      # subset data to single feature
      mtb <- dtst[!is.na(dtst[mtb_e_col_ref]), c("id","treat",mtb_name_e)]
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
    
    names(fc_interim_results_list[[k]])[2:18] <- paste0(names(dataset_list)[datasets_to_run[j]], "_", names(fc_interim_results_list[[k]])[2:18])
  }
  fc_results_list[[j]] <- fc_interim_results_list
}
    

####################################################################################
####################################################################################
### merge results and save out
####################################################################################
####################################################################################

# merge metab
lm_metab <- merge(lm_results_list$rnt$metab[1:1259,],lm_results_list$natlog$metab[1:1259,],by="feature_id")
fc_metab <- merge(fc_results_list$raw$metab[1:1259,],fc_results_list$win$metab[1:1259,],by="feature_id")
metab_results1 <- merge(clean_feature_list$metab,ds_results_list$metab,by.x="compid_to_match",by.y="feature_id",all=T)
names(metab_results1)[1] <- "feature_id"
metab_results2 <- merge(metab_results1,lm_metab,by="feature_id",all=T) 
metab_results3 <- merge(metab_results2,logit_results_list$pa$metab[1:1259,],by="feature_id",all=T)
metab_results  <- merge(metab_results3,fc_metab,by="feature_id")

# merge nmr
lm_nmr <- merge(lm_results_list$rnt$nmr[1:230,],lm_results_list$natlog$nmr[1:230,],by="feature_id")
fc_nmr <- merge(fc_results_list$raw$nmr[1:230,],fc_results_list$win$nmr[1:230,],by="feature_id")
nmr_results1 <- merge(clean_feature_list$nmr,ds_results_list$nmr,by.x="shortname",by.y="feature_id",all=T)
names(nmr_results1)[1] <- "feature_id"
nmr_results2 <- merge(nmr_results1,lm_nmr,by="feature_id",all=T) 
nmr_results3 <- merge(nmr_results2,logit_results_list$pa$nmr[1:230,],by="feature_id",all=T)
nmr_results  <- merge(nmr_results3,fc_nmr,by="feature_id")

# write out 
write.table(metab_results,file=paste0(results_dir,"metabolon_results.txt"),row.names=F,quote=F,sep="\t")
saveRDS(metab_results,file=paste0(results_dir,"metabolon_results.rds"))
write.table(nmr_results,file=paste0(results_dir,"nmr_results.txt"),row.names=F,quote=F,sep="\t")
saveRDS(nmr_results,file=paste0(results_dir,"nmr_results.rds"))

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()

