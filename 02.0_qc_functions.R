## functions for metabolomics QC ##

## load libraries
library("dplyr")
library("knitr")
library("pwr")
library("data.table")
library("ggplot2")
library("psych")
library('pcaMethods')


###########
## POWER ##
###########

eval.power.binary = function(N, effect, alpha) {
  halfN = ceiling(N/2)
  N_case = N_control = halfN
  power_calc = pwr.t.test(n = halfN, d = effect, sig.level = alpha, power = NULL, type="two.sample")
  power <- round(power_calc$power, d=3)
  tmp <- cbind(N, N_case, N_control, effect, alpha, power)
  return(tmp)
}

eval.power.cont = function(N, n_coeff, effect, alpha) {
  power_calc = pwr.f2.test(u = n_coeff, v = N-n_coeff-1, f2 = effect, sig.level = alpha, power = NULL)
  power <- round(power_calc$power, d=3)
  tmp <- cbind(N, effect, alpha, power, n_coeff)
  return(tmp)
}

#################
## MISSINGNESS ##
#################

################### functions to calculate missingness across samples and features ################### 
## take as input a data.frame (dtst) containing metabolite data with samples in columns and metabolite features in rows (ids as column and row names)

## BY SAMPLE
sample.missing = function(dtst) {
  out = colMeans(is.na(dtst))
  return(out)
}

## BY FEATURE
feature.missing = function(dtst) {
  out = rowMeans(is.na(dtst))
  return(out)
}

## function to identify features with lowest/highest missingness
least.missing = function(feature_missing,number){
  d_sorted <- sort(feature_missing,decreasing = FALSE)
  least_missing <- names(d_sorted[1:number])
  return(least_missing)
}

most.missing = function(feature_missing,number){
  d_sorted <- sort(feature_missing,decreasing = TRUE)
  most_missing <- names(d_sorted[1:number])
  return(most_missing)
}

################### functions to check total peak area across samples ################### 

## calculate total peak area (sum across all features)
## takes as input a data.frame (dtst) containing metabolite data with samples in columns and metabolite features in rows (ids as column and row names)
sample.tpa = function(dtst) {
  out = colSums(dtst,na.rm=T)
  return(out)
}

## identify outliers based on stdev from mean
## takes as input the output from sample.tpa (x)
outliers = function(x, nsd = 3){
  # x is a vector of values
  x = as.numeric(x)
  m = mean(na.omit(x))
  s = sd(na.omit(x))
  bottom = m-(nsd*s)
  top = m+(nsd*s)
  cuttoffs = c(top, bottom)
  out = c(which(x > cuttoffs[1]), which(x < cuttoffs[2]))
  output = list(out,bottom,top)
  return(output)
}

################### functions to identify samples and features with high missingness (for exclusion) ################### 
## BY SAMPLE
sample.missing.exclusions = function (dtst, threshold){
  # calculate by sample missingness
  missing = colMeans(is.na(dtst))
  # add as row in data.frame
  new = rbind(dtst,missing)
  # identify those that are greater than the threshold
  exclude = which(new[nrow(new),] > threshold)
  exclude_id = colnames(dtst)[exclude] 
  return(exclude_id)
}

## BY FEATURE (including group)
feature.missing.exclusions.bygroup = function (dtst, threshold, bygroup, phenofile){
  if (bygroup == 0) {
    # calculate by feature missingness in entire dataset
    missing = rowMeans(is.na(dtst))
    # add as row in data.frame
    new = cbind(dtst,missing)
    # identify those that are greater than the threshold
    exclude = which(new[,ncol(new)] > threshold)
    exclude_id = rownames(dtst)[exclude] 
    return(exclude_id)
  } else {
    # identify groups
    phenofile$numeric_group <- as.numeric(phenofile$group)
    n_group <- length(unique(phenofile$numeric_group))
    # make matrix to hold results from loop
    #mat <- matrix(, nrow = nrow(dtst), ncol = n_group)
    total_exclusions <- NULL
    # calculate by feature missingness in group
    for (i in 1:n_group) {
      data_subset_ids <- as.list(subset(phenofile, numeric_group == i,sample_id))
      #print(data_subset_ids)
      data_subset <- dtst[, names(dtst) %in% data_subset_ids$sample_id]
      print(paste0("Number in group ", i ))
      print(ncol(data_subset))
      missing = rowMeans(is.na(data_subset))
      # add as row in data.frame
      new = cbind(data_subset,missing)
      # identify those that are greater than the threshold
      exclude = which(new[,ncol(new)] > threshold)
      exclude_id_temp = rownames(new)[exclude] 
      total_exclusions <- c(total_exclusions,exclude_id_temp)
    }
    ## bring together to work out exclusions - only those that fall below threshold in all groups
    exclusions_tab <- as.data.frame(table(total_exclusions))
    exclude_id <- as.character(exclusions_tab[which(exclusions_tab$Freq == n_group),c("total_exclusions")])
    return(exclude_id)
  }
}

## BY FEATURE (not including group)
feature.missing.exclusions = function (dtst, threshold){
    # calculate by feature missingness in entire dataset
    missing = rowMeans(is.na(dtst))
    # add as row in data.frame
    new = cbind(dtst,missing)
    # identify those that are greater than the threshold
    exclude = which(new[,ncol(new)] > threshold)
    exclude_id = rownames(dtst)[exclude] 
}

################### functions to remove samples and features with high missingness ################### 
## BY SAMPLE
sample.exclude = function(dtst, exclusions){
  new = dtst[, !names(dtst) %in% exclusions]
  return(new)
}

## BY FEATURE
feature.exclude = function(dtst, exclusions){
  new = dtst[!rownames(dtst) %in% exclusions,]
  return(new)
}

######################################
## SUMMARISE OUTLIERS AND WINSORISE ##
######################################

outlier.summary = function(dtst, scaling_flag=0, pdf_filename, nsd=5, cut=0.01, method="sd"){
  
  ## transform dataset
  data_t <- as.data.frame(t(dtst))
  
  ## set up results tables
  outlier_results = as.data.frame(matrix(data = 0, nrow = nrow(data_t), ncol = ncol(data_t)))
  winsorised_data = data_t
  
  # rescale if not imputed
  if (scaling_flag == 1) {
    data_t <- data_t/1000000
  } else {
    data_t <- data_t
  }
  
  # save .pdf with summaries
  pdf(pdf_filename)
  par(mfrow=c(3,2))
  
  ## loop over features
  for (i in 1:ncol(data_t)){
    # continue only if more than one observation
    if (length(which(!is.na(data_t[,i]))) > 1) {
      mtb_name <- names(data_t[i])
      # id outliers
      # by sd
      if (method == "sd") {
        threshold_plus <- mean(data_t[,i],na.rm=T) + (nsd*(sd(data_t[,i],na.rm=T)))
        threshold_minus <- mean(data_t[,i],na.rm=T) - (nsd*(sd(data_t[,i],na.rm=T)))
      }
      #by percentile
      if (method == "percentile") {
        threshold_plus <- quantile(data_t[,i],1-cut,na.rm=T)
        threshold_minus <- quantile(data_t[,i],cut,na.rm=T)
      }
      upperoutliers <- which(data_t[,i] >= threshold_plus)
      loweroutliers <- which(data_t[,i] <= threshold_minus)
      outliers <- c(upperoutliers,loweroutliers)
      # Plot concentration/proportion of metabolite vs observation number 
      y_lim <- max(c(max(data_t[,i],na.rm=T),threshold_plus),na.rm=T)
      plot(seq(1,nrow(data_t),1),(data_t[,i]),ylim=c(0,y_lim),main=(names(data_t)[i]),xlab="sample index", ylab="peak area")
      abline(h=threshold_plus,col="red")   
      text(nrow(data_t)+50,threshold_plus, paste0("upper_threshold"), col="red",cex=0.6, pos = 4,xpd = NA)
      if (threshold_minus > 0) {
        abline(h=threshold_minus,col="red")   
        text(nrow(data_t)+50,threshold_minus, paste0("lower threshold"), col="red",cex=0.6, pos = 4,xpd = NA)
      }
      # Plot histogram with distribution statistics
      if (scaling_flag == 1) {
        summary_stats <- describe(data_t[,i]*1000000)
      } else {
        summary_stats <- describe(data_t[,i])
      }
      #meanvar<-mean(data_t[,i],na.rm=TRUE)
      #medianvar<-median(data_t[,i], na.rm = TRUE)  
      #minvar<-min(data_t[,i], na.rm = TRUE)
      #maxvar<-max(data_t[,i], na.rm = TRUE)
      #kurtosisvar<-kurtosis(data_t[,i])
      #skewnessvar<-skewness(data_t[,i])
      #N<- nrow(data_t) - sum(is.na(data_t[,i]))
      #missingness <- (sum(is.na(data_t[,i]))/nrow(data_t))*100
      
      a<-density(data_t[,i], na.rm=T)
      #thresholdx<-(maxvar+(maxvar/100))
      #thresholdy<-min(a$y)+(max(a$y)/4)
      
      hist(data_t[,i], col="red",main="",prob=TRUE,xlab="peak area") 
      lines(density(data_t[,i], na.rm = TRUE),col="blue", lwd=2)
      #text(thresholdx,thresholdy, cex=0.6, paste("N=", N, "\npercent missing=", signif(missingness, 3), "\nmin=", signif(minvar, 3), " \nmax=",signif(maxvar, 3), 
       #                                          "\nmean=", signif(meanvar, 3), " \nmedian=", signif(medianvar, 3), 
        #                                         sep = ''), pos = 3,xpd = NA)
      
                                                 #"\nkurt=", signif(kurtosisvar, 3), " \nskew=", signif(skewnessvar, 3), sep = ''), pos = 3,xpd = NA)
      text(cex=0.6,paste(summary_stats))
      # output flag info
      outlier_results[c(outliers),i] <- 1
      # replace outliers with limit
      winsorised_data[c(upperoutliers),i] <- threshold_plus*1000000
      winsorised_data[c(loweroutliers),i] <- threshold_minus*1000000
    }
    else {
      outlier_results[,i] <- NA
    }
  }
  dev.off()
  names(outlier_results) <- names(data_t)
  row.names(outlier_results) <- row.names(data_t)
  
  names(winsorised_data) <- names(data_t)
  row.names(winsorised_data) <- row.names(data_t)
  
  return(list(outlier_results,winsorised_data))
}


######################################
## CHECK NORMALITY OF DISTRIBUTIONS ##
######################################

## function to apply shapiro test to each row (feature) of the data
normality.test = function(x,max_missing, min_unique){
  if ( sum(is.na(x)) < max_missing & length(unique(x)) > min_unique) {
    test_result = shapiro.test(x)  
    out = c(test_result$statistic, test_result$p.value)
    } else {
    out = c(NA,NA)
    }
}


## function to calculate proportion normal based on w stat or p-value
proportion.normal.w = function(normality_test_result,threshold){
  sum(normality_test_result[,1]>threshold, na.rm=T)/nrow(normality_test_result)
}
proportion.normal.p = function(normality_test_result,threshold){
  sum(normality_test_result[,2]>threshold, na.rm=T)/nrow(normality_test_result)
}


## function to identify least normal distributions
least.normal = function(normality_test,number){
  d_sorted <- normality_test[order(normality_test$shapiro_w,decreasing = FALSE),]
  least_normal <- row.names(d_sorted[1:number,])
  return(least_normal)
}

most.normal = function(normality_test,number){
  d_sorted <- normality_test[order(normality_test$shapiro_w,decreasing = TRUE),]
  most_normal <- row.names(d_sorted[1:number,])
  return(most_normal)
}

###################################
## CHECK PCA - LOOK FOR OUTLIERS ##
###################################
## make tree to identify independent features
## takes as input a data.frame (dtst) containing metabolite data with samples in columns and metabolite features in rows (ids as column and row names)
make.tree = function(dtst, cor_method, hclust_method){
  cor_matrix <- cor(t(dtst),method=cor_method, use="pairwise.complete.obs")
  dist_matrix <- as.dist(1-abs(cor_matrix))
  SpearTree <- hclust(dist_matrix, method = hclust_method)
  return(SpearTree)
}


## identify PC outliers
## input is a pc vector
pc.outliers = function(x, nsd = 3){
  # x is a vector of values
  x = as.numeric(x)
  m = mean(na.omit(x))
  s = sd(na.omit(x))
  bottom = m-(nsd*s)
  top = m+(nsd*s)
  cuttoffs = c(top, bottom)
  out = c(which(x > cuttoffs[1]), which(x < cuttoffs[2]))
  output = list(out,bottom,top)
  return(output)
}



