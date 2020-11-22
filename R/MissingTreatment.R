#' Multivariate Imputation by Chained Equations
#'
#' The mice package implements a method to deal with missing data. The package creates multiple imputations (replacement values) for multivariate missing data. The method is based on Fully Conditional Specification, where each incomplete variable is imputed by a separate model. The MICE algorithm can impute mixes of continuous, binary, unordered categorical and ordered categorical data. In addition, MICE can impute continuous two-level data, and maintain consistency between imputations by means of passive imputation. Many diagnostic plots are implemented to inspect the quality of the imputations. Generates Multivariate Imputations by Chained Equations (MICE)
#' @param data Any data frame tat needs to be treated for Missing data
#' @param dv The dependant variable that is to be ignored from the given data while treating missing values
#' @param method Can be either a single string, or a vector of strings with length length(blocks), specifying the imputation method to be used for each column in data. If specified as a single string, the same method will be used for all blocks.
#' @return Returns a dataframe treated for missing values
#' @export
#' @usage mice_treatment(data, dv, method)

mice_treatment <- function(data, dv, method)
{
  # require(mice)
  print("Mice Imputation Running.. It runs till 5 5.")
  miceMod <- mice(data[, !names(data) %in% dv], method=method)  # perform mice imputation, based on random forests.
  miceOutput <- complete(miceMod)  # generate the completed data.
  return (miceOutput)
}

#' Treats missing values in a given dataset using the best possible technique
#'
#' Takes in a data frame and performs the best possible missing value treatment to each of the columns in the data frame
#' @param data Any data frame that has to treated for missing data
#' @param dv Dependent variable in the given dataset in order to eliminate that field while treating missing data
#' @description missing_treatments initially eliminates all the columns having more than 70% missing data. Different missing value treatment techiniques considered  are as follows:
#' \itemize{
#' \item Zero Imputation
#' \item Min Imputation
#' \item Max Imputation
#' \item Mean Imputation
#' \item Median Imputation
#' \item Mode Imputation
#' \item Multiple Imputation - Random Forest
#' \item Multiple Imputation - Predictive Mean Matching
#' }
#' Of all these techniques, Best technique is chosed based on the performance of different datasets on a linear/logistic model.
#' @return
#' Returns a list of 3 objects:
#' \describe{
#' \item{missing_treated}{Missing value treated Dataset}
#' \item{fit_model}{A list of model fit files for each of the techniques used}
#' \item{model_perf_metrics}{Performance metrics for eac of the intermediate datasets created}
#' }
#'
#' @export
#' @usage missing_treatments(data, dv)
missing_treatments <- function(data, dv)
{
  # library(BestTransform)
  # library(data.table)

  dvcol <- data[, dv]

  per_missing <- sapply(data, function(x) {round(sum(is.na(x))/length(x)*100, 0)})
  del_cols <- names(per_missing[per_missing>=70])
  data <- data[, !colnames(data) %in% del_cols]

  dist <- data_distribution(data, dv)

  cont_cols <- dist[(dist$distribution == "Continous") & (dist$is_dv == FALSE),]$names
  cat_cols <- dist[(dist$distribution != "Continous") & (dist$is_dv == FALSE),]$names

  cat_data <- data[, cat_cols]
  for (i in cat_cols)
  {
    cat_data[, i] <- as.factor(cat_data[, i])
    zz <- factor(cat_data[,i], levels=c(levels(cat_data[, i]), "Missing_OFE"))
    zz[is.na(zz)] <- "Missing_OFE"
    cat_data[, i] <- zz
  }

  zero_impute <- data.table(data)[,lapply(.SD, function(x) {x[is.na(x)] <- 0; x}), .SDcols=cont_cols]
  zero_impute <- cbind(zero_impute, cat_data, dvcol)

  min_impute <- data.table(data)[,lapply(.SD, function(x) {x[is.na(x)] <- min(x, na.rm = T); x}), .SDcols=cont_cols]
  min_impute <- cbind(min_impute, cat_data, dvcol)
  min_fit <- data.frame(t(data.table(data)[,lapply(.SD, function(x) {min(x, na.rm = T)}), .SDcols=cont_cols]))
  min_fit$columns <- rownames(min_fit)
  colnames(min_fit)[1] <- "min_value"

  max_impute <- data.table(data)[,lapply(.SD, function(x) {x[is.na(x)] <- max(x, na.rm = T); x}), .SDcols=cont_cols]
  max_impute <- cbind(max_impute, cat_data, dvcol)
  max_fit <- data.frame(t(data.table(data)[,lapply(.SD, function(x) {max(x, na.rm = T)}), .SDcols=cont_cols]))
  max_fit$columns <- rownames(max_fit)
  colnames(max_fit)[1] <- "max_value"

  # library(e1071)
  mean_impute <- data.frame(round(impute(data[, cont_cols], what = "mean"),0))
  mean_impute <- cbind(mean_impute, cat_data, dvcol)
  mean_fit <- data.frame(t(data.table(data)[,lapply(.SD, function(x) {round(mean(x, na.rm = T), 0)}), .SDcols=cont_cols]))
  mean_fit$columns <- rownames(mean_fit)
  colnames(mean_fit)[1] <- "mean_value"


  median_impute <- data.frame(round(impute(data[, cont_cols], what = "median"),0))
  median_impute <- cbind(median_impute, cat_data, dvcol)
  median_fit <- data.frame(t(data.table(data)[,lapply(.SD, function(x) {round(median(x, na.rm = T), 0)}), .SDcols=cont_cols]))
  median_fit$columns <- rownames(median_fit)
  colnames(median_fit)[1] <- "median_value"

  mode_impute <- data.table(data)[,lapply(.SD, function(x) {x[is.na(x)] <- names(table(x))[table(x)==max(table(x))][[1]]; x}), .SDcols=c(cont_cols, cat_cols)]
  mode_impute <- cbind(mode_impute, dvcol)
  mode_fit <- data.frame(t(data.table(data)[,lapply(.SD, function(x) {names(table(x))[table(x)==max(table(x))][[1]]}), .SDcols=cont_cols]))
  mode_fit$columns <- rownames(mode_fit)
  colnames(mode_fit)[1] <- "mode_value"

  miceOutput_rf <- mice_treatment(data = data, dv, "rf")
  miceOutput_rf <- cbind(miceOutput_rf,dvcol)

  miceOutput_pmm <- mice_treatment(data = data, dv, "pmm")
  miceOutput_pmm <- cbind(miceOutput_pmm,dvcol)

  total_matrix_all <- list("Zero Imputation" = zero_impute, "Min Imputation" = min_impute, "Max Imputation" = max_impute, "Mean Imputation" = mean_impute, "Median Imputation" = median_impute, "Mode Imputation" = mode_impute, "Multiple Imputation - Random Forest" = miceOutput_rf, "Multiple Imputation - Predictive Mean Matching" = miceOutput_pmm)
  # library(pbmcapply)
  # library(pbmcapply)
  output<-do.call(rbind,pbmclapply(seq(1:length(total_matrix_all)), model_my_data,
                                   data = total_matrix_all,
                                   mc.cores = 1))
  if (dist[dist$is_dv == T, ]$distribution == "Continous")
  {
    output <- data.frame(output[, c("Rsquared", "MAE", "RMSE")])
    var <- c("Rsquared", "MAE", "RMSE")
  } else
  {
    output <- data.frame(output[, c("Mean_F1", "Mean_Precision", "Mean_Recall")])
    var <- c("Mean_F1", "Mean_Precision", "Mean_Recall")
  }

  # var<-nearZeroVar(output, saveMetrics = TRUE)
  # var<-rownames(var[(var$zeroVar=="FALSE"),])
  # out<-data.frame(output)[,var]
  output$Method <- names(total_matrix_all)
  best_missing_treat <- names(total_matrix_all[which(output[,1] == max(output[,1]))[1]])

  missing_treated <- total_matrix_all[[best_missing_treat]]
  # library(Hmisc)
  # form_in <- formula(paste('~', paste(c(cont_cols, cat_cols), collapse = '+')))
  # hmisc_impute <- aregImpute(form_in, data = data[, c(cont_cols, cat_cols)], n.impute = 5)

  output <- output[,c("Method", var)]

  possible_fits <- list("Min Imputation" = min_fit, "Max Imputation" = max_fit, "Mean Imputation" = mean_fit, "Median Imputation" = median_fit, "Mode Imputation" = mode_fit)

  return(list(missing_treated = missing_treated, fit_model = possible_fits, model_perf_metrics = output))


}
