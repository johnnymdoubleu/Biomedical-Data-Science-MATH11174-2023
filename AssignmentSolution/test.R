library(data.table)
library(tidyverse)


# Loading in data
gdm.dt <- fread("data_assignment2/GDM.raw.txt")
# Viewing top of data
#head(gdm.dt)
# Checking for NAs in data
table(is.na(gdm.dt)) # 649 NAs in the dataset\

# Check if any of these are in ID, Sex or Pheno
table(is.na(gdm.dt[,c("ID","sex","pheno")])) # 0 NAs of these

# Function to impute to median
median_imputation <- function(x) {
    # Only apply function to numeric or integer columns
    if (is.numeric(x) || is.integer(x)){
    # Find index of missing values
    na.idx <- is.na(x)
    # Replace missing values with median of the observed values
    x[na.idx] <- median(x, na.rm=TRUE)
    }
return(x)}
# Impute all missing values in SNP cols to median (no NAs in ID, Sex or Pheno cols)
gdm.dt.imputed <- gdm.dt %>% copy() %>% .[, lapply(.SD, median_imputation)]
# Check no NAs exist
table(is.na(gdm.dt.imputed)) # 0 NAs in the dataset

univ.glm.test <- function(x, y, order=FALSE){
    column_names <- c("SNPName","InterceptCoef",
    "SNPCoef","OddsRatio",
    "StandardError","pvalue")
    model_info <- matrix(ncol=6,nrow=0, dimnames = list(c(),column_names))
    # First convert data table to dataframe
    x_df <- data.frame(x)
    # Full list of column names
    col_names <- colnames(x)
    for(col in col_names){
    # Logistic Model
    log_model <- glm(y ~ x_df[,col], family = binomial(link="logit"))
    # Storing summary of model
    model_summary <- summary(log_model)
    # Extracting information from model summary
    reg_intercept <- model_summary$coefficients["(Intercept)","Estimate"]
    reg_coef <- model_summary$coefficients[2,"Estimate"]
    odd_ratios <- exp(reg_coef)
    reg_se <- model_summary$coefficients[2,"Std. Error"]
    reg_pvalue <- model_summary$coefficients[2,"Pr(>|z|)"]
    # Appending information in model_info list
    model_info <- rbind(model_info,
    c(col, reg_intercept, reg_coef, odd_ratios,
    reg_se, reg_pvalue))
    }
    # Converting to data.table
    model_info_dt <- as.data.table(model_info)
    # Converting all apart from the first column of model_info_dt to numeric
    numeric_cols = column_names[-1]
    model_info_dt[,(numeric_cols):= lapply(.SD, as.numeric),
    .SDcols = numeric_cols]
    # If order set to true, reorder data.table by p-value (increasing)
    if(order){
    model_info_dt <- setorderv(model_info_dt, c("pvalue"), c(1))
    }
    return(model_info_dt)
}


gdm_imp_pheno <- gdm.dt.imputed$pheno
# Want to extract SNP names
SNP_names <- colnames(gdm.dt.imputed[,-c("ID","sex","pheno")])
# Extracting all SNPs
gdm_imp_SNPs <- gdm.dt.imputed[,..SNP_names]
# Running a binomial GLM regression for each SNP
# storing results and ordering by most significant
SNP_model_results <- univ.glm.test(gdm_imp_SNPs, gdm_imp_pheno, order=TRUE)
# Viewing most significant SNPs
head(SNP_model_results)

# Defining threshold as 1e-4 for most strong associated SNPs
most_associated <- SNP_model_results[pvalue <= 1e-04,]
# most strongly associated with increased risk of gestational diabetes
most_increased_SNP <- most_associated[which.max(most_associated$SNPCoef),]
# most strongly associated with reduced risk of gestational diabetes
most_reduced_SNP <- most_associated[which.min(most_associated$SNPCoef),]
cat("Most strongly associated with increased risk of gestational diabetes is",
most_increased_SNP$SNPName,"\n")

cat("Most strongly associated with reduced risk of gestational diabetes is",
most_reduced_SNP$SNPName)

# GLM model for SNP rs12243326_A
rs12243326_A_model <- glm(gdm_imp_pheno ~ gdm_imp_SNPs[,c(rs12243326_A)],
family = binomial(link="logit"))
# Summary of rs12243326_A model
summary(rs12243326_A_model)

# GLM model for SNP rs2237897_T
rs2237897_T_model <- glm(gdm_imp_pheno ~ gdm_imp_SNPs[,c(rs2237897_T)],
family = binomial(link="logit"))
# Summary of rs2237897_T model
summary(rs2237897_T_model)

# Confidence intervals for rs12243326_A
# 99% CI (converting to odds-ratios)
rs12243326_A_99CI <- exp(confint(rs12243326_A_model, level = 0.99))
rownames(rs12243326_A_99CI) <- c("Intercept","SNP Coefficient")
# 95% CI
rs12243326_A_95CI <- exp(confint(rs12243326_A_model, level = 0.95))
rownames(rs12243326_A_95CI) <- c("Intercept","SNP Coefficient")
knitr::kable(rs12243326_A_99CI, escape = FALSE, digits = 3,
caption = "99% Confidence Intervals for SNP rs12243326_A")

knitr::kable(rs12243326_A_95CI, escape = FALSE, digits = 3,
caption = "95% Confidence Intervals for SNP rs12243326_A")

# Confidence intervals for rs2237897_T
# 99% CI (converting to odds-ratios)
rs2237897_T_99CI <- exp(confint(rs2237897_T_model, level = 0.99))
rownames(rs2237897_T_99CI) <- c("Intercept","SNP Coefficient")
# 95% CI
rs2237897_T_95CI <- exp(confint(rs2237897_T_model, level = 0.95))
rownames(rs2237897_T_95CI) <- c("Intercept","SNP Coefficient")
knitr::kable(rs2237897_T_99CI, escape = FALSE, digits = 3,
caption = "99% Confidence Intervals for SNP rs2237897_T")

knitr::kable(rs2237897_T_95CI, escape = FALSE, digits = 3,
caption = "95% Confidence Intervals for SNP rs2237897_T")

# Importing table of gene names
GDM_annot <- fread("data_assignment2/GDM.annot.txt")
# Splitting SNP name and effect allele for SNP model results
SNP_model_results[, c("SNP", "allele") := tstrsplit(SNPName, "_", fixed=TRUE)]
# Subsetting SNP_model_results table to only include SNPs
# that have p-value of < 1e-0.4
SNP_hit <- SNP_model_results[pvalue < 1e-04,]
# Joining model results table with gene names table
SNP_w_genes <- GDM_annot[SNP_hit, on = .(snp = SNP)]
# Only reporting SNP name, effect allele, chromosome number and corresponding
# gene name
SNP_report <- SNP_w_genes[,c("snp","allele","chrom","gene")]
knitr::kable(SNP_report, escape = FALSE, digits = 3,
caption = "SNPs with a p-value < 0.0001")


# Hit SNP rs12243326
# Looking for genes that are within 1,000,000 of this SNP
rs12243326_loc <- SNP_w_genes[snp == "rs12243326",pos] # location of snp
# Finding indices of genes that are within 1,000,000 of this position
rs12243326_genes_idx <- which(abs(rs12243326_loc - GDM_annot$pos) <= 1e6)
# name(s) of the genes that are within a 1Mb window from the SNP rs12243326
rs12243326_within1Mb <- unique(GDM_annot[rs12243326_genes_idx,gene])
rs12243326_1MB_window <- matrix(rs12243326_within1Mb,
dimnames = list(c(),
c("Gene Name")))
knitr::kable(rs12243326_1MB_window, escape = FALSE,
caption = "Genes within a 1Mb window of SNP rs12243326")
# Hit SNP rs2237897
# Looking for genes that are within 1,000,000 of this
rs2237897_loc <- SNP_w_genes[snp == "rs2237897",pos] # location of snp
# Finding indices of genes that are within 1,000,000 of this position
rs2237897_genes_idx <- which(abs(rs2237897_loc - GDM_annot$pos) <= 1e6)
# names of the genes that are within a 1Mb window from the SNP rs12243326
rs2237897_within1Mb <- unique(GDM_annot[rs2237897_genes_idx,gene])
rs2237897_1MB_window <- matrix(rs2237897_within1Mb,
dimnames = list(c(),
c("Gene Name")))
knitr::kable(rs2237897_1MB_window, escape = FALSE,
caption = "Genes within a 1Mb window of SNP rs2237897")

#Names of SNPs with pvalues < 10ˆ{-4}
SNP_lessthan_0.0001 <- SNP_model_results[pvalue < 1e-04,SNPName]
# Names of SNPs with pvalues < 10ˆ{-3}
SNP_lessthan_0.001 <- SNP_model_results[pvalue < 1e-03,SNPName]
# Find index of SNPs on FTO gene which match original data
SNP_FTO_gene.idx <- SNP_model_results$SNP %in% GDM_annot[gene == "FTO", snp]

# Obtaining SNPs with allene on FTO gene
SNP_FTO_gene <- SNP_model_results$SNPName[SNP_FTO_gene.idx]
# Weighted score for SNPs with pvalues < 10ˆ{-4}
# calculated as the data matrix of SNPs multiplied by SNP coefs from log regs)
gdm.dt.imputed$weighted_score_0.0001 <-
as.matrix(gdm_imp_SNPs[,..SNP_lessthan_0.0001]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_lessthan_0.0001,"SNPCoef"])
# Weighted score for SNPs with pvalues < 10ˆ{-3}
# calculated as the data matrix of SNPs multiplied by SNP coefs from log regs)
gdm.dt.imputed$weighted_score_0.001 <-
as.matrix(gdm_imp_SNPs[,..SNP_lessthan_0.001]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_lessthan_0.001,"SNPCoef"])
# Weighted score for SNPs on the FTO gene
# calculated as the data matrix of SNPs multiplied by SNP coefs from log regs)
gdm.dt.imputed$weighted_score_FTO <-
as.matrix(gdm_imp_SNPs[,..SNP_FTO_gene]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_FTO_gene,"SNPCoef"])

# Function to report odds-ratio, 95% conf interval and p-value
logreg_report <- function(glm_model, rname){
# Desired column names of report
column_names <- c("Odds-Ratio", "Lower 95% CI","Upper 95% CI", "p-value")
model_summary <- summary(glm_model)
# Extracting information from model summary
odds_ratio <- round(exp(model_summary$coefficients[2,"Estimate"]),3)
conf_int <- round(exp(confint(glm_model)[2,]),3) #conf interval of odds-ratio
reg_pvalue <- round(model_summary$coefficients[2,"Pr(>|z|)"],5)
model_info <- matrix(c(odds_ratio, conf_int[1],conf_int[2], reg_pvalue),
ncol=4, dimnames = list(rname,column_names))
return(model_info)
}
# GLM regressions for each weighted scores
# GLM for Weighted score for SNPs with pvalues < 10ˆ{-4}
weighted_0.0001_glm <- glm(pheno ~ weighted_score_0.0001,
family = binomial(link="logit"), data = gdm.dt.imputed)
# GLM for Weighted score for SNPs with pvalues < 10ˆ{-3}
weighted_0.001_glm <- glm(pheno ~ weighted_score_0.001,
family = binomial(link="logit"), data = gdm.dt.imputed)

# GLM for Weighted score for SNPs on FTO gene
weighted_FTO_glm <- glm(pheno ~ weighted_score_FTO,
family = binomial(link="logit"), data = gdm.dt.imputed)
weighted_models_report <- rbind(
logreg_report(weighted_0.0001_glm, rname= "p-values < $10ˆ{-4}$"),
logreg_report(weighted_0.001_glm, rname= "p-values < $10ˆ{-3}$"),
logreg_report(weighted_FTO_glm, rname= "SNPs on FTO gene"))

# Converting to data.table
weighted_models_report.dt <- data.table(weighted_models_report,
keep.rownames="Group")
knitr::kable(weighted_models_report.dt, escape = FALSE,
caption = "Table of results for each Weight Genetic Risk Score")

# Loading in data
gdm.test <- fread("data_assignment2/GDM.test.txt")
# Renaming columns to match GDM.raw.txt colnames
colnames(gdm.test) <- colnames(gdm.dt)
# Checking for any NAs
table(is.na(gdm.test)) # 0 NAs in the dataset

## Calculating weighted scores for test data, adding these as cols to gdm.test
# SNPs with pvalues < 10ˆ{-4}
gdm.test$weighted_score_0.0001 <-
as.matrix(gdm.test[,..SNP_lessthan_0.0001]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_lessthan_0.0001,"SNPCoef"])
# SNPs with pvalues < 10ˆ{-3}
gdm.test$weighted_score_0.001 <-
as.matrix(gdm.test[,..SNP_lessthan_0.001]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_lessthan_0.001,"SNPCoef"])
# SNPs on the FTO gene
gdm.test$weighted_score_FTO <-
as.matrix(gdm.test[,..SNP_FTO_gene]) %*%
as.matrix(SNP_model_results[SNPName %in% SNP_FTO_gene,"SNPCoef"])
# Showing top of data for these new weighted scores:
head(gdm.test[,c("ID","sex","pheno","weighted_score_0.0001",
"weighted_score_0.001", "weighted_score_FTO")])

# Using models from problem 2.e to predict gestational diabetes status
# Prediction using model with SNPs with pvalues < 10ˆ{-4}
test_0.0001_predicted <- predict(weighted_0.0001_glm,
newdata = data.frame(weighted_score_0.0001 =
gdm.test$weighted_score_0.0001),
type="response")
# Prediction using model with SNPs with pvalues < 10ˆ{-3}
test_0.001_predicted <- predict(weighted_0.001_glm,
newdata = data.frame(weighted_score_0.001 =
gdm.test$weighted_score_0.001),
type="response")
# Prediction using model with SNPs on FTO gene
test_FTO_predicted <- predict(weighted_FTO_glm,
newdata = data.frame(weighted_score_FTO =
gdm.test$weighted_score_FTO),
type="response")

# Computing log-likelihood for the predicted probabilities from the
# three genetic risk score models
# Giving we have a binary outcome (have gestational diabetes or not)
# know the predicted values should follow a Binomial/Bernoulli distribution
loglik_0.0001 <- sum(dbinom(gdm.test$pheno,
prob=test_0.0001_predicted, size=1, log=TRUE))
loglik_0.001 <- sum(dbinom(gdm.test$pheno,
prob=test_0.001_predicted, size=1, log=TRUE))
loglik_FTO <- sum(dbinom(gdm.test$pheno,
prob=test_FTO_predicted, size=1, log=TRUE))
loglik_summary <- matrix(c(loglik_0.0001,loglik_0.001,loglik_FTO),
ncol=3,
dimnames = list(c("Log-likelihood"),
c("p < $10ˆ{-4}$",
"p < $10ˆ{-3}$",
"FTO gene")))
loglik_summary.dt <- data.table(loglik_summary, keep.rownames = "Group")
knitr::kable(loglik_summary.dt, escape = FALSE, digits=3,
caption = "Log-likelihood values for each Weight Genetic Risk Score")

# Loading in data
gdm2.dt <- fread("data_assignment2/GDM.study2.
txt")
# Check if different studies contain all the same SNP identifiers
all(gdm2.dt$snp %in% SNP_model_results$SNP)
all(SNP_model_results$SNP %in% gdm2.dt$snp)

# all SNP identifiers are the same
# Creating new columns in gdm2.dt which is snp with effect allele
gdm2.dt[, SNPName := paste(snp, effect.allele, sep = "_")]
# Check if both studies contain the same SNP with alleles (only effect)
match_SNPNames_gdm <- gdm2.dt$SNPName %in% SNP_model_results$SNPName
match_SNPNames_orig <- SNP_model_results$SNPName %in% gdm2.dt$SNPName
sum_match_gdm <- sum(match_SNPNames_gdm)
sum_match_orig <- sum(match_SNPNames_orig)
cat("There are",nrow(gdm2.dt),"SNPs and alleles in gdm2.dt. \nOf these",
sum_match_gdm, "match both SNPs and alleles (effect allele).","\n")

cat("There are",nrow(SNP_model_results),
"SNPs and alleles in SNP Model Results. \nOf these",
sum_match_orig, "match with the gdm dataset (effect allele).")

# Only keep SNPs which match across both studies (effect allele)
gdm2.dt_meta <- gdm2.dt[match_SNPNames_gdm]
SNP_model_results_meta <- SNP_model_results[match_SNPNames_orig]
# Add new column with study names
SNP_model_results_meta$Study <- "Original Study"
gdm2.dt_meta$Study <- "GDM Study"
# Calculate p-values for gdm.dt data
gdm2.dt_meta[, pvalue := pnorm(-abs(beta)/se)*2]
# Only obtain SNPs with p-values < 10ˆ{-4} for original model results
high_sig_SNP_results <- SNP_model_results_meta[pvalue < 1e-04,]
# Only wanting to keep Study, SNPName, SNPCoef, StandardError, and pvalue
high_sig_SNP_results <- high_sig_SNP_results[,c("Study","SNPName",
"SNPCoef","StandardError",
"pvalue")]
# Only obtain SNPs with p-values < 10ˆ{-4} for original model results
high_sig_gdm_results <- gdm2.dt_meta[pvalue < 1e-04,]
# Only wanting to keep Study, SNPName, beta, se, and pvalue
high_sig_gdm_results <- high_sig_gdm_results[,c("Study","SNPName",
"beta","se",
"pvalue")]
# Columns are in all the same position, binding significant results together
joint_sig_results <- rbind(high_sig_SNP_results, high_sig_gdm_results,
use.names=FALSE)
# Making sure results for SNPs are easily comparable
setorder(joint_sig_results, SNPName)
# Printing results
knitr::kable(joint_sig_results, escape = FALSE,
caption = "Summary of meta-analysis results for SNPs with p-values $< 10ˆ{-4}$")