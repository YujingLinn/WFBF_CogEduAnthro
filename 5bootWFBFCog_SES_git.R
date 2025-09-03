#!/usr/bin/Rscript
#.libPaths(c("/users/k21170717/.conda/envs/r_env/lib/R/library", .libPaths()))

# title: Regression_adjSES with Bootstrap & Weighted Averages 
# author: Yujing Lin
# date: 14th October, 2024
# codes corresponding to publication: https://icajournal.scholasticahq.com/article/140654-polygenic-score-prediction-within-and-between-sibling-pairs-for-intelligence-cognitive-abilities-and-educational-traits-from-childhood-to-early-adul

# Family SES 
#"ases" # "SES composite scale (1st Contact), standardised"
#"gses" # "SES composite variable (7 year parent), standardised"
#"pses" # "SES composite at 16 from parent web questionnaire, standardised"
#"u1pses" # "Parent SES composite (TEDS21 phase 1 parent qnr), standardised"
#"u1pses5in" # "Parental SES item 5, household income (TEDS21 phase 1 parent qnr), see value labels"

nboot <- 10 # change it to 1000 on HPC
ncpus <- 4 # change it to 32 on HPC

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/Data/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'
outFileStemRDS <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/RDS/'
#sourceFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/data_160425/'
#outFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/results_160425/'
#outFileStemRDS <- '/scratch/users/k21170717/WFBF/Cognitive/resultsRDS_160425/'

library(dplyr)
library(scales)
library(tidyr)
library(Hmisc) 
library(data.table)
library(boot)
library(metafor)
library(openxlsx)
library(MASS) # ordered logistic or probit regression
library(pscl)

WFBF_selectunpaired <- read.csv(paste(sourceFileStem, "WFBF_selectunpaired.csv", sep=""))
dat_pair <- read.csv(paste(sourceFileStem, "dat_pair.csv", sep=""))

# add the SES variables to the dat_pair
WFBF_full <- read.csv(paste(sourceFileStem, "WFBF_full.csv", sep=""))
SES_var <- subset(WFBF_full, select=c("id_twin", "id_fam", "ases", "gses", "pses", "u1pses"))
dat_pair <- full_join(dat_pair, SES_var)

# just the main analyses, we are not doing sex-stratified, zygosity-stratified population effects or sex-zygosity stratified within-family effects
dat_pair_selectunpaired <- subset(dat_pair, selectunpaired==1)
dat_DZpair_selectunpaired <- subset(dat_pair_selectunpaired, dat_pair_selectunpaired$zygos==2)

PCs <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# variable list ####
VarPairs_SES <- list( # it has to be a list, otherwise, each row will have auto-refill to match the maximum length of 4, i.e., repeating the continuous variables twice
  c("gt2ac", "gses", "EA4_no23andme_Okbay2022"),
  c("gcg", "gses", "IQ_Savage2018_FRCT1"),
  c("gcl", "gses", "IQ_Savage2018_FRCT1"),
  c("gcn", "gses", "IQ_Savage2018_FRCT1"),
  c("gbmi", "gses", "BMI_Giant2018"), # because if adjusting for PGS, in the context of SES it would be parental PGS, so no need to use childhood BMI
  c("ghtcm", "gses", "Height_Yengo2022"), # 7
  
  c("it3ac", "gses", "EA4_no23andme_Okbay2022"),
  c("icg", "gses", "IQ_Savage2018_FRCT1"),
  c("icvb", "gses", "IQ_Savage2018_FRCT1"),
  c("icnv", "gses", "IQ_Savage2018_FRCT1"), # 9
  
  c("jt3ac", "gses", "EA4_no23andme_Okbay2022"),
  c("jcg", "gses", "IQ_Savage2018_FRCT1"),
  c("jcvb", "gses", "IQ_Savage2018_FRCT1"),
  c("jcnv", "gses", "IQ_Savage2018_FRCT1"), # 10
  
  c("lt3ac", "gses", "EA4_no23andme_Okbay2022"),
  c("lcg", "gses", "IQ_Savage2018_FRCT1"),
  c("lverbal12T", "gses", "IQ_Savage2018_FRCT1"), 
  c("lnonverbal12T", "gses", "IQ_Savage2018_FRCT1"),
  c("lbmi", "gses", "BMI_Giant2018"), 
  c("lchtcm", "gses", "Height_Yengo2022"), # 12
  
  c("ncbmi", "pses", "BMI_Giant2018"),
  c("nchtcm", "pses", "Height_Yengo2022"), # 14
  
  c("pcexgcsecoregrdm", "pses", "EA4_no23andme_Okbay2022"),
  c("pcg", "pses", "IQ_Savage2018_FRCT1"),
  c("pcvctota", "pses", "IQ_Savage2018_FRCT1"),
  c("pcrvtota", "pses", "IQ_Savage2018_FRCT1"),
  c("pcbmi", "pses", "BMI_Giant2018"),
  c("pcqdhtcm", "pses", "Height_Yengo2022"), # 16
  
  c("ra_level_enrol", "pses", "EA4_no23andme_Okbay2022", "rcqalage1"), # binary
  c("rcqalgrdm", "pses", "EA4_no23andme_Okbay2022"),
  c("rcqhe", "pses", "EA4_no23andme_Okbay2022", "rcqalage1"), # 18 # binary
  
  c("u1cdegr1", "u1pses", "EA4_no23andme_Okbay2022"), # degree description at age 21
  c("u1cbmi", "u1pses", "BMI_Giant2018"),
  c("u1chtcm", "u1pses", "Height_Yengo2022"),
  c("ucgt", "u1pses", "IQ_Savage2018_FRCT1"),
  c("ucgvbt", "u1pses", "IQ_Savage2018_FRCT1"),
  c("ucgnvt", "u1pses", "IQ_Savage2018_FRCT1"), # 21
  
  c("zEA", "u1pses", "EA4_no23andme_Okbay2022"),
  c("zmhbmi", "u1pses", "BMI_Giant2018"),
  c("zmhheight", "u1pses", "Height_Yengo2022") # 26
)

length(VarPairs_SES) # 40



# the functions ####
# basically the same as the function from 4bootWFBF, only changing some column # to match the additional covariate of SES
# yea maybe I can just add an extra condition to the original function to easily incorporate SES in or not in the model, next time 
std_fit_regression_boot <- function(dat, indices, x, y, covar) {
  # datx <- dat_std
  datx <- dat[indices, ] # Subset data using bootstrapped indices
  formula <- reformulate(covar, response = y) # Dynamically create the formula
  fit_model <- lm(formula, data = datx) # Fit the model
  #print(formula)
  
  # Calculate slope, standardized beta, p-value, and R²
  x_slope <- coef(fit_model)[x] # Slope 
  x_beta <- x_slope * sd(datx[[x]], na.rm = TRUE) / sd(datx[[y]], na.rm = TRUE) # Standardized beta
  
  se_x <- sqrt(diag(vcov(fit_model))[x])  # SE for the slope
  t_x <- x_beta / se_x  # t-value
  p_x_beta <- 2 * (1 - pt(abs(t_x), df.residual(fit_model)))  # p-value
  R.squared <- summary(fit_model)$r.squared  # R²
  
  return(c(x_beta, p_x_beta, R.squared)) # Return beta, p-value, and R²
  # though it would be nice to include the PGS and trait here, the output from boot cannot be characters; otherwise, everything will be turned into characters
}

std_fit_logistic_boot <- function(dat, indices, x, y, covar) {
  # datx <- dat_std
  datx <- dat[indices, ]
  unique_vals <- unique(datx[[y]][!is.na(datx[[y]])])
  datx[[y]] <- factor(datx[[y]], levels = sort(unique_vals), ordered = TRUE) 
  
  formula <- reformulate(covar, response = y)
  #print(formula)
  
  # when ordered = TRUE, it would apply to meaningful rankings like 1=enter uni, 0=not enter uni; it won't apply to 0=female 1=male
  # in this case, we are ordering the outcome trait, so the order is meaningful
  if (length(unique_vals) == 2) {
    # binary logistic regression
    fit_model <- glm(formula, family = binomial(link = "logit"), data = datx)
  } else if (length(unique_vals) > 2) {
    # ordered logistic regression
    # when there are 3 ordinal levels or more, use polr, which fits a logistic regression model to an ordered factor response 
    # this could not be used for binary variables 
    fit_model <- polr(formula, data = datx, Hess = TRUE)
  }
  #fit_model
  
  logistic_coef <- coef(fit_model)[x]
  odds_ratio <- exp(logistic_coef) # logistic coefficient to odds ratio
  cohens_d <- logistic_coef * sqrt(3) / pi # logistic coefficient to Cohen's d
  pearsons_r <- logistic_coef / sqrt(logistic_coef^2 + (pi^2 / 3)) # logistic coefficient to Pearson's r; this is basically the point-biserial correlation and assumes the predictors are std
  # the above conversion becomes less accurate with extreme probabilities
  # since we have 50% of our participants going to A-level, and 50% of which going to uni, we have a very balanced sample, which is robust to the formula 
  x_beta <- pearsons_r
  
  se_x <- sqrt(diag(vcov(fit_model))[x])  # SE for the slope
  
  z_value <- logistic_coef / se_x
  p_x_beta <- 2 * (1 - pnorm(abs(z_value)))
  
  invisible(capture.output(pseudo_r2_result <- pR2(fit_model)))
  # suppress the output: "fitting null model for pseudo-r2"; otherwise, too many output for bootstrapping
  # Extract the r2ML value (Maximum likelihood pseudo R²)
  R.squared <- pseudo_r2_result["r2ML"] # Maximum likelihood pseudo r-squared--this is supposed to be the closest R² to a linear R²
  return(c(x_beta, p_x_beta, R.squared))
}

# Perform bootstrapping 
beta_linear_results_fun_SES_boot <- function(dat, paste_item_PGS, paste_item_trait, sex_var, nboot, ncpus, rds_name) {
  #dat <- WFBF_selectunpaired # WFBF_selectunpaired   dat_DZpair_selectunpaired   dat_DZpair_selectunpaired_FSS
  
  beta_linear_raw <- lapply(1:length(VarPairs_SES), function(i) {
    #i <- 30 # 29: binary
    cat("Processing item", i, ":", VarPairs_SES[[i]], "\n")
    is_binary <- length(VarPairs_SES[[i]]) == 4 # different from the original function, which is 3
    
    #paste_item_trait <- "1" # "_trait_diff", "1", "_trait_sum"
    #paste_item_PGS <- "" # "_PGS_diff", "", "_PGS_sum"
    #sex_var <- "sex1_label" # "sex1_label_for_pairdiff", "sex1_label", "sex1_label_for_pairmeansum"
    
    # to std the continuous covar (have the std for the outcome trait below)
    if (is_binary) { # binary
      continuous_vars <- c(
        VarPairs_SES[[i]][2],                                   # SES 
        paste(VarPairs_SES[[i]][3], paste_item_PGS, sep = ""),  # PGS variable
        VarPairs_SES[[i]][4],                                   # Age variable
        paste0("PC", 1:10)                                      # PCs 
      )
    } else {
      continuous_vars <- c(
        VarPairs_SES[[i]][2],                                   # SES 
        paste(VarPairs_SES[[i]][3], paste_item_PGS, sep = ""),  # PGS variable
        paste0("PC", 1:10)                                      # PCs
      )
    }
    
    # Standardize them all at once (skip binary/categorical variables like sex)
    # for some reasons, this function sometimes may not work probably due to conflicts among the packages
    # if that happens, just close Rstudio and restart it
    dat_std <- dat %>%
      mutate(across(
        .cols = all_of(continuous_vars),
        .fns = ~ as.numeric(scale(.)),  # Standardize (mean=0, SD=1)
        .names = "{.col}_std"            # Append "_std" to standardized vars
      ))
    
    y <- paste(VarPairs_SES[[i]][1], paste_item_trait, sep = "")  # Trait (will std later)
    x <- paste0(paste(VarPairs_SES[[i]][3], paste_item_PGS, sep = ""), "_std")  # PGS (standardized)
    SES <- paste0(VarPairs_SES[[i]][2], "_std") # SES
    if (is_binary) {
      age <- paste0(VarPairs_SES[[i]][4], "_std")  # Use standardized version
    } else {
      age <- NULL  # Or handle missing age appropriately
    }
    PCs <- paste0("PC", 1:10, "_std")        # Standardized PCs
    
    # Covariates (now using standardized versions)
    covar <- c(x, SES, age, sex_var, "chiptype", PCs)  # Sex and chiptype remains unstandardized
    covar <- covar[!is.na(covar)]
    # variable_name_for_label <- paste(VarPairs_SES[[i]][1], 1, sep="")
    
    # remove sex_var if necessary & std y for continuous outcome
    if (length(unique(dat_std[[sex_var]])) == 1) { # Remove sex if only one level exists
      covar <- covar[covar != sex_var]  # Drop sex if only 1 level (sex-stratified)
      
      if (!is_binary) {
        # don't forget to std the continuous outcome trait y in the sex-stratified analyses
        dat_std[[y]] <- scale(dat_std[[y]])
      }
      
    } else {
      # Use VarPairs_SES length to determine if y is binary (3) or continuous (2)
      
      
      if (!is_binary) {  # Continuous trait (VarPairs_SES[[i]] length = 2)
        # Correct y for sex and remove sex_var from covariates
        formula_sex <- reformulate(sex_var, response = y)
        dat_std[[y]] <- rstandard(lm(formula_sex, data = dat_std, na.action = na.exclude))
        covar <- covar[covar != sex_var]
      }
      # Binary trait (VarPairs_SES[[i]] length = 3): sex_var remains in covar & no std for the outcome trait 
    }
    
    boot_results <- boot(data = dat_std, 
                         if (is_binary) {
                           statistic = std_fit_logistic_boot
                         } else {
                           statistic = std_fit_regression_boot
                         },
                         R = nboot, 
                         parallel = "multicore", 
                         ncpus = ncpus, 
                         x = x, 
                         y = y, 
                         covar = covar)
    
    boot_beta_dist <- boot_results$t[,1]  # just the x_beta from boot
    
    return(boot_results)
  })
  
  boot_dist_list <- lapply(beta_linear_raw, function(boot) boot$t[,1])
  names(boot_dist_list) <- sapply(VarPairs_SES, function(v) paste0(v[1]))
  
  saveRDS(
    boot_dist_list,
    file = paste0(outFileStemRDS, "boot_dist_", rds_name, ".rds")
  )
  
  results_list <- lapply(1:length(beta_linear_raw), function(i) {
    # i <- 30
    boot_results <- beta_linear_raw[[i]]
    
    trait_y <- VarPairs_SES[[i]][1]
    PGS_x <- VarPairs_SES[[i]][2]
    
    x_beta <- boot_results$t0[1] # this is the observed beta (an input), so would be the same as the beta generated by the regression model
    # x_beta
    # there are 3 values in total for my boot_results, b/c I specify for 3 outputs using my std_fit_regression_boot function above
    beta_se_boot <- sd(boot_results$t[,1]) # unlike t0 corresponds with the observed values, t corresponds with the bootstrapped values, so we will use t for se
    #the standard deviation of the bootstrap estimates is the standard error of the sample estimates
    beta_ci <- boot.ci(boot_results, type = "perc")  # Percentile CI
    # beta_ci <- boot.ci(boot_results, type = "bca") # the bias-adjusted option doesn't work: Error in bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  : estimated adjustment 'a' is NA
    lower.x_beta_95CI_boot <- beta_ci$percent[4]
    upper.x_beta_95CI_boot <- beta_ci$percent[5]
    
    p_x_beta <- mean(boot_results$t0[2])
    R.squared <- mean(boot_results$t0[3])
    
    # Return the results as a named list
    list(
      PGS_x = PGS_x,
      trait_y = trait_y,
      x_beta = x_beta,
      beta_se_boot = beta_se_boot,
      lower.x_beta_95CI_boot = lower.x_beta_95CI_boot,
      upper.x_beta_95CI_boot = upper.x_beta_95CI_boot,
      p_x_beta = p_x_beta,
      R.squared = R.squared
    )
  })
  
  beta_linear_results <- do.call(rbind, lapply(results_list, as.data.frame))
  # colnames(beta_linear_results)<-c("PGS_x", "trait_y", "x_beta", "se_beta", "p_x_beta", "R.squared")
  
  beta_linear_results <- beta_linear_results %>%
    mutate(x_beta_se = paste0(round(x_beta, 2), ' (', round(beta_se_boot, 2), ')'))
  
  return(beta_linear_results)
}





# calculation of pop and WF estimates ####
cat("ppl_unrelated_selectunpaired_boot_adjSES \n")
ppl_unrelated_selectunpaired_boot_adjSES <- beta_linear_results_fun_SES_boot(WFBF_selectunpaired, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_boot_adjSES") # suffix for PGS, suffix for trait
ppl_unrelated_selectunpaired_boot_adjSES$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_boot_adjSES$p_x_beta, method="fdr")

cat("WF_DZpairdiff_boot_adjSES \n")
WF_DZpairdiff_boot_adjSES <- beta_linear_results_fun_SES_boot(dat_DZpair_selectunpaired, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZpairdiff_boot_adjSES")
WF_DZpairdiff_boot_adjSES$p_x_beta <- p.adjust(WF_DZpairdiff_boot_adjSES$p_x_beta, method="fdr")



# aggregate the simple regression results ####
ppl_unrelated_selectunpaired_boot_adjSES$Type <- "ppl_unrelated_adjSES"
WF_DZpairdiff_boot_adjSES$Type <- "WF_DZpairdiff_adjSES"

raw_results_boot_adjSES <- full_join(ppl_unrelated_selectunpaired_boot_adjSES, WF_DZpairdiff_boot_adjSES)

# process the results ####
Trait_Domain <- as.data.frame(rbind(
  # g (IQ3 PGS)
  c("gcg", "childhood g", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icg", "childhood g", "Cognitive Abilities"),
  c("jcg", "childhood g", "Cognitive Abilities"),
  c("lcg", "adolescence g", "Cognitive Abilities"), # Adolescence (12+16)
  c("pcg", "adolescence g", "Cognitive Abilities"),
  c("ucgt", "adulthood g", "Cognitive Abilities"), # Adulthood (25)
  # verbal g
  c("gcl", "childhood verbal g", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icvb", "childhood verbal g", "Cognitive Abilities"),
  c("jcvb", "childhood verbal g", "Cognitive Abilities"),
  c("lverbal12T", "adolescence verbal g", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcvctota", "adolescence verbal g", "Cognitive Abilities"),
  c("ucgvbt", "adulthood verbal g", "Cognitive Abilities"), # Adulthood (25)
  # nonverbal g
  c("gcn", "childhood nonverbal g", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icnv", "childhood nonverbal g", "Cognitive Abilities"),
  c("jcnv", "childhood nonverbal g", "Cognitive Abilities"),
  c("lnonverbal12T", "adolescence nonverbal g", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcrvtota", "adolescence nonverbal g", "Cognitive Abilities"),
  c("ucgnvt", "adulthood nonverbal g", "Cognitive Abilities"), # Adulthood (25)
  # Educational achievement
  c("gt2ac", "primary school grades", "Education Achievement"), # Primary (7+9+10+12)
  c("it3ac", "primary school grades", "Education Achievement"), 
  c("jt3ac", "primary school grades", "Education Achievement"), 
  c("lt3ac", "primary school grades", "Education Achievement"), 
  c("pcexgcsecoregrdm", "GCSE grades", "Education Achievement"), # GCSE (16)
  c("rcqalgrdm", "A-level grades", "Education Achievement"), # A-level grades (18) # use A-level grade only, NOT A- & AS-levels
  c("u1cdegr1", "university grades", "Education Achievement"), # university grades	 (21)
  # Educational attainment
  c("ra_level_enrol", "A-level enrollment", "Education Attainment"), # A-level enrolment (16), data obtained at 18
  c("rcqhe", "university enrollment", "Education Attainment"), # university enrolment (18)
  c("zEA", "years of schooling", "Education Attainment"), # Years of schooling (21 & 26)
  # BMI
  c("gbmi", "childhood BMI", "Height & BMI"), # ages 7 childhood
  c("lbmi", "adolescence BMI", "Height & BMI"), # 12, 14, 16 adolescence
  c("ncbmi", "adolescence BMI", "Height & BMI"),
  c("pcbmi", "adolescence BMI", "Height & BMI"),
  c("u1cbmi", "adulthood BMI", "Height & BMI"), # 22, 26 adulthood
  c("zmhbmi", "adulthood BMI", "Height & BMI" ),
  # height
  c("ghtcm", "childhood height", "Height & BMI"), # ages 7 childhood
  c("lchtcm", "adolescence height", "Height & BMI"), # 12, 14, 16 adolescence
  c("nchtcm", "adolescence height", "Height & BMI"), 
  c("pcqdhtcm", "adolescence height", "Height & BMI"),
  c("u1chtcm", "adulthood height", "Height & BMI"), # 22, 26 adulthood
  c("zmhheight", "adulthood height", "Height & BMI")
  
  # removed-these are kept in case the nrow doesn't match
  # c("lcnv", "rm"), # nonverbal g at age 12
  # c("lcvb", "rm"), # verbal g at age 12
  # c("rcqalsgrdm", "rm"), # A- & AS-level grades
  # c("u1cedat", "rm") # degree at age 21
))
dim(Trait_Domain) # 40*3
colnames(Trait_Domain) <- c("trait_y", "Trait_Domain", "V3") # make sure the column name for traits in Trait_Domain is the same as the result column name

raw_results_boot_adjSES <- full_join(raw_results_boot_adjSES, Trait_Domain)



# save the pre-meta-analysed results ####
raw_results_boot_adjSES <- raw_results_boot_adjSES %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("ppl_unrelated_adjSES","WF_DZpairdiff_adjSES"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, V3, Trait_Domain)

write.csv(raw_results_boot_adjSES, paste(outFileStem, "SESadj_WFBF_Cog_Main_boot_raw.csv", sep=""), row.names = FALSE)





cat("compare bootstrapped distributions \n")
# note that we are comparing the bootstrapped distributions of the SES adjusted pop and WF
# as well as the bootstrapped distributions of the SES adjusted pop and the unadjusted pop, SES adjusted WF and the unadjusted WF
boot_pop_adjSES <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_boot_adjSES.rds"))
boot_WF_adjSES <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_boot_adjSES.rds"))

boot_pop <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_boot.rds"))
boot_WF <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_boot.rds"))

fun_compare_boot_dists <- function(boot_dataset1, boot_dataset2, fdr_correction = TRUE) {
  shared_keys <- intersect(names(boot_dataset1), names(boot_dataset2))
  
  results <- lapply(shared_keys, function(key) {
    diffs <- boot_dataset1[[key]] - boot_dataset2[[key]]
    ci_diff <- quantile(diffs, c(0.025, 0.975))
    mean_diff <- mean(diffs)
    se_diff <- sd(diffs)
    z <- mean_diff / se_diff
    p <- 2 * pnorm(-abs(z))
    data.frame(
      trait_y = key,
      diff_mean = mean_diff,
      diff_ci_lower = ci_diff[1],
      diff_ci_upper = ci_diff[2],
      diff_z = z,
      diff_p = p
    )
  })
  
  result_df <- do.call(rbind, results)
  
  if (fdr_correction) {
    result_df$diff_p_fdr <- p.adjust(result_df$diff_p, method = "fdr")
  }
  
  return(result_df)
}

boot_compare_SimReg_popWF_adjSES <- fun_compare_boot_dists(boot_pop_adjSES, boot_WF_adjSES)
boot_compare_SimReg_popWF_adjSES <- full_join(boot_compare_SimReg_popWF_adjSES, Trait_Domain)

boot_compare_SimReg_pop_wwoadjSES <- fun_compare_boot_dists(boot_pop_adjSES, boot_pop) # wwo = with/without
boot_compare_SimReg_pop_wwoadjSES <- full_join(boot_compare_SimReg_pop_wwoadjSES, Trait_Domain)

boot_compare_SimReg_WF_wwoadjSES <- fun_compare_boot_dists(boot_WF_adjSES, boot_WF) # wwo = with/without
boot_compare_SimReg_WF_wwoadjSES <- full_join(boot_compare_SimReg_WF_wwoadjSES, Trait_Domain)

write.csv(boot_compare_SimReg_popWF_adjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_pop_wwoadjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main_POP_wwoadjSES.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_WF_wwoadjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main_WF_wwoadjSES.csv", sep=""), row.names = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="popWF_adjSES")
addWorksheet(wb, sheetName="pop_wwoadjSES")
addWorksheet(wb, sheetName="WF_wwoadjSES")

# Write each table to its corresponding worksheet
writeData(wb, sheet="popWF_adjSES", x=boot_compare_SimReg_popWF_adjSES, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="pop_wwoadjSES", x=boot_compare_SimReg_pop_wwoadjSES, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_wwoadjSES", x=boot_compare_SimReg_WF_wwoadjSES, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_raw_All.xlsx", sep=""), overwrite=TRUE)







cat("meta-analyses for regression_adjSES \n")
run_meta_analysis <- function(result) {
  results_list <- list()
  # result <- raw_results_boot
  #result <- raw_sensitivity_boot
  unique_combinations <- unique(result[, c("Trait_Domain", "Type", "V3")])
  
  # Loop through each unique combination
  for (i in seq_len(nrow(unique_combinations))) {
    # i <- 1
    current_data <- filter(result, 
                           Trait_Domain == unique_combinations$Trait_Domain[i],
                           Type == unique_combinations$Type[i],
                           V3 == unique_combinations$V3[i])
    
    rma_result <- rma(yi = current_data$x_beta, 
                      sei = current_data$beta_se_boot, 
                      weights = 1 / (current_data$beta_se_boot^2), 
                      method = "REML") # random-effect
    
    results_list[[i]] <- list(
      Trait_Domain = unique_combinations$Trait_Domain[i],
      Type = unique_combinations$Type[i],
      V3 = unique_combinations$V3[i],
      beta = rma_result$beta,
      se = rma_result$se,
      ci_lower = rma_result$ci.lb,
      ci_upper = rma_result$ci.ub,
      beta_se = paste0(round(rma_result$beta, 2), ' (', round(rma_result$se, 2), ')')
    )
  }
  result_df <- do.call(rbind, lapply(results_list, as.data.frame))
  
  result_df <- result_df %>%
    mutate(Trait_Domain = factor(Trait_Domain, levels = c(
      "childhood g", "adolescence g", "adulthood g", 
      "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
      "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
      "primary school grades", "GCSE grades", "A-level grades", "university grades", 
      "A-level enrollment", "university enrollment", "years of schooling", 
      "childhood BMI", "adolescence BMI", "adulthood BMI", 
      "childhood height", "adolescence height", "adulthood height"))) %>%
    arrange(Type, V3, Trait_Domain)
  
  return(result_df)
}

MA_results_boot_adjSES <- run_meta_analysis(raw_results_boot_adjSES)



cat("meta-analysis_adjSES done \n")





# tidy the results ####
MA_results_boot_adjSES <- MA_results_boot_adjSES %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("WF_DZpairdiff_adjSES","ppl_unrelated_adjSES"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, V3, Trait_Domain)



# save the results ####
write.csv(MA_results_boot_adjSES , paste(outFileStem, "SESadj_MA_WFBF_Cog_Main_sensi_boot.csv", sep=""), row.names = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="")
writeData(wb, sheet=1, x=MA_results_boot_adjSES, colNames=TRUE, rowNames=FALSE)
saveWorkbook(wb, paste(outFileStem, "SESadj_MA_WFBF_Cog_Main_sensi_boot.xlsx", sep=""), overwrite = TRUE) 

cat("main analysis adj. for SES raw results and meta-analysed results done \n")

cat("meta-analyses for estimate comparisons \n")
run_meta_diff_analysis <- function(diff_result) {
  results_list <- list()
  unique_combinations <- unique(diff_result[, c("Trait_Domain", "V3")])
  
  for (i in seq_len(nrow(unique_combinations))) {
    domain <- unique_combinations$Trait_Domain[i]
    v3 <- unique_combinations$V3[i]
    
    domain_data <- filter(diff_result, 
                          Trait_Domain == domain,
                          V3 == v3)
    
    # Estimate SE from CI
    domain_data <- domain_data %>%
      mutate(diff_se = (diff_ci_upper - diff_ci_lower) / (2 * 1.96))
    
    # Run meta-analysis
    rma_result <- rma(
      yi = diff_mean,
      sei = diff_se,
      data = domain_data,
      method = "REML"
    )
    
    results_list[[i]] <- data.frame(
      Trait_Domain = domain,
      V3 = v3,
      diff_beta = rma_result$b,
      diff_se = rma_result$se,
      diff_ci_lower = rma_result$ci.lb,
      diff_ci_upper = rma_result$ci.ub,
      diff_z = rma_result$zval,
      diff_p = rma_result$pval,
      diff_beta_se = paste0(round(rma_result$b, 2), " (", round(rma_result$se, 2), ")")
    )
  }
  
  result_df <- do.call(rbind, results_list)
  
  # Apply FDR correction
  result_df$diff_p_fdr <- p.adjust(result_df$diff_p, method = "fdr")
  
  result_df <- result_df %>%
    mutate(Trait_Domain = factor(Trait_Domain, levels = c(
      "childhood g", "adolescence g", "adulthood g", 
      "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
      "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
      "primary school grades", "GCSE grades", "A-level grades", "university grades", 
      "A-level enrollment", "university enrollment", "years of schooling", 
      "childhood BMI", "adolescence BMI", "adulthood BMI", 
      "childhood height", "adolescence height", "adulthood height"))) %>%
    arrange(V3, Trait_Domain)
  
  return(result_df)
}

write.csv(boot_compare_SimReg_popWF_adjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_pop_wwoadjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main_POP_wwoadjSES.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_WF_wwoadjSES, paste(outFileStem, "bootSESadj_compare_WFBF_Cog_Main_WF_wwoadjSES.csv", sep=""), row.names = FALSE)


bootMA_compare_SimReg_popWF_adjSES <- run_meta_diff_analysis(boot_compare_SimReg_popWF_adjSES)
bootMA_compare_SimReg_pop_wwoadjSES <- run_meta_diff_analysis(boot_compare_SimReg_pop_wwoadjSES)
bootMA_compare_SimReg_WF_wwoadjSES <- run_meta_diff_analysis(boot_compare_SimReg_WF_wwoadjSES)

# save the results ####
write.csv(bootMA_compare_SimReg_popWF_adjSES, paste(outFileStem, "bootSESadjMA_compare_WFBF_Cog_Main.csv", sep=""), row.names = FALSE)
write.csv(bootMA_compare_SimReg_pop_wwoadjSES, paste(outFileStem, "bootSESadjMA_compare_WFBF_Cog_Main_POP_wwoadjSES.csv", sep=""), row.names = FALSE)
write.csv(bootMA_compare_SimReg_WF_wwoadjSES, paste(outFileStem, "bootSESadjMA_compare_WFBF_Cog_Main_WF_wwoadjSES.csv", sep=""), row.names = FALSE)

wb <- createWorkbook()

# Add a worksheet for each table with descriptive names
addWorksheet(wb, sheetName="popWF_adjSES")
addWorksheet(wb, sheetName="pop_wwoadjSES")
addWorksheet(wb, sheetName="WF_wwoadjSES")

# Write each table to its corresponding worksheet
writeData(wb, sheet="popWF_adjSES", x=bootMA_compare_SimReg_popWF_adjSES, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="pop_wwoadjSES", x=bootMA_compare_SimReg_pop_wwoadjSES, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_wwoadjSES", x=bootMA_compare_SimReg_WF_wwoadjSES, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "bootSESadjMA_compare_WFBF_Cog_All.xlsx", sep=""), overwrite=TRUE)

cat("done \n")








# no more bootstrapping below, just run directly
# modification: 20 April, 2025: add correlations between SES at available ages and PGS in additional to trait 
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#rm(list=ls())
sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/Data/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'
WFBF_selectunpaired <- read.csv(paste(sourceFileStem, "WFBF_selectunpaired.csv", sep=""))

# correlation with SES only ####
cat("correlation with SES only \n")

# Simple regression function for SES with trait and PGS--no bootstrapping
# since SES will always be te outcome, we don't need to do logistic regressions 
std_fit_regression_noboot_SES <- function(dat, x, y, covar) {
  formula <- reformulate(covar, response = y) # Dynamically create the formula
  fit_model <- lm(formula, data = dat) # Fit the model
  #print(formula)
  
  # Calculate slope, standardized beta, p-value, and R²
  x_slope <- coef(fit_model)[x] # Slope 
  x_beta <- x_slope * sd(dat[[x]], na.rm = TRUE) / sd(dat[[y]], na.rm = TRUE) # Standardized beta
  
  se_x <- sqrt(diag(vcov(fit_model))[x])  # SE for the slope
  t_x <- x_beta / se_x  # t-value
  p_x_beta <- 2 * (1 - pt(abs(t_x), df.residual(fit_model)))  # p-value
  R.squared <- summary(fit_model)$r.squared  # R²
  
  CI_lower <- x_beta - 1.96 * se_x
  CI_upper <- x_beta + 1.96 * se_x
  
  return(c(x_beta, se_x, p_x_beta, R.squared, CI_lower, CI_upper)) # Return beta, p-value, and R²
  # though it would be nice to include the PGS and trait here, the output from boot cannot be characters; otherwise, everything will be turned into characters
}


trait_SES_corr_fun <- function(dat) {
  results <- list()
  
  for (i in seq_along(VarPairs_SES)) {
    cat("Processing item", i, ":", VarPairs_SES[[i]], "\n")
    pair <- VarPairs_SES[[i]]
    is_binary <- length(pair) == 4 && !is.na(pair[4])
    
    # Extract variables
    trait <- paste0(pair[1], "1")  # trait with suffix 1
    ses <- pair[2]                 # SES variable (no suffix)
    
    # Standardize SES variable
    dat[[paste0(ses, "_std")]] <- as.numeric(scale(dat[[ses]]))
    
    # Handle trait differently based on whether it's binary
    if (is_binary) {
      age <- pair[4]  # age variable for binary traits
      
      # Standardize age if needed
      if (!is.null(age) && !is.na(age)) {
        dat[[paste0(age, "_std")]] <- as.numeric(scale(dat[[age]]))
        age_var <- paste0(age, "_std")
      } else {
        age_var <- NULL
      }
      
      # Set up covariates - only include age if available
      if (!is.null(age_var)) {
        covar_trait <- c(trait, "sex1", age_var)
      } else {
        covar_trait <- c(trait, "sex1")
      }
      
      # Run regression with binary trait and covariates
      res_trait <- std_fit_regression_noboot_SES(
        dat = dat, 
        x = trait,  # Binary predictor
        y = paste0(ses, "_std"),  # Standardized SES
        covar = covar_trait  # Includes predictor and covariates
      )
    } else {
      # For continuous traits: Pre-adjust for sex
      formula_trait <- reformulate("sex1", response = trait)
      dat[[paste0(trait, "_std")]] <- rstandard(lm(formula_trait, data = dat, na.action = na.exclude))
      
      # Run regression with residualized trait
      trait_var <- paste0(trait, "_std")
      res_trait <- std_fit_regression_noboot_SES(
        dat = dat, 
        x = trait_var,  # Residualized trait
        y = paste0(ses, "_std"),  # Standardized SES
        covar = trait_var  # Only include the trait variable as predictor
      )
    }
    
    # Store results
    results[[i]] <- data.frame(
      SES = ses,
      trait_y = pair[1],
      x_beta = res_trait[1], 
      beta_se_boot = res_trait[2], # I have to use the same name for the rma function above to run; this is no bootstrapping 
      ci_lower = res_trait[5],
      ci_upper = res_trait[6],
      P_value = res_trait[3],
      R2 = res_trait[4],
      Type = "SES_Trait_correlation", # same, I have to have the Type for rma function to run
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  
  return(final_results)
}



Full_List_PGS <- c("EA4_no23andme_Okbay2022", "IQ_Savage2018_FRCT1", 
                   "ChildhoodBMI_Vogelezang2020", "BMI_Giant2018", "Height_Yengo2022")

length(Full_List_PGS) # 5

# Function to calculate PGS-SES correlations for all specified PGS and SES variables
PGS_SES_corr_fun <- function(dat) {
  results <- list()
  
  # Pre-standardize all PCs once (outside the loop)
  for (i in 1:10) {
    pc_name <- paste0("PC", i)
    dat[[paste0(pc_name, "_std")]] <- as.numeric(scale(dat[[pc_name]]))
  }
  
  # Define all SES variables to analyze
  ses_vars <- c("ases", "gses", "pses", "u1pses")
  
  # Use the specified PGS list
  pgs_vars <- Full_List_PGS
  
  # For each PGS and each SES variable, calculate correlation
  for(pgs in pgs_vars) {
    for(ses in ses_vars) {
      cat("Processing item: ", pgs, ses, "\n")
      # Standardize SES variable
      dat[[paste0(ses, "_std")]] <- as.numeric(scale(dat[[ses]]))
      
      # For PGS: Pre-adjust for chiptype and standardized PCs
      std_pcs <- paste0("PC", 1:10, "_std")
      covar_pgs <- c("chiptype", std_pcs)
      formula_PGS <- reformulate(covar_pgs, response = pgs)
      dat[[paste0(pgs, "_std")]] <- rstandard(lm(formula_PGS, data = dat, na.action = na.exclude))
      
      # Run regression with residualized PGS
      pgs_var <- paste0(pgs, "_std")
      res_pgs <- std_fit_regression_noboot_SES(
        dat = dat, 
        x = pgs_var,  # Residualized PGS
        y = paste0(ses, "_std"),  # Standardized SES
        covar = pgs_var  # Only include the PGS variable as predictor
      )
      
      # Store results
      results[[paste(pgs, ses, sep="_")]] <- data.frame(
        SES = ses,
        trait_y = pgs,
        x_beta = res_pgs[1],
        beta_se_boot = res_pgs[2],
        ci_lower = res_pgs[5],
        ci_upper = res_pgs[6],
        P_value = res_pgs[3],
        R2 = res_pgs[4],
        Type = "SES_PGS_correlation",
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all results
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  
  return(final_results)
}

SES_trait_corr <- trait_SES_corr_fun(WFBF_selectunpaired)
SES_trait_corr$P_value <- p.adjust(SES_trait_corr$P_value, method = "fdr")
SES_trait_corr <- full_join(SES_trait_corr, Trait_Domain)
write.csv(SES_trait_corr, paste(outFileStem, "SEScorr_wTrait.csv", sep=""), row.names = FALSE)

# create a PGS_Domain
PGS_Domain <- as.data.frame(rbind(
  c("EA4_no23andme_Okbay2022", "ases", "EA4 at birth"),
  c("IQ_Savage2018_FRCT1", "ases", "IQ3 at birth"),
  c("ChildhoodBMI_Vogelezang2020", "ases", "Childhood BMI at birth"),
  c("BMI_Giant2018", "ases", "BMI at birth"),
  c("Height_Yengo2022", "ases", "Height at birth"),
  c("EA4_no23andme_Okbay2022", "gses", "EA4 at age 7"),
  c("IQ_Savage2018_FRCT1", "gses", "IQ3 at age 7"),
  c("ChildhoodBMI_Vogelezang2020", "gses", "Childhood BMI at age 7"),
  c("BMI_Giant2018", "gses", "BMI at age 7"),
  c("Height_Yengo2022", "gses", "Height at age 7"),
  c("EA4_no23andme_Okbay2022", "pses", "EA4 at age 16"),
  c("IQ_Savage2018_FRCT1", "pses", "IQ3 at age 16"),
  c("ChildhoodBMI_Vogelezang2020", "pses", "Childhood BMI at age 16"),
  c("BMI_Giant2018", "pses", "BMI at age 16"),
  c("Height_Yengo2022", "pses", "Height at age 16"),
  c("EA4_no23andme_Okbay2022", "u1pses", "EA4 at age 21"),
  c("IQ_Savage2018_FRCT1", "u1pses", "IQ3 at age 21"),
  c("ChildhoodBMI_Vogelezang2020", "u1pses", "Childhood BMI at age 21"),
  c("BMI_Giant2018", "u1pses", "BMI at age 21"),
  c("Height_Yengo2022", "u1pses", "Height at age 21")
))
dim(PGS_Domain) # 40*3
colnames(PGS_Domain) <- c("trait_y", "SES", "PGS_Domain")

SES_PGS_corr <- PGS_SES_corr_fun(WFBF_selectunpaired)
SES_PGS_corr$P_value <- p.adjust(SES_PGS_corr$P_value, method = "fdr")
SES_PGS_corr <- full_join(SES_PGS_corr, PGS_Domain)

SES_PGS_corr <- SES_PGS_corr %>%
  mutate(PGS_Domain = factor(PGS_Domain, levels = c("EA4 at birth", "EA4 at age 7", "EA4 at age 16", "EA4 at age 21", 
                                                    "IQ3 at birth", "IQ3 at age 7", "IQ3 at age 16", "IQ3 at age 21", 
                                                    "Childhood BMI at birth", "Childhood BMI at age 7", "Childhood BMI at age 16", "Childhood BMI at age 21", 
                                                    "BMI at birth", "BMI at age 7", "BMI at age 16", "BMI at age 21", 
                                                    "Height at birth", "Height at age 7", "Height at age 16", "Height at age 21")))

write.csv(SES_PGS_corr, paste(outFileStem, "SEScorr_wPGS.csv", sep=""), row.names = FALSE)





# meta-analysis
SEScorr_MA_wTrait <- run_meta_analysis(SES_trait_corr)

SEScorr_MA_wTrait <- SEScorr_MA_wTrait %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI")))%>%
  arrange(Type, V3, Trait_Domain)

# save the results 
write.csv(SEScorr_MA_wTrait, paste(outFileStem, "SEScorr_MA_wTrait.csv", sep=""), row.names = FALSE)

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, sheetName="")
writeData(wb, sheet=1, x=SEScorr_MA_wTrait, colNames=TRUE, rowNames=FALSE)
saveWorkbook(wb, paste(outFileStem, "SEScorr_MA_wTrait.xlsx", sep=""), overwrite = TRUE) 



# plot: trait 
# SEScorr_MA_wTrait <- read.csv(paste(outFileStem, "SEScorr_MA_wTrait.csv", sep=""))
MA_SES_corr_coghtBMI <- SEScorr_MA_wTrait %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI")))

MA_SES_corr_coghtBMI_plot <- ggplot(MA_SES_corr_coghtBMI, aes(y = beta, x = Trait_Domain)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.85, fill = "#0170c0") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.6, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.6, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.6, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.6, label = "BMI & Height", fontface = "bold", size = 7) +
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Standardized Beta Estimates", fill = "") +
  scale_y_continuous(limits=c(-0.2, 0.6), breaks=seq(-0.2, 0.5, 0.1)) + 
  #scale_color_manual(values=c("Category1" = "#FF4040", "Category2" = "#009ACD", "Category3" = "#00FA9A")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

MA_SES_corr_coghtBMI_plot
ggsave("MA_SES_corr_coghtBMI_plot.jpg", plot = MA_SES_corr_coghtBMI_plot, width = 15, height = 10)



# plot: PGS
# SES_PGS_corr <- read.csv(paste(outFileStem, "SEScorr_wPGS.csv", sep=""))

# SES_PGS_corr <- SES_PGS_corr %>%
#   mutate(PGS_Domain = factor(PGS_Domain, levels = c("EA4 at birth", "EA4 at age 7", "EA4 at age 16", "EA4 at age 21", 
#                                                    "IQ3 at birth", "IQ3 at age 7", "IQ3 at age 16", "IQ3 at age 21", 
#                                                    "Childhood BMI at birth", "Childhood BMI at age 7", "Childhood BMI at age 16", "Childhood BMI at age 21", 
#                                                    "BMI at birth", "BMI at age 7", "BMI at age 16", "BMI at age 21", 
#                                                    "Height at birth", "Height at age 7", "Height at age 16", "Height at age 21")))

asterisk_SES_PGS_corr <- SES_PGS_corr %>%
  group_by(PGS_Domain) %>%
  mutate(
    asterisk = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

position_data <- SES_PGS_corr %>%
  group_by(PGS_Domain) %>%
  summarise(
    max_ci_upper = max(ci_upper),
    min_ci_lower = min(ci_lower),
    beta_value = unique(x_beta)  # Assuming x_beta is the same within each PGS_Domain
  )

asterisk_SES_PGS_corr <- asterisk_SES_PGS_corr %>%
  left_join(position_data, by = "PGS_Domain") %>%
  # Determine asterisk position based on whether beta is positive or negative
  mutate(
    asterisk_y_position = case_when(
      beta_value >= 0 ~ max_ci_upper + 0.05,  # Above error bar for positive values
      beta_value < 0 ~ min_ci_lower - 0.05    # Below error bar for negative values
    )
  )

MA_SES_PGS_plot <- ggplot(SES_PGS_corr, aes(y = x_beta, x = PGS_Domain)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.85, fill = "#0170c0") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  
  geom_text(
    data = asterisk_SES_PGS_corr,
    aes(
      x = PGS_Domain, 
      y = asterisk_y_position,  # Position depends on whether beta is positive or negative
      label = asterisk
    ),
    inherit.aes = FALSE,
    size = 8,
    color = "black"
  ) +
  
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Standardized Beta Estimates", fill = "") +
  scale_y_continuous(
    # Adjusted limits to accommodate both positive and negative asterisks
    limits = c(min(SES_PGS_corr$ci_lower) - 0.1, max(SES_PGS_corr$ci_upper) + 0.1), 
    breaks = seq(-0.2, 0.5, 0.1)
  ) + 
  #scale_color_manual(values=c("Category1" = "#FF4040", "Category2" = "#009ACD", "Category3" = "#00FA9A")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

MA_SES_PGS_plot
ggsave("MA_SES_PGS_plot.jpg", plot = MA_SES_PGS_plot, width = 15, height = 10)


