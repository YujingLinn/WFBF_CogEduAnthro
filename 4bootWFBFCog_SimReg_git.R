#!/usr/bin/Rscript
#.libPaths(c("/users/k21170717/.conda/envs/r_env/lib/R/library", .libPaths()))

# title: Simple Regression Approach with Bootstrap & Meta-Analysis
# author: Yujing Lin
# date: 6th October, 2024
# codes corresponding to publication: https://icajournal.scholasticahq.com/article/140654-polygenic-score-prediction-within-and-between-sibling-pairs-for-intelligence-cognitive-abilities-and-educational-traits-from-childhood-to-early-adul

# instruction: this is the simple regression approach with bootstrap and meta-analysis

nboot <- 10 # change it to 1000 on HPC
ncpus <- 4 # change it to 32 on HPC

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/Data/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'
outFileStemRDS <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/RDS/'
#sourceFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/data_160425/'
#outFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/results_160425/'
#outFileStemRDS <- '/scratch/users/k21170717/WFBF/Cognitive/resultsRDS_160425/'

# libraries ####
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

# datasets ####
WFBF_selectunpaired <- read.csv(paste(sourceFileStem, "WFBF_selectunpaired.csv", sep=""))
dat_pair <- read.csv(paste(sourceFileStem, "dat_pair.csv", sep=""))

WFBF_selectunpaired_F <- WFBF_selectunpaired %>% 
  filter(sex1==0) # unrelated females
WFBF_selectunpaired_M <- WFBF_selectunpaired %>% 
  filter(sex1==1) # unrelated males
WFBF_selectunpaired_MZ <- WFBF_selectunpaired %>% 
  filter(zygos==1) # MZ twins
WFBF_selectunpaired_DZ <- WFBF_selectunpaired %>% 
  filter(zygos==2) # DZ twins

dat_pair_selectunpaired <- subset(dat_pair, selectunpaired==1)
dat_DZpair_selectunpaired <- subset(dat_pair_selectunpaired, dat_pair_selectunpaired$zygos==2)
dat_DZpair_selectunpaired_FSS <- dat_DZpair_selectunpaired %>% 
  filter(x3zygos==2 & sex1==0) # same-sex female DZ twins
dat_DZpair_selectunpaired_MSS <- dat_DZpair_selectunpaired %>% 
  filter(x3zygos==2 & sex1==1) # same-sex male DZ twins
dat_DZpair_selectunpaired_OS <- dat_DZpair_selectunpaired %>% 
  filter(x3zygos==3) # opposite-sex DZ twins

PCs <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")



# variable list ####
VarPairs <- list( # it has to be a list, otherwise, each row will have auto-refill to match the maximum length of 4, i.e., repeating the continuous variables twice
  c("gt2ac", "EA4_no23andme_Okbay2022"),
  c("gcg", "IQ_Savage2018_FRCT1"),
  c("gcl", "IQ_Savage2018_FRCT1"),
  c("gcn", "IQ_Savage2018_FRCT1"),
  c("gbmi", "ChildhoodBMI_Vogelezang2020"), # The childhood GWAS participants were aged 2 to 10
  c("ghtcm", "Height_Yengo2022"), # 7
  
  c("it3ac", "EA4_no23andme_Okbay2022"),
  c("icg", "IQ_Savage2018_FRCT1"),
  c("icvb", "IQ_Savage2018_FRCT1"),
  c("icnv", "IQ_Savage2018_FRCT1"), # 9
  
  c("jt3ac", "EA4_no23andme_Okbay2022"),
  c("jcg", "IQ_Savage2018_FRCT1"),
  c("jcvb", "IQ_Savage2018_FRCT1"),
  c("jcnv", "IQ_Savage2018_FRCT1"), # 10
  
  c("lt3ac", "EA4_no23andme_Okbay2022"),
  c("lcg", "IQ_Savage2018_FRCT1"),
  # c("lcvb", "IQ_Savage2018_FRCT1"),
  # c("lcnv", "IQ_Savage2018_FRCT1"),
  c("lverbal12T", "IQ_Savage2018_FRCT1"), 
  c("lnonverbal12T", "IQ_Savage2018_FRCT1"),
  c("lbmi", "BMI_Giant2018"), 
  c("lchtcm", "Height_Yengo2022"), # 12
  
  c("ncbmi", "BMI_Giant2018"),
  c("nchtcm", "Height_Yengo2022"), # 14
  
  c("pcexgcsecoregrdm", "EA4_no23andme_Okbay2022"),
  c("pcg", "IQ_Savage2018_FRCT1"),
  c("pcvctota", "IQ_Savage2018_FRCT1"),
  c("pcrvtota", "IQ_Savage2018_FRCT1"),
  c("pcbmi", "BMI_Giant2018"),
  c("pcqdhtcm", "Height_Yengo2022"), # 16
  
  c("ra_level_enrol", "EA4_no23andme_Okbay2022", "rcqalage1"), # binary
  c("rcqalgrdm", "EA4_no23andme_Okbay2022"),
  # c("rcqalsgrdm", "EA4_no23andme_Okbay2022"), # A- & AS-level
  c("rcqhe", "EA4_no23andme_Okbay2022", "rcqalage1"), # 18 # binary
  
  c("u1cdegr1", "EA4_no23andme_Okbay2022"), # degree description at age 21
  # c("u1cedat", "EA4_no23andme_Okbay2022"), # twin EA composite 
  c("u1cbmi", "BMI_Giant2018"),
  c("u1chtcm", "Height_Yengo2022"),
  c("ucgt", "IQ_Savage2018_FRCT1"),
  c("ucgvbt", "IQ_Savage2018_FRCT1"),
  c("ucgnvt", "IQ_Savage2018_FRCT1"), # 21
  
  c("zEA", "EA4_no23andme_Okbay2022"),
  c("zmhbmi", "BMI_Giant2018"),
  c("zmhheight", "Height_Yengo2022") # 26
)

length(VarPairs) # 40




# the functions ####
# applicable for both continuous and binary variables, no need for an extra step of std beta for binary variables 
# the results for binary variables are comparable to biserial.cor
std_fit_regression_boot <- function(dat, indices, x, y, covar) {
  # datx <- dat
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
  # datx <- dat
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
beta_linear_results_fun_boot <- function(dat, paste_item_PGS, paste_item_trait, sex_var, nboot, ncpus, rds_name) {
  #dat <- dat_DZpair_selectunpaired # WFBF_selectunpaired   dat_DZpair_selectunpaired   dat_DZpair_selectunpaired_FSS
  
  beta_linear_raw <- lapply(1:length(VarPairs), function(i) {
    cat("Processing item", i, ":", VarPairs[[i]], "\n")
    #i <- 29 # 29: binary
    #paste_item_trait <- "_trait_diff" # "_trait_diff", "1", "_trait_sum"
    #paste_item_PGS <- "_PGS_diff" # "_PGS_diff", "", "_PGS_sum"
    #sex_var <- "sex1_label_for_pairdiff" # "sex1_label_for_pairdiff", "sex1_label", "sex1_label_for_pairmeansum"
    
    # to std the continuous covar (have the std for the outcome trait below)
    if (length(VarPairs[[i]]) == 3) { # binary
      continuous_vars <- c(
        paste(VarPairs[[i]][2], paste_item_PGS, sep = ""),  # PGS variable
        VarPairs[[i]][3],                                   # Age variable
        paste0("PC", 1:10)                                  # PCs 
      )
    } else {
      continuous_vars <- c(
        paste(VarPairs[[i]][2], paste_item_PGS, sep = ""),  # PGS variable
        paste0("PC", 1:10)                                  # PCs
      )
    }
    
    # Standardize them all at once (skip binary/categorical variables like sex)
    # for some reasons, this function sometimes may not work probably due to conflicts among the packages
    # if that happens, just close Rstudio and restart it
    dat <- dat %>%
      mutate(across(
        .cols = all_of(continuous_vars),
        .fns = ~ as.numeric(scale(.)),  # Standardize (mean=0, SD=1)
        .names = "{.col}_std"            # Append "_std" to standardized vars
      ))
    
    y <- paste(VarPairs[[i]][1], paste_item_trait, sep = "")  # Trait
    x <- paste0(paste(VarPairs[[i]][2], paste_item_PGS, sep = ""), "_std")  # PGS (standardized)
    if (!is.na(VarPairs[[i]][3])) {
      age <- paste0(VarPairs[[i]][3], "_std")  # Use standardized version
    } else {
      age <- NULL  # Or handle missing age appropriately
    }
    PCs <- paste0("PC", 1:10, "_std")        # Standardized PCs
    
    # Covariates (now using standardized versions)
    covar <- c(x, age, sex_var, "chiptype", PCs)  # Sex and chiptype remains unstandardized
    covar <- covar[!is.na(covar)]
    # variable_name_for_label <- paste(VarPairs[[i]][1], 1, sep="")
    
    is_binary <- length(VarPairs[[i]]) == 3
    if (length(unique(dat[[sex_var]])) == 1) { # Remove sex if only one level exists
      covar <- covar[covar != sex_var]  # Drop sex if only 1 level (sex-stratified)
      
      if (!is_binary) {
        # don't forget to std the continuous outcome trait y in the sex-stratified analyses
        dat[[y]] <- scale(dat[[y]])
      }
      
    } else {
      # Use VarPairs length to determine if y is binary (3) or continuous (2)
      
      
      if (!is_binary) {  # Continuous trait (VarPairs[[i]] length = 2)
        # Correct y for sex and remove sex_var from covariates
        formula_sex <- reformulate(sex_var, response = y)
        dat[[y]] <- rstandard(lm(formula_sex, data = dat, na.action = na.exclude))
        covar <- covar[covar != sex_var]  
      }
      # Binary trait (VarPairs[[i]] length = 3): sex_var remains in covar & no std for the outcome trait 
    }
    
    boot_results <- boot(data = dat, 
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
  names(boot_dist_list) <- sapply(VarPairs, function(v) paste0(v[1]))
  
  saveRDS(
    boot_dist_list,
    file = paste0(outFileStemRDS, "boot_dist_", rds_name, ".rds")
  )
  
  results_list <- lapply(1:length(beta_linear_raw), function(i) {
    # i <- 30
    boot_results <- beta_linear_raw[[i]]
    
    trait_y <- VarPairs[[i]][1]
    PGS_x <- VarPairs[[i]][2]
    
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
# 1. Population-Based Correlations--unrelated samples: trait~PGS w. std. beta est. ####
# [WFBF_selectunpaired] N = 6973
cat("ppl_unrelated_selectunpaired_boot \n")
ppl_unrelated_selectunpaired_boot <- beta_linear_results_fun_boot(WFBF_selectunpaired, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_boot") # suffix for PGS, suffix for trait
ppl_unrelated_selectunpaired_boot$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_boot$p_x_beta, method="fdr")

# 2. Population-Based Correlations--unrelated samples: trait pair sum ~ PGS pair sum ####
# pair sum, conceptually = pair mean 
# even though we could use DZ and MZ, to align with pair diff, we only used DZ twins
cat("ppl_DZpairsum_selectunpaired_boot \n")
ppl_DZpairsum_selectunpaired_boot <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired, "_PGS_sum", "_trait_sum", "sex1_label_for_pairmeansum", nboot = nboot, ncpus = ncpus, "ppl_DZpairsum_selectunpaired_boot") # suffix for PGS, suffix for trait
ppl_DZpairsum_selectunpaired_boot$p_x_beta <- p.adjust(ppl_DZpairsum_selectunpaired_boot$p_x_beta, method="fdr")

cat("ppl_DZpairmean_selectunpaired_boot \n")
ppl_DZpairmean_selectunpaired_boot <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired, "_PGS_mean", "_trait_mean", "sex1_label_for_pairmeansum", nboot = nboot, ncpus = ncpus, "ppl_DZpairmean_selectunpaired_boot") # suffix for PGS, suffix for trait
ppl_DZpairmean_selectunpaired_boot$p_x_beta <- p.adjust(ppl_DZpairmean_selectunpaired_boot$p_x_beta, method="fdr")

# 3. Within-Family Correlations--DZ pairs: trait pair diff ~ PGS pair diff ####
# pair diff, conceptually = deviation from pair mean
# caution: technically, the variation of deviation from pair mean is half of the variation of pair diff
# this will lead to the slope estimates to be twice of pair diff
# but since we transform the PGS and trait simultaneously, the slope estimates are unaffected this way
# plue, we standardise the results
cat("WF_DZpairdiff_boot \n")
WF_DZpairdiff_boot <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZpairdiff_boot") # suffix for PGS, suffix for trait
WF_DZpairdiff_boot$p_x_beta <- p.adjust(WF_DZpairdiff_boot$p_x_beta, method="fdr")

cat("WF_DZmeanpairdiff_boot \n")
WF_DZmeanpairdiff_boot <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired, "_PGS_mean_diff", "_trait_mean_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZmeanpairdiff_boot") # suffix for PGS, suffix for trait
WF_DZmeanpairdiff_boot$p_x_beta <- p.adjust(WF_DZmeanpairdiff_boot$p_x_beta, method="fdr")



# aggregate the simple regression results ####
ppl_unrelated_selectunpaired_boot$Type <- "ppl_unrelated"
ppl_DZpairsum_selectunpaired_boot$Type <- "ppl_DZpairsum"
ppl_DZpairmean_selectunpaired_boot$Type <- "ppl_DZpairmean"
WF_DZpairdiff_boot$Type <- "WF_DZpairdiff"
WF_DZmeanpairdiff_boot$Type <- "WF_DZmeanpairdiff"

raw_results_boot <- full_join(ppl_unrelated_selectunpaired_boot, ppl_DZpairsum_selectunpaired_boot) %>% 
  full_join(ppl_DZpairmean_selectunpaired_boot) %>% 
  full_join(WF_DZpairdiff_boot) %>% 
  full_join(WF_DZmeanpairdiff_boot)





# compare bootstrapped distributions ####
cat("compare bootstrapped distributions \n")
boot_PPL <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_boot.rds"))
names(boot_PPL)
boot_WF <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_boot.rds"))
names(boot_WF)

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

boot_compare_SimReg_popWF <- fun_compare_boot_dists(boot_PPL, boot_WF)



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

raw_results_boot <- full_join(raw_results_boot, Trait_Domain)
boot_compare_SimReg_popWF <- full_join(boot_compare_SimReg_popWF, Trait_Domain)



# save the pre-meta-analysed results ####
raw_results_boot <- raw_results_boot %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("WF_DZpairdiff","WF_DZmeanpairdiff", "ppl_unrelated", "ppl_DZpairsum", "ppl_DZpairmean"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, V3, Trait_Domain)

boot_compare_SimReg_popWF <- boot_compare_SimReg_popWF %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(V3, Trait_Domain)

write.csv(raw_results_boot, paste(outFileStem, "WFBF_Cog_Main_boot_raw.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_popWF, paste(outFileStem, "boot_compare_WFBF_Cog_Main_raw.csv", sep=""), row.names = FALSE)



# raw_results_boot <- read.csv('/Users/yujinglin/Desktop/WFBF Results 140624/4out_raw_results_boot.csv')
# meta-analyses ####
cat("meta-analyses for simple regression \n")
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
                      #r2i = current_data$R.squared,####
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
      #R.square = rma_result$R2,####
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

results_boot_MA <- run_meta_analysis(raw_results_boot)



# save the results ####
write.csv(results_boot_MA, paste(outFileStem, "MA_WFBF_Cog_Main_boot.csv", sep=""), row.names = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="")
writeData(wb, sheet=1, x=results_boot_MA, colNames=TRUE, rowNames=FALSE)
saveWorkbook(wb, paste(outFileStem, "MA_WFBF_Cog_Main_boot.xlsx", sep=""), overwrite = TRUE) 



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

bootMA_compare_SimReg_popWF <- run_meta_diff_analysis(boot_compare_SimReg_popWF)

# save the results ####
write.csv(bootMA_compare_SimReg_popWF, paste(outFileStem, "bootMA_compare_WFBF_Cog_Main.csv", sep=""), row.names = FALSE)








# sensitivity analyses ####
cat("sensitivity analyses \n")
ppl_unrelated_selectunpaired_F_boot <- beta_linear_results_fun_boot(WFBF_selectunpaired_F, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_F_boot") # suffix for PGS, suffix for trait
ppl_unrelated_selectunpaired_F_boot$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_F_boot$p_x_beta, method="fdr")

ppl_unrelated_selectunpaired_M_boot <- beta_linear_results_fun_boot(WFBF_selectunpaired_M, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_M_boot") 
ppl_unrelated_selectunpaired_M_boot$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_M_boot$p_x_beta, method="fdr")

WF_DZpairdiff_FSS <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired_FSS, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZpairdiff_FSS") 
WF_DZpairdiff_FSS$p_x_beta <- p.adjust(WF_DZpairdiff_FSS$p_x_beta, method="fdr")

WF_DZpairdiff_MSS <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired_MSS, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZpairdiff_MSS") 
WF_DZpairdiff_MSS$p_x_beta <- p.adjust(WF_DZpairdiff_MSS$p_x_beta, method="fdr")

WF_DZpairdiff_OS <- beta_linear_results_fun_boot(dat_DZpair_selectunpaired_OS, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff", nboot = nboot, ncpus = ncpus, "WF_DZpairdiff_OS") 
WF_DZpairdiff_OS$p_x_beta <- p.adjust(WF_DZpairdiff_OS$p_x_beta, method="fdr")

ppl_unrelated_selectunpaired_F_boot$Type <- "ppl_unrelated"
ppl_unrelated_selectunpaired_M_boot$Type <- "ppl_unrelated"
WF_DZpairdiff_FSS$Type <- "WF_DZpairdiff"
WF_DZpairdiff_MSS$Type <- "WF_DZpairdiff"
WF_DZpairdiff_OS$Type <- "WF_DZpairdiff"

ppl_unrelated_selectunpaired_F_boot$Subsample <- "Female"
ppl_unrelated_selectunpaired_M_boot$Subsample <- "Male"
WF_DZpairdiff_FSS$Subsample <- "Same-Sex Female"
WF_DZpairdiff_MSS$Subsample <- "Same-Sex Male"
WF_DZpairdiff_OS$Subsample <- "Opposite-Sex"

raw_sensitivity_boot <- full_join(ppl_unrelated_selectunpaired_F_boot, ppl_unrelated_selectunpaired_M_boot) %>% 
  full_join(WF_DZpairdiff_FSS) %>% 
  full_join(WF_DZpairdiff_MSS) %>% 
  full_join(WF_DZpairdiff_OS)

raw_sensitivity_boot <- full_join(raw_sensitivity_boot, Trait_Domain)
write.csv(raw_sensitivity_boot , paste(outFileStem, "WFBF_Cog_Main_sensi_boot_raw.csv", sep=""), row.names = FALSE)



cat("compare bootstrapped distributions for sensitivity analysis \n")
boot_PPL_sensi_F <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_F_boot.rds"))
boot_PPL_sensi_M <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_M_boot.rds"))
boot_WF_sensi_FSS <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_FSS.rds"))
boot_WF_sensi_MSS <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_MSS.rds"))
boot_WF_sensi_OS <- readRDS(paste0(outFileStemRDS, "boot_dist_WF_DZpairdiff_OS.rds"))

boot_compare_SimReg_sensi_popFM <- fun_compare_boot_dists(boot_PPL_sensi_F, boot_PPL_sensi_M)
boot_compare_SimReg_sensi_popFM <- full_join(boot_compare_SimReg_sensi_popFM, Trait_Domain)

boot_compare_SimReg_sensi_WF_FMSS <- fun_compare_boot_dists(boot_WF_sensi_FSS, boot_WF_sensi_MSS)
boot_compare_SimReg_sensi_WF_FMSS <- full_join(boot_compare_SimReg_sensi_WF_FMSS, Trait_Domain)

boot_compare_SimReg_sensi_WF_FSSOS <- fun_compare_boot_dists(boot_WF_sensi_FSS, boot_WF_sensi_OS)
boot_compare_SimReg_sensi_WF_FSSOS <- full_join(boot_compare_SimReg_sensi_WF_FSSOS, Trait_Domain)

boot_compare_SimReg_sensi_WF_MSSOS <- fun_compare_boot_dists(boot_WF_sensi_MSS, boot_WF_sensi_OS)
boot_compare_SimReg_sensi_WF_MSSOS <- full_join(boot_compare_SimReg_sensi_WF_MSSOS, Trait_Domain)

write.csv(boot_compare_SimReg_sensi_popFM, paste(outFileStem, "boot_compare_WFBF_Cog_sensi_popFM_raw.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_sensi_WF_FMSS, paste(outFileStem, "boot_compare_WFBF_Cog_sensi_WF_FMSS_raw.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_sensi_WF_FSSOS, paste(outFileStem, "boot_compare_WFBF_Cog_sensi_WF_FSSOS_raw.csv", sep=""), row.names = FALSE)
write.csv(boot_compare_SimReg_sensi_WF_MSSOS, paste(outFileStem, "boot_compare_WFBF_Cog_sensi_WF_MSSOS_raw.csv", sep=""), row.names = FALSE)








cat("sensitivity analyses meta-analysis \n")
#raw_sensitivity_boot <- read.csv('/Users/yujinglin/Desktop/4out_raw_sensitivity_boot.csv')
raw_sensitivity_bootF <- subset(raw_sensitivity_boot, Subsample == "Female")
sensitivity_bootF_MA <- run_meta_analysis(raw_sensitivity_bootF)
sensitivity_bootF_MA$Subsample <- "Female"

raw_sensitivity_bootM <- subset(raw_sensitivity_boot, Subsample == "Male")
sensitivity_bootM_MA <- run_meta_analysis(raw_sensitivity_bootM)
sensitivity_bootM_MA$Subsample <- "Male"

raw_sensitivity_bootFSS <- subset(raw_sensitivity_boot, Subsample == "Same-Sex Female")
sensitivity_bootFSS_MA <- run_meta_analysis(raw_sensitivity_bootFSS)
sensitivity_bootFSS_MA$Subsample <- "Same-Sex Female"

raw_sensitivity_bootMSS <- subset(raw_sensitivity_boot, Subsample == "Same-Sex Male")
sensitivity_bootMSS_MA <- run_meta_analysis(raw_sensitivity_bootMSS)
sensitivity_bootMSS_MA$Subsample <- "Same-Sex Male"

raw_sensitivity_bootOS <- subset(raw_sensitivity_boot, Subsample == "Opposite-Sex")
sensitivity_bootOS_MA <- run_meta_analysis(raw_sensitivity_bootOS)
sensitivity_bootOS_MA$Subsample <- "Opposite-Sex"

sensitivity_boot_MA <- full_join(sensitivity_bootF_MA, sensitivity_bootM_MA) %>% 
  full_join(sensitivity_bootFSS_MA) %>% 
  full_join(sensitivity_bootMSS_MA) %>% 
  full_join(sensitivity_bootOS_MA)

write.csv(sensitivity_boot_MA , paste(outFileStem, "MA_WFBF_Cog_Main_sensi_boot.csv", sep=""), row.names = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="")
writeData(wb, sheet=1, x=sensitivity_boot_MA, colNames=TRUE, rowNames=FALSE)
saveWorkbook(wb, paste(outFileStem, "MA_WFBF_Cog_Main_sensi_boot.xlsx", sep=""), overwrite = TRUE) 



bootMA_compare_SimReg_sensi_popFM <- run_meta_diff_analysis(boot_compare_SimReg_sensi_popFM)
bootMA_compare_SimReg_sensi_WF_FMSS <- run_meta_diff_analysis(boot_compare_SimReg_sensi_WF_FMSS)
bootMA_compare_SimReg_sensi_WF_FSSOS <- run_meta_diff_analysis(boot_compare_SimReg_sensi_WF_FSSOS)
bootMA_compare_SimReg_sensi_WF_MSSOS <- run_meta_diff_analysis(boot_compare_SimReg_sensi_WF_MSSOS)

# save the results ####
write.csv(bootMA_compare_SimReg_sensi_popFM, paste(outFileStem, "bootMA_compare_WFBF_Cog_sensi_popFM.csv", sep=""), row.names = FALSE)
write.csv(bootMA_compare_SimReg_sensi_WF_FMSS, paste(outFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_FMSS.csv", sep=""), row.names = FALSE)
write.csv(bootMA_compare_SimReg_sensi_WF_FSSOS, paste(outFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_FSSOS.csv", sep=""), row.names = FALSE)
write.csv(bootMA_compare_SimReg_sensi_WF_MSSOS, paste(outFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_MSSOS.csv", sep=""), row.names = FALSE)









# additional sensitivity analysis based on zygosity: POP est ####
cat("additional sensitivity analysis based on zygosity: POP est \n")

ppl_unrelated_selectunpaired_MZ_boot <- beta_linear_results_fun_boot(WFBF_selectunpaired_MZ, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_MZ_boot") # suffix for PGS, suffix for trait
ppl_unrelated_selectunpaired_MZ_boot$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_MZ_boot$p_x_beta, method="fdr")
ppl_unrelated_selectunpaired_MZ_boot <- full_join(ppl_unrelated_selectunpaired_MZ_boot, Trait_Domain)
ppl_unrelated_selectunpaired_MZ_boot$Type <- "ppl_unrelated_MZ"

ppl_unrelated_selectunpaired_MZ_boot_MA <- run_meta_analysis(ppl_unrelated_selectunpaired_MZ_boot)

ppl_unrelated_selectunpaired_DZ_boot <- beta_linear_results_fun_boot(WFBF_selectunpaired_DZ, "", "1", "sex1_label", nboot = nboot, ncpus = ncpus, "ppl_unrelated_selectunpaired_DZ_boot") 
ppl_unrelated_selectunpaired_DZ_boot$p_x_beta <- p.adjust(ppl_unrelated_selectunpaired_DZ_boot$p_x_beta, method="fdr")
ppl_unrelated_selectunpaired_DZ_boot <- full_join(ppl_unrelated_selectunpaired_DZ_boot, Trait_Domain)
ppl_unrelated_selectunpaired_DZ_boot$Type <- "ppl_unrelated_DZ"

ppl_unrelated_selectunpaired_DZ_boot_MA <- run_meta_analysis(ppl_unrelated_selectunpaired_DZ_boot)



ppl_unrelated_selectunpaired_MZDZ_boot <- rbind(ppl_unrelated_selectunpaired_MZ_boot, ppl_unrelated_selectunpaired_DZ_boot)
ppl_unrelated_selectunpaired_MZDZ_boot_MA <- rbind(ppl_unrelated_selectunpaired_MZ_boot_MA, ppl_unrelated_selectunpaired_DZ_boot_MA)

write.csv(ppl_unrelated_selectunpaired_MZDZ_boot , paste(outFileStem, "WFBF_Cog_Main_MZDZsensi_boot_raw.csv", sep=""), row.names = FALSE)
write.csv(ppl_unrelated_selectunpaired_MZDZ_boot_MA , paste(outFileStem, "MA_WFBF_Cog_Main_MZDZsensi_boot.csv", sep=""), row.names = FALSE)


wb <- createWorkbook()
addWorksheet(wb, sheetName="")
writeData(wb, sheet=1, x=ppl_unrelated_selectunpaired_MZDZ_boot_MA, colNames=TRUE, rowNames=FALSE)
saveWorkbook(wb, paste(outFileStem, "MA_WFBF_Cog_Main_MZDZsensi_boot.xlsx", sep=""), overwrite = TRUE) 



cat("compare bootstrapped distributions for zygosity sensitivity \n")
boot_PPL_MZ <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_MZ_boot.rds"))
boot_PPL_DZ <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_DZ_boot.rds"))

boot_compare_SimReg_popMZDZ <- fun_compare_boot_dists(boot_PPL_MZ, boot_PPL_DZ)

boot_compare_SimReg_popMZDZ <- full_join(boot_compare_SimReg_popMZDZ, Trait_Domain)
write.csv(boot_compare_SimReg_popMZDZ, paste(outFileStem, "boot_compare_WFBF_Cog_MZDZsensi_raw.csv", sep=""), row.names = FALSE)

bootMA_compare_SimReg_popMZDZ <- run_meta_diff_analysis(boot_compare_SimReg_popMZDZ)
write.csv(bootMA_compare_SimReg_popMZDZ, paste(outFileStem, "bootMA_compare_WFBF_Cog_MZDZsensi.csv", sep=""), row.names = FALSE)


# save all boot compare raw results to one excel file ####
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="main_PopWF")
addWorksheet(wb, sheetName="PopFM")
addWorksheet(wb, sheetName="PopMZDZ")
addWorksheet(wb, sheetName="WF_FMSS")
addWorksheet(wb, sheetName="WF_FSSOS")
addWorksheet(wb, sheetName="WF_MSSOS")

# Write each table to its corresponding worksheet
writeData(wb, sheet="main_PopWF", x=boot_compare_SimReg_popWF, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="PopFM", x=boot_compare_SimReg_sensi_popFM, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="PopMZDZ", x=boot_compare_SimReg_popMZDZ, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_FMSS", x=boot_compare_SimReg_sensi_WF_FMSS, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_FSSOS", x=boot_compare_SimReg_sensi_WF_FSSOS, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_MSSOS", x=boot_compare_SimReg_sensi_WF_MSSOS, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "boot_compare_WFBF_Cog_raw_All.xlsx", sep=""), overwrite=TRUE)


# save all boot MA compare raw results to one excel file ####
wb <- createWorkbook()

# Add a worksheet for each table with descriptive names
addWorksheet(wb, sheetName="main_PopWF")
addWorksheet(wb, sheetName="PopFM")
addWorksheet(wb, sheetName="PopMZDZ")
addWorksheet(wb, sheetName="WF_FMSS")
addWorksheet(wb, sheetName="WF_FSSOS")
addWorksheet(wb, sheetName="WF_MSSOS")

# Write each table to its corresponding worksheet
writeData(wb, sheet="main_PopWF", x=bootMA_compare_SimReg_popWF, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="PopFM", x=bootMA_compare_SimReg_sensi_popFM, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="PopMZDZ", x=bootMA_compare_SimReg_popMZDZ, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_FMSS", x=bootMA_compare_SimReg_sensi_WF_FMSS, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_FSSOS", x=bootMA_compare_SimReg_sensi_WF_FSSOS, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_MSSOS", x=bootMA_compare_SimReg_sensi_WF_MSSOS, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "bootMA_compare_WFBF_Cog_All.xlsx", sep=""), overwrite=TRUE)



cat("done \n")




# additional: compare population estimates using unrelated population vs pair sums ####
# didn't incorporate in the HPC analyses, ran separately afterwards 
boot_PPL <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_unrelated_selectunpaired_boot.rds"))
boot_PPL_pairsum <- readRDS(paste0(outFileStemRDS, "boot_dist_ppl_DZpairsum_selectunpaired_boot.rds"))

boot_compare_SimReg_Methods_unrelated_pairsum <- fun_compare_boot_dists(boot_PPL, boot_PPL_pairsum)

boot_compare_SimReg_Methods_unrelated_pairsum <- full_join(boot_compare_SimReg_Methods_unrelated_pairsum, Trait_Domain)

boot_compare_SimReg_Methods_unrelated_pairsum <- boot_compare_SimReg_Methods_unrelated_pairsum %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(V3, Trait_Domain)

write.csv(boot_compare_SimReg_Methods_unrelated_pairsum, paste(outFileStem, "boot_compare_SimReg_Methods_unrelated_pairsum.csv", sep=""), row.names = FALSE)




