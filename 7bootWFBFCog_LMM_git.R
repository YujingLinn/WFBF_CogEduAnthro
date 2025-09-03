#!/usr/bin/Rscript
#.libPaths(c("/scratch/users/k21170717/PGS_pipeline_LDPred2/miniconda3/envs/myenv/lib/R/library", .libPaths()))

# title: Mixed-Effects Modeling using full DZ Twins Pairs with Bootstrap & Meta-Analysis
# author: Yujing Lin
# date: 11th October, 2024
# codes corresponding to publication: https://icajournal.scholasticahq.com/article/140654-polygenic-score-prediction-within-and-between-sibling-pairs-for-intelligence-cognitive-abilities-and-educational-traits-from-childhood-to-early-adul

# mixed-effects model: trait ~ pop-PGS + random family effect + error
# mixed-effects model: trait ~ WF-pair mean diff + BF-pair mean + random family effect + error

nboot <- 3 # change it to 1000 on HPC
ncpus <- 4 # change it to 64 on HPC

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/Data/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'
#sourceFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/data_160425/'
#outFileStem <- '/scratch/users/k21170717/WFBF/Cognitive/results_160425/'

dat_pair <- read.csv(paste(sourceFileStem, "dat_pair.csv", sep=""))

library(lme4) # lmer: fit linear mixed-effects models, via REML or maximum likelihood 
library(nlme) # lme: linear mixed-effects models (mainly for ICC), based on the paper (Lindstrom & Bates, 1990)
library(dplyr)
library(tidyverse)
library(MuMIn) # r.squaredGLMM
library(lmerTest)
library(performance) # R2
library(metafor) # meta-analysis
library(boot)
library(openxlsx)
library(glmmTMB)

df <- dat_pair %>% filter(zygos == 2) 
dim(df) # 17293*322

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



# LMM function
# to check (4 conditions): cont variable for within-family analysis, population analysis, binary variable for within-family analysis, population analysis
process_LLM_fun_boot <- function(df, paste_item_trait, paste_item_PGS1, paste_item_PGS2, sex_var) {
  results <- NULL
  
  for (i in 1:length(VarPairs)) {
    #i <- 29 # 29: binary
    cat("Processing item", i, ":", VarPairs[[i]], "\n")
    is_binary <- length(VarPairs[[i]]) == 3
    
    #paste_item_trait <- "1"
    #paste_item_PGS1 <- "_PGS_mean_diff" # "_PGS_mean_diff", "1"
    #paste_item_PGS2 <- "_PGS_mean" # "_PGS_mean", NA
    #sex_var <- "sex1_label"
    
    cont_vars <- c(
      if (!is_binary) paste0(VarPairs[[i]][1], paste_item_trait), # for continuous phenotypes
      paste0(VarPairs[[i]][2], paste_item_PGS1), # PGS
      if (!is.na(paste_item_PGS2)) paste0(VarPairs[[i]][2], paste_item_PGS2), # a second PGS if available
      if (is_binary) VarPairs[[i]][3], # the age variable if available
      paste0("PC", 1:10)) # PCs
    # for continuous phenotypes: pheno, + 1-2 PGS, + 10 PCs = 12-13 components + chiptype later on 
    # for binary phenotypes: 1-2 PGS (no pheno), + age, + 10 PCs = 12-13 components + chiptype later on
    
    # make sure we don't scale() twin pairwise deviation from the pairwise mean & pairwise mean
    # and since we aim to compare pairwise deviation coef to population coef, it is better to not scale the population level predictors either 
    no_scale_vars <- cont_vars[grepl("_PGS_mean_diff|_PGS_mean|EA4_no23andme_Okbay20221|IQ_Savage2018_FRCT11|BMI_Giant20181|Height_Yengo20221|ChildhoodBMI_Vogelezang20201", cont_vars)]
    scale_vars <- cont_vars[!grepl("_PGS_mean_diff|_PGS_mean|EA4_no23andme_Okbay20221|IQ_Savage2018_FRCT11|BMI_Giant20181|Height_Yengo20221|ChildhoodBMI_Vogelezang20201", cont_vars)]
    
    # standardise the inputs: rather than just std the dependent variable for simple regressions, we also std the independent variables here for mixed-effects model
    df_std <- df %>%
      mutate(
        # Scale variables that don't contain _PGS_mean_diff
        across(
          .cols = all_of(scale_vars),
          .fns = ~ as.numeric(scale(.)),   # Standardize (mean=0, SD=1)
          .names = "{.col}_std"            # Append "_std" to standardized vars
        ),
        # Don't scale variables that contain _PGS_mean_diff - just copy them with _std suffix
        across(
          .cols = all_of(no_scale_vars),
          .fns = ~ .,                      # Keep original values
          .names = "{.col}_std"            # Append "_std" for consistency
        )
      )
    
    depvar <- if (is_binary) { # if binary
      paste0(VarPairs[[i]][1], paste_item_trait)
    } else { # use the std version for the continuous phenotypes 
      paste0(VarPairs[[i]][1], paste_item_trait, "_std")
    }
    
    x_PGS1 <- paste0(VarPairs[[i]][2], paste_item_PGS1, "_std")
    x_PGS2 <- if (!is.na(paste_item_PGS2)) paste0(VarPairs[[i]][2], paste_item_PGS2, "_std") else NULL
    
    if (is_binary) {
      age <- paste0(VarPairs[[i]][3], "_std")  # Use standardized version
    } else {
      age <- NULL  # Or handle missing age appropriately
    }
    
    PCs <- paste0("PC", 1:10, "_std")        # Standardized PCs
    
    covar <- c(x_PGS1, x_PGS2, if(is_binary) age, sex_var, "chiptype", PCs)
    covar <- covar[!is.na(covar)]
    #covar

    # this is for the sex-stratified analysis, if there is only one sex in the sample, we won't correct for sex 
    if (length(unique(df_std[[sex_var]])) == 1) {
      covar <- covar[covar != sex_var]  # even though it is possible that NA is involved, when we subset for sex1, we subset for only one value 
    } else {

      if (!is_binary) {  # Continuous trait (VarPairs[[i]] length = 2)
        # Correct y for sex and remove sex_var from covariates
        formula_sex <- reformulate(sex_var, response = depvar)
        df_std[[depvar]] <- rstandard(lm(formula_sex, data = df_std, na.action = na.exclude))
        covar <- covar[covar != sex_var]  
      }
      # Binary trait (VarPairs[[i]] length = 3): sex_var remains in covar
    }

    f <- as.formula(paste(depvar ,"~", paste(covar, collapse = "+"), "+ (1 | id_fam)"))
    #print(f)
    
    if (is_binary) {
      # For binary outcomes
      
      model <- glmer(f,
                     family = binomial(link = "logit"),
                     data = df_std,
                     na.action = na.omit,
                     control = glmerControl(
                       optimizer = "bobyqa",  # More stable default
                       optCtrl = list(maxfun = 2e5))  # Increase iterations if needed
                     )
      
      boot_fun <- function(data, indices) {
        #data <- df_std
        d <- data[indices,]
        model_boot <- glmer(f, data = d, family = binomial(link = "logit"), na.action = na.omit, 
                            control = glmerControl(
                              optimizer = "bobyqa",  # More stable default
                              optCtrl = list(maxfun = 2e5)  # Increase iterations if needed
                            ))
        return(fixef(model_boot))
      }
      # to check the model:
      # summary(model_bobyqa)$optinfo$conv
    } else {
      # For continuous outcomes
      model <- lmer(f, data = df_std, na.action = na.omit)
      #model
      boot_fun <- function(data, indices) {
        d <- data[indices,]
        model_boot <- lmer(f, data = d, na.action = na.omit)
        return(fixef(model_boot))
      }
    }
    
    # Perform bootstrapping
    boot_results <- boot(data = df_std, 
                         statistic = boot_fun, 
                         R = nboot, 
                         parallel = "multicore", 
                         ncpus = ncpus)
    
    coefs <- summary(model)$coefficients
    coefs <- as.data.frame(coefs)
    
    coefs$Estimate <- boot_results$t0
    coefs$`Std. Error` <- apply(boot_results$t, 2, sd)
    
    # Calculate confidence intervals
    #ci_results <- lapply(1:length(boot_results$t0), function(j) {
    #  ci <- boot.ci(boot_results, type = "perc", index = j)
    #  c(ci$percent[4], ci$percent[5])
    #})
    ci_results <- lapply(1:ncol(boot_results$t), function(j) {
      ci <- boot.ci(boot_results, type = "perc", index = j)
      c(ci$percent[4], ci$percent[5])
    })
    
    ci_matrix <- do.call(rbind, ci_results)
    coefs$lower.CI <- ci_matrix[,1]
    coefs$upper.CI <- ci_matrix[,2]
    
    if (!is_binary) { # if continuous
      coefs$`Pr(>|z|)` <- 2 * (1 - pnorm(abs(coefs$Estimate / coefs$`Std. Error`))) # GLMM generates z-value instead of t, and have the p-value automatically; but for simplicity, it doesn't hurt to calculate the p-value again
      coefs <- dplyr::select(coefs, c("Estimate", "Std. Error", "t value", "Pr(>|z|)", "lower.CI", "upper.CI")) 
    } else {
      coefs <- dplyr::select(coefs, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "lower.CI", "upper.CI")) 
    }
    
    R2 <- r2(model)
    
    if (!is.na(paste_item_PGS2)) { # for the within model with 2 PGS
      # fixed effects
      WF_result <- coefs[2, ]
      names(WF_result) <- c("Est_WF", "SE_WF", "z_WF", "p.z_WF", "lower.Est_WF", "upper.Est_WF")
      
      BF_result <- coefs[3, ]
      names(BF_result) <- c("Est_BF", "SE_BF", "z_BF", "p.z_BF", "lower.Est_BF", "upper.Est_BF")
      
      #if (is_binary) { # if binary
      # random effects
      #  vc <- glmmTMB::VarCorr(model)$cond$id_fam
      #  ran_effect <- data.frame(
      #    Variance_random_within = vc[1,1],  # Extract the variance value
      #    SD_random_within = attr(vc, "stddev")[1]  # Extract the SD
      #  )
      #  names(ran_effect) <- c("Variance_random_within", "SD_random_within")
      
      # generate proper pseudo-R²
      #  R2 <- performance::r2(model)
      #} else { # if continuous
      ran_effect <- as.data.frame(VarCorr(model), comp = c("Variance", "Std.Dev."))
      ran_effect <- ran_effect[1,c(4,5)] # variance & sd of id_fam random effect
      names(ran_effect) <- c("Variance_random_within", "SD_random_within")
      #}
      
      # variance explained
      R2_df <- data.frame(
        fixed_R2_within = R2$R2_marginal,
        fixed_and_random_R2_within = R2$R2_conditional
      )
      
      WFBF_result <- cbind(WF_result, BF_result)
      result <- cbind(WFBF_result, ran_effect)
      result <- cbind(result, R2_df)
      
    } else { # for the pop model or other models with 1 PGS
      fix_effect <- coefs[2,]
      
      names(fix_effect) <- c("Est_fixed", "SE_fixed", "z_fixed", "p.z_fixed", "lower.Est_fixed", "upper.Est_fixed")
      
      #if (is_binary) { # if binary
      # random effects
      #  vc <- glmmTMB::VarCorr(model)$cond$id_fam
      #  ran_effect <- data.frame(
      #    Variance_random_within = vc[1,1],  # Extract the variance value
      #    SD_random_within = attr(vc, "stddev")[1]  # Extract the SD
      #  )
      #  names(ran_effect) <- c("Variance_random_within", "SD_random_within")
      
      # generate proper pseudo-R²
      # R2 <- performance::r2(model)
      #} else { # if continuous
      ran_effect <- as.data.frame(VarCorr(model), comp = c("Variance", "Std.Dev."))
      ran_effect <- ran_effect[1,c(4,5)] # variance & sd of id_fam random effect
      names(ran_effect) <- c("Variance_random_within", "SD_random_within")
      
      R2 <- performance::r2(model)
      #}
      
      R2_df <- data.frame(
        fixed_R2 = R2$R2_marginal,
        fixed_and_random_R2 = R2$R2_conditional
      )
      
      result <- cbind(fix_effect, ran_effect)
      result <- cbind(result, R2_df)
    }
    
    result$trait_y <- VarPairs[[i]][1]
    result$PGS_x <- VarPairs[[i]][2]
    results <- rbind(results, result)
  }
  
  rownames(results) <- NULL
  
  return(results)
}



cat("LMM within \n")
within_results <- process_LLM_fun_boot(df, "1", "_PGS_mean_diff", "_PGS_mean", "sex1_label") # paste_item_trait, paste_item_PGS1, paste_item_PGS2
within_results$p.z_WF <- p.adjust(within_results$p.z_WF, method="fdr")
within_results$p.z_BF <- p.adjust(within_results$p.z_BF, method="fdr")

cat("LMM pop \n")
pop_results <- process_LLM_fun_boot(df, "1", "1", NA, "sex1_label")
pop_results$p.z_fixed <- p.adjust(pop_results$p.z_fixed, method="fdr")



within_results$Type <- "LMM_within"
pop_results$Type <- "LMM_pop"

dim(within_results) # 40 * 19
dim(pop_results) # 40 * 13



cat("process the results \n")
Trait_Domain <- as.data.frame(rbind(
  # g (IQ3 PGS)
  c("gcg", "childhood g", "Cognitive Abilities"),
  c("icg", "childhood g", "Cognitive Abilities"),
  c("jcg", "childhood g", "Cognitive Abilities"),
  c("lcg", "adolescence g", "Cognitive Abilities"),
  c("pcg", "adolescence g", "Cognitive Abilities"),
  c("ucgt", "adulthood g", "Cognitive Abilities"),
  # verbal g
  c("gcl", "childhood verbal g", "Cognitive Abilities"),
  c("icvb", "childhood verbal g", "Cognitive Abilities"),
  c("jcvb", "childhood verbal g", "Cognitive Abilities"),
  c("lverbal12T", "adolescence verbal g", "Cognitive Abilities"),
  c("pcvctota", "adolescence verbal g", "Cognitive Abilities"),
  c("ucgvbt", "adulthood verbal g", "Cognitive Abilities"),
  # nonverbal g
  c("gcn", "childhood nonverbal g", "Cognitive Abilities"),
  c("icnv", "childhood nonverbal g", "Cognitive Abilities"),
  c("jcnv", "childhood nonverbal g", "Cognitive Abilities"),
  c("lnonverbal12T", "adolescence nonverbal g", "Cognitive Abilities"),
  c("pcrvtota", "adolescence nonverbal g", "Cognitive Abilities"),
  c("ucgnvt", "adulthood nonverbal g", "Cognitive Abilities"),
  # Educational achievement
  c("gt2ac", "primary school grades", "Education Achievement"),
  c("it3ac", "primary school grades", "Education Achievement"),
  c("jt3ac", "primary school grades", "Education Achievement"),
  c("lt3ac", "primary school grades", "Education Achievement"),
  c("pcexgcsecoregrdm", "GCSE grades", "Education Achievement"),
  c("rcqalgrdm", "A-level grades", "Education Achievement"),
  c("u1cdegr1", "university grades", "Education Achievement"),
  # Educational attainment
  c("ra_level_enrol", "A-level enrollment", "Education Attainment"),
  c("rcqhe", "university enrollment", "Education Attainment"),
  c("zEA", "years of schooling", "Education Attainment"),
  # BMI
  c("gbmi", "childhood BMI", "Height & BMI"),
  c("lbmi", "adolescence BMI", "Height & BMI"),
  c("ncbmi", "adolescence BMI", "Height & BMI"),
  c("pcbmi", "adolescence BMI", "Height & BMI"),
  c("u1cbmi", "adulthood BMI", "Height & BMI"),
  c("zmhbmi", "adulthood BMI", "Height & BMI"),
  # height
  c("ghtcm", "childhood height", "Height & BMI"),
  c("lchtcm", "adolescence height", "Height & BMI"),
  c("nchtcm", "adolescence height", "Height & BMI"),
  c("pcqdhtcm", "adolescence height", "Height & BMI"),
  c("u1chtcm", "adulthood height", "Height & BMI"),
  c("zmhheight", "adulthood height", "Height & BMI")
))

colnames(Trait_Domain) <- c("trait_y", "Trait_Domain", "V3")

raw_LMM_within_boot <- full_join(within_results, Trait_Domain)
raw_LMM_pop_boot <- full_join(pop_results, Trait_Domain)

# save the pre-meta-analysed results ####
raw_LMM_within_boot <- raw_LMM_within_boot %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("LMM_within","LMM_pop"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, V3, Trait_Domain)

raw_LMM_pop_boot <- raw_LMM_pop_boot %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("LMM_within","LMM_pop"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, V3, Trait_Domain)

write.csv(raw_LMM_within_boot, paste(outFileStem, "WFBF_Cog_LMM_within_boot_raw.csv", sep=""), row.names = FALSE)
write.csv(raw_LMM_pop_boot, paste(outFileStem, "WFBF_Cog_LMM_pop_boot_raw.csv", sep=""), row.names = FALSE)



# now, we can do meta-analysis for each of the LMM results
cat("meta-analyses for LMM \n")
run_meta_analysis <- function(result) {
  #result <- raw_LMM_within_boot
  results_list <- list()
  unique_combinations <- unique(result[, c("Trait_Domain", "Type", "V3")])
  
  for (i in seq_len(nrow(unique_combinations))) {
    current_data <- filter(result, 
                           Trait_Domain == unique_combinations$Trait_Domain[i],
                           Type == unique_combinations$Type[i],
                           V3 == unique_combinations$V3[i])
    
    if (unique_combinations$Type[i] == "LMM_within") {
      # For within-family analysis
      rma_result_WF <- rma(yi = current_data$Est_WF, 
                           sei = current_data$SE_WF, 
                           weights = 1 / (current_data$SE_WF^2), 
                           method = "REML")
      
      rma_result_BF <- rma(yi = current_data$Est_BF, 
                           sei = current_data$SE_BF, 
                           weights = 1 / (current_data$SE_BF^2), 
                           method = "REML")
      
      results_list[[i]] <- list(
        Trait_Domain = unique_combinations$Trait_Domain[i],
        Type = unique_combinations$Type[i],
        V3 = unique_combinations$V3[i],
        beta_WF = rma_result_WF$beta,
        se_WF = rma_result_WF$se,
        ci_lower_WF = rma_result_WF$ci.lb,
        ci_upper_WF = rma_result_WF$ci.ub,
        beta_se_WF = paste0(round(rma_result_WF$beta, 2), ' (', round(rma_result_WF$se, 2), ')'),
        beta_BF = rma_result_BF$beta,
        se_BF = rma_result_BF$se,
        ci_lower_BF = rma_result_BF$ci.lb,
        ci_upper_BF = rma_result_BF$ci.ub,
        beta_se_BF = paste0(round(rma_result_BF$beta, 2), ' (', round(rma_result_BF$se, 2), ')')
      )
    } else {
      # For population-level analysis
      rma_result <- rma(yi = current_data$Est_fixed, 
                        sei = current_data$SE_fixed, 
                        weights = 1 / (current_data$SE_fixed^2), 
                        method = "REML")
      
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
  }
  result_df <- do.call(rbind, lapply(results_list, as.data.frame))
  
  result_df <- result_df %>%
    mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
    mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                          "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                          "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                          "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                          "A-level enrollment", "university enrollment", "years of schooling", 
                                                          "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                          "childhood height", "adolescence height", "adulthood height")))%>%
    arrange(V3, Trait_Domain)
  
  return(result_df)
}

MA_LMM_within_boot <- run_meta_analysis(raw_LMM_within_boot)
MA_LMM_pop_boot <- run_meta_analysis(raw_LMM_pop_boot)

write.csv(MA_LMM_within_boot, paste(outFileStem, "bootMA_WFBF_Cog_LMM_within.csv", sep=""), row.names = FALSE)
write.csv(MA_LMM_pop_boot, paste(outFileStem, "bootMA_WFBF_Cog_LMM_pop.csv", sep=""), row.names = FALSE)


#raw_LMM_within_boot <- read.csv(paste(outFileStem, "WFBF_Cog_LMM_within_boot_raw.csv", sep=""))
#MA_LMM_within_boot <- read.csv(paste(outFileStem, "bootMA_WFBF_Cog_LMM_within.csv", sep=""))

wb <- createWorkbook()
addWorksheet(wb, sheetName="raw_LMM_within_boot")
addWorksheet(wb, sheetName="raw_LMM_pop_boot")
addWorksheet(wb, sheetName="MA_LMM_within_boot")
addWorksheet(wb, sheetName="MA_LMM_pop_boot")

writeData(wb, sheet="raw_LMM_within_boot", x=raw_LMM_within_boot, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="raw_LMM_pop_boot", x=raw_LMM_pop_boot, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="MA_LMM_within_boot", x=MA_LMM_within_boot, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="MA_LMM_pop_boot", x=MA_LMM_pop_boot, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "WFBF_Cog_LMM_boot_within_pop_All.xlsx", sep=""), overwrite=TRUE)

cat("main analyses done \n")









