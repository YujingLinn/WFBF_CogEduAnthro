# title: Extreme Analyses: nonlinearity test and stepwise regressions
# author: Yujing Lin
# date: 4 Sept, 2024

# excluding dichotomous traits: squared term for 1 will be 1; if pair diff, for example, squared term for -1 will be 1, which is very strange
# no bootstrapping needed

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(Hmisc)
library(scales)
library(pscl)
library(MASS)
library(patchwork)

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/Data/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'

WFBF_selectunpaired <- read.csv(paste(sourceFileStem, "WFBF_selectunpaired.csv", sep=""))
dat_pair <- read.csv(paste(sourceFileStem, "dat_pair.csv", sep=""))
dat_pair_selectunpaired <- dat_pair %>% filter(selectunpaired == 1)
dim(dat_pair_selectunpaired) # 6973 * 322

dat_DZpair_selectunpaired <- subset(dat_pair_selectunpaired, dat_pair_selectunpaired$zygos==2)
dim(dat_DZpair_selectunpaired) # 4315 * 322

PCs <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")



# 1. nonlinear analyses ####
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
  
  #c("ra_level_enrol", "EA4_no23andme_Okbay2022", "rcqalage1"), # binary
  c("rcqalgrdm", "EA4_no23andme_Okbay2022"),
  # c("rcqalsgrdm", "EA4_no23andme_Okbay2022"), # A- & AS-level
  #c("rcqhe", "EA4_no23andme_Okbay2022", "rcqalage1"), # 18 # binary
  
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

length(VarPairs) # 38, excluding dichotomous traits *2



# the functions ####
# remove dichotomous traits and remove bootstrapping 
# Function to fit linear regression with PGS squared term for continuous traits
std_fit_regression <- function(dat, x, y, covar, x_squared) {
  # Create formula with all covariates including squared term
  formula <- reformulate(c(covar, x, x_squared), response = y)
  
  # Fit the model
  fit_model <- lm(formula, data = dat)
  
  # Calculate metrics for main effect (x)
  x_slope <- coef(fit_model)[x]  # Slope
  x_beta <- x_slope * sd(dat[[x]], na.rm = TRUE) / sd(dat[[y]], na.rm = TRUE) # Standardized beta
  se_x <- sqrt(diag(vcov(fit_model))[x])  # SE for the slope
  t_x <- x_slope / se_x  # t-value
  p_x_beta <- 2 * (1 - pt(abs(t_x), df.residual(fit_model)))  # p-value
  
  # Calculate metrics for squared term (x_squared)
  x_squared_slope <- coef(fit_model)[x_squared]  # Slope for squared term
  x_squared_beta <- x_squared_slope * sd(dat[[x_squared]], na.rm = TRUE) / sd(dat[[y]], na.rm = TRUE) # Standardized beta for squared term
  se_x_squared <- sqrt(diag(vcov(fit_model))[x_squared])  # SE for the squared term
  t_x_squared <- x_squared_slope / se_x_squared  # t-value
  p_x_squared_beta <- 2 * (1 - pt(abs(t_x_squared), df.residual(fit_model)))  # p-value
  
  # Use adjusted R² instead of regular R²
  adj_R_squared <- summary(fit_model)$adj.r.squared
  
  list(
    x_beta = x_beta,
    beta_se = se_x,
    p_x_beta = p_x_beta,
    x_squared_beta = x_squared_beta,
    beta_squared_se = se_x_squared,
    p_x_squared_beta = p_x_squared_beta,
    adj_R_squared = adj_R_squared
  )
}

# Main function for running the analysis on multiple variable pairs
beta_linear_results_fun <- function(dat, paste_item_PGS, paste_item_trait, sex_var) {
  #dat <- df # WFBF_selectunpaired
  
  # Process each variable pair
  results_list <- lapply(1:length(VarPairs), function(i) {
    cat("Processing item", i, ":", VarPairs[[i]], "\n")
    # i <- 29 # no binary
    # paste_item_trait <- "_trait_diff" # "_trait_diff", "1", "_trait_sum"
    # paste_item_PGS <- "_PGS_diff" # "_PGS_diff", "", "_PGS_sum"
    # sex_var <- "sex1_label_for_pairdiff" # "sex1_label_for_pairdiff", "sex1_label", "sex1_label_for_pairmeansum"
    
    # y ~ x + x^2 + chiptype + 10PCs
    
    # Define variables
    trait_y <- VarPairs[[i]][1]
    PGS_x <- VarPairs[[i]][2]
    
    # Create variable names
    y <- paste(trait_y, paste_item_trait, sep = "")  # Trait
    x <- paste(PGS_x, paste_item_PGS, sep = "")  # PGS
    
    # Define principal components
    PCs <- paste0("PC", 1:10)
    
    # First standardize the original variables
    continuous_vars <- c(y, x, PCs)
    
    dat <- dat %>%
      mutate(across(
        .cols = all_of(continuous_vars),
        .fns = ~ as.numeric(scale(.)),  # Standardize (mean=0, SD=1)
        .names = "{.col}_std"           # Append "_std" to standardized vars
      ))
    
    # Update variable names to use standardized versions
    y_std <- paste0(y, "_std")
    x_std <- paste0(x, "_std")
    PCs_std <- paste0(PCs, "_std")
    
    # Create squared term AFTER standardizing x
    x_squared_std <- paste0(x_std, "_squared") # This follows the recommended approach for polynomial terms (standardize first, then create higher-order terms)
    dat[[x_squared_std]] <- dat[[x_std]]^2
    
    # Set up covariates (using standardized variables)
    covar_std <- c("chiptype", PCs_std)
    
    # Handle sex-stratified analysis
    if (length(unique(dat[[sex_var]])) > 1) {
      # If multiple sex values exist, adjust y for sex
      formula_sex <- reformulate(sex_var, response = y_std)
      dat[[y_std]] <- rstandard(lm(formula_sex, data = dat, na.action = na.exclude))
    }
    
    # Fit the model
    model_results <- std_fit_regression(
      dat = dat,
      x = x_std,
      y = y_std,
      covar = covar_std,
      x_squared = x_squared_std
    )
    
    list(
      PGS_x = PGS_x,
      trait_y = trait_y,
      x_beta = model_results$x_beta,
      beta_se = model_results$beta_se,
      p_x_beta = model_results$p_x_beta,
      x_squared_beta = model_results$x_squared_beta,
      beta_squared_se = model_results$beta_squared_se,
      p_x_squared_beta = model_results$p_x_squared_beta,
      adj_R_squared = model_results$adj_R_squared
    )
  })
  
  beta_linear_results <- do.call(rbind, lapply(results_list, as.data.frame))
  
  beta_linear_results <- beta_linear_results %>%
    mutate(
      x_beta_se = paste0(round(x_beta, 2), " (", round(beta_se, 2), ")"),
      x_squared_beta_se = paste0(round(x_squared_beta, 2), " (", round(beta_squared_se, 2), ")")
    )
  
  return(beta_linear_results)
}


ppl_unrelated_PGSnonlinear <- beta_linear_results_fun(WFBF_selectunpaired, "", "1", "sex1_label") # suffix for PGS, suffix for trait
ppl_unrelated_PGSnonlinear$p_x_beta <- p.adjust(ppl_unrelated_PGSnonlinear$p_x_beta, method="fdr")
ppl_unrelated_PGSnonlinear$p_x_squared_beta <- p.adjust(ppl_unrelated_PGSnonlinear$p_x_squared_beta, method="fdr")

WF_DZpairdiff_PGSnonlinear <- beta_linear_results_fun(dat_DZpair_selectunpaired, "_PGS_diff", "_trait_diff", "sex1_label_for_pairdiff")
WF_DZpairdiff_PGSnonlinear$p_x_beta <- p.adjust(WF_DZpairdiff_PGSnonlinear$p_x_beta, method="fdr")
WF_DZpairdiff_PGSnonlinear$p_x_squared_beta <- p.adjust(WF_DZpairdiff_PGSnonlinear$p_x_squared_beta, method="fdr")

raw_Trait_Domain <- as.data.frame(rbind(
  # g (IQ3 PGS)
  c("gcg", "childhood g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icg", "childhood g age 9", "Cognitive Abilities"),
  c("jcg", "childhood g age 10", "Cognitive Abilities"),
  c("lcg", "adolescence g age 12", "Cognitive Abilities"), # Adolescence (12+16)
  c("pcg", "adolescence g age 16", "Cognitive Abilities"),
  c("ucgt", "adulthood g age 25", "Cognitive Abilities"), # Adulthood (25)
  # verbal g
  c("gcl", "childhood verbal g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icvb", "childhood verbal g age 9", "Cognitive Abilities"),
  c("jcvb", "childhood verbal g age 10", "Cognitive Abilities"),
  c("lverbal12T", "adolescence verbal g age 12", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcvctota", "adolescence verbal g age 16", "Cognitive Abilities"),
  c("ucgvbt", "adulthood verbal g age 25", "Cognitive Abilities"), # Adulthood (25)
  # nonverbal g
  c("gcn", "childhood nonverbal g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icnv", "childhood nonverbal g age 9", "Cognitive Abilities"),
  c("jcnv", "childhood nonverbal g age 10", "Cognitive Abilities"),
  c("lnonverbal12T", "adolescence nonverbal g age 12", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcrvtota", "adolescence nonverbal g age 16", "Cognitive Abilities"),
  c("ucgnvt", "adulthood nonverbal g age 25", "Cognitive Abilities"), # Adulthood (25)
  # Educational achievement
  c("gt2ac", "primary school grades age 7", "Education Achievement"), # Primary (7+9+10+12)
  c("it3ac", "primary school grades age 9", "Education Achievement"), 
  c("jt3ac", "primary school grades age 10", "Education Achievement"), 
  c("lt3ac", "primary school grades age 12", "Education Achievement"), 
  c("pcexgcsecoregrdm", "GCSE grades age 16", "Education Achievement"), # GCSE (16)
  c("rcqalgrdm", "A-level grades age 18", "Education Achievement"), # A-level grades (18) # use A-level grade only, NOT A- & AS-levels
  c("u1cdegr1", "university grades age 21", "Education Achievement"), # university grades	 (21)
  # Educational attainment
  #c("ra_level_enrol", "A-level enrollment age 16", "Education Attainment"), # A-level enrolment (16), data obtained at 18
  #c("rcqhe", "university enrollment age 18", "Education Attainment"), # university enrolment (18)
  c("zEA", "years of schooling age 26", "Education Attainment"), # Years of schooling (21 & 26)
  # BMI
  c("gbmi", "childhood BMI age 7", "Height & BMI"), # ages 7 childhood
  c("lbmi", "adolescence BMI age 12", "Height & BMI"), # 12, 14, 16 adolescence
  c("ncbmi", "adolescence BMI age 14", "Height & BMI"),
  c("pcbmi", "adolescence BMI age 16", "Height & BMI"),
  c("u1cbmi", "adulthood BMI age 22", "Height & BMI"), # 22, 26 adulthood
  c("zmhbmi", "adulthood BMI age 26", "Height & BMI" ),
  # height
  c("ghtcm", "childhood height age 7", "Height & BMI"), # ages 7 childhood
  c("lchtcm", "adolescence height age 12", "Height & BMI"), # 12, 14, 16 adolescence
  c("nchtcm", "adolescence height age 14", "Height & BMI"), 
  c("pcqdhtcm", "adolescence height age 16", "Height & BMI"),
  c("u1chtcm", "adulthood height age 22", "Height & BMI"), # 22, 26 adulthood
  c("zmhheight", "adulthood height age 26", "Height & BMI")
))

dim(raw_Trait_Domain) # 40*3
colnames(raw_Trait_Domain) <- c("trait_y", "raw_Trait_Domain", "V3") # make sure the column name for traits in raw_Trait_Domain is the same as the result column name

ppl_unrelated_PGSnonlinear <- full_join(raw_Trait_Domain, ppl_unrelated_PGSnonlinear)
WF_DZpairdiff_PGSnonlinear <- full_join(raw_Trait_Domain, WF_DZpairdiff_PGSnonlinear)

# save the nonlinear analysis results ####
write.csv(ppl_unrelated_PGSnonlinear, paste(outFileStem, "nonlinear_ppl_unrelated_raw.csv", sep=""), row.names = FALSE)
write.csv(WF_DZpairdiff_PGSnonlinear, paste(outFileStem, "nonlinear_WF_DZpairdiff_raw.csv", sep=""), row.names = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName="ppl_unrelated")
addWorksheet(wb, sheetName="WF_pairdiff")

# Write each table to its corresponding worksheet
writeData(wb, sheet="ppl_unrelated", x=ppl_unrelated_PGSnonlinear, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="WF_pairdiff", x=WF_DZpairdiff_PGSnonlinear, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "nonlinear_raw_All.xlsx", sep=""), overwrite=TRUE)



# do the dichotomous traits separately ####
cat("Conduct nonlinear vs linear analyses separately for dichotomous traits")
# a-level enrolment
fit_model <- glm(ra_level_enrol1 ~ scale(EA4_no23andme_Okbay2022) + 
                   scale(rcqalage1) + sex1_label + chiptype + 
                   scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                   scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10),
                 data = WFBF_selectunpaired, family = binomial(link = "logit"))

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the population level (linear): ", R.squared))
# 10.2% -- 0.101686433058046

fit_model <- glm(ra_level_enrol1 ~ scale(EA4_no23andme_Okbay2022) + I(scale(EA4_no23andme_Okbay2022)^2) + 
                   scale(rcqalage1) + sex1_label + chiptype + 
                   scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                   scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10),
                 data = WFBF_selectunpaired, family = binomial(link = "logit"))

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the population level (nonlinear): ", R.squared))
# 10.2% -- 0.101819280283504

dat_DZpair_selectunpaired$ra_level_enrol_trait_diff <- as.factor(dat_DZpair_selectunpaired$ra_level_enrol_trait_diff)
fit_model <- polr(ra_level_enrol_trait_diff ~ scale(EA4_no23andme_Okbay2022_PGS_diff) + 
                    scale(rcqalage1) + sex1_label_for_pairdiff + chiptype + 
                    scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                    scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10), 
                  data = dat_DZpair_selectunpaired, Hess = TRUE)

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the within-family level (linear): ", R.squared))
# 4.1% -- 0.0413279718891087

fit_model <- polr(ra_level_enrol_trait_diff ~ scale(EA4_no23andme_Okbay2022_PGS_diff) + I(scale(EA4_no23andme_Okbay2022_PGS_diff)^2) + 
                    scale(rcqalage1) + sex1_label_for_pairdiff + chiptype + 
                    scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                    scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10), 
                  data = dat_DZpair_selectunpaired, Hess = TRUE)

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the within-family level (nonlinear): ", R.squared))
# 4.2% -- 0.0415347119723466



# university enrolment
fit_model <- glm(rcqhe1 ~ scale(EA4_no23andme_Okbay2022) + 
                   scale(rcqalage1) + sex1_label + chiptype + 
                   scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                   scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10),
                 data = WFBF_selectunpaired, family = binomial(link = "logit"))

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the population level (linear): ", R.squared))
# 8.6% -- 0.0859940141844108

fit_model <- glm(rcqhe1 ~ scale(EA4_no23andme_Okbay2022) + I(scale(EA4_no23andme_Okbay2022)^2) + 
                   scale(rcqalage1) + sex1_label + chiptype + 
                   scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                   scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10),
                 data = WFBF_selectunpaired, family = binomial(link = "logit"))

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the population level (nonlinear): ", R.squared))
# 8.6% -- 0.0863741222614441

dat_DZpair_selectunpaired$rcqhe_trait_diff <- as.factor(dat_DZpair_selectunpaired$rcqhe_trait_diff)
fit_model <- polr(rcqhe_trait_diff ~ scale(EA4_no23andme_Okbay2022_PGS_diff) + 
                    scale(rcqalage1) + sex1_label_for_pairdiff + chiptype + 
                    scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                    scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10), 
                  data = dat_DZpair_selectunpaired, Hess = TRUE)

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the within-family level (linear): ", R.squared))
# 3.4% -- 0.0338668140066877

fit_model <- polr(rcqhe_trait_diff ~ scale(EA4_no23andme_Okbay2022_PGS_diff) + I(scale(EA4_no23andme_Okbay2022_PGS_diff)^2) + 
                    scale(rcqalage1) + sex1_label_for_pairdiff + chiptype + 
                    scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + 
                    scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10), 
                  data = dat_DZpair_selectunpaired, Hess = TRUE)

pseudo_r2_result <- pR2(fit_model)
R.squared <- pseudo_r2_result["r2ML"]
print(paste("R-squared at the within-family level (nonlinear): ", R.squared))
# 3.4% -- 0.0341308591721041










# 2. decile plots ####
# create deciles ####
decile <- 10
# unique_PGS <- unique(unlist(lapply(VarPairs, function(x) x[2]))) # even though I removed the dichotomous traits here, it's fine only extracting the unique PGS
# length(unique_PGS) # 5

# decile plot function ####
show_col(gradient_n_pal(brewer_pal(palette = "Blues")(5))(seq(0.5, 1, length.out = 10)))
blues <- gradient_n_pal(brewer_pal(palette = "Blues")(5))(seq(0.5, 1, length.out = 10))
greens <- gradient_n_pal(brewer_pal(palette = "Greens")(5))(seq(0.5, 1, length.out = 10))
# here, I used brewer_pal(palette = "Greens")(5) to create sequential colors, then within this sequence [0,1] I define a range [0.5,1] to further create 10 gradient within a continuous gradient

create_decile_plot <- function(df,
                               trait_var,
                               pgs_var,
                               num_quantiles = decile,  # Default to 10 deciles
                               sex_var = NULL,      # Make sex_var optional with NULL default
                               x_label = "PGS Deciles",
                               y_label = "Trait Scores",
                               plot_title = "",
                               color,      # Default color
                               rescale_y = FALSE) {
  
  # Filter out rows with NA in key variables
  df <- df[!is.na(df[[trait_var]]) & !is.na(df[[pgs_var]]), ]
  
  # Create quantile column
  pgs_quantile_col <- paste0(pgs_var, "_quantile")
  df[[pgs_quantile_col]] <- cut(df[[pgs_var]], 
                                breaks = quantile(df[[pgs_var]], 
                                                  probs = seq(0, 1, length.out = num_quantiles + 1), 
                                                  na.rm = TRUE),
                                include.lowest = TRUE, 
                                labels = FALSE)
  
  # Check if trait is binary/ordinal based on variable name pattern
  is_binary <- grepl("ra_level_enrol|rcqhe", trait_var)
  
  # Get unique values to determine trait type more precisely
  unique_vals <- unique(df[[trait_var]][!is.na(df[[trait_var]])])
  
  # Create a new trait variable for sex adjustment if needed
  trait_var_adjusted <- trait_var
  
  # Adjust for sex if a sex variable is provided
  if (!is.null(sex_var) && sex_var %in% names(df)) {
    
    # Different approaches based on trait type
    if (is_binary && length(unique_vals) == 2) {
      # For binary traits, we'll handle this later with sex-stratified means
      # Nothing to do at this step - we'll keep original values
    } else if (is_binary && length(unique_vals) > 2) {
      # For ordinal variables with more than 2 levels
      # Convert to factor if not already
      if (!is.factor(df[[trait_var]])) {
        df[[trait_var]] <- factor(df[[trait_var]], levels = sort(unique_vals), ordered = TRUE)
      }
      
      # Use polr for ordinal regression from MASS package if available
      if (requireNamespace("MASS", quietly = TRUE)) {
        complete_cases <- complete.cases(df[, c(trait_var, sex_var)])
        temp_df <- df[complete_cases, ]
        
        sex_model <- tryCatch({
          MASS::polr(as.formula(paste(trait_var, "~", sex_var)), data = temp_df, Hess = TRUE)
        }, error = function(e) {
          warning("Error in ordinal regression: ", e$message, 
                  " - proceeding without sex adjustment for this trait")
          NULL
        })
        
        # If model succeeded, create adjusted trait
        if (!is.null(sex_model)) {
          # This is complex for ordinal models - simplify by using residuals
          # as an approximation (not statistically perfect but practical)
          trait_var_adjusted <- paste0(trait_var, "_sex_adj")
          df[[trait_var_adjusted]] <- df[[trait_var]]  # Initialize with original values
        }
      } else {
        warning("MASS package not available for ordinal regression")
      }
    } else {
      # For continuous traits, use linear regression and get residuals
      complete_cases <- complete.cases(df[, c(trait_var, sex_var)])
      temp_df <- df[complete_cases, ]
      
      # Fit the model on complete cases
      sex_model <- tryCatch({
        lm(as.formula(paste(trait_var, "~", sex_var)), data = temp_df)
      }, error = function(e) {
        warning("Error in linear regression: ", e$message, 
                " - proceeding without sex adjustment for this trait")
        NULL
      })
      
      if (!is.null(sex_model)) {
        # Create a new column with sex-corrected trait values
        trait_var_adjusted <- paste0(trait_var, "_sex_adj")
        
        # Initialize with NA
        df[[trait_var_adjusted]] <- NA
        
        # Calculate and assign residuals + mean for complete cases
        df[complete_cases, trait_var_adjusted] <- 
          residuals(sex_model) + mean(temp_df[[trait_var]], na.rm = TRUE)
      }
    }
  }
  
  # Helper function to calculate summary statistics
  summarize_trait <- function(data, trait_column) {
    data %>%
      summarise(
        n = n(),
        mean_trait = mean(.data[[trait_column]], na.rm = TRUE),
        sd_trait = sd(.data[[trait_column]], na.rm = TRUE),
        se_trait = sd_trait / sqrt(sum(!is.na(.data[[trait_column]]))),
        ci_lower = mean_trait - 1.96 * se_trait,
        ci_upper = mean_trait + 1.96 * se_trait,
        .groups = 'drop'
      )
  }
  
  # Special handling for 3-level ordinal (-1, 0, 1) traits
  if (is_binary && length(unique_vals) == 3 && 
      all(sort(unique_vals) == c(-1, 0, 1))) {
    
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package \"MASS\" needed for ordinal regression. Please install it.")
    }
    
    # Factor the 3-level variable
    df[[trait_var_adjusted]] <- factor(df[[trait_var_adjusted]], 
                                       levels = c(-1, 0, 1), 
                                       ordered = TRUE)
    
    # Generate summary statistics for each quantile
    decile_summary <- df %>%
      filter(!is.na(.data[[pgs_quantile_col]])) %>%
      group_by(.data[[pgs_quantile_col]]) %>%
      summarise(
        n = n(),
        # Calculate proportions for each level
        prop_minus1 = sum(.data[[trait_var_adjusted]] == -1, na.rm = TRUE) / n,
        prop_zero = sum(.data[[trait_var_adjusted]] == 0, na.rm = TRUE) / n,
        prop_plus1 = sum(.data[[trait_var_adjusted]] == 1, na.rm = TRUE) / n,
        # Use a weighted score as the "mean" (-1 * prop_minus1 + 0 * prop_zero + 1 * prop_plus1)
        mean_trait = prop_plus1 - prop_minus1,
        # Approximate SE for this difference of proportions
        se_trait = sqrt((prop_plus1 * (1 - prop_plus1) + prop_minus1 * (1 - prop_minus1)) / n),
        ci_lower = mean_trait - 1.96 * se_trait,
        ci_upper = mean_trait + 1.96 * se_trait,
        .groups = 'drop'
      )
  } else {
    # Standard summarization for binary/continuous traits
    decile_summary <- df %>%
      filter(!is.na(.data[[pgs_quantile_col]])) %>%
      group_by(.data[[pgs_quantile_col]]) %>%
      summarize_trait(trait_var_adjusted)
  }
  
  # Additional sex stratification for binary traits (not already handled by residuals)
  if (!is.null(sex_var) && sex_var %in% names(df) && 
      is_binary && length(unique_vals) == 2 && trait_var_adjusted == trait_var) {
    
    # Calculate sex-stratified means per decile
    sex_stratified <- df %>%
      filter(!is.na(.data[[pgs_quantile_col]])) %>%
      group_by(.data[[pgs_quantile_col]], .data[[sex_var]]) %>%
      summarize_trait(trait_var) %>%
      ungroup()
    
    # Get overall sex proportions in the dataset
    sex_props <- df %>%
      filter(!is.na(.data[[sex_var]])) %>%
      count(.data[[sex_var]]) %>%
      mutate(prop = n / sum(n))
    
    # Compute weighted average based on population sex proportions
    decile_summary <- sex_stratified %>%
      inner_join(sex_props, by = sex_var) %>%
      group_by(.data[[pgs_quantile_col]]) %>%
      summarise(
        mean_trait = sum(mean_trait * prop),
        # Approximate SE for the weighted mean
        se_trait = sqrt(sum((se_trait * prop)^2)),
        ci_lower = mean_trait - 1.96 * se_trait,
        ci_upper = mean_trait + 1.96 * se_trait,
        .groups = 'drop'
      )
  }
  
  # Optionally rescale to mean 100, SD 15 (common in IQ/cognitive testing)
  if (rescale_y) {
    trait_mean <- mean(decile_summary$mean_trait, na.rm = TRUE)
    trait_sd <- sd(decile_summary$mean_trait, na.rm = TRUE)
    
    decile_summary <- decile_summary %>%
      mutate(
        mean_trait = ((mean_trait - trait_mean) / trait_sd) * 15 + 100,
        se_trait = se_trait * (15 / trait_sd),
        ci_lower = mean_trait - 1.96 * se_trait,
        ci_upper = mean_trait + 1.96 * se_trait
      )
  }
  
  # Create the plot
  ggplot(decile_summary, aes(x = factor(.data[[pgs_quantile_col]]), y = mean_trait)) +
    geom_point(size = 1.5, color = color) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = color) +
    geom_line(aes(group = 1), color = color) +
    labs(
      x = x_label,
      y = y_label,
      title = plot_title
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      strip.text.x = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
}


# apply the function
create_decile_plot(
  df = dat_pair_selectunpaired, 
  trait_var = "gcg1",  # Binary 0/1 variable
  pgs_var = "IQ_Savage2018_FRCT11",
  sex_var = "sex1_label",
  x_label = "IQ3 PGS Deciles", 
  y_label = "General cognitive ability trait scores", 
  plot_title = "General cognitive ability",
  color = blues,
  rescale_y = TRUE
)

EA_list <- list(
  c("gt2ac", "EA4_no23andme_Okbay2022", "School Grades Age 7", T),
  c("it3ac", "EA4_no23andme_Okbay2022", "School Grades Age 9", T),
  c("jt3ac", "EA4_no23andme_Okbay2022", "School Grades Age 10", T),
  c("lt3ac", "EA4_no23andme_Okbay2022", "School Grades Age 12", T),
  c("pcexgcsecoregrdm", "EA4_no23andme_Okbay2022", "GCSE Grades Age 16", T),
  c("rcqalgrdm", "EA4_no23andme_Okbay2022", "A-level Grades Age 18", T),
  c("u1cdegr1", "EA4_no23andme_Okbay2022", "University Grades Age 21", T),
  c("ra_level_enrol", "EA4_no23andme_Okbay2022", "A-level Enrolment Proportions", F),
  c("rcqhe", "EA4_no23andme_Okbay2022", "University Enrolment Proportions", F),
  c("zEA", "EA4_no23andme_Okbay2022", "Educational Attainment \n (Years of Schooling) Age 26", F)
)

EA_plot_list <- lapply(1:length(EA_list), function(i) {
  #i <- 1
  trait_var <- paste(EA_list[[i]][1], "1", sep="")
  pgs_var <- paste(EA_list[[i]][2], "1", sep="")
  rescale_y = EA_list[[i]][4]
  plot_title <- EA_list[[i]][3]
  x_label <- "EA4 PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = rescale_y
  )
})

EA_plots_list <- wrap_plots(EA_plot_list, ncol = 4)
print(EA_plots_list)
ggsave("EA_plots_list.png", EA_plots_list, width = 16, height = 12)

EA_plot_list_WF <- lapply(1:length(EA_list), function(i) {
  #i <- 1
  trait_var <- paste(EA_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(EA_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = EA_list[[i]][4]
  plot_title <- EA_list[[i]][3]
  x_label <- "Twin Pairwise EA4 PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

EA_plots_list_WF <- wrap_plots(EA_plot_list_WF, ncol = 4)
print(EA_plots_list_WF)
ggsave("EA_plots_list_WF.png", EA_plots_list_WF, width = 16, height = 12)





g_list <- list(
  c("gcg", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 7", T), 
  c("icg", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 9", T), 
  c("jcg", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 10", T),
  c("lcg", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 12", T), 
  c("pcg", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 16", T), 
  c("ucgt", "IQ_Savage2018_FRCT1", "General Cognitive Ability Age 25", T) 
)

g_plot_list <- lapply(1:length(g_list), function(i) {
  #i <- 1
  trait_var <- paste(g_list[[i]][1], "1", sep="")
  pgs_var <- paste(g_list[[i]][2], "1", sep="")
  rescale_y = g_list[[i]][4]
  plot_title <- g_list[[i]][3]
  x_label <- "IQ3 PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = rescale_y
  )
})

g_plots_list <- wrap_plots(g_plot_list, ncol = 3)
print(g_plots_list)
ggsave("g_plots_list.png", g_plots_list, width = 16, height = 12)

g_plot_list_WF <- lapply(1:length(g_list), function(i) {
  #i <- 1
  trait_var <- paste(g_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(g_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = g_list[[i]][4]
  plot_title <- g_list[[i]][3]
  x_label <- "Twin Pairwise IQ3 PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

g_plots_list_WF <- wrap_plots(g_plot_list_WF, ncol = 3)
print(g_plots_list_WF)
ggsave("g_plots_list_WF.png", g_plots_list_WF, width = 16, height = 12)





verbal_g_list <- list(
  c("gcl", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 7", T),
  c("icvb", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 9", T), 
  c("jcvb", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 10", T), 
  c("lverbal12T", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 12", T), 
  c("pcvctota", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 16", T), 
  c("ucgvbt", "IQ_Savage2018_FRCT1", "Verbal Cognitive Ability Age 25", T) 
)

verbal_g_plot_list <- lapply(1:length(verbal_g_list), function(i) {
  #i <- 1
  trait_var <- paste(verbal_g_list[[i]][1], "1", sep="")
  pgs_var <- paste(verbal_g_list[[i]][2], "1", sep="")
  rescale_y = verbal_g_list[[i]][4]
  plot_title <- verbal_g_list[[i]][3]
  x_label <- "IQ3 PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = rescale_y
  )
})

verbal_g_plots_list <- wrap_plots(verbal_g_plot_list, ncol = 3)
print(verbal_g_plots_list)
ggsave("verbal_g_plots_list.png", verbal_g_plots_list, width = 16, height = 12)



verbal_g_plot_list_WF <- lapply(1:length(verbal_g_list), function(i) {
  #i <- 1
  trait_var <- paste(verbal_g_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(verbal_g_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = verbal_g_list[[i]][4]
  plot_title <- verbal_g_list[[i]][3]
  x_label <- "Twin Pairwise IQ3 PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

verbal_g_plots_list_WF <- wrap_plots(verbal_g_plot_list_WF, ncol = 3)
print(verbal_g_plots_list_WF)
ggsave("verbal_g_plots_list_WF.png", verbal_g_plots_list_WF, width = 16, height = 12)




nonverbal_g_list <- list(
  c("gcn", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 7", T), 
  c("icnv", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 9", T), 
  c("jcnv", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 10", T), 
  c("lnonverbal12T", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 12", T), 
  c("pcrvtota", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 16", T),
  c("ucgnvt", "IQ_Savage2018_FRCT1", "Nonverbal Cognitive Ability Age 25", T) 
)

nonverbal_g_plot_list <- lapply(1:length(nonverbal_g_list), function(i) {
  #i <- 1
  trait_var <- paste(nonverbal_g_list[[i]][1], "1", sep="")
  pgs_var <- paste(nonverbal_g_list[[i]][2], "1", sep="")
  rescale_y = nonverbal_g_list[[i]][4]
  plot_title <- nonverbal_g_list[[i]][3]
  x_label <- "IQ3 PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = rescale_y
  )
})

nonverbal_g_plots_list <- wrap_plots(nonverbal_g_plot_list, ncol = 3)
print(nonverbal_g_plots_list)
ggsave("nonverbal_g_plots_list.png", nonverbal_g_plots_list, width = 16, height = 12)



nonverbal_g_plot_list_WF <- lapply(1:length(nonverbal_g_list), function(i) {
  #i <- 1
  trait_var <- paste(nonverbal_g_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(nonverbal_g_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = nonverbal_g_list[[i]][4]
  plot_title <- nonverbal_g_list[[i]][3]
  x_label <- "Twin Pairwise IQ3 PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

nonverbal_g_plots_list_WF <- wrap_plots(nonverbal_g_plot_list_WF, ncol = 3)
print(nonverbal_g_plots_list_WF)
ggsave("nonverbal_g_plots_list_WF.png", nonverbal_g_plots_list_WF, width = 16, height = 12)





BMI_list <- list(
  c("gbmi", "ChildhoodBMI_Vogelezang2020", "BMI Age 7", F),
  c("lbmi", "BMI_Giant2018", "BMI Age 12", F),
  c("ncbmi", "BMI_Giant2018", "BMI Age 14", F),
  c("pcbmi", "BMI_Giant2018", "BMI Age 16", F),
  c("u1cbmi", "BMI_Giant2018", "BMI Age 22", F),
  c("zmhbmi", "BMI_Giant2018", "BMI Age 26", F)
)

BMI_plot_list <- lapply(1:length(BMI_list), function(i) {
  #i <- 1
  trait_var <- paste(BMI_list[[i]][1], "1", sep="")
  pgs_var <- paste(BMI_list[[i]][2], "1", sep="")
  rescale_y = BMI_list[[i]][4]
  plot_title <- BMI_list[[i]][3]
  x_label <- "BMI PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = FALSE
  )
})

BMI_plots_list <- wrap_plots(BMI_plot_list, ncol = 3)
print(BMI_plots_list)
ggsave("BMI_plots_list.png", BMI_plots_list, width = 16, height = 12)



BMI_plot_list_WF <- lapply(1:length(BMI_list), function(i) {
  #i <- 1
  trait_var <- paste(BMI_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(BMI_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = BMI_list[[i]][4]
  plot_title <- BMI_list[[i]][3]
  x_label <- "Twin Pairwise BMI PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

BMI_plots_list_WF <- wrap_plots(BMI_plot_list_WF, ncol = 3)
print(BMI_plots_list_WF)
ggsave("BMI_plots_list_WF.png", BMI_plots_list_WF, width = 16, height = 12)





Height_list <- list(
  c("ghtcm", "Height_Yengo2022", "Height Age 7", F),
  c("lchtcm", "Height_Yengo2022", "Height Age 12", F),
  c("nchtcm", "Height_Yengo2022", "Height Age 14", F),
  c("pcqdhtcm", "Height_Yengo2022", "Height Age 16", F),
  c("u1chtcm", "Height_Yengo2022", "Height Age 22", F),
  c("zmhheight", "Height_Yengo2022", "Height Age 26", F)
)

Height_plot_list <- lapply(1:length(Height_list), function(i) {
  #i <- 1
  trait_var <- paste(Height_list[[i]][1], "1", sep="")
  pgs_var <- paste(Height_list[[i]][2], "1", sep="")
  rescale_y = Height_list[[i]][4]
  plot_title <- Height_list[[i]][3]
  x_label <- "Height PGS Deciles"
  y_label <- "Within-Decile Phenotypic Means"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_pair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = blues,
    rescale_y = rescale_y
  )
})

Height_plots_list <- wrap_plots(Height_plot_list, ncol = 3)
print(Height_plots_list)
ggsave("Height_plots_list.png", Height_plots_list, width = 16, height = 12)



Height_plot_list_WF <- lapply(1:length(Height_list), function(i) {
  #i <- 1
  trait_var <- paste(Height_list[[i]][1], "_trait_diff", sep="")
  pgs_var <- paste(Height_list[[i]][2], "_PGS_diff", sep="")
  rescale_y = Height_list[[i]][4]
  plot_title <- Height_list[[i]][3]
  x_label <- "Twin Pairwise Height PGS Difference Deciles"
  y_label <- "Within-Decile Twin Pairwise \n Phenotypic Mean Differences"
  
  # Explicitly pass the num_quantiles parameter
  create_decile_plot(
    df = dat_DZpair_selectunpaired, 
    trait_var = trait_var, 
    pgs_var = pgs_var, 
    sex_var = "sex1_label_for_pairdiff",
    x_label = x_label, 
    y_label = y_label, 
    plot_title = plot_title,
    num_quantiles = 10,
    color = greens,
    rescale_y = FALSE
  )
})

Height_plots_list_WF <- wrap_plots(Height_plot_list_WF, ncol = 3)
print(Height_plots_list_WF)
ggsave("Height_plots_list_WF.png", Height_plots_list_WF, width = 16, height = 12)










# 3. quadratic vs linear plots ####
nonlinearity_plots <- function(trait, pgs, names) {
  plot <- ggplot() + 
    # First dataset: POP
    geom_smooth(data=WFBF_selectunpaired, aes(x = get(pgs), y = get(paste0(trait, "1")), color = "Quadratic"), 
                method = "lm", formula = y ~ x + I(x^2), se = T) +
    geom_smooth(data=WFBF_selectunpaired, aes(x = get(pgs), y = get(paste0(trait, "1")), color = "Linear"), 
                method = "lm", formula = y ~ x, se = T) +
    labs(x = "BMI PGS", y = names, title = "", color = "Model") +
    theme(
      axis.text = element_text(size = 16), 
      strip.text.x = element_text(size = 12), 
      legend.text = element_text(size = 16), 
      axis.title = element_text(size = 16)
    ) +
    theme_minimal() +
    # Specify custom colors for each model using scale_color_manual
    # scale_color_manual(values = c("Quadratic" = "brown", "Linear" = "#ea73a4")) + # for PS paper
    scale_color_manual(values = c("Quadratic" = "violet", "Linear" = "#E85C5B")) + # for cog paper
    theme(legend.position = "none")
  return(plot)
}

BMI_list <- list(
  c("pcbmi", "BMI_Giant2018", "BMI at 16"),
  c("u1cbmi", "BMI_Giant2018", "BMI at 22"),
  c("zmhbmi", "BMI_Giant2018", "BMI at 26")
  # c("zeating_disorder_diag", "AnorexiaNervosa_Watson2019", "eating disorder")
)
plots <- lapply(BMI_list, function(x) nonlinearity_plots(x[1], x[2], x[3]))
one_plot_with_legend <- nonlinearity_plots(BMI_list[[3]][1], BMI_list[[3]][2], BMI_list[[3]][3])
legend <- get_legend(one_plot_with_legend + theme(legend.position = "right"))
combined_plots <- wrap_plots(plotlist = plots, nrow = 1, ncol = 3)
final_plot <- combined_plots / legend + plot_layout(heights = c(10, 1)) +
  plot_annotation(title = "Population-Level Linear vs Quadratic Prediction Models")
final_plot
ggsave("nonlinear_linear_plots.jpg", plot = final_plot, width = 10, height = 6) 










# 4. stepwise analyses ####
# Note. some of the numbers are written down before the analyses was run properly, so they may not match the actual results 
df_gt2ac_clean <- na.omit(WFBF_selectunpaired[, c("gt2ac1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(gt2ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_gt2ac_clean), direction = "both")
# nonlinear AIC = -493.46
# linear AIC = -495.26 

fit_gt2ac_nonlinear <- lm(gt2ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_gt2ac_clean)
fit_gt2ac_linear <- lm(gt2ac1 ~ EA4_no23andme_Okbay2022, data = df_gt2ac_clean)
anova(fit_gt2ac_nonlinear, fit_gt2ac_linear)
# p=0.6546



df_gcg_clean <- na.omit(WFBF_selectunpaired[, c("gcg1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(gcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcg_clean), direction = "both")
# nonlinear AIC = -182.2
# linear AIC = -183.23 

fit_gcg_nonlinear <- lm(gcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcg_clean)
fit_gcg_linear <- lm(gcg1 ~ IQ_Savage2018_FRCT1, data = df_gcg_clean)
anova(fit_gcg_nonlinear, fit_gcg_linear)
# p=0.3269



df_gcl_clean <- na.omit(WFBF_selectunpaired[, c("gcl1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(gcl1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcl_clean), direction = "both")
# nonlinear AIC = -169.22
# linear AIC = -169.43 

fit_gcl_nonlinear <- lm(gcl1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcl_clean)
fit_gcl_linear <- lm(gcl1 ~ IQ_Savage2018_FRCT1, data = df_gcl_clean)
anova(fit_gcl_nonlinear, fit_gcl_linear)
# p=0.1813



df_gcn_clean <- na.omit(WFBF_selectunpaired[, c("gcn1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(gcn1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcn_clean), direction = "both")
# nonlinear AIC = -52.12
# linear AIC = -54.09 

fit_gcn_nonlinear <- lm(gcn1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_gcn_clean)
fit_gcn_linear <- lm(gcn1 ~ IQ_Savage2018_FRCT1, data = df_gcn_clean)
anova(fit_gcn_nonlinear, fit_gcn_linear)
# p=0.8627



df_it3ac_clean <- na.omit(WFBF_selectunpaired[, c("it3ac1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(it3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_it3ac_clean), direction = "both")
# nonlinear AIC = -226.52
# linear AIC = -227.6 

fit_it3ac_nonlinear <- lm(it3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_it3ac_clean)
fit_it3ac_linear <- lm(it3ac1 ~ EA4_no23andme_Okbay2022, data = df_it3ac_clean)
anova(fit_it3ac_nonlinear, fit_it3ac_linear)
# p=0.3394



df_icg_clean <- na.omit(WFBF_selectunpaired[, c("icg1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(icg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icg_clean), direction = "both")
# nonlinear AIC = -134.44
# linear AIC = -136.01 

fit_icg_nonlinear <- lm(icg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icg_clean)
fit_icg_linear <- lm(icg1 ~ IQ_Savage2018_FRCT1, data = df_icg_clean)
anova(fit_icg_nonlinear, fit_icg_linear)
# p=0.5106



df_icvb_clean <- na.omit(WFBF_selectunpaired[, c("icvb1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(icvb1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icvb_clean), direction = "both")
# nonlinear AIC = -50.05
# linear AIC = -51.64 

fit_icvb_nonlinear <- lm(icvb1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icvb_clean)
fit_icvb_linear <- lm(icvb1 ~ IQ_Savage2018_FRCT1, data = df_icvb_clean)
anova(fit_icvb_nonlinear, fit_icvb_linear)
# p=0.5245



df_icnv_clean <- na.omit(WFBF_selectunpaired[, c("icnv1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(icnv1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icnv_clean), direction = "both")
# nonlinear AIC = -112.17
# linear AIC = -114.05 

fit_icnv_nonlinear <- lm(icnv1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_icnv_clean)
fit_icnv_linear <- lm(icnv1 ~ IQ_Savage2018_FRCT1, data = df_icnv_clean)
anova(fit_icnv_nonlinear, fit_icnv_linear)
# p=0.7286



df_jt3ac_clean <- na.omit(WFBF_selectunpaired[, c("jt3ac1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(jt3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_jt3ac_clean), direction = "both")
# nonlinear AIC = -203.8
# linear AIC = -205.43 

fit_jt3ac_nonlinear <- lm(jt3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_jt3ac_clean)
fit_jt3ac_linear <- lm(jt3ac1 ~ EA4_no23andme_Okbay2022, data = df_jt3ac_clean)
anova(fit_jt3ac_nonlinear, fit_jt3ac_linear)
# p=0.5431



df_jcg_clean <- na.omit(WFBF_selectunpaired[, c("jcg1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(jcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcg_clean), direction = "both")
# nonlinear AIC = -65.23
# linear AIC = -67.19 

fit_jcg_nonlinear <- lm(jcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcg_clean)
fit_jcg_linear <- lm(jcg1 ~ IQ_Savage2018_FRCT1, data = df_jcg_clean)
anova(fit_jcg_nonlinear, fit_jcg_linear)
# p=0.8455



df_jcvb_clean <- na.omit(WFBF_selectunpaired[, c("jcvb1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(jcvb1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcvb_clean), direction = "both")
# nonlinear AIC = -48.81
# linear AIC = -50.77 

fit_jcvb_nonlinear <- lm(jcvb1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcvb_clean)
fit_jcvb_linear <- lm(jcvb1 ~ IQ_Savage2018_FRCT1, data = df_jcvb_clean)
anova(fit_jcvb_nonlinear, fit_jcvb_linear)
# p=0.8595



df_jcnv_clean <- na.omit(WFBF_selectunpaired[, c("jcnv1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(jcnv1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcnv_clean), direction = "both")
# nonlinear AIC = -47.88
# linear AIC = -49.82 

fit_jcnv_nonlinear <- lm(jcnv1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_jcnv_clean)
fit_jcnv_linear <- lm(jcnv1 ~ IQ_Savage2018_FRCT1, data = df_jcnv_clean)
anova(fit_jcnv_nonlinear, fit_jcnv_linear)
# p=0.81



df_lt3ac_clean <- na.omit(WFBF_selectunpaired[, c("lt3ac1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(lt3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_lt3ac_clean), direction = "both")
# nonlinear AIC = -382.58 
# linear AIC = -382.48 

fit_lt3ac_nonlinear <- lm(lt3ac1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_lt3ac_clean)
fit_lt3ac_linear <- lm(lt3ac1 ~ EA4_no23andme_Okbay2022, data = df_lt3ac_clean)
anova(fit_lt3ac_nonlinear, fit_lt3ac_linear)
# p=0.1468



df_lcg_clean <- na.omit(WFBF_selectunpaired[, c("lcg1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(lcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lcg_clean), direction = "both")
# nonlinear AIC = -544.52 
# linear AIC = -543.39

fit_lcg_nonlinear <- lm(lcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lcg_clean)
fit_lcg_linear <- lm(lcg1 ~ IQ_Savage2018_FRCT1, data = df_lcg_clean)
anova(fit_lcg_nonlinear, fit_lcg_linear)
# p=0.07738



df_lverbal12T_clean <- na.omit(WFBF_selectunpaired[, c("lverbal12T1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(lverbal12T1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lverbal12T_clean), direction = "both")
# nonlinear AIC = -261.2
# linear AIC = -263.2 

fit_lverbal12T_nonlinear <- lm(lverbal12T1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lverbal12T_clean)
fit_lverbal12T_linear <- lm(lverbal12T1 ~ IQ_Savage2018_FRCT1, data = df_lverbal12T_clean)
anova(fit_lverbal12T_nonlinear, fit_lverbal12T_linear)
# p=0.9794



df_lnonverbal12T_clean <- na.omit(WFBF_selectunpaired[, c("lnonverbal12T1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(lnonverbal12T1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lnonverbal12T_clean), direction = "both")
# nonlinear AIC = -393.75
# linear AIC = -395.75 

fit_lnonverbal12T_nonlinear <- lm(lnonverbal12T1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_lnonverbal12T_clean)
fit_lnonverbal12T_linear <- lm(lnonverbal12T1 ~ IQ_Savage2018_FRCT1, data = df_lnonverbal12T_clean)
anova(fit_lnonverbal12T_nonlinear, fit_lnonverbal12T_linear)
# p=0.9746



df_pcexgcsecoregrdm_clean <- na.omit(WFBF_selectunpaired[, c("pcexgcsecoregrdm1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(pcexgcsecoregrdm1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_pcexgcsecoregrdm_clean), direction = "both")
# nonlinear AIC = -919.31
# linear AIC = -920.68 

fit_pcexgcsecoregrdm_nonlinear <- lm(pcexgcsecoregrdm1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_pcexgcsecoregrdm_clean)
fit_pcexgcsecoregrdm_linear <- lm(pcexgcsecoregrdm1 ~ EA4_no23andme_Okbay2022, data = df_pcexgcsecoregrdm_clean)
anova(fit_pcexgcsecoregrdm_nonlinear, fit_pcexgcsecoregrdm_linear)
# p=0.4272



df_pcg_clean <- na.omit(WFBF_selectunpaired[, c("pcg1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(pcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcg_clean), direction = "both")
# nonlinear AIC = -282.52
# linear AIC = -283.59 

fit_pcg_nonlinear <- lm(pcg1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcg_clean)
fit_pcg_linear <- lm(pcg1 ~ IQ_Savage2018_FRCT1, data = df_pcg_clean)
anova(fit_pcg_nonlinear, fit_pcg_linear)
# p=0.3335



df_pcvctota_clean <- na.omit(WFBF_selectunpaired[, c("pcvctota1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(pcvctota1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcvctota_clean), direction = "both")
# nonlinear AIC = -195.77
# linear AIC = -195.97 

fit_pcvctota_nonlinear <- lm(pcvctota1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcvctota_clean)
fit_pcvctota_linear <- lm(pcvctota1 ~ IQ_Savage2018_FRCT1, data = df_pcvctota_clean)
anova(fit_pcvctota_nonlinear, fit_pcvctota_linear)
# p=0.1798



df_pcrvtota_clean <- na.omit(WFBF_selectunpaired[, c("pcrvtota1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(pcrvtota1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcrvtota_clean), direction = "both")
# nonlinear AIC = -219.52
# linear AIC = -221.51 

fit_pcrvtota_nonlinear <- lm(pcrvtota1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_pcrvtota_clean)
fit_pcrvtota_linear <- lm(pcrvtota1 ~ IQ_Savage2018_FRCT1, data = df_pcrvtota_clean)
anova(fit_pcrvtota_nonlinear, fit_pcrvtota_linear)
# p=0.9258



df_ra_level_enrol_clean <- na.omit(WFBF_selectunpaired[, c("ra_level_enrol1", "EA4_no23andme_Okbay2022")])
reg3 <- step(glm(ra_level_enrol1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_ra_level_enrol_clean), family = binomial, direction = "both")
# nonlinear AIC = 6676.82 
# linear AIC = 7152.6 

fit_ra_level_enrol_nonlinear <- glm(ra_level_enrol1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_ra_level_enrol_clean, family = binomial)
fit_ra_level_enrol_linear <- glm(ra_level_enrol1 ~ EA4_no23andme_Okbay2022, data = df_ra_level_enrol_clean, family = binomial)
# because this is glm, I cannot use anova directly to compare the two models
anova_result <- anova(fit_ra_level_enrol_linear, fit_ra_level_enrol_nonlinear)
dev_diff <- anova_result$Deviance[2]  
df_diff <- anova_result$Df[2]       
p_value <- 1 - pchisq(dev_diff, df = df_diff)
print(p_value) 
# p=0.5899346



df_rcqalgrdm_clean <- na.omit(WFBF_selectunpaired[, c("rcqalgrdm1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(rcqalgrdm1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_rcqalgrdm_clean), direction = "both")
# nonlinear AIC = -223.66
# linear AIC = -223.3 

fit_rcqalgrdm_nonlinear <- lm(rcqalgrdm1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_rcqalgrdm_clean)
fit_rcqalgrdm_linear <- lm(rcqalgrdm1 ~ EA4_no23andme_Okbay2022, data = df_rcqalgrdm_clean)
anova(fit_rcqalgrdm_nonlinear, fit_rcqalgrdm_linear)
# p=0.1241



df_rcqhe_clean <- na.omit(WFBF_selectunpaired[, c("rcqhe1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(rcqhe1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_rcqhe_clean), family = binomial, direction = "both")
# nonlinear AIC = -7367.79 
# linear AIC = -6981.2 

fit_rcqhe_nonlinear <- glm(rcqhe1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_rcqhe_clean, family = binomial)
fit_rcqhe_linear <- glm(rcqhe1 ~ EA4_no23andme_Okbay2022, data = df_rcqhe_clean, family = binomial)
# because this is glm, I cannot use anova directly to compare the two models
anova_result <- anova(fit_rcqhe_linear, fit_rcqhe_nonlinear)
dev_diff <- anova_result$Deviance[2]
df_diff <- anova_result$Df[2]       
p_value <- 1 - pchisq(dev_diff, df = df_diff)
print(p_value) 
# p=0.1973849



df_u1cdegr1_clean <- na.omit(WFBF_selectunpaired[, c("u1cdegr11", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(u1cdegr11 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_u1cdegr1_clean), direction = "both")
# nonlinear AIC = 5.2
# linear AIC = 3.22 

fit_u1cdegr1_nonlinear <- lm(u1cdegr11 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_u1cdegr1_clean)
fit_u1cdegr1_linear <- lm(u1cdegr11 ~ EA4_no23andme_Okbay2022, data = df_u1cdegr1_clean)
anova(fit_u1cdegr1_nonlinear, fit_u1cdegr1_linear)
# p=0.868



df_ucgt_clean <- na.omit(WFBF_selectunpaired[, c("ucgt1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(ucgt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgt_clean), direction = "both")
# nonlinear AIC = -278.39
# linear AIC = -280.26   

fit_ucgt_nonlinear <- lm(ucgt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgt_clean)
fit_ucgt_linear <- lm(ucgt1 ~ IQ_Savage2018_FRCT1, data = df_ucgt_clean)
anova(fit_ucgt_nonlinear, fit_ucgt_linear)
# p=0.7172



df_ucgvbt_clean <- na.omit(WFBF_selectunpaired[, c("ucgvbt1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(ucgvbt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgvbt_clean), direction = "both")
# nonlinear AIC = -284.47
# linear AIC = -286.45 

fit_ucgvbt_nonlinear <- lm(ucgvbt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgvbt_clean)
fit_ucgvbt_linear <- lm(ucgvbt1 ~ IQ_Savage2018_FRCT1, data = df_ucgvbt_clean)
anova(fit_ucgvbt_nonlinear, fit_ucgvbt_linear)
# p=0.8845



df_ucgnvt_clean <- na.omit(WFBF_selectunpaired[, c("ucgnvt1", "IQ_Savage2018_FRCT1")])
reg3 <- step(lm(ucgnvt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgnvt_clean), direction = "both")
# nonlinear AIC = -198.71
# linear AIC = -200.52 

fit_ucgnvt_nonlinear <- lm(ucgnvt1 ~ IQ_Savage2018_FRCT1 + I(IQ_Savage2018_FRCT1^2), data = df_ucgnvt_clean)
fit_ucgnvt_linear <- lm(ucgnvt1 ~ IQ_Savage2018_FRCT1, data = df_ucgnvt_clean)
anova(fit_ucgnvt_nonlinear, fit_ucgnvt_linear)
# p=0.6595



df_zEA_clean <- na.omit(WFBF_selectunpaired[, c("zEA1", "EA4_no23andme_Okbay2022")])
reg3 <- step(lm(zEA1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_zEA_clean), direction = "both")
# nonlinear AIC = -160.18 
# linear AIC = -155.83 

fit_zEA_nonlinear <- lm(zEA1 ~ EA4_no23andme_Okbay2022 + I(EA4_no23andme_Okbay2022^2), data = df_zEA_clean)
fit_zEA_linear <- lm(zEA1 ~ EA4_no23andme_Okbay2022, data = df_zEA_clean)
anova(fit_zEA_nonlinear, fit_zEA_linear)
# p=0.01176*



# anthropometric traits 
df_gbmi_clean <- na.omit(WFBF_selectunpaired[, c("gbmi1", "ChildhoodBMI_Vogelezang2020")])
reg3 <- step(lm(gbmi1 ~ ChildhoodBMI_Vogelezang2020 + I(ChildhoodBMI_Vogelezang2020^2), data = df_gbmi_clean), direction = "both")
# nonlinear AIC = -163.24
# linear AIC = -165.07 

fit_gbmi_nonlinear <- lm(gbmi1 ~ ChildhoodBMI_Vogelezang2020 + I(ChildhoodBMI_Vogelezang2020^2), data = df_gbmi_clean)
fit_gbmi_linear <- lm(gbmi1 ~ ChildhoodBMI_Vogelezang2020, data = df_gbmi_clean)
anova(fit_gbmi_nonlinear, fit_gbmi_linear)
# p=0.6815



df_ghtcm_clean <- na.omit(WFBF_selectunpaired[, c("ghtcm1", "Height_Yengo2022")])
reg3 <- step(lm(ghtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_ghtcm_clean), direction = "both")
# nonlinear AIC = -602.66
# linear AIC = -604.46 

fit_ghtcm_nonlinear <- lm(ghtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_ghtcm_clean)
fit_ghtcm_linear <- lm(ghtcm1 ~ Height_Yengo2022, data = df_ghtcm_clean)
anova(fit_ghtcm_nonlinear, fit_ghtcm_linear)
# p=0.6521



df_lbmi_clean <- na.omit(WFBF_selectunpaired[, c("lbmi1", "BMI_Giant2018")])
reg3 <- step(lm(lbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_lbmi_clean), direction = "both")
# nonlinear AIC = -277.34 
# linear AIC = -273.172 

fit_lbmi_nonlinear <- lm(lbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_lbmi_clean)
fit_lbmi_linear <- lm(lbmi1 ~ BMI_Giant2018, data = df_lbmi_clean)
anova(fit_lbmi_nonlinear, fit_lbmi_linear)
# p=0.01305*



df_lchtcm_clean <- na.omit(WFBF_selectunpaired[, c("lchtcm1", "Height_Yengo2022")])
reg3 <- step(lm(lchtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_lchtcm_clean), direction = "both")
# nonlinear AIC = -679.45
# linear AIC = -681.43 

fit_lchtcm_nonlinear <- lm(lchtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_lchtcm_clean)
fit_lchtcm_linear <- lm(lchtcm1 ~ Height_Yengo2022, data = df_lchtcm_clean)
anova(fit_lchtcm_nonlinear, fit_lchtcm_linear)
# p=0.8663



df_ncbmi_clean <- na.omit(WFBF_selectunpaired[, c("ncbmi1", "BMI_Giant2018")])
reg3 <- step(lm(ncbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_ncbmi_clean), direction = "both")
# nonlinear AIC = -230.7 
# linear AIC = 22.956 

fit_ncbmi_nonlinear <- lm(ncbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_ncbmi_clean)
fit_ncbmi_linear <- lm(ncbmi1 ~ BMI_Giant2018, data = df_ncbmi_clean)
anova(fit_ncbmi_nonlinear, fit_ncbmi_linear)
# p=0.01354*
summary(fit_ncbmi_nonlinear) # Adjusted R-squared:  0.1108 
summary(fit_ncbmi_linear) # Adjusted R-squared:  0.1087 


df_nchtcm_clean <- na.omit(WFBF_selectunpaired[, c("nchtcm1", "Height_Yengo2022")])
reg3 <- step(lm(nchtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_nchtcm_clean), direction = "both")
# nonlinear AIC = -497.92
# linear AIC = -499.75 

fit_nchtcm_nonlinear <- lm(nchtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_nchtcm_clean)
fit_nchtcm_linear <- lm(nchtcm1 ~ Height_Yengo2022, data = df_nchtcm_clean)
anova(fit_nchtcm_nonlinear, fit_nchtcm_linear)
# p=0.6771



df_pcbmi_clean <- na.omit(WFBF_selectunpaired[, c("pcbmi1", "BMI_Giant2018")])
reg3 <- step(lm(pcbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_pcbmi_clean), direction = "both")
# nonlinear AIC = -228.888 
# linear AIC = 4.201 

fit_pcbmi_nonlinear <- lm(pcbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_pcbmi_clean)
fit_pcbmi_linear <- lm(pcbmi1 ~ BMI_Giant2018, data = df_pcbmi_clean)
anova(fit_pcbmi_nonlinear, fit_pcbmi_linear)
# p=0.0001643***



df_pcqdhtcm_clean <- na.omit(WFBF_selectunpaired[, c("pcqdhtcm1", "Height_Yengo2022")])
reg3 <- step(lm(pcqdhtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_pcqdhtcm_clean), direction = "both")
# nonlinear AIC = -647.18
# linear AIC = -648.33 

fit_pcqdhtcm_nonlinear <- lm(pcqdhtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_pcqdhtcm_clean)
fit_pcqdhtcm_linear <- lm(pcqdhtcm1 ~ Height_Yengo2022, data = df_pcqdhtcm_clean)
anova(fit_pcqdhtcm_nonlinear, fit_pcqdhtcm_linear)
# p=0.3586



df_u1cbmi_clean <- na.omit(WFBF_selectunpaired[, c("u1cbmi1", "BMI_Giant2018")])
reg3 <- step(lm(u1cbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_u1cbmi_clean), direction = "both")
# nonlinear AIC = -583.52 
# linear AIC = -26.42 

fit_u1cbmi_nonlinear <- lm(u1cbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_u1cbmi_clean)
fit_u1cbmi_linear <- lm(u1cbmi1 ~ BMI_Giant2018, data = df_u1cbmi_clean)
anova(fit_u1cbmi_nonlinear, fit_u1cbmi_linear)
# p=7.408e-05***



df_u1chtcm_clean <- na.omit(WFBF_selectunpaired[, c("u1chtcm1", "Height_Yengo2022")])
reg3 <- step(lm(u1chtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_u1chtcm_clean), direction = "both")
# nonlinear AIC = -1647.46
# linear AIC = -1649.33 

fit_u1chtcm_nonlinear <- lm(u1chtcm1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_u1chtcm_clean)
fit_u1chtcm_linear <- lm(u1chtcm1 ~ Height_Yengo2022, data = df_u1chtcm_clean)
anova(fit_u1chtcm_nonlinear, fit_u1chtcm_linear)
# p=0.7132



df_zmhbmi_clean <- na.omit(WFBF_selectunpaired[, c("zmhbmi1", "BMI_Giant2018")])
reg3 <- step(lm(zmhbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_zmhbmi_clean), direction = "both")
# nonlinear AIC = -580.58 
# linear AIC = -40.11 

fit_zmhbmi_nonlinear <- lm(zmhbmi1 ~ BMI_Giant2018 + I(BMI_Giant2018^2), data = df_zmhbmi_clean)
fit_zmhbmi_linear <- lm(zmhbmi1 ~ BMI_Giant2018, data = df_zmhbmi_clean)
anova(fit_zmhbmi_nonlinear, fit_zmhbmi_linear)
# p=0.2919



df_zmhheight_clean <- na.omit(WFBF_selectunpaired[, c("zmhheight1", "Height_Yengo2022")])
reg3 <- step(lm(zmhheight1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_zmhheight_clean), direction = "both")
# nonlinear AIC = -1571.69
# linear AIC = -1572.51 

fit_zmhheight_nonlinear <- lm(zmhheight1 ~ Height_Yengo2022 + I(Height_Yengo2022^2), data = df_zmhheight_clean)
fit_zmhheight_linear <- lm(zmhheight1 ~ Height_Yengo2022, data = df_zmhheight_clean)
anova(fit_zmhheight_nonlinear, fit_zmhheight_linear)
# p=0.2763



# within-family stepwise ####
dfWF_gt2ac_clean <- na.omit(dat_DZpair_selectunpaired[, c("gt2ac_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(gt2ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_gt2ac_clean), direction = "both")
# nonlinear AIC = -176.19 
# linear AIC = -174.59 

fitWF_gt2ac_nonlinear <- lm(gt2ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_gt2ac_clean)
fitWF_gt2ac_linear <- lm(gt2ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_gt2ac_clean)
anova(fitWF_gt2ac_nonlinear, fitWF_gt2ac_linear)
# p=0.05798



dfWF_gcg_clean <- na.omit(dat_DZpair_selectunpaired[, c("gcg_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(gcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcg_clean), direction = "both")
# nonlinear AIC = 6.94
# linear AIC = 5.16 

fitWF_gcg_nonlinear <- lm(gcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcg_clean)
fitWF_gcg_linear <- lm(gcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_gcg_clean)
anova(fitWF_gcg_nonlinear, fitWF_gcg_linear)
# p=0.6411



dfWF_gcl_clean <- na.omit(dat_DZpair_selectunpaired[, c("gcl_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(gcl_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcl_clean), direction = "both")
# nonlinear AIC = -98.44
# linear AIC = -100.35 

fitWF_gcl_nonlinear <- lm(gcl_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcl_clean)
fitWF_gcl_linear <- lm(gcl_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_gcl_clean)
anova(fitWF_gcl_nonlinear, fitWF_gcl_linear)
# p=0.7649



dfWF_gcn_clean <- na.omit(dat_DZpair_selectunpaired[, c("gcn_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(gcn_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcn_clean), direction = "both")
# nonlinear AIC = 366.32
# linear AIC = 365.42 

fitWF_gcn_nonlinear <- lm(gcn_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_gcn_clean)
fitWF_gcn_linear <- lm(gcn_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_gcn_clean)
anova(fitWF_gcn_nonlinear, fitWF_gcn_linear)
# p=0.2935



dfWF_it3ac_clean <- na.omit(dat_DZpair_selectunpaired[, c("it3ac_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(it3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_it3ac_clean), direction = "both")
# nonlinear AIC = -35.56
# linear AIC = -37.38 

fitWF_it3ac_nonlinear <- lm(it3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_it3ac_clean)
fitWF_it3ac_linear <- lm(it3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_it3ac_clean)
anova(fitWF_it3ac_nonlinear, fitWF_it3ac_linear)
# p=0.6713



dfWF_icg_clean <- na.omit(dat_DZpair_selectunpaired[, c("icg_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(icg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icg_clean), direction = "both")
# nonlinear AIC = -169.28
# linear AIC = -170.8 

fitWF_icg_nonlinear <- lm(icg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icg_clean)
fitWF_icg_linear <- lm(icg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_icg_clean)
anova(fitWF_icg_nonlinear, fitWF_icg_linear)
# p=0.4901



dfWF_icvb_clean <- na.omit(dat_DZpair_selectunpaired[, c("icvb_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(icvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icvb_clean), direction = "both")
# nonlinear AIC = -112.03
# linear AIC = -113.85 

fitWF_icvb_nonlinear <- lm(icvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icvb_clean)
fitWF_icvb_linear <- lm(icvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_icvb_clean)
anova(fitWF_icvb_nonlinear, fitWF_icvb_linear)
# p=0.672



dfWF_icnv_clean <- na.omit(dat_DZpair_selectunpaired[, c("icnv_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(icnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icnv_clean), direction = "both")
# nonlinear AIC = -31.79
# linear AIC = -33.78 

fitWF_icnv_nonlinear <- lm(icnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_icnv_clean)
fitWF_icnv_linear <- lm(icnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_icnv_clean)
anova(fitWF_icnv_nonlinear, fitWF_icnv_linear)
# p=0.9092



dfWF_jt3ac_clean <- na.omit(dat_DZpair_selectunpaired[, c("jt3ac_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(jt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_jt3ac_clean), direction = "both")
# nonlinear AIC = -91.7
# linear AIC = -93.68 

fitWF_jt3ac_nonlinear <- lm(jt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_jt3ac_clean)
fitWF_jt3ac_linear <- lm(jt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_jt3ac_clean)
anova(fitWF_jt3ac_nonlinear, fitWF_jt3ac_linear)
# p=0.8749



dfWF_jcg_clean <- na.omit(dat_DZpair_selectunpaired[, c("jcg_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(jcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcg_clean), direction = "both")
# nonlinear AIC = -12.49 
# linear AIC = -11.799 

fitWF_jcg_nonlinear <- lm(jcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcg_clean)
fitWF_jcg_linear <- lm(jcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_jcg_clean)
anova(fitWF_jcg_nonlinear, fitWF_jcg_linear)
# p=0.1018



dfWF_jcvb_clean <- na.omit(dat_DZpair_selectunpaired[, c("jcvb_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(jcvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcvb_clean), direction = "both")
# nonlinear AIC = 29.48 
# linear AIC = 30.107

fitWF_jcvb_nonlinear <- lm(jcvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcvb_clean)
fitWF_jcvb_linear <- lm(jcvb_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_jcvb_clean)
anova(fitWF_jcvb_nonlinear, fitWF_jcvb_linear)
# p=0.1058



dfWF_jcnv_clean <- na.omit(dat_DZpair_selectunpaired[, c("jcnv_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(jcnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcnv_clean), direction = "both")
# nonlinear AIC = 105.2
# linear AIC = 105.01 

fitWF_jcnv_nonlinear <- lm(jcnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_jcnv_clean)
fitWF_jcnv_linear <- lm(jcnv_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_jcnv_clean)
anova(fitWF_jcnv_nonlinear, fitWF_jcnv_linear)
# p=0.1798



dfWF_lt3ac_clean <- na.omit(dat_DZpair_selectunpaired[, c("lt3ac_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(lt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_lt3ac_clean), direction = "both")
# nonlinear AIC = -196.37
# linear AIC = -198.14 

fitWF_lt3ac_nonlinear <- lm(lt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_lt3ac_clean)
fitWF_lt3ac_linear <- lm(lt3ac_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_lt3ac_clean)
anova(fitWF_lt3ac_nonlinear, fitWF_lt3ac_linear)
# p=0.6333



dfWF_lcg_clean <- na.omit(dat_DZpair_selectunpaired[, c("lcg_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(lcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lcg_clean), direction = "both")
# nonlinear AIC = 29.44
# linear AIC = 27.65 

fitWF_lcg_nonlinear <- lm(lcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lcg_clean)
fitWF_lcg_linear <- lm(lcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_lcg_clean)
anova(fitWF_lcg_nonlinear, fitWF_lcg_linear)
# p=0.6497



dfWF_lverbal12T_clean <- na.omit(dat_DZpair_selectunpaired[, c("lverbal12T_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(lverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lverbal12T_clean), direction = "both")
# nonlinear AIC = -255.8
# linear AIC = -257.69 

fitWF_lverbal12T_nonlinear <- lm(lverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lverbal12T_clean)
fitWF_lverbal12T_linear <- lm(lverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_lverbal12T_clean)
anova(fitWF_lverbal12T_nonlinear, fitWF_lverbal12T_linear)
# p=0.7407



dfWF_lnonverbal12T_clean <- na.omit(dat_DZpair_selectunpaired[, c("lnonverbal12T_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(lnonverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lnonverbal12T_clean), direction = "both")
# nonlinear AIC = 34.36
# linear AIC = 33.32 

fitWF_lnonverbal12T_nonlinear <- lm(lnonverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_lnonverbal12T_clean)
fitWF_lnonverbal12T_linear <- lm(lnonverbal12T_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_lnonverbal12T_clean)
anova(fitWF_lnonverbal12T_nonlinear, fitWF_lnonverbal12T_linear)
# p=0.3276



dfWF_pcexgcsecoregrdm_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcexgcsecoregrdm_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(pcexgcsecoregrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_pcexgcsecoregrdm_clean), direction = "both")
# nonlinear AIC = -589.71
# linear AIC = -591.59 

fitWF_pcexgcsecoregrdm_nonlinear <- lm(pcexgcsecoregrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_pcexgcsecoregrdm_clean)
fitWF_pcexgcsecoregrdm_linear <- lm(pcexgcsecoregrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_pcexgcsecoregrdm_clean)
anova(fitWF_pcexgcsecoregrdm_nonlinear, fitWF_pcexgcsecoregrdm_linear)
# p=0.7304



dfWF_pcg_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcg_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(pcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcg_clean), direction = "both")
# nonlinear AIC = 145.8
# linear AIC = 143.8 

fitWF_pcg_nonlinear <- lm(pcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcg_clean)
fitWF_pcg_linear <- lm(pcg_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_pcg_clean)
anova(fitWF_pcg_nonlinear, fitWF_pcg_linear)
# p=0.9328



dfWF_pcvctota_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcvctota_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(pcvctota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcvctota_clean), direction = "both")
# nonlinear AIC = 293.66
# linear AIC = 291.66 

fitWF_pcvctota_nonlinear <- lm(pcvctota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcvctota_clean)
fitWF_pcvctota_linear <- lm(pcvctota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_pcvctota_clean)
anova(fitWF_pcvctota_nonlinear, fitWF_pcvctota_linear)
# p=0.9685



dfWF_pcrvtota_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcrvtota_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(pcrvtota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcrvtota_clean), direction = "both")
# nonlinear AIC = 216.47
# linear AIC = 214.54 

fitWF_pcrvtota_nonlinear <- lm(pcrvtota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_pcrvtota_clean)
fitWF_pcrvtota_linear <- lm(pcrvtota_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_pcrvtota_clean)
anova(fitWF_pcrvtota_nonlinear, fitWF_pcrvtota_linear)
# p=0.7972



dfWF_ra_level_enrol_clean <- na.omit(dat_DZpair_selectunpaired[, c("ra_level_enrol_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(ra_level_enrol_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_ra_level_enrol_clean), direction = "both")
# nonlinear AIC = -3305.12
# linear AIC = -3306.68 

fitWF_ra_level_enrol_nonlinear <- lm(ra_level_enrol_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_ra_level_enrol_clean)
fitWF_ra_level_enrol_linear <- lm(ra_level_enrol_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_ra_level_enrol_clean)
anova(fitWF_ra_level_enrol_nonlinear, fitWF_ra_level_enrol_linear)
# p=0.5085



dfWF_rcqalgrdm_clean <- na.omit(dat_DZpair_selectunpaired[, c("rcqalgrdm_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(rcqalgrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_rcqalgrdm_clean), direction = "both")
# nonlinear AIC = 177.31
# linear AIC = 175.51 

fitWF_rcqalgrdm_nonlinear <- lm(rcqalgrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_rcqalgrdm_clean)
fitWF_rcqalgrdm_linear <- lm(rcqalgrdm_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_rcqalgrdm_clean)
anova(fitWF_rcqalgrdm_nonlinear, fitWF_rcqalgrdm_linear)
# p=0.6532



dfWF_rcqhe_clean <- na.omit(dat_DZpair_selectunpaired[, c("rcqhe_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(rcqhe_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_rcqhe_clean), direction = "both")
# nonlinear AIC = -2813.11 #非线性略好，虽然是binary variable，但没有显著区别
# linear AIC = -2814.49 

fitWF_rcqhe_nonlinear <- lm(rcqhe_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_rcqhe_clean)
fitWF_rcqhe_linear <- lm(rcqhe_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_rcqhe_clean)
anova(fitWF_rcqhe_nonlinear, fitWF_rcqhe_linear)
# p=0.43



dfWF_u1cdegr1_clean <- na.omit(dat_DZpair_selectunpaired[, c("u1cdegr1_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(u1cdegr1_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_u1cdegr1_clean), direction = "both")
# nonlinear AIC = 232.01
# linear AIC = 231.94 

fitWF_u1cdegr1_nonlinear <- lm(u1cdegr1_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_u1cdegr1_clean)
fitWF_u1cdegr1_linear <- lm(u1cdegr1_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_u1cdegr1_clean)
anova(fitWF_u1cdegr1_nonlinear, fitWF_u1cdegr1_linear)
# p=0.4678



dfWF_ucgt_clean <- na.omit(dat_DZpair_selectunpaired[, c("ucgt_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(ucgt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgt_clean), direction = "both")
# nonlinear AIC = 57.97
# linear AIC = 57.01   

fitWF_ucgt_nonlinear <- lm(ucgt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgt_clean)
fitWF_ucgt_linear <- lm(ucgt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_ucgt_clean)
anova(fitWF_ucgt_nonlinear, fitWF_ucgt_linear)
# p=0.3081



dfWF_ucgvbt_clean <- na.omit(dat_DZpair_selectunpaired[, c("ucgvbt_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(ucgvbt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgvbt_clean), direction = "both")
# nonlinear AIC = 74.37 
# linear AIC = 75.966 

fitWF_ucgvbt_nonlinear <- lm(ucgvbt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgvbt_clean)
fitWF_ucgvbt_linear <- lm(ucgvbt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_ucgvbt_clean)
anova(fitWF_ucgvbt_nonlinear, fitWF_ucgvbt_linear)
# p=0.05886



dfWF_ucgnvt_clean <- na.omit(dat_DZpair_selectunpaired[, c("ucgnvt_trait_diff", "IQ_Savage2018_FRCT1_PGS_diff")])
reg3WF <- step(lm(ucgnvt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgnvt_clean), direction = "both")
# nonlinear AIC = 146.4
# linear AIC = 144.51 

fitWF_ucgnvt_nonlinear <- lm(ucgnvt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff + I(IQ_Savage2018_FRCT1_PGS_diff^2), data = dfWF_ucgnvt_clean)
fitWF_ucgnvt_linear <- lm(ucgnvt_trait_diff ~ IQ_Savage2018_FRCT1_PGS_diff, data = dfWF_ucgnvt_clean)
anova(fitWF_ucgnvt_nonlinear, fitWF_ucgnvt_linear)
# p=0.7414



dfWF_zEA_clean <- na.omit(dat_DZpair_selectunpaired[, c("zEA_trait_diff", "EA4_no23andme_Okbay2022_PGS_diff")])
reg3WF <- step(lm(zEA_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_zEA_clean), direction = "both")
# nonlinear AIC = 105.16 
# linear AIC = 105.01 

fitWF_zEA_nonlinear <- lm(zEA_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff + I(EA4_no23andme_Okbay2022_PGS_diff^2), data = dfWF_zEA_clean)
fitWF_zEA_linear <- lm(zEA_trait_diff ~ EA4_no23andme_Okbay2022_PGS_diff, data = dfWF_zEA_clean)
anova(fitWF_zEA_nonlinear, fitWF_zEA_linear)
# p=0.1745



# anthropometric traits 
dfWF_gbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("gbmi_trait_diff", "ChildhoodBMI_Vogelezang2020_PGS_diff")])
reg3WF <- step(lm(gbmi_trait_diff ~ ChildhoodBMI_Vogelezang2020_PGS_diff + I(ChildhoodBMI_Vogelezang2020_PGS_diff^2), data = dfWF_gbmi_clean), direction = "both")
# nonlinear AIC = 58.27
# linear AIC = 56.27 

fitWF_gbmi_nonlinear <- lm(gbmi_trait_diff ~ ChildhoodBMI_Vogelezang2020_PGS_diff + I(ChildhoodBMI_Vogelezang2020_PGS_diff^2), data = dfWF_gbmi_clean)
fitWF_gbmi_linear <- lm(gbmi_trait_diff ~ ChildhoodBMI_Vogelezang2020_PGS_diff, data = dfWF_gbmi_clean)
anova(fitWF_gbmi_nonlinear, fitWF_gbmi_linear)
# p=0.9831



dfWF_ghtcm_clean <- na.omit(dat_DZpair_selectunpaired[, c("ghtcm_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(ghtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_ghtcm_clean), direction = "both")
# nonlinear AIC = -493.79 
# linear AIC = -493.30 

fitWF_ghtcm_nonlinear <- lm(ghtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_ghtcm_clean)
fitWF_ghtcm_linear <- lm(ghtcm_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_ghtcm_clean)
anova(fitWF_ghtcm_nonlinear, fitWF_ghtcm_linear)
# p=0.1156



dfWF_lbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("lbmi_trait_diff", "BMI_Giant2018_PGS_diff")])
reg3WF <- step(lm(lbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_lbmi_clean), direction = "both")
# nonlinear AIC = 103.31 
# linear AIC = 104.45 

fitWF_lbmi_nonlinear <- lm(lbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_lbmi_clean)
fitWF_lbmi_linear <- lm(lbmi_trait_diff ~ BMI_Giant2018_PGS_diff, data = dfWF_lbmi_clean)
anova(fitWF_lbmi_nonlinear, fitWF_lbmi_linear)
# p=0.07688



dfWF_lchtcm_clean <- na.omit(dat_DZpair_selectunpaired[, c("lchtcm_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(lchtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_lchtcm_clean), direction = "both")
# nonlinear AIC = -160.91
# linear AIC = -162.87 

fitWF_lchtcm_nonlinear <- lm(lchtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_lchtcm_clean)
fitWF_lchtcm_linear <- lm(lchtcm_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_lchtcm_clean)
anova(fitWF_lchtcm_nonlinear, fitWF_lchtcm_linear)
# p=0.8418



dfWF_ncbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("ncbmi_trait_diff", "BMI_Giant2018_PGS_diff")])
reg3WF <- step(lm(ncbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_ncbmi_clean), direction = "both")
# nonlinear AIC = 96.58
# linear AIC = 94.81 

fitWF_ncbmi_nonlinear <- lm(ncbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_ncbmi_clean)
fitWF_ncbmi_linear <- lm(ncbmi_trait_diff ~ BMI_Giant2018_PGS_diff, data = dfWF_ncbmi_clean)
anova(fitWF_ncbmi_nonlinear, fitWF_ncbmi_linear)
# p=0.6327



dfWF_nchtcm_clean <- na.omit(dat_DZpair_selectunpaired[, c("nchtcm_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(nchtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_nchtcm_clean), direction = "both")
# nonlinear AIC = 37.51
# linear AIC = 35.86 

fitWF_nchtcm_nonlinear <- lm(nchtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_nchtcm_clean)
fitWF_nchtcm_linear <- lm(nchtcm_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_nchtcm_clean)
anova(fitWF_nchtcm_nonlinear, fitWF_nchtcm_linear)
# p=0.5574



dfWF_pcbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcbmi_trait_diff", "BMI_Giant2018_PGS_diff")])
reg3WF <- step(lm(pcbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_pcbmi_clean), direction = "both")
# nonlinear AIC = 81.55 
# linear AIC = 82.156 

fitWF_pcbmi_nonlinear <- lm(pcbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_pcbmi_clean)
fitWF_pcbmi_linear <- lm(pcbmi_trait_diff ~ BMI_Giant2018_PGS_diff, data = dfWF_pcbmi_clean)
anova(fitWF_pcbmi_nonlinear, fitWF_pcbmi_linear)
# p=0.1075



dfWF_pcqdhtcm_clean <- na.omit(dat_DZpair_selectunpaired[, c("pcqdhtcm_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(pcqdhtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_pcqdhtcm_clean), direction = "both")
# nonlinear AIC = 1.09
# linear AIC = -0.88 

fitWF_pcqdhtcm_nonlinear <- lm(pcqdhtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_pcqdhtcm_clean)
fitWF_pcqdhtcm_linear <- lm(pcqdhtcm_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_pcqdhtcm_clean)
anova(fitWF_pcqdhtcm_nonlinear, fitWF_pcqdhtcm_linear)
# p=0.8562



dfWF_u1cbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("u1cbmi_trait_diff", "BMI_Giant2018_PGS_diff")])
reg3WF <- step(lm(u1cbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_u1cbmi_clean), direction = "both")
# nonlinear AIC = 293.05
# linear AIC = 291.16 

fitWF_u1cbmi_nonlinear <- lm(u1cbmi_trait_diff ~ BMI_Giant2018_PGS_diff + I(BMI_Giant2018_PGS_diff^2), data = dfWF_u1cbmi_clean)
fitWF_u1cbmi_linear <- lm(u1cbmi_trait_diff ~ BMI_Giant2018_PGS_diff, data = dfWF_u1cbmi_clean)
anova(fitWF_u1cbmi_nonlinear, fitWF_u1cbmi_linear)
# p=0.7365



dfWF_u1chtcm_clean <- na.omit(dat_DZpair_selectunpaired[, c("u1chtcm_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(u1chtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_u1chtcm_clean), direction = "both")
# nonlinear AIC = -313.38 
# linear AIC = 150.39
# summary(lm(u1chtcm_trait_diff ~ Height_Yengo2022_PGS_diff , data = dfWF_u1chtcm_clean)) # 0.274
# summary(lm(u1chtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_u1chtcm_clean)) # 0.276

fitWF_u1chtcm_nonlinear <- lm(u1chtcm_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_u1chtcm_clean)
fitWF_u1chtcm_linear <- lm(u1chtcm_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_u1chtcm_clean)
anova(fitWF_u1chtcm_nonlinear, fitWF_u1chtcm_linear)
# p=0.02388*



dfWF_zmhbmi_clean <- na.omit(dat_DZpair_selectunpaired[, c("zmhbmi_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(zmhbmi_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_zmhbmi_clean), direction = "both")
# nonlinear AIC = 217.56
# linear AIC = 215.68 

fitWF_zmhbmi_nonlinear <- lm(zmhbmi_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_zmhbmi_clean)
fitWF_zmhbmi_linear <- lm(zmhbmi_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_zmhbmi_clean)
anova(fitWF_zmhbmi_nonlinear, fitWF_zmhbmi_linear)
# p=0.7256



dfWF_zmhheight_clean <- na.omit(dat_DZpair_selectunpaired[, c("zmhheight_trait_diff", "Height_Yengo2022_PGS_diff")])
reg3WF <- step(lm(zmhheight_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_zmhheight_clean), direction = "both")
# nonlinear AIC = -344.34 
# linear AIC = -343.38 

fitWF_zmhheight_nonlinear <- lm(zmhheight_trait_diff ~ Height_Yengo2022_PGS_diff + I(Height_Yengo2022_PGS_diff^2), data = dfWF_zmhheight_clean)
fitWF_zmhheight_linear <- lm(zmhheight_trait_diff ~ Height_Yengo2022_PGS_diff, data = dfWF_zmhheight_clean)
anova(fitWF_zmhheight_nonlinear, fitWF_zmhheight_linear)
# p=0.08611