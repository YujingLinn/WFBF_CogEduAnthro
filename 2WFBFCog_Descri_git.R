# title: "Descriptive Statistics"
# author: "Yujing"
# date: "2024-04-02"

# Load libraries and datasets ####
library(dplyr)
library(tibble) # as.tibble
library(purrr) # map
library(moments)
library(magrittr)
library(Hmisc) # label
library(openxlsx)

origin <- haven::read_sav('/Users/yujinglin/Downloads/path')

WFBF_raw_full <- read.csv('/Users/yujinglin/Desktop/WF&BF/Data/WFBF_raw_full.csv')
WFBF_full <- read.csv('/Users/yujinglin/Desktop/WF&BF/Data/WFBF_full.csv')
WFBF_raw_selectunpaired <- read.csv('/Users/yujinglin/Desktop/WF&BF/Data/WFBF_raw_selectunpaired.csv')

dat_pair <- read.csv('/Users/yujinglin/Desktop/WF&BF/Data/dat_pair.csv')

# from 1WFBF_Clean.Rmd, the class of the output data is [1] "tbl_df"     "tbl"        "data.frame"
# in order for count(!is.na(variable)) to work, the class of the data need to be align with the above
# the approach described above is extremely inflexible therefore replaced by summary / summarise
origin <- as.tibble(origin)
class(origin) # [1] "tbl_df"     "tbl"        "data.frame"
WFBF_raw_full <- as.tibble(WFBF_raw_full)
WFBF_full <- as.tibble(WFBF_full)
WFBF_raw_selectunpaired <- as.tibble(WFBF_raw_selectunpaired)
dat_pair <- as.tibble(dat_pair)

dat_selectunpaired <- dat_pair %>% filter(selectunpaired == 1) # all unrelated genotyped N: 【6972】
dat_DZ <- dat_selectunpaired %>% filter(zygos == 2) # unrelated DZ: 4314
dat_MZ <- dat_selectunpaired %>% filter(zygos == 1) # unrelated MZ: 2658





# Step 2: trait counts ####
CogVars1 <- c( # "EA4_no23andme_Okbay20221", 
               "gcg1", "icg1", "jcg1", "lcg1", "pcg1", "ucgt1", # g
               "gcl1", "icvb1", "jcvb1", "lverbal12T1", "pcvctota1", "ucgvbt1", # g-verbal
               "gcn1",  "icnv1", "jcnv1", "lnonverbal12T1", "pcrvtota1", "ucgnvt1", # g-nonverbal
               "gt2ac1", "it3ac1", "jt3ac1", "lt3ac1", "pcexgcsecoregrdm1", "rcqalgrdm1", "u1cdegr11", # edu achieve
               "ra_level_enrol1", "rcqhe1", "zEA1", # edu attain
               "gbmi1", "lbmi1", "ncbmi1", "pcbmi1", "u1cbmi1", "zmhbmi1", # BMI
               "ghtcm1", "lchtcm1", "nchtcm1", "pcqdhtcm1", "u1chtcm1", "zmhheight1", # height
               
               "zmhmfqt1", "zmhganxt1", "zmhpclt1", "zmhalcoaudit1", "zmhmdq_bipolar1", "zmhconnt1", "zmhraadst1", "zeating_disorder_diag1" # psychopathology
) # length = 48

calculate_N <- function(var, data) {
  # var <- 'gcg1'
  # data = dat_pair
  N_POPunrelated <- data %>%
    filter(selectunpaired == 1) %>%
    summarise(sum = sum(!is.na(get(var)))) %>%
    pull(sum)
  Female_Percentage_POPunrelated <- data %>%
    filter(selectunpaired == 1, sex1 == 0) %>%
    summarise(Percentage = sum(!is.na(get(var))) / N_POPunrelated) %>%
    pull(Percentage)
  
  N_DZunrelated <- data %>%
    filter(selectunpaired == 1, zygos == 2) %>%
    summarise(sum = sum(!is.na(get(var)))) %>%
    pull(sum)
  Female_Percentage_DZunrelated <- data %>%
    filter(selectunpaired == 1, zygos == 2, sex1 == 0) %>%
    summarise(Percentage = sum(!is.na(get(var))) / N_DZunrelated) %>%
    pull(Percentage)
  
  N_DZfullpair <- data %>% # the unit for N counts is "pair", so the filter[selectunpaired == 1] is applied; but should rm for descriptive stats
    filter(selectunpaired == 1, zygos == 2, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(sum = sum(!is.na(get(var)))) %>%
    pull(sum)
  Female_Percentage_DZfullpair <- data %>%
    filter(selectunpaired == 1, zygos == 2, sex1 == 0, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(Percentage = sum(!is.na(get(var))) / N_DZfullpair) %>%
    pull(Percentage)
  
  N_DZSS <- data %>%
    filter(selectunpaired == 1, x3zygos == 2, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(sum = sum(!is.na(get(var)))) %>%
    pull(sum)
  Female_Percentage_DZSS <- data %>%
    filter(selectunpaired == 1, x3zygos == 2, sex1 == 0, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(Percentage = sum(!is.na(get(var))) / N_DZSS) %>%
    pull(Percentage)
  
  N_DZOS <- data %>%
    filter(selectunpaired == 1, x3zygos == 3, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(sum = sum(!is.na(get(var)))) %>%
    pull(sum)
  Female_Percentage_DZOS <- data %>%
    filter(selectunpaired == 1, x3zygos == 3, sex1 == 0, !is.na(EA4_no23andme_Okbay2022_PGS_mean)) %>%
    summarise(Percentage = sum(!is.na(get(var))) / N_DZOS) %>%
    pull(Percentage)
  
  return(c(N_POPunrelated, Female_Percentage_POPunrelated, N_DZunrelated, Female_Percentage_DZunrelated, 
           N_DZfullpair, Female_Percentage_DZfullpair, N_DZSS, Female_Percentage_DZSS, N_DZOS, Female_Percentage_DZOS))
}
 
N_list <- sapply(CogVars1, calculate_N, data = dat_pair, simplify = "array")

N_df <- data.frame(
  variable = CogVars1,
  N_POPunrelated = N_list[1, ],
  Female_Percentage_POPunrelated = N_list[2, ],
  N_DZunrelated = N_list[3, ],
  Female_Percentage_DZunrelated = N_list[4, ],
  N_DZfullpair = N_list[5, ],
  Female_Percentage_DZfullpair = N_list[6, ],
  N_DZSS = N_list[7, ],
  Female_Percentage_DZSS = N_list[8, ],
  N_DZOS = N_list[9, ],
  Female_Percentage_DZOS = N_list[10, ]
)

N_df_original <- N_df

# N_DZunrelated should NOT be used

N_df[c(3,5,7,9, 11)] <- round(N_df[c(3,5,7,9, 11)], digits = 3)
N_df$N_POPunrelated <- paste0(N_df$N_POPunrelated, " (", N_df$Female_Percentage_POPunrelated, ")")
N_df$N_DZunrelated <- paste0(N_df$N_DZunrelated, " (", N_df$Female_Percentage_DZunrelated, ")")
N_df$N_DZfullpair <- paste0(N_df$N_DZfullpair, " (", N_df$Female_Percentage_DZfullpair, ")")
N_df$N_DZSS <- paste0(N_df$N_DZSS, " (", N_df$Female_Percentage_DZSS, ")")
N_df$N_DZOS <- paste0(N_df$N_DZOS, " (", N_df$Female_Percentage_DZOS, ")")

N_df <- N_df[-c(3,5,7,9, 11)]

print(N_df)



# add labels
trait_labels <- list()
index <- 1

for (Trait in N_df$variable) {
  trait_label <- label(origin[[Trait]])
  
  trait_labels[[index]] <- list(variable = Trait, trait_label = trait_label)
  index <- index + 1
}

trait_labels <- do.call(rbind, lapply(trait_labels, data.frame, stringsAsFactors = FALSE))

N_df <- N_df %>% full_join(trait_labels)

# outFileStem <- '/Users/yujinglin/Desktop/'
# write.csv(N_df, paste(outFileStem, "N_df.csv", sep=""), row.names = FALSE)







# Step 3: Descriptives for traits ####

### var-age pairs ####
# 7: 
"gciage1" # IQ-child
"gtqage1" # edu achieve-teacher
"gpbage" # bmi, ht-parent
# 9: 
"icpage" # IQ, bmi, ht-child
"itage1" # edu achieve-teacher
# 10: 
"jcstage1" # child 
"jtqage1" # teacher
# 12: 
"lctwage1" # IQ-TOWER test, child
"lcwage1" # bmi, ht & IQ-vb, nv, same instruments as other ages, not as good as TOWER
"ltqage1" # edu achieve-teacher
# 14:
"ncqage1" # bmi, ht-child
# 16:
"pcwebage1" # IQ, bmi, ht-child web study
"pcexgcseage1" # gcse
# 18:  
"rcqalage1" # only one age variable
# 21: 
"u1cage1" # variables started with u1c is 1st wave child data, including bmi, ht, EA, current degree
"ucgage1" # variables started with ucg is the g-game, including all IQ measures
# 26:  
"zmhage1" # only one age variable

Var_Age_Pairs <- c( # "EA4_no23andme_Okbay20221", 
              # g
              c("gcg1", "gciage1"),
              c("icg1", "icpage"),
              c("jcg1", "jcstage1"),
              c("lcg1", "lctwage1"),
              c("pcg1", "pcwebage1"),
              c("ucgt1", "ucgage1"),
              # g-verbal
              c("gcl1", "gciage1"),
              c("icvb1", "icpage"),
              c("jcvb1", "jcstage1"),
              # c("lcvb1", "lcwage1"),
              c("lverbal12T1", "lctwage1"),
              c("pcvctota1", "pcwebage1"),
              c("ucgvbt1", "ucgage1"),
              # g-nonverbal
              c( "gcn1", "gciage1"),
              c("icnv1", "icpage"),
              c("jcnv1", "jcstage1"),
              # c("lcnv1", "lcwage1"),
              c("lnonverbal12T1", "lctwage1"),
              c("pcrvtota1", "pcwebage1"),
              c("ucgnvt1", "ucgage1"),
              # edu achieve
              c( "gt2ac1", "gtqage1"),
              c("it3ac1", "itage1"),
              c("jt3ac1", "jtqage1"),
              c("lt3ac1", "ltqage1"),
              c( "pcexgcsecoregrdm1", "pcexgcseage1"),
              c("rcqalgrdm1", "rcqalage1"),
              c("u1cdegr11", "u1cage1"),
              # edu attain
               c("ra_level_enrol1", "rcqalage1"),
              c("rcqhe1", "rcqalage1"),
              c("zEA1", "zmhage1"),
              # BMI
              c("gbmi1", "gpbage"),
              c("lbmi1", "lcwage1"),
              c("ncbmi1", "ncqage1"),
              c("pcbmi1", "pcwebage1"),
              c("u1cbmi1", "u1cage1"),
              c("zmhbmi1", "zmhage1" ),
              # height
              c("ghtcm1", "gpbage"),
              c("lchtcm1", "lcwage1"),
              c("nchtcm1", "ncqage1"),
              c("pcqdhtcm1", "pcwebage1"),
              c("u1chtcm1", "u1cage1"),
              c("zmhheight1", "zmhage1")
)
length(Var_Age_Pairs) # 96 variables (48 pairs)

### FUN descriptive_stat ####
descriptive_stat <- function(x) {
  stat<-c(
    mean=mean(x, na.rm=TRUE),
    sd=sd(x, na.rm=TRUE),
    min=min(x, na.rm=TRUE),
    max=max(x, na.rm=TRUE),
    skew=skewness(x, na.rm=TRUE),
    kurtosis=kurtosis(x, na.rm=TRUE)
  )
  
  stat<-round(stat, digits = 3)
  # combined<-paste(stat["mean"], "(", stat["sd"], ") [", stat["min"], ", ", stat["max"], "]", sep = "") # with [min, max]
  combined<-paste(stat["mean"], " (", stat["sd"], ")", sep = "")
  
  result<-as.data.frame(t(stat))
  colnames(result)<-c("Mean", "SD", "Min", "Max", "Skewness", "Kurtosis")
  result$Combined<-combined
  
  return(result)
}

### FUN descriptive_stat_Var_Age_Pairs ####
descriptive_stat_Var_Age_Pairs <- function(data, Var_Age_Pairs) {
  # data <- data_POPunrelated <- WFBF_raw_full %>% filter(selectunpaired == 1)
  x <- data[, Var_Age_Pairs]
  
  descriptive_stat_results <- t(sapply(x, descriptive_stat)) # use FUN descriptive_stat from above
  descriptive_stat_results <- as.data.frame(descriptive_stat_results)
  
  # since the age is every other row (2n, where n=[1,50]), create a new column named "Age" and then take the "descriptive_stat_age_POPunrelated$Combined" up
  # 1. take out the age 
  Age <- t(as.data.frame(descriptive_stat_results$Combined[seq(2, nrow(descriptive_stat_results), by = 2)]))
  # 2. remove all 2n rows starting from the 2nd row
  descriptive_stat_results <- descriptive_stat_results[-seq(2, nrow(descriptive_stat_results), by = 2),]
  # 3. cbind
  descriptive_stat_results <- cbind(descriptive_stat_results, Age)
  
  category_vector <- c(
    rep("g1_childhood", 3),
    rep("g2_adolescence", 2), # here age 12 was put in adolescence, easily change the category by changing the rep from 2 to 1, then for the above line from 3 to 4; remember to change the BMI and height categorisation below to be consistent; and also how I pair up the trait and PGS in the analyses
    rep("g3_adulthood", 1),
    
    rep("gverbal1_childhood", 3),
    rep("gverbal2_adolescence", 2),
    rep("gverbal3_adulthood", 1),
    
    rep("gnonverbal1_childhood", 3),
    rep("gnonverbal2_adolescence", 2),
    rep("gnonverbal3_adulthood", 1),

    rep("edu_achieve1_primary", 4), "edu_achieve2_GCSE", "edu_achieve3_AlvGrade", "edu_achieve4_UniGrade",
    "edu_attain1_AlvTaken", "edu_attain2_EnterUni", "edu_attain3_EA",
    
    rep("BMI1_childhood", 1), rep("BMI2_adolescence", 3), rep("BMI3_adulthood", 2), # age 7 childhood; 12, 14, 16 adolescence; 22, 26 adulthood
    rep("height1_childhood", 1), rep("height2_adolescence", 3), rep("height3_adulthood", 2)
  ) # length=48
  descriptive_stat_results$Category <- category_vector
  
  descriptive_stat_results$variable <- rownames(descriptive_stat_results) # this is added to be easily combined with N_df
  rownames(descriptive_stat_results) <- NULL
  
  names(descriptive_stat_results)[7] <- "Mean(SD)" # the 7th column should be the "Combined" one
  
  return(descriptive_stat_results)
}



### apply the function "descriptive_stat_Var_Age_Pairs" above to POPunrelated subset, DZunrelated subset, and DZfullpair subset ####
# first, create the dataset; then, apply the function
# 1. POPunrelated descriptives
data_POPunrelated <- WFBF_raw_full %>% # WFBF_raw_full is used here to align with the DZfullpair subset below
  filter(selectunpaired == 1)
# 6972 * 643

descriptive_stat_age_POPunrelated <- descriptive_stat_Var_Age_Pairs(data_POPunrelated, Var_Age_Pairs)



# 2. DZunrelated descriptives
data_DZunrelated <- WFBF_raw_full %>%
  filter(selectunpaired == 1, zygos == 2)
# 4314 * 643

descriptive_stat_age_DZunrelated <- descriptive_stat_Var_Age_Pairs(data_DZunrelated, Var_Age_Pairs)



# 3. DZfullpair
# to create a dataset contains DZ full pairs only
subset_for_filter <- dat_pair %>%
  select(id_twin, EA4_no23andme_Okbay2022_PGS_mean)  # Adjust the column names as necessary
WFBF_raw_full_merged <- WFBF_raw_full %>%
  inner_join(subset_for_filter, by = "id_twin")
data_DZfullpair <- WFBF_raw_full_merged %>%
  filter(# selectunpaired == 1, # we have all the full DZ twin pairs here, so the unit is "individual"
    zygos == 2, !is.na(EA4_no23andme_Okbay2022_PGS_mean))
# 6612 * 644 (the extra column is EA4_no23andme_Okbay2022_PGS_mean to filter full DZ twins)
# 6612 individual DZ twins from 3306 full DZ twin pairs

### very important/confusing note: the unit for "N_DZfullpair" is pair, while the unit for "data_DZfullpair" is individual ####

descriptive_stat_age_DZfullpair <- descriptive_stat_Var_Age_Pairs(data_DZfullpair, Var_Age_Pairs)



# 4. Same-Sex DZ full pairs
data_DZSSfullpair <- WFBF_raw_full_merged %>%
  filter(# selectunpaired == 1, # we have all the full DZ twin pairs here, so the unit is "individual"
    x3zygos == 2, !is.na(EA4_no23andme_Okbay2022_PGS_mean))
# x3zygos == 1: MZ, x3zygos == 2: DZSS, x3zygos == 3: DZOS

# 3464 same-sex DZ individuals, 1732 pairs

descriptive_stat_age_DZSSfullpair <- descriptive_stat_Var_Age_Pairs(data_DZSSfullpair, Var_Age_Pairs)



# 4. Opposite-Sex DZ full pairs
data_DZOSfullpair <- WFBF_raw_full_merged %>%
  filter(# selectunpaired == 1, # we have all the full DZ twin pairs here, so the unit is "individual"
    x3zygos == 3, !is.na(EA4_no23andme_Okbay2022_PGS_mean))

# 3148 opposite-sex DZ individuals, 1574 pairs

descriptive_stat_age_DZOSfullpair <- descriptive_stat_Var_Age_Pairs(data_DZOSfullpair, Var_Age_Pairs)



# combine "descriptive_stat_age_POPunrelated" & "descriptive_stat_age_DZfullpair", 
# then DZSS and DZOS
# then with "N_df"
descriptive_stat_age_POPunrelated <- descriptive_stat_age_POPunrelated[,-c(1,2,6)]

descriptive_stat_age_DZfullpair <- descriptive_stat_age_DZfullpair[,-c(1,2,6)]

descriptive_stat_age_POP_DZfullpair <- full_join(descriptive_stat_age_POPunrelated, 
                                                 descriptive_stat_age_DZfullpair, 
                                                 by = c("Category", "variable"), 
                                                 suffix = c(".pop", ".fullDZ"))

descriptive_stat_age_DZSSfullpair <- descriptive_stat_age_DZSSfullpair[,-c(1,2,6)]

descriptive_stat_age_DZOSfullpair <- descriptive_stat_age_DZOSfullpair[,-c(1,2,6)]

descriptive_stat_age_DZSS_DZOS <- full_join(descriptive_stat_age_DZSSfullpair, 
                                            descriptive_stat_age_DZOSfullpair, 
                                            by = c("Category", "variable"), 
                                            suffix = c(".DZSS", ".DZOS"))

descriptive_stat_age_POP_DZfull_DZSS_DZOS <- full_join(descriptive_stat_age_POP_DZfullpair, descriptive_stat_age_DZSS_DZOS)



N_df <- N_df[,-3]

N_descriptive_stat_age_POP_DZfull_DZSS_DZOS <- full_join(descriptive_stat_age_POP_DZfull_DZSS_DZOS, N_df)

# reorder columns
# colnames(N_descriptive_stat_age_POP_DZfull_DZSS_DZOS)

desired_order <- c("variable", "trait_label", "Category", 
                   "N_POPunrelated", "Age.pop", "Mean(SD).pop", "Min.pop", "Max.pop", "Skewness.pop", 
                   "N_DZfullpair", "Age.fullDZ", "Mean(SD).fullDZ", "Min.fullDZ", "Max.fullDZ", "Skewness.fullDZ",
                   "N_DZSS", "Age.DZSS", "Mean(SD).DZSS", "Min.DZSS", "Max.DZSS", "Skewness.DZSS",
                   "N_DZOS", "Age.DZOS", "Mean(SD).DZOS", "Min.DZOS", "Max.DZOS", "Skewness.DZOS")

order_indices <- match(desired_order, colnames(N_descriptive_stat_age_POP_DZfull_DZSS_DZOS))

N_descriptive_stat_age_POP_DZfull_DZSS_DZOS <- N_descriptive_stat_age_POP_DZfull_DZSS_DZOS[, order_indices]

dim(N_descriptive_stat_age_POP_DZfull_DZSS_DZOS) # 48*27

###【N_descriptive_stat_age_POP_DZfull_DZSS_DZOS】IS THE DESCRIPTIVE RESULT ####







# Step 4: Output the descriptive result file ####
setwd("/Users/yujinglin/Desktop/")

wb <- createWorkbook()
addWorksheet(wb, sheetName="descriptive_table")

writeData(wb, sheet=1, x=N_descriptive_stat_age_POP_DZfull_DZSS_DZOS, colNames=TRUE, rowNames=FALSE)

saveWorkbook(wb, "WFBF_descriptives.xlsx") 







# Step 5: Averaging the N by category and developmental periods ####
Category <- select(N_descriptive_stat_age_POP_DZfull_DZSS_DZOS, variable, Category)
N_df_original <- full_join(N_df_original, Category)

grouped_N <- N_df_original %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(
    Mean_N_POPunrelated = mean(N_POPunrelated, na.rm = TRUE),
    Mean_Female_POPunrelated = mean(Female_Percentage_POPunrelated, na.rm = TRUE),
    Mean_N_DZfullpair = mean(N_DZfullpair, na.rm = TRUE),
    Mean_Female_DZfullpair = mean(Female_Percentage_DZfullpair, na.rm = TRUE),
    Mean_N_DZSS = mean(N_DZSS, na.rm = TRUE),
    Mean_Female_DZSS = mean(Female_Percentage_DZSS, na.rm = TRUE),
    Mean_N_DZOS = mean(N_DZOS, na.rm = TRUE),
    Mean_Female_DZOS = mean(Female_Percentage_DZOS, na.rm = TRUE)
  )

outFileStem <- '/Users/yujinglin/Desktop/'
write.csv(grouped_N, paste(outFileStem, "grouped_N.csv", sep=""), row.names = FALSE)










# SOM: sensitivity analysis--the effect of age, sex, same-sex, opposite-sex, birth order on twin phenotypes ####
t_test_results <- function(data, vars, group_var) {
  results <- data.frame(
    Variable = character(),
    F_Statistic = numeric(),
    P_Value = numeric(),
    R_Squared = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in vars) {
    ttest <- t.test(data[[var]] ~ data[[group_var]], data = data)
    
    mean_diff <- ttest$estimate[1] - ttest$estimate[2]
    pooled_sd <- sqrt((ttest$stderr)^2 * 2)
    R_squared <- (mean_diff / pooled_sd)^2 / 
      ((mean_diff / pooled_sd)^2 + 1)
    
    results <- rbind(
      results,
      data.frame(
        Variable = var,
        F_Statistic = round(ttest$statistic, 3),
        P_Value = round(ttest$p.value, 3),
        R_Squared = round(R_squared, 3)
      )
    )
  }
  
  return(results)
}

# POP sex effect: yes
pheno_sensi_sex_MZDZ <- t_test_results(data = WFBF_raw_selectunpaired, vars = CogVars1, group_var = "sex1")
pheno_sensi_sex_MZDZ$p_adjusted <- p.adjust(pheno_sensi_sex_MZDZ$P_Value, method="fdr")
# POP birth order effect: no (only 4 variables display birth order effect)
pheno_sensi_birthorder_MZDZ <- t_test_results(data = WFBF_raw_selectunpaired, vars = CogVars1, group_var = "twin")
pheno_sensi_birthorder_MZDZ$p_adjusted <- p.adjust(pheno_sensi_birthorder_MZDZ$P_Value, method="fdr")
# POP zygosity effect: yes, for about half of the variables, but we still decide to include both DZ and MZ twins, but we provide zygosity estimates in SOM
pheno_sensi_zygos_MZDZ <- t_test_results(data = WFBF_raw_selectunpaired, vars = CogVars1, group_var = "zygos")
pheno_sensi_zygos_MZDZ$p_adjusted <- p.adjust(pheno_sensi_zygos_MZDZ$P_Value, method="fdr")



# WF sex effect (aka SS v OS): not at all
pheno_sensi_SSOS_DZ <- t_test_results(data = dat_DZ, vars = CogVars1, group_var = "x3zygos")
pheno_sensi_SSOS_DZ$p_adjusted <- p.adjust(pheno_sensi_SSOS_DZ$P_Value, method="fdr")



# Age effect is a bit tricky, b/c age is dependent on the variable, but luckily I have created a Var_Age_Pairs variable above that I can use here
# POP age effect: yes
regression_results_paired <- function(data, var_age_pairs) {
  results <- data.frame(
    Variable = character(),
    Age_Variable = character(),
    F_Statistic = numeric(),
    P_Value = numeric(),
    R_Squared = numeric(),
    stringsAsFactors = FALSE
  )
  # var_age_pairs <- Var_Age_Pairs
  i <- 2
  for (i in seq(1, length(var_age_pairs), by = 2)) {
    dep_var <- var_age_pairs[i]
    age_var <- var_age_pairs[i + 1]
    
    if (!all(c(dep_var, age_var) %in% colnames(data))) {
      warning(paste("Skipping pair:", dep_var, "and", age_var, "- not found in the dataset"))
      next
    }
    
    formula <- as.formula(paste(dep_var, "~", age_var))
    model <- lm(formula, data = data)
    
    model_summary <- summary(model)
    F_stat <- model_summary$fstatistic[1]
    P_val <- pf(F_stat, model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
    R_squared <- model_summary$r.squared
    
    results <- rbind(
      results,
      data.frame(
        Variable = dep_var,
        Age_Variable = age_var,
        F_Statistic = round(F_stat, 3),
        P_Value = round(P_val, 3),
        R_Squared = round(R_squared, 3)
      )
    )
  }
  
  return(results)
}

pheno_sensi_age_MZDZ <- regression_results_paired(data = WFBF_raw_selectunpaired, var_age_pairs = Var_Age_Pairs)
pheno_sensi_age_MZDZ$p_adjusted <- p.adjust(pheno_sensi_age_MZDZ$P_Value, method="fdr")
# if "Error: object 'Var_Age_Pairs' not found", don't panic
# run the Var_Age_Pairs above 

# add variable names to the results
raw_Trait_Domain <- as.data.frame(rbind(
  # g (IQ3 PGS)
  c("gcg1", "childhood g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icg1", "childhood g age 9", "Cognitive Abilities"),
  c("jcg1", "childhood g age 10", "Cognitive Abilities"),
  c("lcg1", "adolescence g age 12", "Cognitive Abilities"), # Adolescence (12+16)
  c("pcg1", "adolescence g age 16", "Cognitive Abilities"),
  c("ucgt1", "adulthood g age 25", "Cognitive Abilities"), # Adulthood (25)
  # verbal g
  c("gcl1", "childhood verbal g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icvb1", "childhood verbal g age 9", "Cognitive Abilities"),
  c("jcvb1", "childhood verbal g age 10", "Cognitive Abilities"),
  c("lverbal12T1", "adolescence verbal g age 12", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcvctota1", "adolescence verbal g age 16", "Cognitive Abilities"),
  c("ucgvbt1", "adulthood verbal g age 25", "Cognitive Abilities"), # Adulthood (25)
  # nonverbal g
  c("gcn1", "childhood nonverbal g age 7", "Cognitive Abilities"), # Childhood (7+9+10)
  c("icnv1", "childhood nonverbal g age 9", "Cognitive Abilities"),
  c("jcnv1", "childhood nonverbal g age 10", "Cognitive Abilities"),
  c("lnonverbal12T1", "adolescence nonverbal g age 12", "Cognitive Abilities"), # Adolescence (12+16): regular + TOWER test for age 12
  c("pcrvtota1", "adolescence nonverbal g age 16", "Cognitive Abilities"),
  c("ucgnvt1", "adulthood nonverbal g age 25", "Cognitive Abilities"), # Adulthood (25)
  # Educational achievement
  c("gt2ac1", "primary school grades age 7", "Education Achievement"), # Primary (7+9+10+12)
  c("it3ac1", "primary school grades age 9", "Education Achievement"), 
  c("jt3ac1", "primary school grades age 10", "Education Achievement"), 
  c("lt3ac1", "primary school grades age 12", "Education Achievement"), 
  c("pcexgcsecoregrdm1", "GCSE grades age 16", "Education Achievement"), # GCSE (16)
  c("rcqalgrdm1", "A-level grades age 18", "Education Achievement"), # A-level grades (18) # use A-level grade only, NOT A- & AS-levels
  c("u1cdegr11", "university grades age 21", "Education Achievement"), # university grades	 (21)
  # Educational attainment
  c("ra_level_enrol1", "A-level enrollment age 16", "Education Attainment"), # A-level enrolment (16), data obtained at 18
  c("rcqhe1", "university enrollment age 18", "Education Attainment"), # university enrolment (18)
  c("zEA1", "years of schooling age 26", "Education Attainment"), # Years of schooling (21 & 26)
  # BMI
  c("gbmi1", "childhood BMI age 7", "Height & BMI"), # ages 7 childhood
  c("lbmi1", "adolescence BMI age 12", "Height & BMI"), # 12, 14, 16 adolescence
  c("ncbmi1", "adolescence BMI age 14", "Height & BMI"),
  c("pcbmi1", "adolescence BMI age 16", "Height & BMI"),
  c("u1cbmi1", "adulthood BMI age 22", "Height & BMI"), # 22, 26 adulthood
  c("zmhbmi1", "adulthood BMI age 26", "Height & BMI" ),
  # height
  c("ghtcm1", "childhood height age 7", "Height & BMI"), # ages 7 childhood
  c("lchtcm1", "adolescence height age 12", "Height & BMI"), # 12, 14, 16 adolescence
  c("nchtcm1", "adolescence height age 14", "Height & BMI"), 
  c("pcqdhtcm1", "adolescence height age 16", "Height & BMI"),
  c("u1chtcm1", "adulthood height age 22", "Height & BMI"), # 22, 26 adulthood
  c("zmhheight1", "adulthood height age 26", "Height & BMI")
))

dim(raw_Trait_Domain) # 48*3
colnames(raw_Trait_Domain) <- c("Variable", "Trait_Domain", "V3") # make sure the column name for traits in Trait_Domain is the same as the result column name

pheno_sensi_sex_MZDZ <- full_join(pheno_sensi_sex_MZDZ, raw_Trait_Domain)
pheno_sensi_birthorder_MZDZ <- full_join(pheno_sensi_birthorder_MZDZ, raw_Trait_Domain)
pheno_sensi_zygos_MZDZ <- full_join(pheno_sensi_zygos_MZDZ, raw_Trait_Domain)
pheno_sensi_SSOS_DZ <- full_join(pheno_sensi_SSOS_DZ, raw_Trait_Domain)
pheno_sensi_age_MZDZ <- full_join(pheno_sensi_age_MZDZ, raw_Trait_Domain)



# Save the results to an Excel file
library(openxlsx)
wb <- createWorkbook()
output_list <- list(
  pheno_sensi_sex_MZDZ = pheno_sensi_sex_MZDZ,
  pheno_sensi_birthorder_MZDZ = pheno_sensi_birthorder_MZDZ,
  pheno_sensi_zygos_MZDZ = pheno_sensi_zygos_MZDZ,
  pheno_sensi_SSOS_DZ = pheno_sensi_SSOS_DZ,
  pheno_sensi_age_MZDZ = pheno_sensi_age_MZDZ
)

for (name in names(output_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = output_list[[name]], colNames = TRUE, rowNames = FALSE)
}

# Save the workbook
saveWorkbook(wb, "pheno_sensi.xlsx", overwrite = TRUE)

