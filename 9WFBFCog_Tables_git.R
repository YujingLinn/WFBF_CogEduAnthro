# title: Tables
# author: Yujing Lin
# date: 10 Oct, 2024

# to reproduce raw result tables: because originally I have to keep the trait domain the same for each category for weighted averages
# but when publishing the results, I have to add the annotation for age to differentiate among the variables within the same trait domain

library(dplyr)
library(openxlsx)

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/'
outFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results_SupTables/'

WFBF_Cog_Main_raw <- read.csv(paste(sourceFileStem, "WFBF_Cog_Main_boot_raw.csv", sep=""))
SEScorr_wTrait <- read.csv(paste(sourceFileStem, "SEScorr_wTrait.csv", sep=""))
SEScorr_wPGS <- read.csv(paste(sourceFileStem, "SEScorr_wPGS.csv", sep=""))
SESadj_WFBF_Cog_Main_boot_raw <- read.csv(paste(sourceFileStem, "SESadj_WFBF_Cog_Main_boot_raw.csv", sep=""))
nonlinear_ppl_unrelated_raw <- read.csv(paste(sourceFileStem, "nonlinear_ppl_unrelated_raw.csv", sep=""))
nonlinear_WF_DZpairdiff_raw <- read.csv(paste(sourceFileStem, "nonlinear_WF_DZpairdiff_raw.csv", sep=""))
LMM_raw_pop <- read.csv(paste(sourceFileStem, "WFBF_Cog_LMM_pop_boot_raw.csv", sep=""))
LMM_raw_within <- read.csv(paste(sourceFileStem, "WFBF_Cog_LMM_within_boot_raw.csv", sep=""))

boot_compare_WFBF_Cog_Main_raw <- read.csv(paste(sourceFileStem, "boot_compare_WFBF_Cog_Main_raw.csv", sep=""))
bootSESadj_compare_WFBF_Cog_Main <- read.csv(paste(sourceFileStem, "bootSESadj_compare_WFBF_Cog_Main.csv", sep=""))

WFBF_Cog_Main_raw_pop <- WFBF_Cog_Main_raw %>%
  filter(Type == "ppl_unrelated")
WFBF_Cog_Main_raw_WF <- WFBF_Cog_Main_raw %>%
  filter(Type == "WF_DZpairdiff")
SESadj_WFBF_Cog_Main_boot_raw_pop <- SESadj_WFBF_Cog_Main_boot_raw %>%
  filter(Type == "ppl_unrelated_adjSES")
SESadj_WFBF_Cog_Main_boot_raw_WF <- SESadj_WFBF_Cog_Main_boot_raw %>%
  filter(Type == "WF_DZpairdiff_adjSES")

# process the results ####
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
  c("ra_level_enrol", "A-level enrollment age 16", "Education Attainment"), # A-level enrolment (16), data obtained at 18
  c("rcqhe", "university enrollment age 18", "Education Attainment"), # university enrolment (18)
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

WFBF_Cog_Main_raw_pop <- full_join(raw_Trait_Domain, WFBF_Cog_Main_raw_pop)
WFBF_Cog_Main_raw_WF <- full_join(raw_Trait_Domain, WFBF_Cog_Main_raw_WF)
SEScorr_wTrait <- full_join(raw_Trait_Domain, SEScorr_wTrait)
SESadj_WFBF_Cog_Main_boot_raw_pop <- full_join(raw_Trait_Domain, SESadj_WFBF_Cog_Main_boot_raw_pop)
SESadj_WFBF_Cog_Main_boot_raw_WF <- full_join(raw_Trait_Domain, SESadj_WFBF_Cog_Main_boot_raw_WF)
LMM_raw_pop <- full_join(raw_Trait_Domain, LMM_raw_pop)
LMM_raw_within <- full_join(raw_Trait_Domain, LMM_raw_within)

boot_compare_WFBF_Cog_Main_raw <- full_join(raw_Trait_Domain, boot_compare_WFBF_Cog_Main_raw)
boot_compare_WFBF_Cog_Main_raw <- boot_compare_WFBF_Cog_Main_raw %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

bootSESadj_compare_WFBF_Cog_Main <- full_join(raw_Trait_Domain, bootSESadj_compare_WFBF_Cog_Main)
bootSESadj_compare_WFBF_Cog_Main <- bootSESadj_compare_WFBF_Cog_Main %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()



wb <- createWorkbook()
addWorksheet(wb, sheetName="SupT3.Main_raw_pop")
addWorksheet(wb, sheetName="SupT4.Main_raw_WF")
addWorksheet(wb, sheetName="SupT5.SEScorr_wTrait")
addWorksheet(wb, sheetName="SupT6.SEScorr_wPGS")
addWorksheet(wb, sheetName="SupT7.SESadj_raw_pop")
addWorksheet(wb, sheetName="SupT8.SESadj_raw_WF")
addWorksheet(wb, sheetName="SupT9a.nonlinear_pop")
addWorksheet(wb, sheetName="SupT9b.nonlinear_within")
addWorksheet(wb, sheetName="SupT10.LMM_raw_pop")
addWorksheet(wb, sheetName="SupT11.LMM_raw_within")

writeData(wb, sheet="SupT3.Main_raw_pop", x=WFBF_Cog_Main_raw_pop, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT4.Main_raw_WF", x=WFBF_Cog_Main_raw_WF, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT5.SEScorr_wTrait", x=SEScorr_wTrait, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT6.SEScorr_wPGS", x=SEScorr_wPGS, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT7.SESadj_raw_pop", x=SESadj_WFBF_Cog_Main_boot_raw_pop, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT8.SESadj_raw_WF", x=SESadj_WFBF_Cog_Main_boot_raw_WF, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT9a.nonlinear_pop", x=nonlinear_ppl_unrelated_raw, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT9b.nonlinear_within", x=nonlinear_WF_DZpairdiff_raw, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT10.LMM_raw_pop", x=LMM_raw_pop, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="SupT11.LMM_raw_within", x=LMM_raw_within, colNames=TRUE, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, paste(outFileStem, "SupT3_T11.xlsx", sep=""), overwrite=TRUE)


wb <- createWorkbook()
addWorksheet(wb, sheetName="Sup_boot_compare_Main")
addWorksheet(wb, sheetName="Sup_boot_compare_SESadj")
writeData(wb, sheet="Sup_boot_compare_Main", x=boot_compare_WFBF_Cog_Main_raw, colNames=TRUE, rowNames=FALSE)
writeData(wb, sheet="Sup_boot_compare_SESadj", x=bootSESadj_compare_WFBF_Cog_Main, colNames=TRUE, rowNames=FALSE)
# Save the workbook with raw traits
saveWorkbook(wb, paste(outFileStem, "Sup_boot_compare.xlsx", sep=""), overwrite=TRUE)



