# title: Plots
# author: Yujing Lin
# date: 10 Oct, 2024

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats) # reorder variables from largest to smallest for example

sourceFileStem <- '/Users/yujinglin/Desktop/WF&BF/WFBF_Cog_Results/' # the path where the results are saved
setwd('/Users/yujinglin/Desktop/') # the plots will be saved to this path 

bootMA_results <- read.csv(paste(sourceFileStem, "MA_WFBF_Cog_Main_boot.csv", sep=""))
bootMA_sensitivity <- read.csv(paste(sourceFileStem, "MA_WFBF_Cog_Main_sensi_boot.csv", sep=""))
bootMA_sensitivityMZDZ <- read.csv(paste(sourceFileStem, "MA_WFBF_Cog_Main_MZDZsensi_boot.csv", sep=""))
bootMA_SESresults <- read.csv(paste(sourceFileStem, "SESadj_MA_WFBF_Cog_Main_sensi_boot.csv", sep=""))

# compare between pop and WF weighted results: for bootMA_results
bootMA_compare_main_PopWF <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_Main.csv", sep=""))
# comparison for bootMA_sensitivity at pop level: 
bootMA_compare_WFBF_Cog_sensi_popFM <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_sensi_popFM.csv", sep=""))
# comparisons for bootMA_sensitivity at WF level: 3 comparisons due to 3 types of data:
bootMA_compare_WFBF_Cog_sensi_WF_FSSOS <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_FSSOS.csv", sep=""))
bootMA_compare_WFBF_Cog_sensi_WF_MSSOS <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_MSSOS.csv", sep=""))
bootMA_compare_WFBF_Cog_sensi_WF_FMSS <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_sensi_WF_FMSS.csv", sep=""))
# comparison for bootMA_sensitivityMZDZ
bootMA_compare_WFBF_Cog_MZDZsensi <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_MZDZsensi.csv", sep=""))
# comparison for bootMA_SESresults
bootMA_compare_WFBF_Cog_MZDZsensi <- read.csv(paste(sourceFileStem, "bootMA_compare_WFBF_Cog_MZDZsensi.csv", sep=""))
# comparison for the bootMA_SESresults
bootSESadjMA_compare_WFBF_Cog_Main <- read.csv(paste(sourceFileStem, "bootSESadjMA_compare_WFBF_Cog_Main.csv", sep=""))

# method comparison for population est
boot_compare_SimReg_Methods_unrelated_pairsum <- read.csv(paste(sourceFileStem, "boot_compare_SimReg_Methods_unrelated_pairsum.csv", sep=""))

bootMA_results <- bootMA_results %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("WF_DZpairdiff","WF_DZmeanpairdiff", "ppl_unrelated", "ppl_DZpairsum", "ppl_DZpairmean"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI")))

bootMA_sensitivity <- bootMA_sensitivity %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  mutate(Type = factor(Type, levels = c("WF_DZpairdiff","WF_DZmeanpairdiff", "ppl_unrelated", "ppl_DZpairsum", "ppl_DZpairmean"))) %>%
  mutate(Subsample = factor(Subsample, levels = c("Female", "Male", "Same-Sex Female", "Same-Sex Male", "Opposite-Sex")))%>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  arrange(Type, Subsample, V3, Trait_Domain)


bootMA_sensitivityMZDZ <- bootMA_sensitivityMZDZ %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(Type = factor(Type, levels = c("ppl_unrelated_MZ","ppl_unrelated_DZ"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) %>%
  arrange(Type, V3, Trait_Domain)

bootMA_SESresults <- bootMA_SESresults %>%
  mutate(Type = factor(Type, levels = c("WF_DZpairdiff_adjSES","ppl_unrelated_adjSES"))) %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c("childhood g", "adolescence g", "adulthood g", 
                                                        "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                                                        "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                                                        "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                                                        "A-level enrollment", "university enrollment", "years of schooling", 
                                                        "childhood BMI", "adolescence BMI", "adulthood BMI", 
                                                        "childhood height", "adolescence height", "adulthood height"))) 



# Fig. 2 (main text) cognitive domains: reg pop v reg WF ####
# yes, we start with Fig. 2, b/c Fig. 1 is the method illustration
MA_WFBF_SimCorr <- bootMA_results %>%
  subset(Type == "WF_DZpairdiff" | Type == "ppl_unrelated") %>%
  subset(V3 == "Cognitive Abilities" | V3 == "Education Achievement" | V3 == "Education Attainment" | V3 == "Height & BMI") %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  # Arrange Type so that "ppl_unrelated" (Population) comes before "WF_DZpairdiff" (Within-Family)
  arrange(V3, Trait_Domain, factor(Type, levels = c("ppl_unrelated", "WF_DZpairdiff")))

trait_order <- c("childhood g", "adolescence g", "adulthood g", 
                 "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
                 "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
                 "primary school grades", "GCSE grades", "A-level grades", "university grades", 
                 "A-level enrollment", "university enrollment", "years of schooling", 
                 "childhood BMI", "adolescence BMI", "adulthood BMI", 
                 "childhood height", "adolescence height", "adulthood height")

MA_WFBF_SimCorr <- MA_WFBF_SimCorr %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = trait_order))

asterisk_bootMA_compare_main_PopWF <- bootMA_compare_main_PopWF %>%
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

max_y_positions <- MA_WFBF_SimCorr %>%
  group_by(Trait_Domain) %>%
  summarise(
    max_ci_upper = max(ci_upper),
    max_beta = max(beta)
  )

asterisk_bootMA_compare_main_PopWF <- asterisk_bootMA_compare_main_PopWF %>%
  left_join(max_y_positions, by = "Trait_Domain")

MA_WFBF_SimCorr_plot <- ggplot(MA_WFBF_SimCorr, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.2), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.2)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.7, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.7, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.7, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.7, label = "BMI & Height", fontface = "bold", size = 7) +
  # coord_flip version:
  #geom_vline(xintercept = c(6.5, 9.5, 13.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  #annotate("text", x = 17.5, y = 0.8, label = "Cognitive \nAbilities", fontface = "bold", size = 7) + 
  #annotate("text", x = 11.5, y = 0.8, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  #annotate("text", x = 8, y = 0.8, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  #annotate("text", x = 3.5, y = 0.8, label = "BMI & Height", fontface = "bold", size = 7) +
  geom_text(
    data = asterisk_bootMA_compare_main_PopWF,
    aes(
      x = Trait_Domain, 
      y = max_ci_upper + 0.05,  # Slightly above the highest error bar
      label = asterisk
    ),
    position = position_dodge(width = 0.2),
    inherit.aes = FALSE,
    size = 8,
    color = "black"
  ) +
  
  theme_minimal() +
  #theme_pubclean() + 
  #coord_flip() +
  #scale_x_discrete(limits = rev(trait_order)) +
  # theme(legend.position = "top") +
  labs(x = "", y = "Standardized Beta Estimates", fill = "") +
  scale_y_continuous(
    limits = c(-0.05, max(MA_WFBF_SimCorr$ci_upper) + 0.1),  # Extend y-axis to accommodate asterisks
    breaks = seq(0, 0.6, 0.1)
  ) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                    labels = c("Within-Family Estimates (WF)", "Population Estimates (POP)")) +
  #scale_color_manual(values=c("Category1" = "#FF4040", "Category2" = "#009ACD", "Category3" = "#00FA9A")) +
  theme(
    axis.text=element_text(size=20), 
    strip.text.x=element_text(size=20), 
    legend.text=element_text(size=20), 
    axis.title=element_text(size=20),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

MA_WFBF_SimCorr_plot
ggsave("MA_WFBF_SimCorr_plot.png", plot = MA_WFBF_SimCorr_plot, width = 20, height = 10)




# Fig. 3 (main text) combined WFBF and .SES cognitive paper ####
bootMA_ratio <- bootMA_results %>% 
  filter(Type %in% c("WF_DZpairdiff", "ppl_unrelated")) %>%
  dplyr::select(-any_of(c("se", "ci_lower", "ci_upper", "beta_se"))) %>%
  pivot_wider(names_from = Type, values_from = beta, 
              values_fill = NA) %>%
  mutate(Ratio = WF_DZpairdiff / ppl_unrelated)

bootMA_SESratio <- bootMA_SESresults %>% 
  filter(Type %in% c("WF_DZpairdiff_adjSES", "ppl_unrelated_adjSES")) %>%
  dplyr::select(-c(`se`, `ci_lower`, `ci_upper`, `beta_se`)) %>%
  pivot_wider(names_from = Type, values_from = beta, 
              values_fill = NA) %>%
  mutate(Ratio_adjSES = `WF_DZpairdiff_adjSES` / `ppl_unrelated_adjSES`)

ratio_before_after_SES <- bind_rows(
  bootMA_ratio %>% 
    pivot_longer(cols = c(`Ratio`), 
                 names_to = "Type", 
                 values_to = "Ratio") %>%
    mutate(Status = "WF/Pop Ratio"),
  
  bootMA_SESratio %>% 
    pivot_longer(cols = c(`Ratio_adjSES`), 
                 names_to = "Type", 
                 values_to = "Ratio") %>%
    mutate(Status = "Ratio After Corrected for SES")
)

ratio_before_after_SES_Cog <- ratio_before_after_SES %>%
  dplyr::select(any_of(c("Trait_Domain", "V3", "Type", "Status", "Ratio"))) %>%
  filter(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(Status = factor(Status, levels = c("WF/Pop Ratio", "Ratio After Corrected for SES")))

asterisk_bootSESadjMA_compare_WFBF_Cog_Main <- bootSESadjMA_compare_WFBF_Cog_Main %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

asterisk_bootMA_compare_main_PopWF <- bootMA_compare_main_PopWF %>% # this version does not have ns
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

ratio_before_after_SES_Cog_plot <- ggplot(ratio_before_after_SES_Cog, aes(y = Ratio, x = Trait_Domain, shape = Status, fill = Status)) + 
  geom_point(size = 4, alpha = 0.85) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  
  geom_text(
    data = asterisk_bootMA_compare_main_PopWF,
    aes(
      x = Trait_Domain, 
      y = 0.2,  # Slightly above the highest error bar
      label = asterisk
    ),
    inherit.aes = FALSE,
    size = 12,
    color = "#8C4FAE"
  ) +
  
  geom_text(
    data = asterisk_bootSESadjMA_compare_WFBF_Cog_Main,
    aes(
      x = Trait_Domain, 
      y = 1.2,  # Slightly above the highest error bar
      label = asterisk
    ),
    inherit.aes = FALSE,
    size = 12,
    color = "#F5D44B"
  ) +
  
  theme_pubclean() + 
  labs(x = "", y = "Within-Family/Population Ratio", shape = "", fill = "") + 
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.2, 0.1)) + 
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Unadjusted Ratio", "SES-Adjusted Ratio")) + 
  scale_fill_manual(values = c("#8C4FAE", "#F5D44B"), #"#FF7518", "#00796B"
                    labels = c("Unadjusted Ratio", "SES-Adjusted Ratio")) + 
  theme(
    axis.text = element_text(size = 20), 
    strip.text.x = element_text(size = 20), 
    legend.text = element_text(size = 20), 
    axis.title = element_text(size = 20), 
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets 
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"), # Major grid lines
    panel.grid.minor = element_blank()
  ) + 
  coord_flip() + 
  scale_x_discrete(limits = rev(unique(ratio_before_after_SES_Cog$Trait_Domain)))
ratio_before_after_SES_Cog_plot

ggsave("ratio_before_after_SES_Cog_plot.jpg", plot = ratio_before_after_SES_Cog_plot, width = 12, height = 10)




# supplementary material figures ####
# Fig. S1 sensitivity reg pop cognitive: zygos ####
# no significant difference between MZ and DZ, so didn't add annotation for significance levels 
bootMA_sensitivityMZDZ_Plot <- ggplot(bootMA_sensitivityMZDZ, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.7, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.7, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.7, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.7, label = "BMI & Height", fontface = "bold", size = 7) +
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Population Standardized Beta Estimates", fill = "") +
  scale_y_continuous(limits=c(-0.05, 0.7), breaks=seq(0, 0.6, 0.1)) + 
  scale_fill_manual(values = c("#FF7518", "#00796B"),
                    labels = c("Unrelated MZ Twins", "Unrelated DZ Twins")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

bootMA_sensitivityMZDZ_Plot
ggsave("bootMA_sensitivityMZDZ_Plot.jpg", plot = bootMA_sensitivityMZDZ_Plot, width = 15, height = 10)




# Fig. S2 method comparison reg pop cog: unrelated v DZpairsum / DZpairmean ####
SimCorr_Method_Compare_Cog <- bootMA_results %>%
  subset(Type == "ppl_unrelated" | Type == "ppl_DZpairsum" ) %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI")

asterisk_boot_compare_SimReg_Methods_unrelated_pairsum <- boot_compare_SimReg_Methods_unrelated_pairsum %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

max_y_positions <- SimCorr_Method_Compare_Cog %>%
  group_by(Trait_Domain) %>%
  summarise(
    max_ci_upper = max(ci_upper),
    max_beta = max(beta)
  )

asterisk_boot_compare_SimReg_Methods_unrelated_pairsum <- asterisk_boot_compare_SimReg_Methods_unrelated_pairsum %>%
  left_join(max_y_positions, by = "Trait_Domain")

SimCorr_Method_Compare_Cog_Plot <- ggplot(SimCorr_Method_Compare_Cog, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.75, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.75, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.75, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.75, label = "BMI & Height", fontface = "bold", size = 7) +
  
  geom_text(
    data = asterisk_boot_compare_SimReg_Methods_unrelated_pairsum,
    aes(
      x = Trait_Domain, 
      y = max_ci_upper + 0.05,  # Slightly above the highest error bar
      label = asterisk
    ),
    position = position_dodge(width = 0.2),
    inherit.aes = FALSE,
    size = 8,
    color = "black"
  ) +
  
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Population Standardized Beta Estimates", fill = "") +
  scale_y_continuous(
    limits = c(-0.05, max(SimCorr_Method_Compare_Cog$ci_upper) + 0.1),  # Extend y-axis to accommodate asterisks
    breaks = seq(0, 0.7, 0.1)
  ) +
  scale_fill_manual(values = c("#F8766D", "#FFB84D"),
                    labels = c("Unrelated Population", "DZ Pair Sum")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

SimCorr_Method_Compare_Cog_Plot
ggsave("SimCorr_Method_Compare_Cog_Plot.jpg", plot = SimCorr_Method_Compare_Cog_Plot, width = 15, height = 10)




# Fig. S3 method comparison reg WF Cog: unrelated v DZpairsum / DZpairmean ####
SimCorrWF_Method_Compare_Cog <- bootMA_results %>%
  subset(Type == "WF_DZpairdiff" | Type == "WF_DZmeanpairdiff" ) %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI")

SimCorrWF_Method_Compare_Cog_Plot <- ggplot(SimCorrWF_Method_Compare_Cog, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.65, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.65, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.65, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.65, label = "BMI & Height", fontface = "bold", size = 7) +
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Within-Family Standardized Beta Estimates", fill = "") +
  scale_y_continuous(limits=c(-0.05, 0.65), breaks=seq(0, 0.6, 0.1)) + 
  scale_fill_manual(values = c("#00BFC4", "blue"),
                    labels = c("DZ Pair Differences", "DZ Pair Deviations")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

SimCorrWF_Method_Compare_Cog_Plot
ggsave("SimCorrWF_Method_Compare_Cog_Plot.jpg", plot = SimCorrWF_Method_Compare_Cog_Plot, width = 15, height = 10)




# Fig. S4a & S4b SES-trait & SES-PGS correlations: the codes are in the SES script ####




# Fig. S5 .SES cognitive domains: reg pop v reg WF ####
MA_WFBF_SimCorrSES <- bootMA_SESresults %>%
  subset(Type == "WF_DZpairdiff_adjSES" | Type == "ppl_unrelated_adjSES") %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI"))) %>%
  # Arrange Type so that "ppl_unrelated" (Population) comes before "WF_DZpairdiff" (Within-Family)
  arrange(V3, Trait_Domain, factor(Type, levels = c("ppl_unrelated_adjSES", "WF_DZpairdiff_adjSES")))

asterisk_bootSESadjMA_compare_WFBF_Cog_Main <- bootSESadjMA_compare_WFBF_Cog_Main %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

max_y_positions <- MA_WFBF_SimCorrSES %>%
  group_by(Trait_Domain) %>%
  summarise(
    max_ci_upper = max(ci_upper),
    max_beta = max(beta)
  )

asterisk_bootSESadjMA_compare_WFBF_Cog_Main <- asterisk_bootSESadjMA_compare_WFBF_Cog_Main %>%
  left_join(max_y_positions, by = "Trait_Domain")

MA_WFBF_SimCorrSES_plot <- ggplot(MA_WFBF_SimCorrSES, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.2), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.2)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.7, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.7, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.7, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.7, label = "BMI & Height", fontface = "bold", size = 7) +
  
  geom_text(
    data = asterisk_bootSESadjMA_compare_WFBF_Cog_Main,
    aes(
      x = Trait_Domain, 
      y = max_ci_upper + 0.05,  # Slightly above the highest error bar
      label = asterisk
    ),
    position = position_dodge(width = 0.2),
    inherit.aes = FALSE,
    size = 8,
    color = "black"
  ) +
  
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Standardized Beta Estimates", fill = "") +
  scale_y_continuous(
    limits = c(-0.05, max(MA_WFBF_SimCorr$ci_upper) + 0.1),  # Extend y-axis to accommodate asterisks
    breaks = seq(0, 0.6, 0.1)
  ) +
  scale_fill_manual(values = c("#00A9B8", "#E85C5B"),
                    labels = c("Within-Family Estimates (WF)", "Population Estimates (POP)")) +
  #scale_color_manual(values=c("Category1" = "#FF4040", "Category2" = "#009ACD", "Category3" = "#00FA9A")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

MA_WFBF_SimCorrSES_plot
ggsave("MA_WFBF_SimCorrSES_plot.jpg", plot = MA_WFBF_SimCorrSES_plot, width = 15, height = 10)




# Fig. S6-S12 quadratic plot and decile plots: codes in the extreme analyses script ####




# Fig. S13 sensitivity reg pop cognitive: female v male ####
# no significant differences between females and males, do didn't add annotation for significance levels
Sensitivity_Pop <- bootMA_sensitivity %>%
  subset(Subsample == "Female" | Subsample == "Male") %>%
  subset(Type == "ppl_unrelated") %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI")))

Sensitivity_Pop_Plot <- ggplot(Sensitivity_Pop, aes(y = beta, x = Trait_Domain, fill = Subsample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.65, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.65, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.65, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.65, label = "BMI & Height", fontface = "bold", size = 7) +
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Population Standardized Beta Estimates", fill = "") +
  scale_y_continuous(limits=c(-0.05, 0.65), breaks=seq(0, 0.6, 0.1)) + 
  scale_fill_manual(values = c("#FFD700", "#C0C0C0"),
                    labels = c("Unrelated Females", "Unrelated Males")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

Sensitivity_Pop_Plot
ggsave("Sensitivity_Pop_Plot.jpg", plot = Sensitivity_Pop_Plot, width = 15, height = 10)




# Fig. S14 sensitivity reg WF: SS-F, SS-M, OS ####
# only two estimate points significantly differ between SS-M and OS for childhood and adolescence height 
Sensitivity_WF_Cog <- bootMA_sensitivity %>%
  subset(Subsample == "Same-Sex Female" | Subsample == "Same-Sex Male" | Subsample == "Opposite-Sex") %>%
  mutate(Subsample = factor(Subsample, levels = c("Same-Sex Female", "Same-Sex Male", "Opposite-Sex"))) %>%
  subset(Type == "WF_DZpairdiff") %>%
  subset(V3 == "Cognitive Abilities" |V3 == "Education Achievement"|V3 == "Education Attainment"|V3 == "Height & BMI") %>%
  mutate(V3 = factor(V3, levels = c("Cognitive Abilities","Education Achievement", "Education Attainment", "Height & BMI")))

asterisk_bootMA_compare_WFBF_Cog_sensi_WF_MSSOS <- bootMA_compare_WFBF_Cog_sensi_WF_MSSOS %>%
  filter(V3 %in% c("Cognitive Abilities", "Education Achievement", "Education Attainment", "Height & BMI")) %>%
  group_by(Trait_Domain) %>%
  mutate(
    asterisk = case_when(
      diff_p_fdr < 0.001 ~ "***",
      diff_p_fdr < 0.01 ~ "**",
      diff_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

max_y_positions <- Sensitivity_WF_Cog %>%
  group_by(Trait_Domain) %>%
  summarise(
    max_ci_upper = max(ci_upper),
    max_beta = max(beta)
  )

asterisk_bootMA_compare_WFBF_Cog_sensi_WF_MSSOS <- asterisk_bootMA_compare_WFBF_Cog_sensi_WF_MSSOS %>%
  left_join(max_y_positions, by = "Trait_Domain")

Sensitivity_WF_Cog_Plot <- ggplot(Sensitivity_WF_Cog, aes(y = beta, x = Trait_Domain, fill = Subsample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(9.5, 13.5, 16.5), linetype = "dashed", color = "gray") +  # Light vertical lines for separation
  annotate("text", x = 5, y = 0.65, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 11.5, y = 0.65, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 15, y = 0.65, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 19.5, y = 0.65, label = "BMI & Height", fontface = "bold", size = 7) +
  
  geom_text(
    data = asterisk_bootMA_compare_WFBF_Cog_sensi_WF_MSSOS,
    aes(
      x = Trait_Domain, 
      y = max_ci_upper + 0.05,  # Slightly above the highest error bar
      label = asterisk
    ),
    position = position_dodge(width = 0.3),
    inherit.aes = FALSE,
    size = 8,
    color = "black"
  ) +
  
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Within-Family Standardized Beta Estimates", fill = "") +
  scale_y_continuous(limits=c(-0.2, 0.7), breaks=seq(-0.2, 0.6, 0.1)) + 
  scale_fill_manual(values = c("#FFD700", "#C0C0C0", "#228B22"),
                    labels = c("DZ Same-Sex Females", "DZ Same-Sex Males", "DZ Opposite-Sex Twins")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

Sensitivity_WF_Cog_Plot
ggsave("Sensitivity_WF_Cog_Plot.jpg", plot = Sensitivity_WF_Cog_Plot, width = 15, height = 10)










#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

MA_LMM_within_boot <- read.csv(paste(sourceFileStem, "bootMA_WFBF_Cog_LMM_within.csv", sep=""))
MA_LMM_pop_boot <- read.csv(paste(sourceFileStem, "bootMA_WFBF_Cog_LMM_pop.csv", sep=""))

# cognitive domains: LMM pop v LMM WF ####
MA_LMM_within_boot_WF <- MA_LMM_within_boot %>%
  transmute(
    Trait_Domain,
    V3,
    beta = beta_WF,
    ci_lower = ci_lower_WF,
    ci_upper = ci_upper_WF,
    Type = "LMM_WF"
  )

MA_LMM_within_boot_BF <- MA_LMM_within_boot %>%
  transmute(
    Trait_Domain,
    V3,
    beta = beta_BF,
    ci_lower = ci_lower_BF,
    ci_upper = ci_upper_BF,
    Type = "LMM_BF"
  )

MA_LMM_within_boot_long <- bind_rows(MA_LMM_within_boot_WF, MA_LMM_within_boot_BF)
MA_LMM_within_boot_long <- MA_LMM_within_boot_long %>%
  mutate(Type = factor(Type, levels = c("LMM_WF","LMM_BF")))

MA_LMM_pop_modified <- MA_LMM_pop_boot %>% select(Trait_Domain, V3, beta, ci_lower, ci_upper, Type)
MA_LMM_pop_within <- bind_rows(MA_LMM_within_boot_long, MA_LMM_pop_modified) 

MA_LMM_pop_within <- subset(MA_LMM_pop_within, V3 %in% c("Height & BMI", "Cognitive Abilities", "Education Achievement", "Education Attainment"))

MA_LMM_pop_within <- MA_LMM_pop_within %>%
  mutate(Type = factor(Type, levels = c("LMM_pop", "LMM_WF", "LMM_BF")))%>%
  mutate(V3 = factor(V3, levels = c("Height & BMI", "Cognitive Abilities","Education Achievement", "Education Attainment"))) %>%
  mutate(Trait_Domain = factor(Trait_Domain, levels = c(#"anxiety", "depression", "mania", "autism", "ADHD", "PTSD", "alcohol", "eating disorder", 
    "childhood BMI", "adolescence BMI", "adulthood BMI", 
    "childhood height", "adolescence height", "adulthood height",
    "childhood g", "adolescence g", "adulthood g", 
    "childhood verbal g", "adolescence verbal g", "adulthood verbal g", 
    "childhood nonverbal g", "adolescence nonverbal g", "adulthood nonverbal g", 
    "primary school grades", "GCSE grades", "A-level grades", "university grades", 
    "A-level enrollment", "university enrollment", "years of schooling")))%>%
  arrange(Type, V3, Trait_Domain)

MA_LMM_within_pop_boot_coghtBMI_long_plot <- ggplot(MA_LMM_pop_within, aes(y = beta, x = Trait_Domain, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.8)) +
  # facet_wrap(~V3, scales="free_x", nrow = 1) +
  geom_vline(xintercept = c(6.5, 15.5, 19.5), linetype = "dashed", color = "gray") +
  annotate("text", x = 11, y = 5.5, label = "Cognitive Abilities", fontface = "bold", size = 7) + 
  annotate("text", x = 17.5, y = 5.5, label = "Educational \nAchievement", fontface = "bold", size = 7) + 
  annotate("text", x = 21, y = 5.5, label = "Educational\nAttainment", fontface = "bold", size = 7) +
  annotate("text", x = 3.5, y = 5.5, label = "BMI & Height", fontface = "bold", size = 7) +
  theme_pubclean() +
  # theme(legend.position = "top") +
  labs(x = "", y = "Mixed-Effects Model Beta/Log-odds Coefficients", fill = "") +
  scale_y_continuous(limits=c(-0.2, 5.5), breaks=seq(0, 5.5, 0.5)) + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#8B5FFF"),
                   labels = c("Population (Unrelated Twin) Coefficient", "DZ Twin Pairwise Deviation Coefficient", "DZ Twin Pairwise Mean Coefficient")) +
  #scale_color_manual(values=c("Category1" = "#FF4040", "Category2" = "#009ACD", "Category3" = "#00FA9A")) +
  theme(
    axis.text=element_text(size=16), 
    strip.text.x=element_text(size=16), 
    legend.text=element_text(size=16), 
    axis.title=element_text(size=16),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between facets
    axis.text.x = element_text(angle = 35, hjust = 1))

MA_LMM_within_pop_boot_coghtBMI_long_plot
ggsave("MA_LMM_within_boot_coghtBMI_plot.jpg", plot = MA_LMM_within_pop_boot_coghtBMI_long_plot, width = 20, height = 10)

