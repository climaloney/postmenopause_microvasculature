# MENOHRT.R
# Analysis of Menopause and Hormone Replacement Therapy effects on CVC
# Author: [Your Name]
# Date: [Current Date]

setwd("C:/Users/carma/Documents/R/MD")
options(scipen = 9999)

### LIBRARY LOADING ####
library(nlme)
library(car)
library(tidyverse)
library(ggplot2)
library(lme4)
library(emmeans)
library(ggthemes)
library(broom)
library(performance)


### DATA PREPARATION ####

# Load participant data
participant <- read.csv("MENO2_Participant.csv", header = TRUE)
participant$YearsMeno <- participant$Age - participant$MenopauseOnset

# Load and process main datasets
DATA <- read.csv("MENO2_DATA.csv", header = TRUE)
SITE <- read.csv("MENO2_SITES.csv", header = TRUE)
CVC <- read.csv("MENO2_CV.csv", header = TRUE)
rep<-read.csv("hrtrep.csv", header = TRUE)

  # Reshape CVC data from wide to long format
CVC_P<-CVC
CVC_P$POOLED<-rowMeans(CVC_P[, 4:7], na.rm = TRUE)


CVC_P <- CVC_P%>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF, POOLED), names_to = "Region", values_to = "CVCMax") %>%
  mutate(across(where(is.character), as.factor))


CVC <- CVC %>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF), names_to = "Region", values_to = "CVC") %>%
  mutate(across(where(is.character), as.factor)) %>%
  pivot_wider(names_from = "Phase", values_from = "CVC")

# Create pooled measure and reshape DATA
DATA_E <- DATA
DATA_E$POOLED <- rowMeans(DATA_E[, 4:7], na.rm = TRUE)


DATA_E <- DATA_E %>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF, POOLED), names_to = "Region", values_to = "CVCMax") %>%
  mutate(across(where(is.character), as.factor))

DATA <- DATA %>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF), names_to = "Region", values_to = "CVCMax") %>%
  mutate(across(where(is.character), as.factor)) %>%
  pivot_wider(names_from = "Phase", values_from = "CVCMax")

# Process site data
SITE <- SITE %>%
  pivot_longer(
    cols = c(FOREARM, CHEST, ABDOMEN, CALF),
    names_to = "Region",
    values_to = "Site_Order"
  ) %>%
  mutate(Region = toupper(Region))

# Merge data with site information
DATA <- DATA %>%
  left_join(SITE, by = c("ID" = "Participant", "Region", "Group" = "GROUP")) %>%
  mutate(across(where(is.character), as.factor))

CVC <- CVC %>%
  left_join(SITE, by = c("ID" = "Participant", "Region", "Group" = "GROUP")) %>%
  mutate(across(where(is.character), as.factor))

# Set reference levels
DATA$Group <- relevel(DATA$Group, ref = "NON")
DATA$Region <- relevel(DATA$Region, ref = "FOREARM")

### SUMMARY STATISTICS ####

# Participant characteristics summary
S_STAT <- participant %>%
  group_by(GROUP) %>%
  summarise(across(c(Age, YearsMeno, Height, Weight, BF.,Vo2max,MAP,BL.HR,BL.ST),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        min = ~min(., na.rm = TRUE),
                        max = ~max(., na.rm = TRUE)))) %>%
  pivot_longer(cols = -GROUP, names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "statistic"), sep = "_") %>%
  pivot_wider(names_from = statistic, values_from = value)

# CVC summary statistics by phase, group, and region
STAT <- DATA_E %>%
  group_by(Phase, Group, Region) %>%
  summarise_at(vars(CVCMax), list(mean = mean, sd = sd, min = min, max = max),
               na.rm = TRUE)

STAT_C <- CVC_P %>%
  group_by(Phase, Group, Region) %>%
  summarise_at(vars(CVCMax), list(mean = mean, sd = sd, min = min, max = max),
               na.rm = TRUE)
# Statistical tests for group differences in anthropometric measures
t.test(participant$Height ~ participant$GROUP)
t.test(participant$Weight ~ participant$GROUP)
t.test(participant$BF. ~ participant$GROUP)
t.test(participant$Vo2max ~ participant$GROUP)
t.test(participant$Age ~ participant$GROUP)

# Statistical tests for cardiovascular measures
t.test(participant$BL.HR ~ participant$GROUP)
t.test(participant$MAP ~ participant$GROUP)
t.test(participant$End.HR ~ participant$GROUP)

# Statistical tests for other measures
t.test(participant$USG ~ participant$GROUP)
t.test(participant$BL.ST ~ participant$GROUP)
t.test(participant$MenopauseOnset ~ participant$GROUP)
t.test(participant$YearsMeno ~ participant$GROUP)

# List of variables to test
variables <- c("Height", "Weight", "BF.", "Vo2max", "Age", "BL.HR", "MAP", 
               "End.HR", "USG", "BL.ST", "MenopauseOnset", "YearsMeno")

# Create dataframe using map and broom::tidy
t_test_df <- map_dfr(variables, ~{
  t.test(as.formula(paste(.x, "~ GROUP")), data = participant) %>% 
    tidy() %>%
    mutate(variable = .x)
}) %>%
  select(variable, everything())



### LINEAR MIXED EFFECTS MODELS ####

#### CVCMAX MODELS ####

# Baseline model with interaction and variance structure
BL_final <- lme(BL ~ Group * Region,
                random = ~ 1 | ID,
                weights = varIdent(form = ~ 1 | Site_Order),
                data = DATA,
                na.action = na.exclude,
                method = "REML")
plot(BL_final)

# Heating phase model adjusting for baseline
HEAT_final <- lme(HEAT ~  Group * Region,
                  random = ~ 1 | ID,
                  weights = varIdent(form = ~ 1 | Site_Order),
                  data = DATA,
                  na.action = na.exclude)

# L-NAME phase model
LNA_final <- lme(LNAME ~  Group * Region,
                 random = ~ 1 | ID,
                 weights = varIdent(form = ~ 1 | Site_Order),
                 data = DATA,
                 na.action = na.exclude)

# NOS phase model
NOS_final <- lme(NOS ~  Group * Region,
                 random = ~ 1 | ID,
                 weights = varIdent(form = ~ 1 | Site_Order),
                 data = DATA,
                 na.action = na.exclude)

###PEAK
PEAK_final<-lme(PEAK~  Group * Region,
                random = ~ 1 | ID,
                weights = varIdent(form = ~ 1 | Site_Order),
                data = DATA,
                na.action = na.exclude)
summary(PEAK_final)

# Models without weighting for comparison
BL_2 <- lme(BL ~ Group * Region,
            random = ~ 1 | ID,
            data = DATA,
            na.action = na.exclude)
HEAT_2 <- lme(HEAT ~ BL + Region,
              random = ~ 1 | ID,
              data = DATA[which(DATA$Group == "HRT"),],
              na.action = na.exclude)

NOS_2 <- lme(NOS ~ BL + Group * Region,
             random = ~ 1 | ID,
             data = DATA,
             na.action = na.exclude)
HEAT_2
# Model summaries
summary(BL_final)
summary(BL_2)
summary(HEAT_final)
summary(HEAT_2)
summary(NOS_final)
summary(NOS_2)
summary(LNA_final)

###Assumption test
plot(BL_final)
plot(BL_2)
plot(HEAT_final)
plot(HEAT_2)
qqnorm(residuals(BL_final, type = "normalized"))
qqnorm(residuals(BL_2, type = "normalized"))
shapiro.test(residuals(BL_final, type = "normalized"))
shapiro.test(residuals(BL_2, type = "normalized"))
qqnorm(residuals(HEAT_final, type = "normalized"))
qqnorm(residuals(HEAT_2, type = "normalized"))
shapiro.test(residuals(HEAT_final, type = "normalized"))
shapiro.test(residuals(HEAT_2, type = "normalized"))
qqnorm(residuals(NOS_final))
qqnorm(residuals(NOS_2))
shapiro.test(residuals(NOS_final, type = "normalized"))
shapiro.test(residuals(NOS_2, type = "normalized"))

### MARGINAL MEANS AND POST-HOC COMPARISONS ####

# Calculate estimated marginal means for each phase
X1 <- emmeans(BL_final, ~ Region * Group, vcov. = vcov(BL_final))
X2 <- emmeans(HEAT_final, ~ Region * Group, vcov. = vcov(HEAT_final))
X3 <- emmeans(NOS_final, ~ Region * Group, vcov. = vcov(NOS_final))
X4 <- emmeans(LNA_final, ~ Region * Group, vcov. = vcov(LNA_final))
X5 <- emmeans(PEAK_final, ~ Region * Group, vcov. = vcov(PEAK_final))

# Pairwise comparisons within groups
B1 <- as.data.frame(pairs(X1, by = "Group", adjust = "holm"))
H1 <- as.data.frame(pairs(X2, by = "Group", adjust = "holm"))
N1 <- as.data.frame(pairs(X3, by = "Group", adjust = "holm"))
L1 <- as.data.frame(pairs(X4, by = "Group", adjust = "holm"))
P1<-as.data.frame(pairs(X5, by = "Group", adjust = "holm"))
P1
# Pairwise comparisons within regions
B2 <- as.data.frame(pairs(X1, by = "Region", adjust = "holm"))
H2 <- as.data.frame(pairs(X2, by = "Region", adjust = "holm"))
N2 <- as.data.frame(pairs(X3, by = "Region", adjust = "holm"))
L2 <- as.data.frame(pairs(X4, by = "Region", adjust = "holm"))
P2<-as.data.frame(pairs(X5, by = "Region", adjust = "holm"))
P2
# Confidence intervals for comparisons
B3 <- as.data.frame(confint(pairs(X1, by = "Group", adjust = "holm")))
H3 <- as.data.frame(confint(pairs(X2, by = "Group", adjust = "holm")))
N3 <- as.data.frame(confint(pairs(X3, by = "Group", adjust = "holm")))
L3 <- as.data.frame(confint(pairs(X4, by = "Group", adjust = "holm")))
P3<- as.data.frame(confint(pairs(X5, by = "Group", adjust = "holm")))
B4 <- as.data.frame(confint(pairs(X1, by = "Region", adjust = "holm")))
H4 <- as.data.frame(confint(pairs(X2, by = "Region", adjust = "holm")))
N4 <- as.data.frame(confint(pairs(X3, by = "Region", adjust = "holm")))
L4 <- as.data.frame(confint(pairs(X4, by = "Region", adjust = "holm")))

B1$PHASE<-"BL"
B2$PHASE<-"BL"
B3$PHASE<-"BL"
B4$PHASE<-"BL"
H1$PHASE<-"HEAT"
H2$PHASE<-"HEAT"
H3$PHASE<-"HEAT"
H4$PHASE<-"HEAT"
N1$PHASE<-"NOS"
N2$PHASE<-"NOS"
N3$PHASE<-"NOS"
N4$PHASE<-"NOS"
L1$PHASE<-"LNAME"
L2$PHASE<-"LNAME"
L3$PHASE<-"LNAME"
L4$PHASE<-"LNAME"


# Merge and process within-group comparisons
WITHIN1 <- merge(B1, B3, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))
WITHIN2 <- merge(H1, H3, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))
WITHIN3 <- merge(L1, L3, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))
WITHIN4 <- merge(N1, N3, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))

WITH <- rbind(WITHIN1, WITHIN2, WITHIN3, WITHIN4)

# Reverse specific contrasts for consistent interpretation
rows_to_change <- c(3, 4, 7, 8, 15, 16, 19, 20, 27, 28, 31, 32, 39, 40, 43, 44)
WITH$estimate[rows_to_change] <- -WITH$estimate[rows_to_change]
WITH$lower.CL[rows_to_change] <- -WITH$lower.CL[rows_to_change]
WITH$upper.CL[rows_to_change] <- -WITH$upper.CL[rows_to_change]
WITH$t.ratio[rows_to_change] <- -WITH$t.ratio[rows_to_change]

# Recode factors for better visualization
WITH <- WITH %>%
  mutate(PHASE = recode_factor(PHASE,
                               BL = "Baseline",
                               HEAT = "Heating Plateau",
                               LNAME = "L-NAME Plateau",
                               NOS = "Nitric Oxide Contribution")) %>%
  mutate(Region = recode_factor(contrast,
                                "ABDOMEN - CALF" = "Abdomen - Calf",
                                "ABDOMEN - CHEST" = "Abdomen - Chest",
                                "FOREARM - ABDOMEN" = "Abdomen - Forearm",
                                "CALF - CHEST" = "Calf - Chest",
                                "FOREARM - CALF" = "Forearm - Calf",
                                "FOREARM - CHEST" = "Forearm - Chest")) %>%
  mutate(Group = recode_factor(Group,
                               HRT = "Hormone Replacement Therapy",
                               NON = "Control"))

# Process between-group comparisons
TWEEN1 <- merge(B2, B4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN2 <- merge(H2, H4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN3 <- merge(L2, L4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN4 <- merge(N2, N4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))

TWEEN <- rbind(TWEEN1, TWEEN2, TWEEN3, TWEEN4)
B2
# Overall group comparisons
M <- emmeans(BL_final, ~ Group)
M1 <- emmeans(HEAT_final, ~ Group)
M2 <- emmeans(LNA_final, ~ Group)
M3 <- emmeans(NOS_final, ~ Group)
M4<-emmeans(PEAK_final, ~Group)
M4
E1 <- as.data.frame(confint(pairs(M)))
E2 <- as.data.frame(confint(pairs(M1)))
E3 <- as.data.frame(confint(pairs(M2)))
E4 <- as.data.frame(confint(pairs(M3)))
E5 <- as.data.frame(confint(pairs(M4)))

E <- rbind(E1, E2, E3, E4)
E$PHASE <- c("BL", "HEAT", "LNAME", "NOS")

# Add overall comparison to between-group dataframe
E_for_bind <- E %>%
  mutate(Region = "POOLED") %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL)

TWEEN_for_bind <- TWEEN %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL)

TWEEN <- rbind(TWEEN_for_bind, E_for_bind)

X1 <- emmeans(BL_final, ~ Region, vcov. = vcov(BL_final))
X2 <- emmeans(HEAT_final, ~ Region, vcov. = vcov(HEAT_final))
X3 <- emmeans(NOS_final, ~ Region, vcov. = vcov(NOS_final))
X4 <- emmeans(LNA_final, ~ Region, vcov. = vcov(NOS_final))
X5 <- emmeans(PEAK_final, ~ Region, vcov. = vcov(PEAK_final))

B2 <- as.data.frame(confint(contrast(X1, method = "pairwise", adjust = "holm")))
B3 <- as.data.frame(confint(contrast(X2, method = "pairwise", adjust = "holm")))
B4 <- as.data.frame(confint(contrast(X3, method = "pairwise", adjust = "holm")))
B5 <- as.data.frame(confint(contrast(X4, method = "pairwise", adjust = "holm")))
B6<-as.data.frame(confint(contrast(X5, method = "pairwise", adjust = "holm")))
B2
stufF<-rbind(B2,B3,B4,B5)
#### CVC MODELS ####

BLCVC <- lme(BL ~ Group * Region,
             random = ~ 1 | ID,
             weights = varIdent(form = ~ 1 | Site_Order),
             data = CVC,
             na.action = na.exclude)

MAX <- lme(MAX ~  Group * Region,
           random = ~ 1 | ID,
           weights = varIdent(form = ~ 1 | Site_Order),
           data = CVC,
           na.action = na.exclude)

summary(MAX)
summary(BLCVC)


# Calculate estimated marginal means for CVC phases
X1_cvc <- emmeans(BLCVC, ~ Group*Region, vcov. = vcov(BLCVC))
X2_cvc <- emmeans(MAX, ~ Group, vcov. = vcov(MAX))
X1_cvc
# Pairwise comparisons within groups
B1_cvc <- as.data.frame(pairs(X1_cvc, adjust = "holm"))
H1_cvc <- as.data.frame(pairs(X2_cvc, adjust = "holm"))
B1_cvc
H1_cvc
# Pairwise comparisons within regions
B2_cvc <- as.data.frame(pairs(X1_cvc, adjust = "holm"))
H2_cvc <- as.data.frame(pairs(X2_cvc,  adjust = "holm"))
B2_cvc
# Confidence intervals for comparisons
B3_cvc <- as.data.frame(confint(pairs(X1_cvc,adjust = "holm")))
H3_cvc <- as.data.frame(confint(pairs(X2_cvc, adjust = "holm")))

B4_cvc <- as.data.frame(confint(pairs(X1_cvc, by = "Region", adjust = "holm")))
H4_cvc <- as.data.frame(confint(pairs(X2_cvc, by = "Region", adjust = "holm")))
B4_cvc
# Add phase labels
B1_cvc$PHASE <- "BL"
B2_cvc$PHASE <- "BL"
B3_cvc$PHASE <- "BL"
B4_cvc$PHASE <- "BL"
H1_cvc$PHASE <- "MAX"
H2_cvc$PHASE <- "MAX"
H3_cvc$PHASE <- "MAX"
H4_cvc$PHASE <- "MAX"

# Merge and process within-group comparisons
WITHIN1_cvc <- merge(B1_cvc, B3_cvc, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))
WITHIN2_cvc <- merge(H1_cvc, H3_cvc, by = c("contrast", "Group", "PHASE"), suffixes = c("", "_ci"))

WITH_cvc <- rbind(WITHIN1_cvc, WITHIN2_cvc)

# Recode factors for better visualization
WITH_cvc <- WITH_cvc %>%
  mutate(PHASE = recode_factor(PHASE,
                               BL = "Baseline",
                               MAX = "Maximum Vasodilation")) %>%
  mutate(Region = recode_factor(contrast,
                                "ABDOMEN - CALF" = "Abdomen - Calf",
                                "CHEST - ABDOMEN" = "Abdomen - Chest",
                                "FOREARM - ABDOMEN" = "Abdomen - Forearm",
                                "CHEST - CALF" = "Chest - Calf",
                                "FOREARM - CALF" = "Forearm - Calf",
                                "FOREARM - CHEST" = "Forearm - Chest")) %>%
  mutate(Group = recode_factor(Group,
                               HRT = "Hormone Replacement Therapy",
                               NON = "Control"))

# Process between-group comparisons
TWEEN1_cvc <- merge(B2_cvc, B4_cvc, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN2_cvc <- merge(H2_cvc, H4_cvc, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))

TWEEN_cvc <- rbind(TWEEN1_cvc, TWEEN2_cvc)

# Overall group comparisons for CVC
M_cvc <- emmeans(BLCVC, ~ Group)
M1_cvc <- emmeans(MAX, ~ Group)

E1_cvc <- as.data.frame(confint(pairs(M_cvc)))
E2_cvc <- as.data.frame(confint(pairs(M1_cvc)))

E_cvc <- rbind(E1_cvc, E2_cvc)
E_cvc$PHASE <- c("BL", "MAX")

# Add overall comparison to between-group dataframe
E_for_bind_cvc <- E_cvc %>%
  mutate(Region = "POOLED") %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL)

TWEEN_for_bind_cvc <- TWEEN_cvc %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL)

TWEEN_cvc <- rbind(TWEEN_for_bind_cvc, E_for_bind_cvc)

# Final dataframes for CVC models
CVC_WITHIN_COMPARISONS <- WITH_cvc
CVC_BETWEEN_COMPARISONS <- TWEEN_cvc

# View the results
head(CVC_WITHIN_COMPARISONS)
head(CVC_BETWEEN_COMPARISONS)

####Sensitivity analysis#####
DATA<- DATA%>%
  mutate(`Therapy Type` = case_when(
    ID %in% c("MENOMD-H03") ~ "Estrogen Only",
    ID %in% c("MENOMD-H01", "MENOMD-H05") ~ "Estrogen, Progesterone, Tesosterone",
    ID %in% c("MENOMD-H02", "MENOMD-H04", "MENOMD-H06",
              "MENOMD-H07", "MENOMD-H08", "MENOMD-H09",
              "MENOMD-H10") ~ "Estrogen, Progesterone",
    
    TRUE ~ "None"  # For any IDs not specified
  ))
DATAT<-DATA[which(DATA$`Therapy Type` != "Estrogen, Progesterone, Tesosterone"),]

BL_final <- lme(BL ~ Group * Region,
                random = ~ 1 | ID,
                weights = varIdent(form = ~ 1 | Site_Order),
                data = DATAT,
                na.action = na.exclude,
                method = "REML")

# Heating phase model adjusting for baseline
HEAT_final <- lme(HEAT ~  Group * Region+BL,
                  random = ~ 1 | ID,
                  weights = varIdent(form = ~ 1 | Site_Order),
                  data = DATAT,
                  na.action = na.exclude)

# L-NAME phase model
LNA_final <- lme(LNAME ~  Group * Region+BL,
                 random = ~ 1 | ID,
                 weights = varIdent(form = ~ 1 | Site_Order),
                 data = DATAT,
                 na.action = na.exclude)

# NOS phase model
NOS_final <- lme(NOS ~  Group * Region+BL,
                 random = ~ 1 | ID,
                 weights = varIdent(form = ~ 1 | Site_Order),
                 data = DATAT,
                 na.action = na.exclude)

summary(BL_final)

###RE-RUN MARGINAL MEANS SECTION
CVC<- CVC%>%
  mutate(`Therapy Type` = case_when(
    ID %in% c("MENOMD-H03") ~ "Estrogen Only",
    ID %in% c("MENOMD-H01", "MENOMD-H05") ~ "Estrogen, Progesterone, Tesosterone",
    ID %in% c("MENOMD-H02", "MENOMD-H04", "MENOMD-H06",
              "MENOMD-H07", "MENOMD-H08", "MENOMD-H09",
              "MENOMD-H10") ~ "Estrogen, Progesterone",
    
    TRUE ~ "None"  # For any IDs not specified
  ))
CVCT<-CVC[which(CVC$`Therapy Type` != "Estrogen, Progesterone, Tesosterone"),]


BLCVC <- lme(BL ~ Group * Region,
             random = ~ 1 | ID,
             weights = varIdent(form = ~ 1 | Site_Order),
             data = CVC,
             na.action = na.exclude)

MAX <- lme(MAX ~  Group * Region,
           random = ~ 1 | ID,
           weights = varIdent(form = ~ 1 | Site_Order),
           data = CVC,
           na.action = na.exclude)

summary(MAX)
summary(BLCVC)
# Calculate estimated marginal means for each phase
X1 <- emmeans(BLCVC, ~ Region * Group, vcov. = vcov(BLCVC))
X2 <- emmeans(MAX, ~ Region * Group, vcov. = vcov(MAX))


# Pairwise comparisons within groups
B1 <- as.data.frame(pairs(X1, by = "Group", adjust = "holm"))
H1 <- as.data.frame(pairs(X2, by = "Group", adjust = "holm"))
B2 <- as.data.frame(pairs(X1, by = "Region", adjust = "holm"))
H2 <- as.data.frame(pairs(X2, by = "Region", adjust = "holm"))
confint(pairs(X1, by = "Region", adjust = "holm"))

B1### DATA VISUALIZATION ####
B2
H1
H2
# Prepare data for plotting

DATA_E <- DATA_E %>%
  mutate(Phase = recode_factor(Phase,
                               BL = "Baseline",
                               HEAT = "Heating Plateau",
                               LNAME = "L-NAME Plateau",
                               NOS = "Nitric Oxide Contribution",
                               PEAK = "Initial Vasodilator Peak")) %>%
  mutate(Region = recode_factor(Region,
                                POOLED = "All Regions",
                                ABDOMEN = "Abdomen",
                                CALF = "Calf",
                                CHEST = "Chest",
                                FOREARM = "Forearm")) %>%
  mutate(Group = recode_factor(Group,
                               HRT = "Hormone Replacement Therapy",
                               NON = "Control"))

# Main boxplot visualization
PLOTA <- ggplot(data = DATA_E[which (DATA_E$Phase != "Initial Vasodilator Peak",)], aes(x = Region, y = CVCMax, fill = Group)) +
  scale_fill_manual(values = c("#882255","#117733")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
               width = 0.6, lwd = 0.7, color = "black") +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = expression("%CVC"["Max"]), fill = "Group") +
  #facet_wrap(~Phase, ncol = 4) +
  theme_few() +
  theme(legend.position = c(0.90, -0.44),
        legend.justification = c(0.5, 0),
        legend.direction = "vertical",
        legend.box.just = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.margin = margin(t = 10),
        plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0.05, 0.1)))

PLOTA 
DATA_M<-DATA_E



DATA_M<- DATA_M %>%
  mutate(`Therapy Type` = case_when(
    ID %in% c("MENOMD-H03") ~ "Estrogen Only",
    ID %in% c("MENOMD-H01", "MENOMD-H05") ~ "Estrogen, Progesterone, Tesosterone",
    ID %in% c("MENOMD-H02", "MENOMD-H04", "MENOMD-H06",
              "MENOMD-H07", "MENOMD-H08", "MENOMD-H09",
              "MENOMD-H10") ~ "Estrogen, Progesterone",
    
    TRUE ~ "None"  # For any IDs not specified
  ))


PLOTB <- ggplot(data = DATA_M, aes(x = Region, y = CVCMax, fill = `Therapy Type`)) +
  #geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
   #            width = 0.6, lwd = 0.4, color = "black") +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = expression("%CVC"["Max"]), fill = "Group") +
  facet_wrap(~Phase, ncol = 4) +
  theme_few() +
  theme(#legend.position = c(0.90, -0.44),
        #legend.justification = c(0.5, 0),
        #legend.direction = "vertical",
        #legend.box.just = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.margin = margin(t = 10),
        plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0.05, 0.1)))
PLOTB
# Prepare data for between-group comparison plot
TWEEN <- TWEEN %>%
  mutate(PHASE = recode_factor(PHASE,
                               BL = "Baseline",
                               HEAT = "Heating Plateau",
                               LNAME = "L-NAME Plateau",
                               NOS = "Nitric Oxide Contribution")) %>%
  mutate(Region = recode_factor(Region,
                                POOLED = "All Regions",
                                ABDOMEN = "Abdomen",
                                CALF = "Calf",
                                CHEST = "Chest",
                                FOREARM = "Forearm")) %>%
  mutate(color_group = ifelse(Region == "All Regions", "All Regions", "Other Regions"))


# Between-group comparison plot
PLOT_REGION <- ggplot(data = TWEEN, aes(x = estimate, y = Region, color = color_group)) +
  annotate(geom = "rect", xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", alpha = 0.7) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2, show.legend = FALSE) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_text(aes(x = 25, label = "*"), 
            hjust = 0, size = 6, vjust = 0.75, color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("#062A78", "#45B1E8")) +
  facet_wrap(~PHASE, ncol = 1, strip.position = "right") +
  scale_y_discrete(limits = rev(levels(TWEEN$Region))) +
  scale_x_continuous(limits = c(-40, 40), expand = c(0, 0)) +
  labs(x = expression("Mean %CVC"["Max"]*" Difference"), y = "") +
  theme_few() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        plot.margin = margin(t = 10, r = 10, b = 5, l = 5),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside")
PLOT_REGION
# Within-group comparison plot
PLOT_GROUP <- ggplot(data = WITH, aes(x = estimate, y = Region, colour = Group)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2, alpha = 0.5) +
  annotate(geom = "rect", xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "grey50") +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.1, show.legend = FALSE) +
  geom_point(aes(shape = Group), size = 2, show.legend = FALSE) +
  scale_colour_manual(values = c("#882255", "#117733")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = c("#882255", "#117733")) +
  scale_y_discrete(limits = rev(levels(WITH$Region))) +
  coord_cartesian(xlim = c(-40, 40)) +
  labs(x = expression("Mean %CVC"["Max"]*" Difference"), y = "") +
  facet_grid(PHASE ~ Group) +
  geom_text(aes(x = 30, label = "*"), 
            hjust = 0, size = 6, vjust = 0.75,color = "black", show.legend = FALSE) +
  theme_few() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 9),
        strip.text.y = element_text(angle = 270, size = 8),
        panel.spacing = unit(0.8, "lines"))
PLOT_GROUP

###CVC comparision
CVC_P <- CVC_P%>%
  mutate(Phase = recode_factor(Phase,
                               BL = "Baseline",
                               MAX = "Maximum Vasodilation")) %>%
  mutate(Region = recode_factor(Region,
                                POOLED = "All Regions",
                                ABDOMEN = "Abdomen",
                                CALF = "Calf",
                                CHEST = "Chest",
                                FOREARM = "Forearm")) %>%
  mutate(Group = recode_factor(Group,
                               HRT = "Hormone Replacement Therapy",
                               NON = "Control"))


PLOTC <- ggplot(data = CVC_P, aes(x = Region, y = CVCMax, fill =Group)) +
  scale_fill_manual(values = c("#882255","#117733")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.6,
               width = 0.6, lwd = 0.4, color = "black") +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = "CVC", fill = "Group") +
  facet_wrap(~Phase, ncol = 4) +
  theme_few() +
  theme(legend.position = c(0.90, -0.44),
        legend.justification = c(0.5, 0),
        legend.direction = "vertical",
        legend.box.just = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.margin = margin(t = 10),
        plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))

PLOTC <- ggplot(data = CVC_P, aes(x = Region, y = CVCMax, fill = Group)) +
  scale_fill_manual(values = c("#882255","#117733")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.6,
               width = 0.6, lwd = 0.4, color = "black") +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = "CVC (PU·mmHg⁻¹)", fill = "Group") +
  facet_wrap(~Phase, ncol = 4) +
  theme_few() +
  theme(legend.position = c(0.90, -0.44),
        legend.justification = c(0.5, 0),
        legend.direction = "vertical",
        legend.box.just = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.margin = margin(t = 10),
        plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill = guide_legend(reverse = TRUE)) #scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0.05, 0.1)))
PLOTC

### OUTPUT RESULTS ####

write.csv(WITH, "HRT_within_means.csv")
write.csv(TWEEN, "HRT_tween_means.csv")
write.csv(S_STAT, "HRT_S_characteristics.csv")
write.csv(STAT, "HRT_S_CVC.csv")
write.csv(t_test_df,"HRT_test.csv")

ggsave("BIGPLOT.png", PLOTA, dpi = 1000, height = 8, width = 11)
ggsave("HRTTWEENCVC.png", PLOT_REGION, dpi = 700, height = 7, width = 7)
ggsave("HRTWITHINCVC.png", PLOT_GROUP, dpi = 700, height = 7, width = 7)
ggsave("HRTCVC.png",PLOTC,dpi = 700, height = 6, width = 9)

###CVC COMPARISON FIGURE####

REP <- ggplot(data = rep, aes(x = TIME, y = CVCA, colour = GROUP)) +
  geom_point(size = 2) +
  scale_fill_manual(values = c("#882255", "#117733"),
                    labels = c("Hormone Replacement Therapy", "Control")) +
  scale_colour_manual(values = c("#882255", "#117733"),
                      labels = c("Hormone Replacement Therapy", "Control")) +
  # Change y-axis label and scale (25% increments = every 25 units)
  labs(y = expression("%CVC"["max"])) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 25)) +
  # Change x-axis label and add tick marks every 10 minutes
  labs(x = "Time (minutes)",
       colour = "group") +  # lowercase the legend title
  scale_x_continuous(breaks = seq(0, max(rep$TIME, na.rm = TRUE), 10)) +
  # Solid grey background with no grid or box, and bold axis labels
  theme(
    panel.background = element_rect(fill = "grey90"),
    plot.background = element_rect(fill = "grey90"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(face = "bold"),  # Bold x-axis label
    axis.title.y = element_text(face = "bold")   # Bold y-axis label
  )

REP
REP <- ggplot(data = rep, aes(x = TIME, y = CVCA, colour = GROUP)) +
  geom_point(size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#882255", "#117733"),
                    labels = c("Hormone Replacement Therapy", "Control")) +
  scale_colour_manual(values = c("#882255", "#117733"),
                      labels = c("Hormone Replacement Therapy", "Control")) +
  # Change y-axis label and scale (25% increments = every 25 units)
  labs(y = expression(bold("%CVC"["max"]))) +  # Bold the y-axis label
  scale_y_continuous(limits = c(0, 110), 
                     breaks = seq(0, 100, 25)) +
  # Change x-axis label and add tick marks every 10 minutes
  labs(x = "Time (minutes)",
       color = "Group") +
  scale_x_continuous(breaks = seq(0, 120, 20),
                     limits = c(0, 121),  # Anchor the axis to end at 120
                     expand = c(0, 0))+
  # Grey background only in the plot area, no grid, no box
  theme(
    panel.background = element_rect(fill = "grey95"),  # Only the plot area gets grey
    plot.background = element_rect(fill = "white"),    # Outside area stays white
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text = element_text( face = "bold"),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(face = "bold", size = 12),        # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 12)         # Bold y-axis label
  )

REP

REP2 <- ggplot(data = rep, aes(x = TIME, y = PUA, colour = GROUP)) +
  geom_point(size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#882255", "#117733"),
                    labels = c("Hormone Replacement Therapy", "Control")) +
  scale_colour_manual(values = c("#882255", "#117733"),
                      labels = c("Hormone Replacement Therapy", "Control")) +
  # Change y-axis label and scale (25% increments = every 25 units)
  labs(y = "Red blood cell flux (perfusion units)")+
  
  scale_y_continuous(limits = c(0, 250), 
                    breaks = seq(0, 250, 50),
                    expand = c(0,0)) +
  # Change x-axis label and add tick marks every 10 minutes
  labs(x = "Time (minutes)",
       color = "Group") +
  scale_x_continuous(breaks = seq(0, 120, 10),
                     limits = c(0, 121),  # Anchor the axis to end at 120
                     expand = c(0, 0))+
  # Grey background only in the plot area, no grid, no box
  theme(
    panel.background = element_rect(fill = "grey95"),  # Only the plot area gets grey
    plot.background = element_rect(fill = "white"),    # Outside area stays white
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(face = "bold"),        # Bold x-axis label
    axis.title.y = element_text(face = "bold")         # Bold y-axis label
  )
REP2

ggsave("REP1.png", REP, dpi = 700, height = 7, width = 10)
ggsave("REP2.png", REP2, dpi = 700, height = 7, width = 10)
