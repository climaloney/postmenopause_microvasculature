
####LIBRARY####

library(lmerTest) ###mixed model--prefer nlme
library(lmtest) ###LRT
library(lme4) ###mixed model
library(ggplot2) ###plots
library(ggdist) ###plots
library(tidyverse) ###DATa prep
library(car)
library(ggrain)
library(effects)
library(rptR)
library(nlme) ###
library(car)
library(emmeans)
library(parallel)
library(rptR) ##to compare results with nlme
library(boot)
library(effects)
library(ggthemes)
####CVCMAX PREP####

###LOAD CVCMAX SHEETS
CVCMAX<-read.csv("OF_REGIONAL_CVCMAX.csv",header = T)
CVC<-read.csv("OF_REGIONAL_CVC.csv", header = T)
SITES<-read.csv("OF_REGIONAL_SITES.csv", header = T)
rep<-read.csv("flashrep.csv", header = TRUE)


CVCMAX <- CVCMAX %>%
  filter(
    # Keep all rows where Study is "FLASH"
    Study == "FLASH")
CVC <- CVC %>%
  filter(
    Study == "FLASH")


CVC<-CVC[which(CVC$STAGE != "HEAT"),]
CVC<-CVC[which(CVC$STAGE != "NOS"),]
CVCMAX <- CVCMAX %>%
  mutate(GROUP = case_when(
    substr(ID, 2, 2) == "H" ~ "FLASH",
    substr(ID, 2, 2) == "N" ~ "NONE",
    
    TRUE ~ NA_character_  
  ))
CVC <- CVC %>%
  mutate(GROUP = case_when(
    substr(ID, 2, 2) == "H" ~ "FLASH",
    substr(ID, 2, 2) == "N" ~ "NONE",
    TRUE ~ NA_character_  
  ))
write.csv(CVC, "CVC_VMS.csv")
write.csv(CVCMAX,"CVCMAX_VMS.csv")
#change
CVCMAX$YearsMeno<-CVCMAX$Age-CVCMAX$MenoOnset
CVC$YearsMeno<-CVC$Age-CVC$MenoOnset
###CREATE POOLED MEASUREMENT FOR FIGURE (MOSTLY)

CVCMAX$POOLED<-rowMeans(CVCMAX[, 5:8], na.rm = TRUE)
CVC$POOLED<-rowMeans(CVC[, 5:8], na.rm = TRUE)
###RESHAPE CVCMAX FOR USABILITY FOR ANALYSIS

CVCMAX<- CVCMAX%>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF, POOLED), names_to = "Region", values_to = "CVCMax") %>%
  mutate(across(where(is.character), as.factor))

CVC<- CVC%>%
  pivot_longer(c(FOREARM, CHEST, ABDOMEN, CALF,POOLED), names_to = "Region", values_to = "CVC") %>%
  mutate(across(where(is.character), as.factor))

SITES <- SITES %>%
  pivot_longer(
    cols = c(FOREARM, CHEST, ABDOMEN, CALF,POOLED),
    names_to = "Region",
    values_to = "Site_Order"
  ) %>%
  mutate(Region = toupper(Region))

CVCMAX$BSA<-0.007184*CVCMAX$Height^0.725*CVCMAX$Weight^0.425
CVC$BSA<-0.007184*CVC$Height^0.725*CVC$Weight^0.425

###BIND SITES TO CVCMAX
CVCMAX <- CVCMAX %>%
  left_join(SITES, by = c("ID", "Region", "Study")) %>%
  mutate(across(where(is.character), as.factor))

###BIND SITES TO CVC
CVC<- CVC %>%
  left_join(SITES, by = c("ID", "Region", "Study")) %>%
  mutate(across(where(is.character), as.factor))

CVCMAX$Site_Order<-as.factor(CVCMAX$Site_Order)
CVC$Site_Order<-as.factor(CVC$Site_Order)
CVCMAX$Region <- relevel(CVCMAX$Region, ref = "FOREARM")
CVC$Region <- relevel(CVC$Region, ref = "FOREARM") 
##pivot wider

CVCMAX<-CVCMAX %>%
  pivot_wider(names_from = "STAGE", values_from = "CVCMax")
CVC<-CVC %>%
  pivot_wider(names_from = "STAGE", values_from = "CVC")
CVCMAX$LNAME<-as.numeric(CVCMAX$HEAT-CVCMAX$NOS)

##dataframes for plots
CVCPLOT<-CVCMAX %>%
  pivot_longer(c(BL,HEAT,LNAME,NOS), names_to = "STAGE", values_to = "CVCMAX") %>%
  mutate(across(where(is.character), as.factor))
CPLOT<-CVC %>%
  pivot_longer(c(BL,MAX), names_to = "STAGE", values_to = "CVC") %>%
  mutate(across(where(is.character), as.factor))

CVCMAX<-CVCMAX[which(CVCMAX$Region != "POOLED"),]
CVC<-CVC[which(CVC$Region != "POOLED"),]

CVCMAX$GROUP<-relevel(CVCMAX$GROUP, ref = "NONE")
CVCMAX$Region<-relevel(CVCMAX$Region, ref = "FOREARM")

##test something
ONLY<-CVCMAX[which(CVCMAX$Region == "FOREARM" & CVCMAX$GROUP == "FLASH"),]
ONLY2<-CVCMAX[which(CVCMAX$Region == "CHEST" & CVCMAX$GROUP == "NONE"),]
t.test(ONLY$BLHR,ONLY$ENDHR, paired = T)
t.test(ONLY$BLST,ONLY$ENDST,paired=T)
t.test(ONLY2$BLHR,ONLY2$ENDHR, paired = T)
t.test(ONLY2$BLST,ONLY2$ENDST,paired=T)
t.test(ONLY$BLST,ONLY2$BLST)
mean(ONLY$BLST)
ONLY2<-ONLY2[!is.na(ONLY2$BLST),]

ONLY<-CVCMAX[which(CVCMAX$GROUP == "FLASH"),]
mean(ONLY$HEAT[!is.na(ONLY$HEAT)])
sd(ONLY$HEAT[!is.na(ONLY$HEAT)])
ONLY2<-CVCMAX[which(CVCMAX$GROUP == "NONE"),]
mean(ONLY2$HEAT[!is.na(ONLY2$HEAT)])
sd(ONLY2$HEAT[!is.na(ONLY2$HEAT)])
31.97-30.35
mean(ONLY2$BLST)
sd (ONLY2$BLST)

cvc_summary <- CPLOT%>%
  group_by(Region, STAGE,GROUP) %>%
  summarise(
    n_observations = sum(!is.na(CVC)),
    mean_CVC = mean(CVC, na.rm = TRUE),
    sd_CVC = sd(CVC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Region, STAGE)


cvcmax_summary <- CVCPLOT%>%
  group_by(Region, STAGE,GROUP) %>%
  summarise(
    n_observations = sum(!is.na(CVCMAX)),
    mean_CVC = mean(CVCMAX, na.rm = TRUE),
    sd_CVC = sd(CVCMAX, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Region, STAGE)

###LINEAR MIXED MODELS####

#####CVCMAX####

BL<- lme(BL ~ GROUP * Region+Site_Order+SEASON,
                random = ~ 1 | ID,
                data = CVCMAX,
                na.action = na.exclude,
                method = "REML")
leveneTest(residuals(BL) ~ CVCMAX$Region)
leveneTest(residuals(BL) ~ CVCMAX$Site_Order)
leveneTest(residuals(BL) ~ CVCMAX$GROUP)
leveneTest(residuals(BL) ~ CVCMAX$SEASON)

BL1<-lme(BL ~ GROUP * Region,
         random = ~ 1 | ID,
         data = CVCMAX,
         na.action = na.exclude,
         method = "REML")
BL2<-lme(BL ~ GROUP * Region,
         random = ~ 1 | ID,
         #weights = varIdent(form = ~ 1 | Site_Order),
         data = CVCMAX,
         na.action = na.exclude,
         method = "REML")

lrtest(BL1,BL2)
summary(BL2)
anova(BL2)
# Heating phase model adjusting for baseline
HEAT <- lme(HEAT ~ BL + GROUP * Region+Site_Order+SEASON,
                  random = ~ 1 | ID,
                  data = CVCMAX,
                  na.action = na.exclude)

leveneTest(residuals(HEAT) ~ CVCMAX$Region)
leveneTest(residuals(HEAT) ~ CVCMAX$Site_Order)
leveneTest(residuals(HEAT) ~ CVCMAX$GROUP)
leveneTest(residuals(HEAT) ~ CVCMAX$SEASON)

HEAT1 <- lme(HEAT ~ BL + GROUP * Region,
            random = ~ 1 | ID,
            data = CVCMAX,
            na.action = na.exclude)
HEAT2<-lme(HEAT ~ BL + GROUP * Region+Site_Order+SEASON,
           random = ~ 1 | ID,
           #weights = varIdent(form = ~ 1 | SEASON),
           data = CVCMAX,
           na.action = na.exclude)



lrtest(HEAT1,HEAT2)
summary(HEAT1)
summary(HEAT2)
anova(HEAT2)
# L-NAME phase model

LNA<-lme(LNAME ~ BL + GROUP * Region+Site_Order+SEASON,
         random = ~ 1 | ID,
         data = CVCMAX,
         na.action = na.exclude)
leveneTest(residuals(LNA) ~ CVCMAX$Region)
leveneTest(residuals(LNA) ~ CVCMAX$Site_Order)
leveneTest(residuals(LNA) ~ CVCMAX$GROUP) ###hetrogenity
leveneTest(residuals(LNA) ~ CVCMAX$SEASON)

LNA1 <- lme(LNAME ~ BL + GROUP * Region,
                 random = ~ 1 | ID,
                 weights = varIdent(form = ~ 1 | GROUP),
                 data = CVCMAX,
                 na.action = na.exclude)
LNA2 <- lme(LNAME ~ BL + GROUP * Region+Site_Order+SEASON,
            random = ~ 1 | ID,
            weights = varIdent(form = ~ 1 | GROUP),
            data = CVCMAX,
            na.action = na.exclude)

lrtest(LNA1,LNA2) ###last one again 
summary(LNA2)
anova(LNA2)

# NOS phase model
NOS<- lme(NOS ~ BL + GROUP * Region+Site_Order+SEASON,
                 random = ~ 1 | ID,
                 data = CVCMAX,
                 na.action = na.exclude)
leveneTest(residuals(NOS) ~ CVCMAX$Region)
leveneTest(residuals(NOS) ~ CVCMAX$Site_Order)
leveneTest(residuals(NOS) ~ CVCMAX$GROUP)
leveneTest(residuals(NOS) ~ CVCMAX$SEASON)

NOS1<- lme(NOS ~ BL + GROUP * Region,
           random = ~ 1 | ID,
           data = CVCMAX,
           na.action = na.exclude)
NOS2<-lme(NOS ~ BL + GROUP * Region+Site_Order+SEASON,
          random = ~ 1 | ID,
          data = CVCMAX,
          na.action = na.exclude)

lrtest(NOS1,NOS2)
summary(NOS2)
anova(NOS2)
####assumption test

plot(BL2)
plot(HEAT2)
plot(LNA2)
plot(NOS2)




######MARGINAL MEANS AND POST-HOC COMPARISONS #####

# Calculate estimated marginal means for each phase
X1 <- emmeans(BL2, ~ Region * GROUP, vcov. = vcov(BL2))
X2 <- emmeans(HEAT2, ~ Region * GROUP, vcov. = vcov(HEAT2))
X3 <- emmeans(NOS2, ~ Region * GROUP, vcov. = vcov(NOS2))
X4 <- emmeans(LNA2, ~ Region * GROUP, vcov. = vcov(NOS2))

# Pairwise comparisons within groups
B1 <- as.data.frame(pairs(X1, by = "GROUP", adjust = "holm"))
H1 <- as.data.frame(pairs(X2, by = "GROUP", adjust = "holm"))
N1 <- as.data.frame(pairs(X3, by = "GROUP", adjust = "holm"))
L1 <- as.data.frame(pairs(X4, by = "GROUP", adjust = "holm"))

# Pairwise comparisons within regions
B2 <- as.data.frame(pairs(X1, by = "Region", adjust = "holm"))
H2 <- as.data.frame(pairs(X2, by = "Region", adjust = "holm"))
N2 <- as.data.frame(pairs(X3, by = "Region", adjust = "holm"))
L2 <- as.data.frame(pairs(X4, by = "Region", adjust = "holm"))

# Confidence intervals for comparisons
B3 <- as.data.frame(confint(pairs(X1, by = "GROUP", adjust = "holm")))
H3 <- as.data.frame(confint(pairs(X2, by = "GROUP", adjust = "holm")))
N3 <- as.data.frame(confint(pairs(X3, by = "GROUP", adjust = "holm")))
L3 <- as.data.frame(confint(pairs(X4, by = "GROUP", adjust = "holm")))

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
WITHIN1 <- merge(B1, B3, by = c("contrast", "GROUP", "PHASE"), suffixes = c("", "_ci"))
WITHIN2 <- merge(H1, H3, by = c("contrast", "GROUP", "PHASE"), suffixes = c("", "_ci"))
WITHIN3 <- merge(L1, L3, by = c("contrast", "GROUP", "PHASE"), suffixes = c("", "_ci"))
WITHIN4 <- merge(N1, N3, by = c("contrast", "GROUP", "PHASE"), suffixes = c("", "_ci"))

WITH <- rbind(WITHIN1, WITHIN2, WITHIN3, WITHIN4)

# Process between-group comparisons
TWEEN1 <- merge(B2, B4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN2 <- merge(H2, H4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN3 <- merge(L2, L4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))
TWEEN4 <- merge(N2, N4, by = c("contrast", "Region", "PHASE"), suffixes = c("", "_ci"))

TWEEN <- rbind(TWEEN1, TWEEN2, TWEEN3, TWEEN4)

###Overall group comparisons
M <- emmeans(BL2, ~ GROUP)
M1 <- emmeans(HEAT2, ~ GROUP)
M2 <- emmeans(LNA2, ~ GROUP)
M3 <- emmeans(NOS2, ~ GROUP)

# Overall region comparisons
N <- emmeans(BL2, ~ Region)
N1 <- emmeans(HEAT2, ~ Region)
N2 <- emmeans(LNA2, ~ Region)
N3 <- emmeans(NOS2, ~ Region)

E1_pairs <- pairs(M)
E2_pairs <- pairs(M1)
E3_pairs <- pairs(M2)
E4_pairs <- pairs(M3)


F1_pairs <- pairs(N)
F2_pairs <- pairs(N1)
F3_pairs <- pairs(N2)
F4_pairs <- pairs(N3)


E1 <- as.data.frame(confint(E1_pairs))
E2 <- as.data.frame(confint(E2_pairs))
E3 <- as.data.frame(confint(E3_pairs))
E4 <- as.data.frame(confint(E4_pairs))


E1$p.value <- summary(E1_pairs)$p.value
E2$p.value <- summary(E2_pairs)$p.value
E3$p.value <- summary(E3_pairs)$p.value
E4$p.value <- summary(E4_pairs)$p.value


F1 <- as.data.frame(confint(F1_pairs))
F2 <- as.data.frame(confint(F2_pairs))
F3 <- as.data.frame(confint(F3_pairs))
F4 <- as.data.frame(confint(F4_pairs))

F1$p.value <- summary(F1_pairs)$p.value
F2$p.value <- summary(F2_pairs)$p.value
F3$p.value <- summary(F3_pairs)$p.value
F4$p.value <- summary(F4_pairs)$p.value

E <- rbind(E1, E2, E3, E4)
E$PHASE <- rep(c("BL", "HEAT", "LNAME", "NOS"), each = nrow(E1))
E_for_bind <- E %>%
  mutate(Region = "POOLED") %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL, p.value)
TWEEN_for_bind <- TWEEN %>%
  select(contrast, Region, PHASE, estimate, lower.CL, upper.CL, p.value)

TWEEN <- rbind(TWEEN_for_bind, E_for_bind)
E$COMPARISON_TYPE <- "Overall Group Comparison"

Fx <- rbind(F1, F2, F3, F4)
Fx$PHASE <- rep(c("BL", "HEAT", "LNAME", "NOS"), each = nrow(F1))
Fx$COMPARISON_TYPE <- "Overall Region Comparison"

overall_comparisons <- rbind(
  E %>% rename(Comparison = contrast),
  Fx %>% rename(Comparison = contrast)
)


#####CVC####

BLCVC <- lme(BL ~ GROUP * Region+Site_Order,
             random = ~ 1 | ID,
             data = CVC,
             na.action = na.exclude)

MAX <- lme(MAX ~ BL + GROUP * Region+Site_Order,
           random = ~ 1 | ID,
           data = CVC,
           weights = varIdent(form = ~ 1 | Site_Order),
           na.action = na.exclude)

summary(MAX)
summary(BLCVC)



# Calculate estimated marginal means for CVC phases
X1 <- emmeans(BLCVC, ~ Region * GROUP, vcov. = vcov(BLCVC))
X2 <- emmeans(MAX, ~ Region * GROUP, vcov. = vcov(MAX))
X1
B1 <- confint(pairs(X1, by = "GROUP", adjust = "holm"))
H1 <- as.data.frame(pairs(X2, by = "GROUP", adjust = "holm"))

B2 <- confint(pairs(X1, by = "Region", adjust = "holm"))
H2 <- confint(pairs(X2, by = "GROUP", adjust = "holm"))

pairs(emmeans(BLCVC, ~GROUP*Region))
pairs(emmeans(BLCVC,~Region))
pairs(emmeans(MAX, ~Region))

BX<-pairs(emmeans(BLCVC,~Region*GROUP))
HX<-pairs(emmeans(MAX,~Region))
HX1<-pairs(emmeans(MAX,~GROUP))
BX1<-pairs(emmeans(BLCVC,~GROUP))
HX
HX1
BX
BX1
####CVCMAX VISUALIZATION####

CVCPLOT <- CVCPLOT %>%
  mutate(STAGE = recode_factor(STAGE,
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
  mutate(GROUP = recode_factor(GROUP,
                                FLASH = "Vasomotor Symptoms",
                                NONE = "No Vasomotor Symptoms"))
  
PLOTA <- ggplot(data = CVCPLOT, aes(x = Region, y = CVCMAX, fill = GROUP)) +
  scale_fill_manual(values = c("#CC6677","#332288")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
               width = 0.6, lwd = 0.4, color = "black") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0.05, 0.1)))+
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = expression("%CVC"["Max"]), fill = "") +
  facet_wrap(~STAGE, ncol = 4) +
  theme_few() +
  theme(legend.position = c(0.8, -0.3),
        legend.justification = c(0.4, 0),
        legend.direction = "horizontal",
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
  guides(fill = guide_legend(reverse = TRUE)) 

PLOTA 

CPLOT <- CPLOT %>%
  mutate(STAGE = recode_factor(STAGE,
                               BL = "Baseline",
                               MAX = "Maximum Vasodilation")) %>%
  mutate(Region = recode_factor(Region,
                                POOLED = "All Regions",
                                ABDOMEN = "Abdomen",
                                CALF = "Calf",
                                CHEST = "Chest",
                                FOREARM = "Forearm")) %>%
  mutate(GROUP = recode_factor(GROUP,
                               FLASH = "Vasomotor Symptoms",
                               NONE = "No Vasomotor Symptoms"))
PLOTB<- ggplot(data = CPLOT, aes(x = Region, y = CVC, fill = GROUP)) +
  scale_fill_manual(values = c("#CC6677","#332288")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
               width = 0.6, lwd = 0.4, color = "black") +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.7, color = "black", stroke = 0.3) +
  labs(x = "Region", y = expression(CVC~(PU%.%mmHg^{-1})), fill = "") +
  facet_wrap(~STAGE, ncol = 4) +
  theme_few() +
  theme(legend.position = c(0.85, -0.3),
        legend.justification = c(0.5, 0),
        legend.direction = "horizontal",
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
  guides(fill = guide_legend(reverse = TRUE)) 
PLOTB

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
  scale_color_manual(values = c("#542788", "#9972B8")) +
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

WITH<-WITH %>%
  mutate(PHASE = recode_factor(PHASE,BL = "Baseline",
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
  mutate(GROUP = recode_factor(GROUP, FLASH = "Vasomotor Symptoms",
                               NONE = "No Vasomotor Symptoms"))

PLOT_GROUP <- ggplot(data = WITH, aes(x = estimate, y = Region, colour = GROUP)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2, alpha = 0.5) +
  annotate(geom = "rect", xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "grey50") +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.1, show.legend = FALSE) +
  geom_point(aes(shape = GROUP), size = 2, show.legend = FALSE) +
  scale_colour_manual(values = c("#CC6677","#332288")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = c("#CC6677","#332288")) +
  scale_y_discrete(limits = rev(levels(WITH$Region))) +
  coord_cartesian(xlim = c(-40, 40)) +
  labs(x = expression("Mean %CVC"["Max"]*" Difference"), y = "") +
  facet_grid(PHASE ~ GROUP) +
  theme_few() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 9),
    strip.text.y = element_text(angle = 270, size = 8),
    panel.spacing = unit(0.8, "lines")
  )

PLOT_GROUP




REP <- ggplot(data = rep, aes(x = TIME, y = CVC)) +
  geom_point(size = 1, shape = 21, colour = "#CC6677", fill = "#CC6677")+
  labs(
    y = expression(bold("%CVC"["max"])),
    x = "Time (minutes)"
  ) +
  coord_cartesian(
    ylim = c(0, 110),
    xlim = c(0, 135),
    expand = FALSE  
  ) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, max(rep$TIME, na.rm = TRUE), 10)) +
  theme(
    panel.background = element_rect(fill = "grey95"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

REP

ggsave("TWEENCVC.png", PLOT_REGION, dpi = 700, height = 7, width = 7)
ggsave("WITHINCVC.png", PLOT_GROUP, dpi = 700, height = 7, width = 7)
ggsave("CVCMAX.png",PLOTA,dpi = 700, height = 6, width = 9)
ggsave("CVC.png",PLOTB,dpi = 700, height = 6, width = 9)
ggsave("REPFLASH.png", REP, dpi = 700, height = 7, width = 10)
