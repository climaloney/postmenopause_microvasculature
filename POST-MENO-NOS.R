
####LIBRARY####

library(lmerTest) ###mixed model--prefer nlme
library(lmtest) ###LRT
library(lme4) ###mixed model
library(ggplot2) ###plots
library(tidyverse) ###DATa prep
library(car)
library(effects)
library(rptR)
library(nlme) ###
library(emmeans)
library(boot)
library(ggthemes)
library(MuMIn)
library(broom.mixed)
library(patchwork)
library(ggeffects)
####DATA PREP####

###LOAD DATA SHEETS
CVCMAX<-read.csv("OF_REGIONAL_CVCMAX.csv",header = T)

CVC<-read.csv("OF_REGIONAL_CVC.csv", header = T)
SITES<-read.csv("OF_REGIONAL_SITES.csv", header = T)

#change
CVCMAX$YearsMeno<-CVCMAX$Age-CVCMAX$MenoOnset
CVC$YearsMeno<-CVC$Age-CVC$MenoOnset
CVCMAX$BMI <- CVCMAX$Weight / (CVCMAX$Height / 100)^2
###CREATE POOLED MEASUREMENT FOR FIGURE (MOSTLY)

CVCMAX$POOLED<-rowMeans(CVCMAX[, 5:8], na.rm = TRUE)
CVC$POOLED<-rowMeans(CVC[, 5:8], na.rm = TRUE)
###RESHAPE DATA FOR USABILITY FOR ANALYSIS

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



###BIND SITES TO CVCMAX
CVCMAX <- CVCMAX %>%
  left_join(SITES, by = c("ID", "Region", "Study")) %>%
  mutate(across(where(is.character), as.factor))

###BIND SITES TO CVC
CVC<- CVC %>%
  left_join(SITES, by = c("ID", "Region", "Study")) %>%
  mutate(across(where(is.character), as.factor))

###factor and relevel everything

CVCMAX$Site_Order<-as.factor(CVCMAX$Site_Order)
CVC$Site_Order<-as.factor(CVC$Site_Order)
CVCMAX$Region <- relevel(CVCMAX$Region, ref = "FOREARM")
CVC$Region <- relevel(CVC$Region, ref = "FOREARM") 
CVCMAX$Chronic<-relevel(CVCMAX$Chronic, ref = "N")
CVC$Chronic<-relevel(CVC$Chronic, ref = "N")
CVCMAX$SEASON <-relevel(CVCMAX$SEASON, ref = "SUMMER")

##pivot wider

CVCMAX <- CVCMAX %>%
  pivot_wider(names_from = STAGE,
              values_from = c(CVCMax, ST, HR, DBP2)) %>%
  rename(BL   = CVCMax_BL,
         HEAT = CVCMax_HEAT,
         NOS  = CVCMax_NOS)
CVC<-CVC %>%
  pivot_wider(names_from = "STAGE", values_from = "CVC")
CVCMAX$LNAME<-as.numeric(CVCMAX$HEAT-CVCMAX$NOS)
CVC$LNAME<-as.numeric(CVC$HEAT-CVC$NOS)

##dataframes for plots
CVCPLOT<-CVCMAX %>%
  pivot_longer(c(BL,HEAT,LNAME,NOS), names_to = "STAGE", values_to = "CVCMAX") %>%
  mutate(across(where(is.character), as.factor))
CPLOT<-CVC %>%
  pivot_longer(c(BL,HEAT,LNAME,NOS,MAX), names_to = "STAGE", values_to = "CVC") %>%
  mutate(across(where(is.character), as.factor))

##All-in-one model vs. separate for each 
CVCMAX<-CVCMAX[which(CVCMAX$Region != "POOLED"),]
CVC<-CVC[which(CVC$Region != "POOLED"),]
CVCMAX1 <-CVCMAX
CVC1<-CVC
CVCMAX<-CVCMAX %>%
  pivot_longer(c(HEAT,LNAME,NOS), names_to = "STAGE", values_to = "CVCMAX") %>%
  mutate(across(where(is.character), as.factor))
CVC<-CVC %>%
  pivot_longer(c(HEAT,LNAME,NOS,MAX), names_to = "STAGE", values_to = "CVC") %>%
  mutate(across(where(is.character), as.factor))

###PARTICIPANT CHARACTERISTICS####
##total number of participants
length(unique(CVC$ID))

##count for each observation
unique_id_counts <- CVC%>%
  group_by(ID) %>%
  summarise(
    has_Weight = any(!is.na(Weight)),
    has_Height = any(!is.na(Height)),
    has_Age = any(!is.na(Age)),
    has_BF = any(!is.na(BF)),
    has_CVC = any(!is.na(CVC)),
    has_HRT = any(!is.na(HRT)),
    has_Meno = any(!is.na(MenoOnset)),
    has_HYPE = any(!is.na(HYPE)),
    has_T2D = any(!is.na(T2D)),
    has_CH = any(!is.na(Chronic)),
    .groups = 'drop'
  ) %>%
  summarise(across(starts_with("has_"), ~sum(.))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_unique_ids"
  ) %>%
  mutate(variable = gsub("has_", "", variable))

##physical characteristics 

anthro_continous <- CVCMAX %>%
  distinct(ID, Weight, Height, Age,YearsMeno, BF, MAP, VO2Rel,BLHR,BLST,BMI,SBP, DBP) %>%
  pivot_longer(
    cols = c(Weight, Height, Age,YearsMeno, BF, MAP, VO2Rel, BLHR,BLST,BMI,SBP, DBP),
    names_to ="variable",
    values_to = "value"
  )%>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

##chronic conditions/med use

anthropo_categorical <- CVC%>%
  distinct(ID, HRT, HYPE, T2D,MOOD,Chronic) %>%
  pivot_longer(
    cols = c(HRT, HYPE, T2D,MOOD,Chronic),
    names_to = "variable",
    values_to = "value"
  ) %>%
  count(variable, value) %>%
  group_by(variable) %>%
  mutate(
    percent = n / sum(n) * 100,
    display = paste0(n, " (", round(percent, 1), "%)")
  ) %>%
  pivot_wider(
    id_cols = variable,
    names_from = value,
    values_from = display,
    values_fill = "0 (0%)"
  )

test<-CVCMAX1[which(CVCMAX1$Region == "CHEST"),]

###check obesity vs high bf

pp <- CVCMAX1[!duplicated(CVCMAX1$ID), ]

sum(pp$BMI < 30 & pp$BF >= 35, na.rm = TRUE)
round(100 * summary(lm(BF ~ BMI, data = pp))$r.squared, 1)

###Test changes in HR and skin temperature
t.test (test$BLHR,test$ENDHR)
t.test(test$BLST,test$ENDST)

###other quick test:

pp <- CVCMAX1[!duplicated(CVCMAX1$ID), ]   # one row per participant

cor.test(pp$BLST, pp$ENDST) 

pp$dST <- pp$ENDST - pp$BLST
cor.test(pp$BLST, pp$dST)

###CVC CHARACTERISTICS####
plot(CVCMAX$CVCMAX~CVCMAX$STAGE)
cvc_summary <- CPLOT%>%
  group_by(Region, STAGE) %>%
  summarise(
    n_observations = sum(!is.na(CVC)),
    mean_CVC = mean(CVC, na.rm = TRUE),
    sd_CVC = sd(CVC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Region, STAGE)

cvcmax_summary <- CVCPLOT %>%
  group_by(Region, STAGE) %>%
  summarise(
    n_observations = sum(!is.na(CVCMAX)),
    mean_CVC = mean(CVCMAX, na.rm = TRUE),
    sd_CVC = sd(CVCMAX, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Region, STAGE)



###DATA INSPECTION####

hist(CVCMAX1$BL)
hist(CVCMAX1$HEAT)
hist(CVCMAX1$NOS)##nos and lname perfectly distributed

hist(CVC1$BL)
hist(CVC1$HEAT)
hist(CVC1$LNAME)
hist(CVC1$NOS)##nos and lname perfectly distributed
hist(CVC1$MAX)



###LMM ANALYSIS####

#####CVCMAX####
##INITIAL MODEL
###TRY ONE MODEL FOR EVERYTHING:
####including each enviro/measurement factor
###have to analyze BL separately

full_model1 <- lme(CVCMAX ~ Region * STAGE+BL,
                   random = ~ 1|ID,  
                   #weights = varIdent(form = ~ 1 | Site_Order),
                   data = CVCMAX,
                   na.action = na.exclude,
                   method = "REML")

model_residuals <- residuals(full_model1)
clean_residuals <- model_residuals[!is.na(model_residuals)]
clean_data <- CVCMAX[!is.na(CVCMAX$CVCMAX), ]
clean_stage <- clean_data$STAGE
leveneTest(clean_residuals ~ clean_data$STAGE)
leveneTest(clean_residuals ~ clean_data$Site_Order)
leveneTest(clean_residuals ~ clean_data$SEASON) ###season

###significant heterodiacity with STAGE & season, not with other covariates. 
full_model2 <- lme(CVCMAX ~ Region * STAGE+SEASON+BL,
                   random = ~ 1 | ID,  
                   weights = varIdent(form = ~ 1 | STAGE),
                   data = CVCMAX,
                   na.action = na.exclude,
                   method = "REML")
full_model3 <- lme(CVCMAX ~ Region * STAGE+Site_Order+BL,
                   random = ~ 1 | ID,  
                   weights = varIdent(form = ~ 1 | STAGE),
                   data = CVCMAX,
                   na.action = na.exclude,
                   method = "REML")
full_model4 <- lme(CVCMAX ~ Region * STAGE+SEASON+Site_Order+BL,
                   random = ~ 1 | ID,  
                   data = CVCMAX,
                   weights = varIdent(form = ~ 1 | STAGE),
                   na.action = na.exclude,
                   method = "REML")

full_model5<-lme(CVCMAX ~ Region * STAGE+SEASON+Site_Order+BL,
                 random = ~ 1 | ID,
                 data = CVCMAX,
                 weights = varIdent(form = ~ 1 | SEASON),
                 na.action = na.exclude,
                 method = "REML")

lrtest(full_model1,full_model2,full_model3,full_model4,full_model5)

BL_CVCMAX<-lme(BL ~ Region+SEASON+Site_Order,
               random = ~ 1 | ID, 
               data = CVCMAX,
               na.action = na.exclude,
               method = "REML")
anova(BL_CVCMAX)
summary(BL_CVCMAX) ##all sites significant effect, disapears when modelling across stages, same with season
summary(full_model4)
emm_BL<-emmeans(BL_CVCMAX, ~ Region)
confint(pairs(emm_BL))
# Pairwise comparisons of regions within each stage
region_comparisons <- pairs(emm_region_stage)
bl_comparisons<-pairs(emm_BL)
EM_CVCMAX <- cbind(
  as.data.frame(region_comparisons),
  as.data.frame(confint(region_comparisons))[, c("lower.CL", "upper.CL")]
)
EM_BL_CVCMAX <- cbind(
  as.data.frame(bl_comparisons),
  as.data.frame(confint(bl_comparisons))[, c("lower.CL", "upper.CL")]
)


###now each model individually
####Heating to 42degrees
HEAT<- lme(HEAT ~  BL+Region+SEASON+Site_Order,
           random = ~ 1 | ID,
           data = CVCMAX1,
           na.action = na.exclude,
           method = "REML")
anova(HEAT)
summary(HEAT)
emm_HEAT<-emmeans(HEAT, ~ Region)
pairs(emm_HEAT)

###NOS contribution
NOS<-lme(NOS ~BL+Region+Site_Order+SEASON,
         random = ~ 1 | ID,
         data = CVCMAX1,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(NOS)
anova(NOS)
emm_NOS<-emmeans(NOS, ~ Region)
confint(pairs(emm_NOS))
pairs(emm_NOS)

LNAME<-lme(LNAME ~BL+Region+Site_Order+SEASON,
           random = ~ 1 | ID,
           data = CVCMAX1,
           na.action = na.exclude,
           method = "REML")
summary(LNAME)
anova(LNAME)
emm_LNAME<-emmeans(LNAME, ~ Region)
confint(pairs(emm_LNAME))
pairs(emm_LNAME)
table(CVCMAX1$Site_Order, CVCMAX1$Region)


######sensitivity analysis####
SENSE <- CVCMAX[-which(CVCMAX$T2D == "Y" | 
                         CVCMAX$Chronic != "N" | 
                         CVCMAX$HYPE == "Y" | 
                         CVCMAX$MOOD == "Y" | 
                         CVCMAX$HRT == "Y"), ]
SENSEMODEL<-lme(CVCMAX ~ Region * STAGE+SEASON+Site_Order+BL,
                random = ~ 1 | ID, 
                data = SENSE,
                na.action = na.exclude,
                weights = varIdent(form = ~ 1 | STAGE),
                method = "REML")


anova(SENSEMODEL)
summary(SENSEMODEL)

emm_SENSE <- emmeans(SENSEMODEL, ~ Region | STAGE)
# Pairwise comparisons of regions within each stage
pairs(emm_SENSE) 

###no difference.



SENSE_BL<-lme(BL ~ Region+SEASON+Site_Order,
              random = ~ 1 | ID, 
              data = SENSE,
              na.action = na.exclude,
              method = "REML")
anova(SENSE_BL)
summary(SENSE_BL)
confint(pairs(emmeans(SENSE_BL,~Region)))
pairs(emmeans(SENSE_BL,~Region))


#####CVC####
full_modelCVC <- lme(CVC ~ Region * STAGE+SEASON+Site_Order+BL,
                     random = ~ 1 | ID,  
                     data = CVC[which(CVC$STAGE != "NOS"),],
                     na.action = na.exclude,
                     method = "REML")


leveneTest(residuals(full_modelCVC) ~ full_modelCVC$data$STAGE)
leveneTest(residuals(full_modelCVC) ~ full_modelCVC$data$Site_Order)
leveneTest(residuals(full_modelCVC) ~ full_modelCVC$data$SEASON) ##SEASON
leveneTest(residuals(full_modelCVC) ~ full_modelCVC$data$Region)

CVC<-CVC[which(CVC$STAGE != "NOS"),]
###season has a lot of variability, region too but the season is more important 
full_modelCVC1 <- lme(CVC ~ Region * STAGE+SEASON+Site_Order+BL,
                      random = ~ 1 | ID,  
                      data = CVC,
                      na.action = na.exclude,
                      weights = varIdent(form = ~ 1 | SEASON),
                      method = "REML")
full_modelCVC2 <- lme(CVC ~ Region * STAGE+SEASON+Site_Order+BL,
                      random = ~ 1 | ID,  
                      data = CVC,
                      na.action = na.exclude,
                      weights = varIdent(form = ~ 1 | STAGE),
                      method = "REML")
##2 IS SAME AS CVCMAX MODEL


lrtest(full_modelCVC,full_modelCVC1,full_modelCVC2) ###stage right one  
summary(full_modelCVC2)
anova(full_modelCVC2)
BL_CVC<- lme(BL ~  Region+Site_Order+SEASON,
             random = ~ 1 | ID,
             data = CVC1,
             na.action = na.exclude,
             weights = varIdent(form = ~ 1 | SEASON),
             method = "REML")

summary(BL_CVC)
anova(BL_CVC)
emm_region_CVC <- emmeans(full_modelCVC2, ~ Region | STAGE)
emm_BL_CVC<-emmeans(BL_CVC, ~ Region)
region_CVC<- pairs(emm_region_CVC)
bl_CVC<-pairs(emm_BL_CVC)
EM_CVC <- cbind(
  as.data.frame(region_CVC),
  as.data.frame(confint(region_CVC))[, c("lower.CL", "upper.CL")]
)
EM_BL_CVC <- cbind(
  as.data.frame(bl_CVC),
  as.data.frame(confint(bl_CVC))[, c("lower.CL", "upper.CL")]
)

x <-lme(CVC~Region*STAGE+SEASON+Site_Order,
        random = ~1 | ID,
        weights = varIdent(form = ~ 1 | STAGE),
        data =  CPLOT[which(CPLOT$Site_Order != 5),],
        na.action = na.exclude)
confint(pairs(emmeans(x, ~STAGE)))
######Sensitivity analysis####

SENSECVC <- CVC[-which(CVC$T2D == "Y" | 
                         CVC$Chronic != "N" | 
                         CVC$HYPE == "Y" | 
                         CVC$MOOD == "Y" | 
                         CVC$HRT == "Y"), ]
SENSE <- CVCMAX1[-which(CVCMAX1$T2D == "Y" | 
                          CVCMAX1$Chronic != "N" | 
                          CVCMAX1$HYPE == "Y" | 
                          CVCMAX1$MOOD == "Y" | 
                          CVCMAX1$HRT == "Y"), ]

unique_id_counts <- SENSECVC%>%
  group_by(ID) %>%
  summarise(
    has_Weight = any(!is.na(Weight)),
    has_Height = any(!is.na(Height)),
    has_Age = any(!is.na(Age)),
    has_BF = any(!is.na(BF)),
    has_CVC = any(!is.na(CVC)),
    has_HRT = any(!is.na(HRT)),
    has_Meno = any(!is.na(MenoOnset)),
    has_HYPE = any(!is.na(HYPE)),
    has_T2D = any(!is.na(T2D)),
    has_CH = any(!is.na(Chronic)),
    .groups = 'drop'
  ) %>%
  summarise(across(starts_with("has_"), ~sum(.))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_unique_ids"
  ) %>%
  mutate(variable = gsub("has_", "", variable))

SENSEC<-lme(CVC ~ Region * STAGE+SEASON+Site_Order+BL,
            random = ~ 1 | ID, 
            data = SENSECVC,
            na.action = na.exclude,
            weights = varIdent(form = ~ 1 | STAGE),
            method = "REML")

emm_SENSEC <- emmeans(SENSEC, ~ Region | STAGE)
confint(pairs(emm_SENSEC))

SENSECB<-lme(BL ~ Region * STAGE+SEASON+Site_Order,
             random = ~ 1 | ID, 
             data = SENSECVC,
             na.action = na.exclude,
             weights = varIdent(form = ~ 1 | STAGE),
             method = "REML")
emm_SENSECB<-emmeans(SENSECB, ~Region)
pairs(emm_SENSECB)

####distibution check######

# body fat & vo2peak--> centre for possible interaction
CVCMAX1$BF_c     <- as.numeric(scale(CVCMAX1$BF,     scale = FALSE))
CVCMAX1$VO2Rel_c <- as.numeric(scale(CVCMAX1$VO2Rel, scale = FALSE))
CVCMAX1$BMI_c     <- as.numeric(scale(CVCMAX1$BMI,     scale = FALSE))

###PREDICTOR ANALYSIS########

calc_partial_r2 <- function(model_full, predictor_to_remove) {
  model_reduced <- update(model_full, as.formula(paste(". ~ . -", predictor_to_remove)))
  r2_reduced <- r.squaredGLMM(model_reduced)[1]
  partial_r2 <- r2_full - r2_reduced
  return(partial_r2)
}


###baseline CVCMax, participant characteristics 
C1<-lme(BL~YearsMeno+SEASON+
          BF+VO2Rel,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        weights = varIdent(form = ~ 1 | SEASON),
        method = "REML")
summary(C1)

##Age, sensitivity analysis 
CA1<-lme(BL~YearsMeno+Age+
           BF_c*VO2Rel_c,         
         random = ~ 1 | ID,  
         data = CVCMAX1,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(CA1)


r.squaredGLMM(C1)
r2_full<-r.squaredGLMM(C1)[1]
r2_full

CBL_partial_r2_results <- data.frame(
  Predictor = c("YearsMeno","SEASON",
                "BF", "VO2Rel"),
  Partial_R2 = c( calc_partial_r2(C1, "YearsMeno"),
                  calc_partial_r2(C1, "SEASON"),
                  calc_partial_r2(C1, "BF"),
                  calc_partial_r2(C1, "VO2Rel")  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)

print(CBL_partial_r2_results)

####baseline CVCmax, physiological characteristics

P1<-lme(BL~BLST+BLHR+DBP,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        weights = varIdent(form = ~ 1 | SEASON),
        method = "REML")
summary(P1)

r.squaredGLMM(P1)
r2_full<-r.squaredGLMM(P1)[1]
r2_full


PBL_partial_r2_results <- data.frame(
  Predictor = c("BLST","BLHR", "DBP"),
  Partial_R2 = c( calc_partial_r2(P1, "BLST"),
                  calc_partial_r2(P1, "BLHR"),
                  calc_partial_r2(P1, "DBP")
  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)

print(PBL_partial_r2_results)

####Heating!

C2<-lme(HEAT~YearsMeno+SEASON+
          BF+VO2Rel,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        weights = varIdent(form = ~ 1 | SEASON),
        method = "REML")
summary(C2)

##Age, sensitivity analysis 
CA2<-lme(HEAT~YearsMeno+Age+
           BF_c*VO2Rel_c+SEASON,         
         random = ~ 1 | ID,  
         data = CVCMAX1,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(CA2) ###AGE EXPLAINING THE VARIANCE AS YEARSMENO?

r.squaredGLMM(C2)
r2_full<-r.squaredGLMM(C2)[1]
r2_full

CBL_partial_r2_results <- data.frame(
  Predictor = c("YearsMeno", "SEASON",
                "BF", "VO2Rel"),
  Partial_R2 = c( calc_partial_r2(C2, "YearsMeno"),
                  calc_partial_r2(C2, "BF"),
                  calc_partial_r2(C2, "SEASON"),
                  calc_partial_r2(C2, "VO2Rel")
  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)


print(CBL_partial_r2_results)

#### physiological characteristics

P2<-lme(HEAT~ST_HEAT+HR_HEAT+DBP2_HEAT+BL,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        weights = varIdent(form = ~ 1 | SEASON),
        method = "REML")
summary(P2)

r.squaredGLMM(P2)
r2_full<-r.squaredGLMM(P2)[1]
r2_full


PBL_partial_r2_results <- data.frame(
  Predictor = c("ST_HEAT","HR_HEAT", "DBP2_HEAT", "BL"),
  Partial_R2 = c( calc_partial_r2(P2, "ST_HEAT"),
                  calc_partial_r2(P2, "HR_HEAT"),
                  calc_partial_r2(P2, "DBP2_HEAT"),
                  calc_partial_r2(P2, "BL")
  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)

print(PBL_partial_r2_results)
P2_bl<-lme(HEAT~BLST+BLHR+DBP+BL,         
           random = ~ 1 | ID,  
           data = CVCMAX1,
           na.action = na.exclude,
           method = "REML")
summary(P2_bl)

C3<-lme(NOS~YearsMeno+SEASON+
          BF+VO2Rel,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        weights = varIdent(form = ~ 1 | SEASON),
        method = "REML")
summary(C3)

##Age, sensitivity analysis 
CA3<-lme(NOS~YearsMeno+Age+
           BF_c*VO2Rel_c,         
         random = ~ 1 | ID,  
         data = CVCMAX1,
         na.action = na.exclude,
         method = "REML")
summary(CA3) 

r.squaredGLMM(C3)
r2_full<-r.squaredGLMM(C3)[1]
r2_full

CBL_partial_r2_results <- data.frame(
  Predictor = c("YearsMeno", "SEASON",
                "BF", "VO2Rel"),
  Partial_R2 = c( calc_partial_r2(C3, "YearsMeno"),
                  calc_partial_r2(C3, "BF"),
                  calc_partial_r2(C3, "SEASON"),
                  calc_partial_r2(C3, "VO2Rel")
  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)
print(CBL_partial_r2_results)



P3<-lme(NOS~ST_NOS+HR_NOS+DBP2_NOS+BL,         
        random = ~ 1 | ID,  
        data = CVCMAX1,
        na.action = na.exclude,
        method = "REML")
summary(P3)

r.squaredGLMM(P3)
r2_full<-r.squaredGLMM(P3)[1]
r2_full



PBL_partial_r2_results <- data.frame(
  Predictor = c("ST_NOS","HR_NOS", "DBP2_NOS", "BL"),
  Partial_R2 = c( calc_partial_r2(P3, "ST_NOS"),
                  calc_partial_r2(P3, "HR_NOS"),
                  calc_partial_r2(P3, "DBP2_NOS"),
                  calc_partial_r2(P3, "BL")
  )
) %>%
  arrange(desc(Partial_R2)) %>%
  mutate(Percent = Partial_R2 * 100)

print(PBL_partial_r2_results)

P3_bl<-lme(NOS~BLST+BLHR+DBP+BL,         
           random = ~ 1 | ID,  
           data = CVCMAX1,
           na.action = na.exclude,
           method = "REML")
summary(P3_bl)




####sensitivity analysis 

SENSE <- CVCMAX1[-which(CVCMAX1$T2D == "Y" | 
                          CVCMAX1$Chronic != "N" | 
                          CVCMAX1$HYPE == "Y" | 
                          CVCMAX1$MOOD == "Y" | 
                          CVCMAX1$HRT == "Y"), ]
SENSE$BF_c <- scale(SENSE$BF, scale = FALSE)
SENSE$VO2Rel_c <- scale(SENSE$VO2Rel, scale = FALSE)
SC1<-lme(BL~YearsMeno+SEASON+
           BF_c*VO2Rel_c,     
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(SC1)

SC2<-lme(HEAT~YearsMeno+SEASON+
           BF_c*VO2Rel_c,     
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(SC2)

SC3<-lme(NOS~YearsMeno+SEASON+
           BF_c*VO2Rel_c,     
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         method = "REML")
summary(SC3)

SP1<-lme(BL~BLST+BLHR+DBP,
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(SP1)

SP2<-lme(HEAT~ST_HEAT+HR_HEAT+DBP2_HEAT+BL,
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(SP2)

SP3<-lme(NOS~ST_NOS+HR_NOS+DBP2_NOS,
         random = ~ 1 | ID,  
         data = SENSE,
         na.action = na.exclude,
         weights = varIdent(form = ~ 1 | SEASON),
         method = "REML")
summary(SP3)




####effect sizes######

CVCMAX_s <- CVCMAX1 %>%
  mutate(
    BF_VO2_interact = BF_c * VO2Rel_c  # Create interaction first
  ) %>%
  mutate(across(c(YearsMeno, BF_c, VO2Rel_c, BF_VO2_interact, BLST, BLHR, MAP, BL), 
                ~as.numeric(scale(.))))

C1_std <- lme(BL ~ SEASON+YearsMeno+
                BF+VO2Rel,         
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(C1_std)
anova(C1_std)


# Extract standardized coefficients
BL_effects <- tidy(C1_std, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(BL_effects)


P1_std <- lme(BL~BLST+BLHR+DBP,        
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(P1_std)
anova(P1_std)

BL_effects <- tidy(P1_std, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(BL_effects)

C2_std <- lme(HEAT ~ SEASON+YearsMeno+
                BF+VO2Rel,         
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(C2_std)
anova(C2_std)

HEAT_effects <- tidy(C2_std, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(HEAT_effects)

P2_std <- lme(HEAT~ST_HEAT+HR_HEAT+DBP2_HEAT+BL,        
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(P2_std)
anova(P2_std)
# Extract standardized coefficients
HEAT_effects <- tidy(P2_std, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(HEAT_effects)

C3_std <- lme(NOS ~ SEASON+YearsMeno+
                BF+VO2Rel,         
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(C3_std)
anova(C3_std)





P3_std <- lme(NOS~ST_NOS+HR_NOS+DBP2_NOS+BL,        
              random = ~ 1 | ID,  
              data = CVCMAX_s,
              na.action = na.exclude,
              weights = varIdent(form = ~ 1 | SEASON),
              method = "REML")
summary(P3_std)
anova(P3_std)
# Extract standardized coefficients
NOS_effects <- tidy(C3_std, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

print(NOS_effects)



###PLOTS####


#####Significant predictors#####


CVCMAX1$partial_resid_YearsMeno <- residuals(C2) + 
  fixef(C2)["YearsMeno"] * CVCMAX1$YearsMeno
CVCMAX1$partial_resid_YearsMeno <- CVCMAX1$partial_resid_YearsMeno - mean(CVCMAX1$partial_resid_YearsMeno, na.rm = TRUE)

season_coefs <- fixef(C2)[grep("^SEASON", names(fixef(C2)))]
X_season <- model.matrix(~ SEASON, data = CVCMAX1)[, names(season_coefs), drop = FALSE]

CVCMAX1$partial_resid_SEASONHEAT <- residuals(C2) + as.vector(X_season %*% season_coefs)
CVCMAX1$partial_resid_SEASONHEAT <- CVCMAX1$partial_resid_SEASONHEAT -
  mean(CVCMAX1$partial_resid_SEASONHEAT, na.rm = TRUE)

CVCMAX1 <- CVCMAX1 |>
  mutate(Region = factor(Region,
                         levels = c("FOREARM", "CHEST", "ABDOMEN", "CALF"),
                         labels = c("Forearm", "Chest", "Abdomen", "Calf")))

region_cols <- c("Forearm" = "#C6DBEF",
                 "Chest"   = "#6BAED6",
                 "Abdomen" = "#1D9FB0",
                 "Calf"    = "#08306B")
CVCMAX1 <- CVCMAX1 |>
  mutate(
    n_condition = as.integer(Chronic == "T") + as.integer(Chronic == "C") +
      as.integer(HYPE == "Y") + as.integer(T2D == "Y"),
    Condition = factor(
      case_when(
        n_condition >= 2 ~ "Two or more conditions",
        n_condition == 1 ~ "One condition",
        TRUE             ~ "No conditions"
      ),
      levels = c("No conditions", "One condition", "Two or more conditions")
    )
  )
stopifnot(!any(is.na(CVCMAX1$Condition)))
table(CVCMAX1$Condition[!duplicated(CVCMAX1$ID)])

shp <- c("No conditions"           = 18,
         "One condition"           = 15,
         "Two or more conditions"  = 25)

key <- guides(shape = guide_legend(
  ncol = 1,
  override.aes = list(colour = "grey30", size = 2.5)
))

key <- guides(
  shape  = guide_legend(ncol = 1, override.aes = list(colour = "grey30", size = 2.5)),
  colour = guide_legend(ncol = 1)
)


Meno <- ggplot(CVCMAX1, aes(x = YearsMeno, y = partial_resid_YearsMeno)) +
  geom_point(aes(shape = Condition, colour = Region, fill = after_scale(colour)),
             show.legend = F) +
  scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_cols, name = "Region") +
  geom_smooth(method = "lm", se = FALSE, colour = "#253494") +
  geom_rug(color = "grey70", size = 0.5) +
  labs(x = "Menopausal duration (years)",
       y = expression(bold("Partial Residuals %CVC"["Max"]))) +
  key + theme_few() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_blank())
Meno
Season <- ggplot(CVCMAX1, aes(x = SEASON, y = partial_resid_SEASONHEAT)) +
  geom_boxplot(colour = "grey40", fill = "grey80", alpha = 0.7,
               outlier.shape = NA) +
  geom_jitter(aes(shape = Condition, colour = Region, fill = after_scale(colour)),
              size = 2, width = 0.1) +
  scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_cols, name = "Region") +
  labs(x = "Season",
       y = expression(bold("Partial Residuals %CVC"["Max"]))) +
  key + theme_few() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"),
        legend.box = "vertical",
        legend.title = element_blank())
Season

PRED<- (Meno | Season)
PRED
ggsave("PRED.png", PRED,dpi  =1000, height = 4, width = 10)


##### Group Plot#####
CVCPLOT<- CVCPLOT %>%
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
                                FOREARM = "Forearm")) 
CVCPLOT<-CVCPLOT[-which(CVCPLOT$STAGE == "Baseline" & CVCPLOT$CVCMAX >50 & CVCPLOT$Region != "Chest"),]


CVCPLOT <- CVCPLOT |>
  mutate(
    n_morbid = as.integer(Chronic == "T") + as.integer(Chronic == "C") +
      as.integer(HYPE == "Y") + as.integer(T2D == "Y"),
    condition = factor(
      case_when(
        n_morbid >= 2 ~ "Two or more conditions",
        n_morbid == 1 ~ "One condition",
        TRUE          ~ "No conditions"
      ),
      levels = c("No conditions", "One condition", "Two or more conditions")
    )
  )

CPLOT <- CPLOT |>
  mutate(
    n_morbid = as.integer(Chronic == "T") + as.integer(Chronic == "C") +
      as.integer(HYPE == "Y") + as.integer(T2D == "Y"),
    condition = factor(
      case_when(
        n_morbid >= 2 ~ "Two or more conditions",
        n_morbid == 1 ~ "One condition",
        TRUE          ~ "No conditions"
      ),
      levels = c("No conditions", "One condition", "Two or more conditions")
    )
  )

stopifnot(!any(is.na(CVCPLOT$condition)), !any(is.na(CPLOT$condition)))

shp <- c("No conditions"           = 18,
         "One condition"           = 15,
         "Two or more conditions"  = 25)

region_fill_vals <- c("#C6DBEF", "#6BAED6", "#1D9FB0", "#0D6B78", "#08306B")

PLOT1 <- ggplot(data = CVCPLOT, aes(x = Region, y = CVCMAX, fill = Region)) +
  scale_fill_manual(values = region_fill_vals, guide = "none") +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
               width = 0.6, lwd = 0.4) +
  geom_point(aes(shape = condition, colour = Region, fill = after_scale(colour)),
             size = 1.5,
             position = position_jitter(width = 0.1, seed = 1),
             alpha = 0.7, stroke = 0.3) +
  scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_fill_vals, guide = "none") +
  guides(shape = guide_legend(nrow = 1,
                              override.aes = list(colour = "grey30", fill = "grey30",
                                                  size = 2.5, alpha = 1))) +
  labs(x = "Region", y = expression(bold("%CVC"["Max"]))) +
  facet_wrap(~STAGE, ncol = 4) +
  theme_few() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     expand = expansion(mult = c(0.05, 0.1)))

PLOT1

CPLOT<-CPLOT[which(CPLOT$STAGE != "NOS"),]
CPLOT<- CPLOT %>%
  mutate(STAGE = recode_factor(STAGE,
                               BL = "Baseline",
                               HEAT = "Heating Plateau",
                               LNAME = "L-NAME Plateau",
                               MAX = "Maximum Vasodilation")) %>%
  mutate(Region = recode_factor(Region,
                                POOLED = "All Regions",
                                ABDOMEN = "Abdomen",
                                CALF = "Calf",
                                CHEST = "Chest",
                                FOREARM = "Forearm")) 

PLOT2 <- ggplot(data = CPLOT, aes(x = Region, y = CVC, fill = Region)) +
  scale_fill_manual(values = region_fill_vals, guide = "none") +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, coef = 0, alpha = 0.5,
               width = 0.6, lwd = 0.4) +
  geom_point(aes(shape = condition, colour = Region, fill = after_scale(colour)),
             size = 2,
             position = position_jitter(width = 0.1, seed = 1),
             alpha = 0.7, stroke = 0.3) +
  scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_fill_vals, guide = "none") +
  guides(shape = guide_legend(nrow = 1,
                              override.aes = list(colour = "grey30", fill = "grey30",
                                                  size = 2.5, alpha = 1))) +
  labs(x = "Region", y = expression(bold("CVC (PU·mmHg"^{-1}~")"))) +
  facet_wrap(~STAGE, ncol = 4) +
  theme_few() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        strip.text = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm")) +
  scale_y_continuous(breaks = seq(0, 10, by = 1),
                     expand = expansion(mult = c(0.05, 0.1)))

PLOT2
ggsave("CVCMAX_OF.png", PLOT1, dpi = 700, height = 6, width = 9)
ggsave("CVC_OF.png", PLOT2, dpi = 700, height = 6, width = 9)


#####Supplement A1, A2, A3#####




labs_char <- c(
  YearsMeno = "Menopause duration (years)",
  BF        = "Body fat (%)",
  VO2Rel    = "VO2peak (mL/kg/min)"
)

pp <- CVCMAX1[!duplicated(CVCMAX1$ID), ]

d_char <- pp |>
  select(all_of(names(labs_char))) |>
  pivot_longer(everything(), names_to = "var", values_to = "value") |>
  filter(!is.na(value), is.finite(value)) |>
  mutate(var = factor(var, levels = names(labs_char), labels = labs_char))

ann_char <- d_char |> count(var) |> mutate(lab = sprintf("n = %d", n))

A2<-ggplot(d_char, aes(value)) +
  geom_density(fill = "#4C72B0", colour = "#2A4A7F",
               linewidth = 0.4, alpha = 0.45) +
  geom_rug(sides = "b", alpha = 0.4, length = unit(0.03, "npc"),
           colour = "grey40") +
  geom_text(data = ann_char, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.15, vjust = 1.4, size = 2.7,
            colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~ var, scales = "free", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     breaks = scales::breaks_pretty(n = 3)) +
  labs(x = NULL, y = "Density") +
  theme_bw(base_size = 9) +
  theme(
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(size = 7.5, face = "bold",
                                      margin = margin(3,2,3,2)),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(colour = "grey80"),
    axis.text          = element_text(size = 7),
    axis.title.y       = element_text(size = 10, face = "bold")
  )
A2
ggsave("A2.png", A2,dpi  =1000, height = 4, width = 10)

labs_phys <- c(
  BLST      = "Skin temperature (°C)",
  ST_HEAT   = "Skin temperature (°C)",
  ST_NOS    = "Skin temperature (°C)",
  BLHR      = "Heart rate (bpm)",
  HR_HEAT   = "Heart rate (bpm)",
  HR_NOS    = "Heart rate (bpm)",
  DBP       = "Diastolic blood pressure (mmHg)",
  DBP2_HEAT = "Diastolic blood pressure (mmHg)",
  DBP2_NOS  = "Diastolic blood pressure (mmHg)"
)

phase_map <- c(
  BLST = "Baseline", ST_HEAT = "Heating plateau", ST_NOS = "NOS-dependent phase",
  BLHR = "Baseline", HR_HEAT = "Heating plateau", HR_NOS = "NOS-dependent phase",
  DBP  = "Baseline", DBP2_HEAT = "Heating plateau", DBP2_NOS = "NOS-dependent phase"
)

d_phys <- pp |>
  select(all_of(names(labs_phys))) |>
  pivot_longer(everything(), names_to = "col", values_to = "value") |>
  filter(!is.na(value), is.finite(value)) |>
  mutate(
    var   = factor(labs_phys[col], levels = unique(labs_phys)),
    Phase = factor(phase_map[col],
                   levels = c("Baseline", "Heating plateau", "NOS-dependent phase"))
  )

ann_phys <- d_phys |> count(var, Phase) |> mutate(lab = sprintf("n = %d", n))

d_st   <- d_phys[d_phys$var == "Skin temperature (°C)", ]
ann_st <- ann_phys[ann_phys$var == "Skin temperature (°C)", ]

row_st <- ggplot(d_st, aes(value)) +
  geom_density(fill = "#4C72B0", colour = "#4C72B0",
               linewidth = 0.4, alpha = 0.45) +
  geom_rug(sides = "b", alpha = 0.4, length = unit(0.03, "npc"), colour = "grey40") +
  geom_text(data = ann_st, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.15, vjust = 1.4, size = 2.5,
            colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~ Phase, scales = "free_x", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = scales::breaks_pretty(n = 3)) +
  labs(x = NULL, y = "Skin temperature (°C)") +
  theme_bw(base_size = 9) +
  theme(
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(size = 7.5, face = "bold",
                                      margin = margin(3,2,3,2)),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(colour = "grey80"),
    axis.text          = element_text(size = 6.5),
    axis.title.y       = element_text(size = 8, face = "bold")
  )


d_hr   <- d_phys[d_phys$var == "Heart rate (bpm)", ]
ann_hr <- ann_phys[ann_phys$var == "Heart rate (bpm)", ]

row_hr <- ggplot(d_hr, aes(value)) +
  geom_density(fill = "#DD8452", colour = "#DD8452",
               linewidth = 0.4, alpha = 0.45) +
  geom_rug(sides = "b", alpha = 0.4, length = unit(0.03, "npc"), colour = "grey40") +
  geom_text(data = ann_hr, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.15, vjust = 1.4, size = 2.5,
            colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~ Phase, scales = "free_x", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = scales::breaks_pretty(n = 3)) +
  labs(x = NULL, y = "Heart rate (bpm)") +
  theme_bw(base_size = 9) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(colour = "grey80"),
    axis.text          = element_text(size = 6.5),
    axis.title.y       = element_text(size = 8, face = "bold")
  )

d_dbp   <- d_phys[d_phys$var == "Diastolic blood pressure (mmHg)", ]
ann_dbp <- ann_phys[ann_phys$var == "Diastolic blood pressure (mmHg)", ]

row_dbp <- ggplot(d_dbp, aes(value)) +
  geom_density(fill = "#55A868", colour = "#55A868",
               linewidth = 0.4, alpha = 0.45) +
  geom_rug(sides = "b", alpha = 0.4, length = unit(0.03, "npc"), colour = "grey40") +
  geom_text(data = ann_dbp, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.15, vjust = 1.4, size = 2.5,
            colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~ Phase, scales = "free_x", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = scales::breaks_pretty(n = 3)) +
  labs(x = NULL, y = "Diastolic blood pressure (mmHg)") +
  theme_bw(base_size = 9) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(colour = "grey80"),
    axis.text          = element_text(size = 6.5),
    axis.title.y       = element_text(size = 8, face = "bold")
  )

A3 <- row_st / row_hr / row_dbp

A3

ggsave("A3.png",A3, width = 10, height = 10, dpi = 1000)
#####Supplemental Figure A4 & A5####

CVCMAX1 <- CVCMAX1 |>
  mutate(
    n_medication = as.integer(HYPE == "Y") + as.integer(T2D == "Y") +
      as.integer(Chronic == "C") + as.integer(Chronic == "T") +
      as.integer(MOOD == "Y") + as.integer(HRT == "Y"),
    Medication = factor(
      case_when(
        n_medication >= 2  ~ "More than one medication",
        HYPE == "Y"        ~ "Hypertensive medication",
        T2D  == "Y"        ~ "Type 2 diabetes medication",
        Chronic == "C"     ~ "Cholesterol medication",
        Chronic == "T"     ~ "Levothyroxine",
        MOOD == "Y"        ~ "Psychotropic medication",
        HRT  == "Y"        ~ "Hormone replacement therapy",
        TRUE               ~ "No medication"
      ),
      levels = c("No medication", "Hypertensive medication",
                 "Type 2 diabetes medication", "Cholesterol medication",
                 "Levothyroxine", "Psychotropic medication",
                 "Hormone replacement therapy", "More than one medication")
    )
  )
stopifnot(!any(is.na(CVCMAX1$Medication)))
table(CVCMAX1$Medication[!duplicated(CVCMAX1$ID)])

med_shp <- c("No medication"                = 16,
             "Hypertensive medication"       = 1,
             "Type 2 diabetes medication"    = 2,
             "Cholesterol medication"        = 0,
             "Levothyroxine"                 = 5,
             "Psychotropic medication"       = 6,
             "Hormone replacement therapy"   = 3,
             "More than one medication"      = 4)

med_key <- guides(shape = guide_legend(
  ncol = 1,
  override.aes = list(colour = "grey30", size = 2.5)
))

A4 <- ggplot(CVCMAX1, aes(x = YearsMeno, y = partial_resid_YearsMeno)) +
  geom_point(aes(shape = Medication, colour = Region),show.legend = F) +
  scale_shape_manual(values = med_shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_cols, name = "Region") +
  geom_smooth(method = "lm", se = FALSE, colour = "#253494") +
  geom_rug(color = "grey70", size = 0.5) +
  labs(x = "Menopausal duration (years)",
       y = expression(bold("Partial Residuals %CVC"["Max"]))) +
  key + theme_few() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_blank())
A4


A5 <- ggplot(CVCMAX1, aes(x = SEASON, y = partial_resid_SEASONHEAT)) +
  geom_boxplot(colour = "grey40", fill = "grey80", alpha = 0.7,
               outlier.shape = NA) +
  geom_jitter(aes(shape = Medication, colour = Region),
              size = 2, width = 0.05) +
  scale_shape_manual(values = med_shp, name = NULL, drop = FALSE) +
  scale_colour_manual(values = region_cols, name = "Region") +
  labs(x = "Season",
       y = expression(bold("Partial Residuals %CVC"["Max"]))) +
  key + theme_few() +
  theme(axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold"),
        legend.box = "vertical",
        legend.title = element_blank())
A5

PRED2<- (A4 | A5)

PRED2
ggsave("PRED2.png", PRED2,dpi  =1000, height = 4, width = 10)



