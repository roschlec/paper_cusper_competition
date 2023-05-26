####  Effect of MRO in vitro  ####

##    Dependencies
source("~/cusper_libraries.R")

##    DATA
df = readRDS("invitro_roi_full.rds")
df_mean = readRDS("invitro_roi.rds")
df_estimates = df_mean[-1,2:7]
colnames(df_estimates) = c("Profile_all", "Profile_r", "Profile_K", "Growth", "MRO", "PD")

####  Correlations  ####
pearson_r = cor(df_estimates, method = "pearson")
pearson_r_correlationmatrix <- rcorr(as.matrix(df_estimates), type = "pearson")

####  Regression models ####
# Pe299R was excluded from this analysis
df_noPe299R = df[df$competitor!="peuca299",]
df_mean_noPe299R = df_mean[df_mean$competitor!="peuca299",]

# Pe299R and PssB728a were excluded from this analysis
df_noPe299R_noPssB728a = df_noPe299R[df_noPe299R$competitor!="psyrib728",]

##    Carbon profile vs MRO
lmMRO_Cprof = lm(Cprof ~ mro, data = df_mean_noPe299R)
  ncvTest(lmMRO_Cprof)  # Breusch-Pagan test of equal variances
  lmMRO_Cprof %>% resid %>% qqnorm  # Distribution of residuals
  summary(lmMRO_Cprof)  # Summary regression model
  etaSquared(lmMRO_Cprof) # Effect size (partial R2)

lmMRO_Cprof_r = lm(Cprof_r ~ mro, data = df_mean_noPe299R)
  ncvTest(lmMRO_Cprof_r)  # Breusch-Pagan test of equal variances
  lmMRO_Cprof_r %>% resid %>% qqnorm  # Distribution of residuals
  summary(lmMRO_Cprof_r)  # Summary regression model
  etaSquared(lmMRO_Cprof_r) # Effect size (partial R2)

lmMRO_Cprof_k = lm(Cprof_k ~ mro, data = df_mean_noPe299R)
  ncvTest(lmMRO_Cprof_k)  # Breusch-Pagan test of equal variances
  lmMRO_Cprof_k %>% resid %>% qqnorm  # Distribution of residuals
  summary(lmMRO_Cprof_k)  # Summary regression model
  etaSquared(lmMRO_Cprof_k) # Effect size (partial R2)

##    Competition vs MRO
lmMRO_competition1 = lm(comp_scale ~ mro, data = df_noPe299R)
  ncvTest(lmMRO_competition1)  # Breusch-Pagan test of equal variances
  lmMRO_competition1 %>% resid %>% qqnorm  # Distribution of residuals
  summary(lmMRO_competition1)  # Summary regression model
  confint(lmMRO_competition1)  # Confidence intervals
  etaSquared(lmMRO_competition1) # Effect size (partial R2)

lmMRO_competition2 = lm(comp_scale ~ mro, data = df_noPe299R_noPssB728a)
  ncvTest(lmMRO_competition2)  # Breusch-Pagan test of equal variances
  lmMRO_competition2 %>% resid %>% qqnorm  # Distribution of residuals
  summary(lmMRO_competition2)  # Summary regression model
  confint(lmMRO_competition2)  # Confidence intervals
  etaSquared(lmMRO_competition2) # Effect size (partial R2)

glmMRO = glm(comp_scale ~ mro*pd, data = df_noPe299R, family=Gamma(link = 'log'))
  glmMRO %>% resid(., type="resp") %>% qqnorm  # Distribution of residuals
  summary(glmMRO)  # Summary regression model
  confint(glmMRO)  # Confidence intervals
  (pseudoR2 <- 1-(glmMRO$deviance/glmMRO$null.deviance))  # Pseudo-R2
  anova(glmMRO, test='F')  # ANOVA on GLM
  emmglmMRO = emmeans(glmMRO, specs = ~ mro * pd, data = df_noPe299R)  # Multiple comparison
  summary(emmglmMRO, type="response")