##  Analysis of Pe299R in competition with a second epiphyte in Arabidopsis thaliana (CFU) ##


##    Dependencies
source("~/cusper_libraries.R")

##    DATA
# The presented files contain CFU gFW^-1 data of Pe299R (peuca299) or total in linear and logarithmic scale in the presence of a competitor (trt) at different sampling points (time/hpi).
# cfu data was used to calculate competition scores (comp_rate and comp_scale). df is a data frame containing the relationship between Pe299R and a second strain. Information about
# competition scores, phylogenetic distances (pd), MRO (complete to leaf5), and phylogroup are provided.
cfu = readRDS("cfu_inplanta.rds")
df = readRDS("arabidopsis_roi.rds")

##    Analysis CFU data
hist(cfu$cfu)
qqnorm(cfu$cfu)
plot(cfu$trt, cfu$cfu)
ncvTest(lm(cfu ~ trt * time, data = cfu)) # Breusch-Pagan test

# ANOVA
aov_cfu = aov(cfu ~ trt * time, data = cfu)
summary.aov(aov_cfu)
plot(resid(aov_cfu))
etaSquared(aov_cfu)

# Multiple Comparisons
emm.cfu = emmeans(aov_cfu, ~ trt | time, adjust = "bonferroni")
contrast(emm.cfu, "trt.vs.ctrl", ref = "mono", adjust = "bonferroni")

##    Regression models
df2_arabidopsis = df %>% filter(competitor != "peuca299")

lm0 = lm(comp_scale ~ pd, data = df2_arabidopsis)
lm1 = lm(comp_scale ~ complete * pd, data = df2_arabidopsis)
lm2 = lm(comp_scale ~ invitro * pd, data = df2_arabidopsis)
lm3 = lm(comp_scale ~ leaf1 * pd, data = df2_arabidopsis)
lm4 = lm(comp_scale ~ leaf2 * pd, data = df2_arabidopsis)
lm5 = lm(comp_scale ~ leaf3 * pd, data = df2_arabidopsis)
lm5.2 = lm(comp_scale ~ leaf3, data = df2_arabidopsis)
lm6 = lm(comp_scale ~ leaf4 * pd, data = df2_arabidopsis)
lm7 = lm(comp_scale ~ leaf5 * pd, data = df2_arabidopsis)

AIC(lm0,lm1,lm2,lm3,lm4,lm5,lm6,lm7)

summary(lm0)
summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)
summary(lm5) # best model
summary(lm6)
summary(lm7)

anova(lm5)
etaSquared(lm5, type = 1, anova = TRUE)
qqnorm(resid(lm5))
confint(lm5)