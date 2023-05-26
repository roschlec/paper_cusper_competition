##  Analysis of Pe299R CUSPER in competition with a second epiphyte in Arabidopsis thaliana (single cell) ##

##    Dependencies
source("~/cusper_libraries.R")
set.seed(123131)

##    DATA
rs = readRDS("cusper_rs.rds") # Relative fraction of CUSPER subpopulations at different time points and in the presence of different competitors
rs.df = readRDS("cusper_rs_df.rds") # Data frame of 'rs' with MRO, PD, phylogroups provided
rs.mat = as.matrix(rs.df[,8:length(rs.df)]) # Matrix of rs.df
cusper_increase = readRDS("cusper_increase.rds") # Population increase (logsum) based on CUSPER cells
cusper_vs_cfu_increase = readRDS("cusper_vs_cfu_increase.rds") # Population increase (logsum) based on CUSPER cells vs population increase by CFU data (increase)

####  Multivariate analysis ####
# NMDS
rs.mds = metaMDS(rs.mat, k = 2)
stressplot(rs.mds)
plot(rs.mds)

# Bray-Curtis distances between samples
rs.dist = vegdist(rs.mat, method = "bray")

# Multivariate dispersions
rs.dispersion = betadisper(rs.dist, rs.df$competitor)
plot(rs.dispersion, hull = FALSE, ellipse = TRUE)

# Permutation test
permutest(rs.dispersion, permutations = 999, pairwise = TRUE)

#  PERMANOVA
adonis2(rs.dist ~ time, data = rs.df, permutations = 999, method = "bray")
adonis2(rs.dist ~ competitor, data = rs.df, permutations = 999, method = "bray")
adonis2(rs.dist ~ time + competitor, data = rs.df, permutations = 999, method = "bray")

#   Pairwise comparisons
pairwise_time = pairwise.adonis2(rs.dist ~ time, data = rs.df, permutations = 999, method = "bray")
pairwise_comp = pairwise.adonis2(rs.dist ~ competitor, data = rs.df, permutations = 999, method = "bray")


####  Population increase ####
#   Histogram
hist(cusper_increase$logsum[cusper_increase$t==0])
hist(cusper_increase$logsum[cusper_increase$t==24])
hist(cusper_increase$logsum[cusper_increase$t==36])

# Brausch-Pagan test
ncvTest(lm(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t0",]))
ncvTest(lm(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t1",]))
ncvTest(lm(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t2",]))

# ANOVA
aov_increase_t0 = aov(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t0",])
aov_increase_t1 = aov(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t1",])
aov_increase_t2 = aov(logsum ~ sample, data = cusper_increase[cusper_increase$time=="t2",])

# ANOVA Table
summary(aov_increase_t0)
summary(aov_increase_t1)
summary(aov_increase_t2)


# Compared to CFU
hist(cusper_vs_cfu_increase$logsum[cusper_vs_cfu_increase$time!="t0"])
hist(cusper_vs_cfu_increase$loginc[cusper_vs_cfu_increase$time!="t0"])

lm_increase = lm(loginc ~ logsum, data = cusper_vs_cfu_increase[cusper_vs_cfu_increase$time!="t0",])
ncvTest(lm_increase)
summary(lm_increase)
qqnorm(resid(lm_increase))
anova(lm_increase)
cor(cusper_vs_cfu_increase$logsum, cusper_vs_cfu_increase$loginc)
