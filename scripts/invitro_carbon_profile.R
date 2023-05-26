#### In vitro carbon profile

##    Dependencies
source("~/cusper_libraries.R")

##    DATA
cProfile = readRDS("carbon_profile_curves.rds")
num_analyses <- length(names(cProfile)) - 1

##    Logistic curve fitting
df_cProfile = data.frame(strain = character(num_analyses),
                     k = numeric(num_analyses),
                     r = numeric(num_analyses),
                     auc_l = numeric(num_analyses),
                     auc_e = numeric(num_analyses),
                     stringsAsFactors = FALSE)
lcurve = list()
n = 1

par(mfcol = c(7,5))
par(mar = c(0.25,0.25,0.25,0.25))

for (col_name in names(cProfile)){
  if (col_name != "time"){
    if(col_name == "malate_meth85" | col_name == "methanol_meth85"){d_loop <- cProfile[, c("time", col_name)]}
    else{d_loop <- cProfile[-6:-5, c("time", col_name)]}
    min_value <- min(d_loop[, col_name])
    d_loop[,col_name] <- d_loop[,col_name] - min_value
    gc_fit <- SummarizeGrowth(data_t = d_loop[, "time"],
                              data_n = d_loop[, col_name],
                              bg_correct = "none")
    df_cProfile$strain[n] <- col_name
    df_cProfile[n,2:length(df_cProfile)] = c(gc_fit$vals$k,
                                     gc_fit$vals$r,
                                     gc_fit$vals$auc_l,
                                     gc_fit$vals$auc_e)
    n = n + 1
    n_obs <- length(gc_fit$data$t)
    idx_to_plot <- 1:100/100*n_obs
    plot(gc_fit$data$t[idx_to_plot], gc_fit$data$N[idx_to_plot],
         pch=20,
         xlim=c(0,130),
         ylim=c(0,1),
         cex=0.6, xaxt="n", yaxt="n")
    text(x=24, y=1, labels=col_name, pos=1)
    if (gc_fit$vals$r > 0) {
      lines(gc_fit$data$t, predict(gc_fit$model), col="red")
    }
  }
}
dev.off()

df_cProfile = df_cProfile %>% separate(strain, into = c("cond", "strain"), sep = "_")

#   Rescaling data
df_cProfile_std = df_cProfile %>% 
  group_by(cond) %>% 
  mutate(k = scale(k, center = TRUE, scale = TRUE),
         r = scale(r, center = TRUE, scale = TRUE))

#   Data frame to Matrix
df_cProfile_mat = df_cProfile_std %>% 
  pivot_longer(k:r) %>% 
  unite("label", c("cond", "name")) %>% 
  pivot_wider(id_cols = strain, names_from = label, values_from = value) %>% 
  column_to_rownames(var="strain") %>% 
  as.matrix()

##    Distance Matrix K&r
dist.mat = dist(df_cProfile_mat, method = "euclidean")
clust = hclust(dist.mat, method = "complete")
clust = dendextend::rotate(clust, c("peuca299", "arth145", "psyrib728", "pkore", "rhod225", "meth85", "smelo"))
plot(clust)
dist.mat.all = melt(as.matrix(dist.mat)) %>% filter(Var2 == "peuca299") %>%dplyr::select(competitor = Var1, Cprof = value)

##    Distance Matrix r
dist.r = dist(df_cProfile_mat[,c(2,4,6,8,10)], method = "euclidean")
clust.r = hclust(dist.r, method = "complete")
plot(clust.r)
clust.r = dendextend::rotate(clust.r, c("peuca299", "pkore", "arth145", "psyrib728", "meth85", "smelo", "rhod225"))
dist.mat.r = melt(as.matrix(dist.r)) %>% filter(Var2 == "peuca299") %>%dplyr::select(competitor = Var1, Cprof_r = value)

##    Distance Matrix K
dist.K = dist(df_cProfile_mat[,c(1,3,5,7,9)], method = "euclidean")
clust.K = hclust(dist.mat, method = "complete")
plot(clust.K)
clust.K = dendextend::rotate(clust.K, c("peuca299", "arth145", "psyrib728", "pkore", "rhod225", "meth85", "smelo"))
dist.mat.k = melt(as.matrix(dist.K)) %>% filter(Var2 == "peuca299") %>%dplyr::select(competitor = Var1, Cprof_k = value)