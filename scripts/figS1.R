###   SUPPLEMENTAL FIGURE S1


##    Dependencies
source("~/cusper_libraries.R")

# DATA
data_full = read.csv("pe_rfu.csv", header = T)
data_full = data_full %>% dplyr::select("time_min":"X50_4") %>% mutate(trim = "full")

# DATA TRIMMING
data_trim1 = data_full[seq(1,nrow(data_full),4),] %>% mutate(trim = "trim1")
data_trim2 = data_full[seq(1,nrow(data_full),8),] %>% mutate(trim = "trim2")
data_trim3 = data_full[seq(1,nrow(data_full),16),] %>% mutate(trim = "trim3")
data_trim4 = data_full[seq(1,nrow(data_full),20),] %>% mutate(trim = "trim4")
data_trim5 = data_full[c(1,15,70,80),] %>% mutate(trim = "trim5")

d = rbind(data_full,data_trim1,data_trim2,data_trim3,data_trim4,data_trim5)
d = d[,c(26,1:5)]
trim_lev = levels(factor(d$trim))

# LOGISTIC CURVE FITTING
num_analyses <- length(names(d)) - 2
d_gc <- data.frame(trim = character(num_analyses),
                   sample = character(num_analyses),
                   k = numeric(num_analyses),
                   n0  = numeric(num_analyses),
                   r = numeric(num_analyses),
                   t_mid = numeric(num_analyses),
                   t_gen = numeric(num_analyses),
                   auc_l = numeric(num_analyses),
                   auc_e = numeric(num_analyses),
                   sigma = numeric(num_analyses),
                   stringsAsFactors = FALSE)
d_gc_full = data.frame()

par(mfcol = c(4,6))
par(mar = c(0.25,0.25,0.25,0.25))
n <- 1    # keeps track of the current row in the output data frame

for (i in 1:length(trim_lev)) {
  d2 = d %>% filter(trim == trim_lev[i]) %>% dplyr::select(-trim)
  n <- 1
  for (col_name in names(d2)) {
    
    # Don't process the column called "time". 
    # It contains time and not absorbance data.
    if (col_name != "time_min") {
      
      # Create a temporary data frame that contains just the time and current col
      d_loop <- d2[, c("time_min", col_name)]
      
      # Do the background correction.
      # Background correction option 1: subtract the minimum value in a column
      #                                 from all measurements in that column
      min_value <- min(d_loop[, col_name])
      d_loop[, col_name] <- d_loop[, col_name] - min_value
      
      # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
      gc_fit <- SummarizeGrowth(data_t = d_loop[, "time_min"], 
                                data_n = d_loop[, col_name],
                                bg_correct = "none")
      
      # Now, add the metrics from this column to the next row (n) in the 
      # output data frame, and increment the row counter (n)
      d_gc$sample[n] <- col_name
      d_gc$trim[n] <- trim_lev[i]
      d_gc[n, 3:10] <- c(gc_fit$vals$k,
                        gc_fit$vals$n0,
                        gc_fit$vals$r,
                        gc_fit$vals$t_mid,
                        gc_fit$vals$t_gen,
                        gc_fit$vals$auc_l,
                        gc_fit$vals$auc_e,
                        gc_fit$vals$sigma)
      
      d_gc_full <- rbind(d_gc_full, d_gc)
      n <- n + 1
      
      # Finally, plot the raw data and the fitted curve
      # Here, I'll just print some of the data points to keep the file size smaller
      n_obs <- length(gc_fit$data$t)
      idx_to_plot <- 1:20 / 20 * n_obs
      plot(gc_fit$data$t[idx_to_plot], gc_fit$data$N[idx_to_plot], 
           pch = 20,
           cex = 0.6, xaxt = "n", yaxt = "n")
      lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
    }
  }
}
dev.off()

# ANOVA
anova_trim = aov(r ~ trim, data = d_gc_full)
summary(anova_trim)
