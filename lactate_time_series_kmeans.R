# Author: Michelle Fang
# Date: 7/4/2025 
# K-means clustering of lactate

source("/studies/cardiac/support/cpb/adequate/analyses/fangm/utils.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/cluster/utils_cluster.R")
library(mgcv)


# Histogram of maximum recorded lactate time ----
max_iv_df <- cldta_onpump %>%
  group_by(surg_id) %>%
  summarize(max_iv_evth = max(iv_evth, na.rm = TRUE), .groups = "drop")

ggplot(max_iv_df, aes(x = max_iv_evth)) +
  geom_histogram(fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "Time of Last Recorded Lactate per Patient",
    x = "Time recorded (hours)",
    y = "Number of Patients"
  ) +
  theme_minimal()

sum(max_iv_df$max_iv_evth > 4)


# Time grid code ----

target_time_grid <- sort(unique(cldta_onpump$iv_evth))
target_time_grid <- target_time_grid[target_time_grid >= 0 & target_time_grid <= 4]
min(target_time_grid)
max(target_time_grid)
length(target_time_grid)

lactate_df <- cldta_onpump_lactate_mean %>%
  rename(lactate = mean_lactate)

valid_lactate_df <- lactate_df %>%
  group_by(surg_id) %>%
  filter(sum(!is.na(lactate)) >= 2) %>%
  ungroup()

all_ids <- unique(lactate_df$surg_id)
valid_ids <- unique(valid_lactate_df$surg_id)
length(setdiff(all_ids, valid_ids))

# LME ----

library(lme4)

set.seed(123)

lmer_lactate <- lmer(
  lactate ~ iv_evth + (1 + iv_evth | surg_id),
  data = valid_lactate_df,
  REML = TRUE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(lmer_lactate)

pred_grid <- expand.grid(
  surg_id = unique(valid_lactate_df$surg_id),
  iv_evth = target_time_grid
)

lactate_interp_df <- pred_grid %>%
  mutate(
    lactate = predict(
      lmer_lactate,
      newdata = pred_grid,
      re.form = NULL   # include random effects (BLUPs)
    )
  )


lactate_mat <- lactate_interp_df %>%
  pivot_wider(names_from = iv_evth, values_from = lactate) %>%
  arrange(surg_id)

surg_ids <- lactate_mat$surg_id
lactate_matrix <- lactate_mat %>% select(-surg_id) %>% as.matrix()


# Elbow Plot ----
# wss <- numeric()
# max_k <- 10
# 
# for (k in 1:max_k) {
#   set.seed(123)
#   print(k)
#   km <- kmeans(lactate_matrix, centers = k, nstart = 20)
#   wss[k] <- km$tot.withinss
# }
# 
# plot(1:max_k, wss, type = "b", pch = 19,
#      xlab = "Number of clusters k",
#      ylab = "Total within-cluster sum of squares",
#      main = "Elbow Method for Choosing k")

# K = 4 Clustering ---- 
k <- 4
cluster_df <- get_kmeans_cluster_df(lactate_matrix, k = k)
# saveRDS(cluster_df, file = file.path(figure_dir, "cluster_df.rds"))

## Compare between clusters ---- 
get_spaghetti_lactate(cluster_df, x_lim = 4, y_lim = 10)
get_tables_comparing_clusters(cluster_df)
