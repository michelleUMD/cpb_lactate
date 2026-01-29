# CPB Lactate Code

Author: Michelle Fang\
Date: 2026-01-12

Note: some code has been removed to abide by patient confidentiality and HIPPA

## File Tree

```         
.
├── xgboost_end_cpb_lactate.R: XGBoost end-CPB lactate modeling 
├── xgboost_peak_lactate.R: XGBoost peak lactate modeling 
├── xgboost_cross_clamp_lactate.R: XGBoost end-cross-clamp lactate modeling 
├── utils_xgboost.R: utility functions xgboost modeling 
├── utils_3d.R: utility functions for 2D and 3D figures 
├── xgb_models
|   └── json objects for end-lactate XGBoost models
|
├── do2_pressor_perturb.R: hypothetical lactate if DO2 increased and pressors decreased 
|
├── lactate_time_series_kmeans.R: k-means clustering of lactate as time series
├── utils_cluster.R: utility functions for k means clustering
|
├── map_ci_modeling.R: model MAP vs CI  
└── utils_map_ci_modeling.R: utility functions for MAP vs CI  
```
