# Fn for grid search RF 

# Construct a large Cartesian hyper-parameter space
ntrees_opts = c(10000)       # early stopping will stop earlier
max_depth_opts = seq(1,20)
min_rows_opts = c(1,5,10,20,50,100)
learn_rate_opts = seq(0.001,0.01,0.001)
sample_rate_opts = seq(0.3,1,0.05)
col_sample_rate_opts = seq(0.3,1,0.05)
col_sample_rate_per_tree_opts = seq(0.3,1,0.05)
mtries_opts = seq(2,20,2)
#nbins_cats_opts = seq(100,10000,100) # no categorical features
# in this dataset

fold_column = "block"
data_in = test22
predictors = predictor_columns
response = "ENS_Plants"#"richness_Plants"


random_grid_h2oRF = function(data_in, predictors, response,fold_column,
                             ntrees_opts = ntrees_opts, ntrees_opts = ntrees_opts,
                             min_rows_opts = min_rows_opts, learn_rate_opts=learn_rate_opts, 
                             sample_rate_opts = sample_rate_opts, 
                             col_sample_rate_per_tree_opts = col_sample_rate_per_tree_opts){
  
  hyper_params = list( ntrees = ntrees_opts,
                       max_depth = max_depth_opts,
                       min_rows = min_rows_opts,
                       # learn_rate = learn_rate_opts,
                       mtries = mtries_opts,
                       sample_rate = sample_rate_opts#,
                       # col_sample_rate = col_sample_rate_opts,
                       # col_sample_rate_per_tree = col_sample_rate_per_tree_opts
                       #,nbins_cats = nbins_cats_opts
  )
  
  search_criteria = list(strategy = "RandomDiscrete",
                         max_runtime_secs = 120,
                         max_models = 100,
                         stopping_metric = "AUTO",
                         stopping_tolerance = 0.00001,
                         stopping_rounds = 5,
                         seed = 123456)
  
  
  data_h20_kbsplots <- as.h2o(x = data_in %>% filter(!is.na(get(response))))
  
  kbsdata_drf_summary <-  h2o.randomForest(x = predictors,#,"Treatment_Category"
                                           y = response ,ntrees = 200,stopping_rounds = 20,
                                           seed = 10000, score_each_iteration = T,
                                           balance_classes = T,
                                           max_depth = max_depth, mtries = mtries,fold_column = fold_column,
                                           keep_cross_validation_predictions= TRUE,
                                           training_frame = data_h20_kbsplots)
  
  
  gbm_grid <- h2o.grid("randomForest",
                       grid_id = "mygrid_ENS",
                       x = predictors,
                       y = response,
                       # faster to use a 80/20 split
                       training_frame = data_h20_kbsplots,
                       # validation_frame = validSplit,
                       fold_column = fold_column,
                       # alternatively, use N-fold cross-validation:
                       # training_frame = train,
                       # nfolds = 5,
                       # Gaussian is best for MSE loss, but can try 
                       # other distributions ("laplace", "quantile"):
                       distribution="gaussian",
                       # stop as soon as mse doesn't improve by 
                       # more than 0.1% on the validation set, 
                       # for 2 consecutive scoring events:
                       stopping_rounds = 10,
                       stopping_tolerance = 1e-4,
                       stopping_metric = "MSE",
                       # how often to score (affects early stopping):
                       score_tree_interval = 100,
                       ## seed to control the sampling of the 
                       ## Cartesian hyper-parameter space:
                       seed = 123456,
                       hyper_params = hyper_params,
                       search_criteria = search_criteria)
  gbm_sorted_grid2 <- h2o.getGrid(grid_id = "mygrid_ENS", sort_by = "mse")
  print(gbm_sorted_grid2)
  best_model <- h2o.getModel(gbm_sorted_grid@model_ids[[1]])
  summary(best_model)
  
  
  
}
