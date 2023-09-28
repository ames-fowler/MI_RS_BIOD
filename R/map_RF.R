cross_validate = function(data_in = kbs_summary_data, 
                          predictors = predictors_all,
                          response =response,
                          model_id_in = paste(response, scope_label),
                          label,
                          scope_label,
                          max_depth= 4,
                          fold_column = "block",
                          mtries = 16){
  
  gc()
  
  data_h20_kbsplots <- as.h2o(x = data_in %>% filter(!is.na(get(response))))
  
  kbsdata_drf_summary <-  h2o.randomForest(x = predictors,#,"Treatment_Category"
                                           y = response,
                                           ntrees = 900,
                                           min_rows = 100,
                                           sample_rate = .95,  
                                           stopping_rounds = 20,
                                           seed = 10000, 
                                           score_each_iteration = T,
                                           balance_classes = T,
                                           max_depth = max_depth, 
                                           mtries = mtries,
                                           model_id = model_id_in,
                                           fold_column = fold_column,
                                           keep_cross_validation_predictions= TRUE,
                                           training_frame = data_h20_kbsplots)
  
  path <- paste0(getwd(),"./processed_data/rf_models/",label)
  if(!dir.exists(path)){dir.create(path, recursive = T)}
  mojo_destination <- h2o.saveModel(kbsdata_drf_summary, path = path,force = T,
                                    export_cross_validation_predictions = T)
  
  ## check KBS fit
  # kbs_plants = kbsdata_drf_summary
  pout_temp = out <- h2o.predict(kbsdata_drf_summary, newdata =  data_h20_kbsplots)
  names(pout_temp) = "prediction"
  pout = cbind( as.data.frame(data_h20_kbsplots), 
                as.data.frame(pout_temp)) %>% mutate(response = .[,response])
  
  cvpreds <- h2o.getFrame(kbsdata_drf_summary@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  
  Nfold_data  = data_h20_kbsplots %>% #data_in_split[[1]] %>% 
    as.data.frame() %>% mutate(response = .[,response])%>% 
    cbind(as.data.frame(cvpreds)) %>% dplyr::rename(prediction = predict)
  Nfold_data$val = "Cross validation \nprediction"
  pout$val = "Full model"
  mod_compare = rbind(pout[,c("Treatment_Category","response","prediction","val")], 
                      Nfold_data[,c("Treatment_Category","response","prediction","val")])
  
  
  write_rds(mod_compare, paste0("./figures/crossvalidation/paper_model_230726/kbs_model_",label,"_",response,".rds"))
  g = ggplot( mod_compare %>% filter(!is.na(prediction)),
             aes(x = response, y = prediction))+
    geom_point(aes(col = val),alpha =.1 )+ #, aes(col =val), shape = Treatment_Category
    geom_abline(slope = 1)+
    theme_bw()+coord_equal()+
    scale_fill_gradient(low = "grey90", high = "grey15")+
    scale_color_manual(values = c( "darkgreen", "darkblue"))+
    labs(x = paste("Obs.", 
                   (response %>% gsub(x = .,pattern = "_", replacement = " "))), 
         y = paste("Predicted", 
                   (response %>% gsub(x = .,pattern = "_", replacement = " "))),
         col = "Model")#+
    # facet_wrap(~(val %>% factor(levels = c("Full model","Cross validation \nprediction"))))
  
  ggsave(paste0("./figures/crossvalidation/paper_model_230726/kbs_model_",label,"_",response,".jpeg"), 
         plot = g, height = 4, width = 6)
  
  
  g = ggplot(Nfold_data %>% filter(!is.na(prediction)),
             aes(x = response, y = prediction))+
    geom_point(aes(),alpha =.1 )+ #, aes(col =val), shape = Treatment_Category
    geom_abline(slope = 1)+
    theme_bw()+coord_equal()+
    # scale_fill_gradient(low = "grey90", high = "grey15")+
    # scale_color_manual(values = c( "darkgreen", "darkblue"))+
    labs(x = paste("Obs.", 
                   (response %>% gsub(x = .,pattern = "_", replacement = " "))), 
         y = paste("Predicted", 
                   (response %>% gsub(x = .,pattern = "_", replacement = " "))),
         col = "Model")#+
  # facet_wrap(~(val %>% factor(levels = c("Full model","Cross validation \nprediction"))))
  
  ggsave(paste0("./figures/crossvalidation/paper_model_230702/kbs_crossvalidation_",label,"_",response,".jpeg"), 
         plot = g, height = 4, width = 4)
  return(list(Nfold_data, kbsdata_drf_summary))
}



get_h2o_rf <- function(response,predictors,dataset,depth , train, valid){
  data_drf <-  h2o.randomForest(x = predictors,
                                y = response ,ntrees = 200,stopping_rounds = 20,
                                seed = 10000, score_each_iteration = T,
                                balance_classes = T,
                                max_depth = depth,
                                training_frame = train,
                                validation_frame = valid,
                                model_id = paste0("RF_", response, 
                                                  "_",dataset,"_", 
                                                  gsub(pattern = " |-|:",
                                                       replacement = "_",Sys.time())))
  path <- paste0(getwd(),"./processed_data/rf_models")
  mojo_destination <- h2o.save_mojo(data_drf, path = path)
  
  return(data_drf)
}


map_h2o_rf_CV <- function(data_frame,predictors, response_list,depth = 4,mtries = 16,label = "paper_230822",
                          scope_label = "BCSE"){
  h2o.init(nthreads = -1)
  print("change flag")
  # data_h20_plots <- as.h2o(x = data_frame)
  # train <- data_in_split[[1]]
  # valid <- data_in_split[[2]]
  # test <- data_in_split[[3]]
  
  models <- purrr::map(response_list, function(H) cross_validate(data_in = data_frame,
                                                            predictors =predictors,
                                                            response = H,
                                                            label =label,
                                                            scope_label = scope_label,
                                                            max_depth = depth,
                                                            fold_column = "block",
                                                            mtries = mtries))
  
  # perf_all <- purrr::map(models, h2o.performance())
  # r2_model <- purrr::map(models, h2o.r2)
  # imp <- purrr::map2(models,response, function(h,g) {
  #   temp <- h2o.varimp(h)
  #   temp$rank <- rev(rank(temp$percentage))
  #   temp$name <- g
  #   return(temp %>% as.data.frame())
  # }) %>% rbindlist()
  # 
  # performance_report <-  list(perf_all, r2_model, imp)
  # save(performance_report, "./processed_data/rf_models")
  
  
  return(models)
}

model_plot = function(label, data_in = test22, Treatment_Category_in, response_in){
  
  
  data_frame = data_in %>% filter(Treatment_Category == Treatment_Category_in)
  
  data_h20_kbsplots = as.h2o(x = data_frame %>% filter(!is.na(get(response_in))))
  
  kbsdata_drf_summary <-  h2o.randomForest(x = predictors,#,"Treatment_Category"
                                           y = response_in,
                                           ntrees = 900,
                                           min_rows = 100,
                                           sample_rate = .95,  
                                           stopping_rounds = 20,
                                           seed = 10000, 
                                           score_each_iteration = T,
                                           balance_classes = T,
                                           max_depth = max_depth, 
                                           mtries = mtries,
                                           model_id = label,
                                           fold_column = fold_column,
                                           keep_cross_validation_predictions= TRUE,
                                           training_frame = data_h20_kbsplots)
  
  path <- paste0(getwd(),"./processed_data/rf_models/paper_230822")
  mojo_destination <- h2o.saveModel(kbsdata_drf_summary, path = path,force = T,
                                    export_cross_validation_predictions = T)

  cvpreds <- h2o.getFrame(kbsdata_drf_summary@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  
  Nfold_data  =  data_h20_kbsplots %>% #data_in_split[[1]] %>% 
    as.data.frame() %>% mutate(response = .[,response_in])%>% 
    cbind(as.data.frame(cvpreds)) %>% dplyr::rename(prediction = predict)
  
  
  g = ggplot(Nfold_data %>% filter(!is.na(prediction)),
             aes(x = response, y = prediction))+
    geom_point(aes(),alpha =.1 )+ #, aes(col =val), shape = Treatment_Category
    geom_abline(slope = 1)+
    theme_bw()+coord_equal()+
    # scale_fill_gradient(low = "grey90", high = "grey15")+
    # scale_color_manual(values = c( "darkgreen", "darkblue"))+
    labs(x = paste("Obs.", 
                   (response_in %>% gsub(x = .,pattern = "_", replacement = " "))), 
         y = paste("Predicted", 
                   (response_in %>% gsub(x = .,pattern = "_", replacement = " "))),
         col = "Model")#+
  # facet_wrap(~(val %>% factor(levels = c("Full model","Cross validation \nprediction"))))
  
  ggsave(paste0("./figures/crossvalidation/paper_model_230702/kbs_crossvalidation_",label,"_",response_in,".jpeg"), 
         plot = g, height = 4, width = 4)
  return(list(Nfold_data, kbsdata_drf_summary))
}


## function for partials.  
plot_pd <- function(model = plant_scaleUP_model_2_rf_planet_10_top10[c(1)],
                    data = as.h2o(scaleUP_data),
                    name = "plants_top10" 
                    ){
  testit <- h2o.explain(object = model, 
                        newdata = data,
                        top_n_features = 10)
  
  figure = ggarrange(testit[[5]]$plots[[1]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[2]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[3]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[4]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[5]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[6]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[7]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[8]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[9]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[10]]+ rremove("ylab")+labs(title = ""), 
                     common.legend = T, 
                     legend = "bottom",
                     nrow = 2, ncol = 5, 
                     align = "hv", labels="auto")
  
  figure_out <- annotate_figure(figure, 
                                left = textGrob("Mean Response\n", rot = 90, 
                                                vjust = 1, gp = gpar(cex = 1)))
  
  ggsave(plot = figure_out, filename = paste0("./figures/crossvalidation/paper_model_230726/partials",name,"_pd.jpeg"), 
         height = 4, width = 8)
  
}

## function for states and table.
# df =model_out$richness_Plants[[1]]
# lab_in =  "plant richness"
# stats_in = all_stats %>% filter(Response=="richness_Plants")

RF_pairs_plot= function(df, lab_in, stats_in, xshift = .75, yshit = .12){
  temp = ggplot(df%>% filter(!is.na(prediction)),
               aes(x = response, y = prediction))+
    geom_point(alpha =.1 )+ #, aes(col =val), shape = Treatment_Category
    geom_smooth(se=F, method = "lm", col = "blue", alpha = .5)+
    geom_abline(slope = 1)+
    # geom_text(x = 10, y = 30, label = lm_eqn(model_out$richness_Plants[[1]],x = response, y = prediction), parse = TRUE)
    theme_bw()+coord_equal()+
    scale_fill_gradient(low = "grey90", high = "grey15")+
    labs(x = paste("Observed", lab_in), 
         y = paste("Predicted", lab_in)) #+ theme(element_text(size = 12))
  
  eqn <- bquote(~~ R^2 == .(stats_in$R2) * "," ~~ RMSE == .(stats_in$RMSE))
  
  temp_annotate1 = annotate("text", 
                           x = max(df$response)*xshift, y = min(df$prediction)+max(df$response)*yshit, 
                           label= eqn)
  temp_out = temp+temp_annotate1 +theme(text = element_text(size=12))
  
}
response_stats= function(model_in_list= model_out[[3]][[1]], response){
  MAPE_temp = ((mean(abs(model_in_list$response-model_in_list$prediction)))/
                 mean(model_in_list$response)*100 )%>% round(digits = 3)
  RRMSE_temp = (sqrt(mean((model_in_list$response-model_in_list$prediction)^2))/
                  mean(model_in_list$response)*100 )%>% round(digits = 3)
  RMSE_temp = (sqrt(mean((model_in_list$response-model_in_list$prediction)^2)))%>%
    round(digits = 2)
  R2_temp = cor(model_in_list$response,model_in_list$prediction)^2 %>%
    round(digits = 2)
  
  repsonce_preformance = list("Response"= response,
                              "MAPE" = MAPE_temp, 
                              "RRMSE" = RRMSE_temp,
                              "RMSE" = RMSE_temp,
                              "R2" = R2_temp) %>% as.data.frame()
  
}

get_scale_prediction= function(model_in, data_in){
  scale_temp = h2o.predict(object = model_in, as.h2o(data_in))
  data_in$prediction = as.vector(scale_temp)
  data_in$response = data_in[[(model_in@model$names[length(model_in@model$names)])]]
  return(data_in)
}


response_stats_scale= function(model_in_list= model_out,response_name){
  MAPE_temp = ((mean(abs(model_in_list$response-model_in_list$prediction)))/
                 mean(model_in_list$response)*100 )%>% round(digits = 3)
  RRMSE_temp = (sqrt(mean((model_in_list$response-model_in_list$prediction)^2))/
                  mean(model_in_list$response)*100 )%>% round(digits = 3)
  RMSE_temp = (sqrt(mean((model_in_list$response-model_in_list$prediction)^2)))%>%
    round(digits = 2)
  R2_temp = cor(model_in_list$response,model_in_list$prediction)^2 %>%
    round(digits = 2)
  
  repsonce_preformance = list("Response"= response_name,
                              "MAPE" = MAPE_temp, 
                              "RRMSE" = RRMSE_temp,
                              "RMSE" = RMSE_temp,
                              "R2" = R2_temp) %>% as.data.frame()
  
}
