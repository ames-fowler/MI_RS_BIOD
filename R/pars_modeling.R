# parsimonious model fitting 


parsmodel_defined = function(model_in, intlist, predictors_in = predictors, 
         pars_label,response_in = response, data_in = test22,
         fold_column= "block", max_depth = 4, mtries = 16){
  temp_model = h2o.loadModel(model_in)
  data_h20_kbsplots= as.h2o(x = data_in %>% filter(!is.na(get(response_in))))
  i= intlist[1]
  j=1
  samplelist=list()
  for(i in intlist){
    print(names(intlist)[j])
    if(length(unlist(i))==1&all(is.na(as.numeric(unlist(i))))){
    predictors_temp = predictors_in[grepl(predictors_in, pattern = paste0("^",i))]
    }else{predictors_temp = predictors_in[grepl(predictors_in, pattern = paste(unlist(i), collapse = "|"))]
    }
    kbsdata_drf_summary <-  h2o.randomForest(x = predictors_temp,#,"Treatment_Category"
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
                                             model_id = paste0("richness_Plants_",names(intlist)[j]),
                                             fold_column = fold_column,
                                             keep_cross_validation_predictions= TRUE,
                                             training_frame = data_h20_kbsplots)
    
    path <- paste0(getwd(),"./processed_data/rf_models/paper_230726/parsimonious")
    mojo_destination <- h2o.saveModel(kbsdata_drf_summary, path = path,force = T,
                                      export_cross_validation_predictions = T)
    
    
    cvpreds <- h2o.getFrame(kbsdata_drf_summary@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
    
    Nfold_data  =  data_h20_kbsplots %>% #data_in_split[[1]] %>% 
      as.data.frame() %>% mutate(response = .[,response])%>% 
      cbind(as.data.frame(cvpreds)) %>% 
      dplyr::rename_at(.vars = vars(predict), .funs = ~paste0(response,"_prediction"))%>%
      dplyr::select(c("x","y","Treatment_Category", 
                            "trtmt_code", "block", "year",
                            eval(response),paste0(response,"_prediction")))
    
    write_csv(x = Nfold_data, file = paste0("./processed_data/Paper_out/parimonious/richness_Plants_",
                                            names(intlist)[j],".csv"))
    samplelist = c(samplelist,kbsdata_drf_summary)
    j =j+1
  }
  temp = lapply(c(temp_model, samplelist), function(H){
    out = H@model$cross_validation_metrics_summary$mean%>% as_tibble()  %>% t()%>% as_tibble()
    names(out) = c("mae","mean_residual_deviance","mse", "r2", "residual_deviance", "rmse", "rmsle")
    return(out)
    
  })
  
  temp_out = temp %>% data.table::rbindlist() %>% mutate(model =  c("All", names(intlist)))
  write_csv(x = temp_out, file = paste0( paste0("./processed_data/Paper_out/parimonious/richness_Plants_",pars_label,".csv")))

  return(temp_out)
}


parsmodel_iterative = function(model_in, intlist, predictors_in = predictors, 
                               pars_label,response_in = response, data_in = test22,
                               fold_column= "block", max_depth = 4, mtries = 16){
  temp_model= temp_model_in = h2o.loadModel(model_in)
  data_h20_kbsplots= as.h2o(x = data_in %>% filter(!is.na(get(response_in))))
  i= intlist[1]
  samplelist=list()
  for(i in intlist){
    
    predictors_temp = h2o.varimp(temp_model)$variable[1:unlist(i)]
    print((i))
    
    temp_model <-  h2o.randomForest(x = predictors_temp,#,"Treatment_Category"
                                    y = response_in,
                                    ntrees = 900,
                                    min_rows = 100,
                                    sample_rate = .95,  
                                    stopping_rounds = 20,
                                    seed = 10000, 
                                    score_each_iteration = T,
                                    balance_classes = T,
                                    max_depth = max_depth, 
                                    mtries = min(mtries, sqrt(length(predictors_temp)) %>% round(0)),
                                    model_id = paste0("richness_Plants_features_",i),
                                    fold_column = fold_column,
                                    keep_cross_validation_predictions= TRUE,
                                    training_frame = data_h20_kbsplots)
    
    path <- paste0(getwd(),"./processed_data/rf_models/paper_230726/parsimonious")
    mojo_destination <- h2o.saveModel(temp_model, path = path,force = T,
                                      export_cross_validation_predictions = T)
    
    
    cvpreds <- h2o.getFrame(temp_model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
    
    Nfold_data  =  data_h20_kbsplots %>% #data_in_split[[1]] %>% 
      as.data.frame() %>% mutate(response = .[,response])%>% 
      cbind(as.data.frame(cvpreds)) %>% 
      dplyr::rename_at(.vars = vars(predict), .funs = ~paste0(response,"_prediction"))%>%
      dplyr::select(c("x","y","Treatment_Category", 
                      "trtmt_code", "block", "year",
                      eval(response),paste0(response,"_prediction")))
    
    write_csv(x = Nfold_data, file = paste0("./processed_data/Paper_out/parimonious/richness_Plants_features_",
                                            (i),".csv"))
    samplelist = c(samplelist,temp_model)
    
  }
  temp = lapply(c(temp_model_in, samplelist), function(H){
    out = H@model$cross_validation_metrics_summary$mean%>% as_tibble()  %>% t()%>% as_tibble()
    names(out) = c("mae","mean_residual_deviance","mse", "r2", "residual_deviance", "rmse", "rmsle")
    return(out)
    
  })
  
  temp_out = temp %>% data.table::rbindlist() %>% mutate(model =  c("All", paste("Features =", intlist)))
  write_csv(x = temp_out, file = paste0( paste0("./processed_data/Paper_out/parimonious/richness_Plants_",pars_label,".csv")))
  
  return(temp_out)
}


plot_parse= function(pars_performace, names_in, fig_dir, pars_label){

  pars_performace = pars_performace %>% mutate(model = names_in)
  
  g = 
    ggplot(pars_performace)+ 
      geom_point(aes(x = model %>% 
                       factor(., levels=unique(model[order(r2,decreasing = T)]), ordered=TRUE),
                     y = r2))+
      # geom_line(aes(x = model, y = r2))+
      theme_bw()+ labs(y = ("R-squared"), x = "", col = "", title = paste(pars_label ,"models"))
  
  g2 = 
    ggplot(pars_performace)+ 
    geom_point(aes(x = model %>% factor(., levels=unique(model[order(r2,decreasing = T)]), ordered=TRUE), 
                   y = rmse ))+
    # geom_line(aes(x = model, y = r2))+
    theme_bw()+ labs(y = ("RMSE"), x = "")
  
  gout = ggpubr::ggarrange(g,g2, common.legend = T, align = "hv")
  
}


