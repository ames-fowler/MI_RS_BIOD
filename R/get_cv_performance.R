# get comn CV stats
# models = model_out %>% sapply(.,"[",2)
# models = plant_rich_by_treatment %>% sapply(.,"[",2)

# names(models) = response_list
# names(models) = TC

H = models$richness_Bee_Groups
get_cv_states = function(models){
  
  temp = lapply(models, function(H){
    out = H@model$cross_validation_metrics_summary$mean%>% as_tibble()  %>% t()%>% as_tibble()
    names(out) = c("mae","mean_residual_deviance","mse", "r2", "residual_deviance", "rmse", "rmsle")
    return(out)
    
    })
  
  temp %>% data.table::rbindlist() %>% mutate(model = names(models))
  
}
