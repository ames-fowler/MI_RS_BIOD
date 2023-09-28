#TD KBS data #### 
  #see site_2_test_2023.R and kbs_all_dates and thermal_time_scratch and targets for DEM  


predictors = predictor_columns[predictor_columns %in% names(lalal)]
# Run model main ###
source("./R/papre functions/map_RF.R")
data_frame = test22
predictors = (test22 %>% names)[c(7:9,35:984)]#982]
response_list = names(test22)[c(10:19, 30:34, 985:992)]
model_out = map_h2o_rf_CV(data_frame = data_frame,
                          predictors = predictors,
                          response_list = response_list)
# model_out_plants = map_h2o_rf_CV(data_frame = data_frame,
#                           predictors = predictors,
#                           response_list = response_list[5])

names(model_out) = response_list

#save out all results? or done in function 
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
# find RRMSE and R2 of cross validated #####
all_stats = map2(model_out, names(model_out), function(G,H){
  temp_stats = response_stats(model_in_list = G[[1]], response = H)
}) %>% rbindlist()

# Model and included model treatment category in the PlantGrowth #########

model_out_treat_annual = cross_validate(data_in = test22 %>% 
                                         mutate(Treatment_Category = Treatment_Category %>% factor())%>%
                                         filter(Treatment_Category== "Annual Monoculture"),
                          predictors = predictors,
                          response = "richness_Plants", 
                          model_id_in = "richness_Plants_AM",
                          label = "paper_072623",
                          scope_label = "BCSE",
                          max_depth = 4,
                          fold_column = "block",
                          mtries = 16)

model_out_treat_low_div = cross_validate(data_in = test22 %>% 
                                         mutate(Treatment_Category = Treatment_Category %>% factor())%>%
                                         filter(Treatment_Category== "Low-diversity Perennial Grass"),
                                         predictors = predictors,
                                         response = "richness_Plants", 
                                         model_id_in = "richness_Plants_LD",
                                         label = "paper_072623",
                                         scope_label = "BCSE",
                                         max_depth = 4,
                                         fold_column = "block",
                                         mtries = 16)

model_out_treat_high_div = cross_validate(data_in = test22 %>% 
                                           mutate(Treatment_Category = Treatment_Category %>% factor())%>%
                                           filter(Treatment_Category== "High-diversity Perennial Polyculture"),
                                          predictors = predictors,
                                          response = "richness_Plants", 
                                          model_id_in = "richness_Plants_HD",
                                          label = "paper_072623",
                                          scope_label = "BCSE",
                                          max_depth = 4,
                                          fold_column = "block",
                                          mtries = 16)
temp_stats_AM = response_stats(model_in_list = model_out_treat_annual[[1]],
                               response = "Plant richness - Annual monoculture")
temp_stats_LD = response_stats(model_in_list = model_out_treat_low_div[[1]],
                               response = "Plant richness - Low-diversity Perennial")
temp_stats_HD = response_stats(model_in_list = model_out_treat_high_div[[1]],
                               response = "Plant richness - High-diversity Perennial")
temp_stats_groups = rbind(temp_stats_AM, temp_stats_LD,temp_stats_HD)

#top ten variables 
AM_top_ten_Vars = h2o.varimp(model_out_treat_annual[[2]]) %>% head(10)
LD_top_ten_Vars = h2o.varimp(model_out_treat_low_div[[2]]) %>% head(10)
HD_top_ten_Vars = h2o.varimp(model_out_treat_high_div[[2]]) %>% head(10)
important_plants = cbind(AM_top_ten_Vars$variable, LD_top_ten_Vars$variable, HD_top_ten_Vars$variable)
###parse modeling ####

model_in= "./processed_data/rf_models/paper_230810//richness_Plants"

richness_model = h2o.loadModel(model_in)#"./processed_data/rf_models/paper_230726/richness_Plants")
ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230810/ENS_Plants")


# data_frame = test22 %>% mutate(Treatment_Category = Treatment_Category %>% factor())
# predictors = predictors
# response = ("richness_Plants")#names(test22)[c(7:16, 27:30, 981:990)]
# data_frame = data_frame %>% filter(Treatment_Category== "Annual Monoculture")
# 

# Temporal 

(h2o.varimp(model_out$richness_Plants[[2]]) %>% head(75))$variable %>% strsplit("_") %>% sapply(., tail, 1) %>% table
(h2o.varimp(model_out$richness_Plants[[2]]) %>% head(75))$percentage %>% sum
temporal_list = list("three_periods" = c(133 ,223 ,343  ),"Spring_Fall"= c(133 ,343 ),
                     "Fall" = c(343), "Spring" = c(133))

# Spatial 
predictors_no_CV =predictors[!grepl(predictors, pattern="cv")]
predictors_1x=predictors_no_CV[!grepl(predictors_no_CV, pattern="mean|DEM")]
predictors_3x=predictors[grepl(predictors, pattern="3mean")]
predictors_5x=predictors[grepl(predictors, pattern="5mean")]
spatial_list = list("predictors_1x" = predictors_1x, 
                    "predictors_3x" = predictors_3x, 
                    "predictors_5x" = predictors_5x)


# Spectral 
bands = predictors_1x %>% strsplit("_") %>% sapply(.,"[",1) %>% unique()
VNIR = bands[1:4]
VI_list = split(bands[5:10],1:6)
spectral_list = list(VNIR, VI_list)
# Minimal feature importance 
list_mostimp = list(75,35,16,8,4,2,1 )

pars_spectral = parsmodel_defined(model_in = model_in, intlist = VI_list,
                                  pars_label = "spectral")

pars_temporal = parsmodel_defined(model_in = model_in, intlist = temporal_list,
                                  pars_label = "temporal")

pars_spatial = parsmodel_defined(model_in = model_in, intlist = spatial_list,
                                 pars_label = "spatial")



pars_importantce = parsmodel_iterative(model_in = model_in, intlist = list_mostimp, pars_label = "importance")

testit = h2o.loadModel("./processed_data/rf_models/paper_230822/parsimonious/richness_Plants_features_2")
h2o.varimp(testit)

pars_performace = pars_importantce

pars_models = list(pars_spectral, pars_temporal, pars_spatial, pars_importantce)
pars_model_names = list(
  c("All", unlist(VI_list)),
  c("All DOY", "3 periods", "Spring,Fall", "Fall","Spring"),
  c("All res", "1x res", "3x res", "5x res"),
  c("All", "75", "35", "16", "8", "4", "2", "1")
)
pars_model_label = list("Spectral", "Temporal", "Spatial", "Importance")
pars_plot_out = purrr::pmap(list(pars_models,pars_model_names, pars_model_label), 
                            function(a,b,c){plot_parse(pars_performace = a,
                                                       names_in = b, 
                                                       fig_dir = NA,
                                                       pars_label = c)})

gout_all = ggpubr::ggarrange(pars_plot_out[[1]],
                             pars_plot_out[[2]],
                             pars_plot_out[[3]],
                             pars_plot_out[[4]],
                             ncol = 1)

ggsave(plot = gout_all, filename =  "./figures/crossvalidation/paper_model_230726/Pars_stats_plant_richness.jpeg",
       width = 7.5, height = 10)


lm_eqn <- function(df, x, y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == #a + b %.% italic(x)*","
                     ~~italic(r)^2~"="~r2, 
                   list(
                   #a = format(unname(coef(m)[1]), digits = 2),
                   #      b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
##### Plant Figure 3 -----
gp1 = RF_pairs_plot(df = model_out$richness_Plants[[1]]%>% filter(!is.na(prediction)),
                    lab_in = "plant richness",stats_in = all_stats %>% filter(Response=="richness_Plants")
                    )
gp2 = RF_pairs_plot(df = model_out$ENS_Plants[[1]]%>% filter(!is.na(prediction)),
                    lab_in = "plant Hill's number",stats_in = all_stats %>% filter(Response=="ENS_Plants")
)

gp_out = ggpubr::ggarrange(gp1, gp2, ncol = 2,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plant_model_predictions.jpeg", 
       plot = gp_out, width = 8.5, height = 4, dpi = 700)


##### Plant pollinator supp figure s2-5  -----
gpr1 = RF_pairs_plot(df = model_out$richness_Bee_Groups[[1]]%>% filter(!is.na(prediction)),
                    lab_in = "bee groups richness",
                    stats_in = all_stats %>% filter(Response=="richness_Bee_Groups")
)
gpr2 = RF_pairs_plot(df = model_out$richness_Bumblebees[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "bumblebees richness",
                     stats_in = all_stats %>% filter(Response=="richness_Bumblebees")
)
gpr3 = RF_pairs_plot(df = model_out$richness_Butterflies[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "butterflies richness",
                     stats_in = all_stats %>% filter(Response=="richness_Butterflies")
)

gpr_out = ggpubr::ggarrange(gpr1, gpr2, gpr3, ncol = 3,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plot_pollinator_richness_s3.jpeg", 
       plot = gpr_out, width =12.5, height = 3.5, dpi = 700)

#hills number
gph1 = RF_pairs_plot(df = model_out$ENS_Bee_Groups[[1]]%>% filter(!is.na(prediction)),
                    lab_in = "bee group Hill's number",stats_in = all_stats %>% 
                      filter(Response=="ENS_Bee_Groups")
)
gph2 = RF_pairs_plot(df = model_out$ENS_Bumblebees[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "bumblebees Hill's number",stats_in = all_stats %>% 
                      filter(Response=="ENS_Bumblebees")
)
gph3 = RF_pairs_plot(df = model_out$ENS_Butterflies[[1]]%>% filter(!is.na(prediction)),
                    lab_in = "butterflies Hill's number",stats_in = all_stats %>% 
                      filter(Response=="ENS_Butterflies")
)

gph_out = ggpubr::ggarrange(gph1, gph2, gph3, ncol = 3,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plot_pollinator_ens_s3.jpeg", 
       plot = gph_out, width =12.5, height = 3.5, dpi = 700)

# Richneass_residuals
gprr1 = RF_pairs_plot(df = model_out$richness_res_Bee_Groups[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "richness bee groups",stats_in = all_stats %>% 
                       filter(Response=="richness_res_Bee_Groups")
)
gprr2 = RF_pairs_plot(df = model_out$richness_res_Bumblebees[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "richness Bumblebees",stats_in = all_stats %>% 
                       filter(Response=="richness_res_Bumblebees")
)
gprr3 = RF_pairs_plot(df = model_out$richness_res_Butterflies[[1]]%>% filter(!is.na(prediction)),
                     lab_in = "richness Butterflies",stats_in = all_stats %>% 
                       filter(Response=="richness_res_Butterflies")
)

gprr1_out = ggpubr::ggarrange(gprr1, gprr2, gprr3, ncol = 1,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plot_pollinator_richness_res_s5.jpeg", 
       plot = gprr1_out, width =7, height = 14, dpi = 700)

# plant treatment residuals
gptr1 = RF_pairs_plot(df = model_out_treat_annual[[1]]%>% filter(!is.na(prediction)),
                      lab_in = "plant richness",stats_in = temp_stats_groups  %>% 
                        filter(Response=="Plant richness - Annual monoculture")
)
gptr2 = RF_pairs_plot(df = model_out_treat_low_div[[1]]%>% filter(!is.na(prediction)),
                      lab_in = "plant richness",stats_in = temp_stats_groups %>% 
                        filter(Response=="Plant richness - Low-diversity Perennial")
)
gptr3 = RF_pairs_plot(df = model_out_treat_high_div[[1]]%>% filter(!is.na(prediction)),
                      lab_in = "plant richness",stats_in = temp_stats_groups %>% 
                        filter(Response=="Plant richness - High-diversity Perennial")
)

gpptr_out = ggpubr::ggarrange(gptr1, gptr2, gptr3, ncol = 1,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plot_plant_treatmentgroup_richness_s5.jpeg", 
       plot = gpptr_out, width =7, height = 14, dpi = 700)




ENS_inplotvar = model_out$ENS_Plants[[1]] %>% group_by(trtmt_code, block) %>% summarize(min = min(prediction),
                                                                        max = max(prediction),
                                                                        range = max-min) %>% 
  group_by(trtmt_code) %>% summarise(range = mean(range))
mean(ENS_inplotvar$range)
richness_inplotvar = model_out$richness_Plants[[1]] %>% group_by(trtmt_code, block) %>% summarize(min = min(prediction),
                                                                                        max = max(prediction),
                                                                                        range = max-min) %>% 
  group_by(trtmt_code) %>% summarise(range = mean(range))
mean(richness_inplotvar$range)

##linear fixed effect model  ##########

lm_plant_richness = lm(richness_Plants ~Treatment_Category, data = test22)
summary(lm_plant_richness)

lm_plant_ENS = lm(ENS_Plants ~Treatment_Category, data = test22)
summary(lm_plant_ENS)

### growth curve and figure 4 ######## 

testitall = test2 %>% 
  group_by(Treatment_Category, trtmt_code, block, bin, year) %>% 
  summarise_all(~mean(.x, na.rm=T)) 

testitall2 = testitall%>% 
  group_by(Treatment_Category, trtmt_code, block,year)%>% 
  mutate(cumsum_GNDVI = cumsum(GNDVI))%>% 
  group_by( bin, year) %>%
  mutate(mean_cumsum_GNDVI = mean(cumsum_GNDVI),
         mean_GNDVI = mean(GNDVI),
         sd_GNDVI = sd(GNDVI))%>% 
  group_by(Treatment_Category, trtmt_code, block,year)%>%
  mutate(norm_cumsum_GNDVI = cumsum_GNDVI/mean_cumsum_GNDVI,
         norm_GNDVI = (GNDVI-mean_GNDVI)/sd_GNDVI)

test5 = testitall %>% group_by(x,y, year) %>% 
  mutate(max_NDVI = max(NDVI, na.rm = T),
         TT_atMax = TT[which.max(NDVI)])

g = ggplot(data = testitall2 %>% filter(!is.na(richness_Plants)))+
  geom_line(aes(x = doy, y= GNDVI,
                col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                  factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                # size = ENS_Plants %>% cut(4) %>% factor,
                group = interaction(x,y)))+# ,alpha = .5
  facet_grid(cols = vars(year))+theme_minimal()+
  labs(x = "Day of the year", col = "Plant richness")+
  scale_color_viridis_d(direction = -1) +
  geom_rect(aes(xmin = 126, xmax = 141, ymin = .35, ymax = .85),
            fill = "transparent",color = "purple",linewidth =1)+
  geom_vline(xintercept = 133, col= "purple", alpha = .5, 
             aes(ymin = .35, ymax = .85))

g

#suplimant by land use 
gs1 = ggplot(data = testitall2 %>% filter(!is.na(richness_Plants)))+
  geom_line(aes(x = doy, y= GNDVI,
                col = interaction(Treatment_Category,trtmt_code%>% 
                  factor(labels = c("Corn", "Prairie", "Sorghum", "Sorghum + covercrop",
                                    "Transitional Switchgrass", "Switchgrass",
                                    "Miscanthus", "Native Grasses", "Poplar" ,
                                    "Successional"))),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                # size = ENS_Plants %>% cut(4) %>% factor,
                linetype = Treatment_Category , group = interaction(x,y)))+# ,alpha = .5
  facet_grid(cols = vars(year))+theme_minimal()+
  labs(x = "Day of the year", col = "Plot Treatment")#+
  # scale_color_viridis_d(direction = -1) +
  # geom_rect(aes(xmin = 126, xmax = 141, ymin = .35, ymax = .85),
  #           fill = "transparent",color = "purple",linewidth =1)+
  # geom_vline(xintercept = 133, col= "purple", alpha = .5, 
  #            aes(ymin = .35, ymax = .85))

gs1

gs2 = ggplot(data = testitall2 %>% filter(!is.na(richness_Plants)) %>% 
               group_by(year, trtmt_code, block) %>% 
               mutate(GNDVI_cumsum_treatment = cumsum(GNDVI)) %>% 
               group_by(doy) %>% 
               mutate(GNDVI_min_cumsum = min(GNDVI_cumsum_treatment), 
                       GNDVI_max_cumsum = max(GNDVI_cumsum_treatment),
                       normalized = (GNDVI_cumsum_treatment-GNDVI_min_cumsum)/
                         (GNDVI_max_cumsum-GNDVI_min_cumsum)))+
  geom_line(aes(x = doy, y= normalized,
                 col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                   factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                # size = ENS_Plants %>% cut(4) %>% factor,
                group = interaction(x,y)))+# ,alpha = .5
  facet_grid(cols = vars(year))+theme_minimal()+
  labs(x = "Day of the year", col = "Plot Treatment")+
  scale_color_viridis_d(direction = -1) 


gs2


g2 = ggplot(data = test2 %>% filter(doy ==133))+
  geom_point( 
             aes(y = richness_Plants, x= red_3mean,
                 col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                   factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                 # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                 # size = ENS_Plants %>% cut(4) %>% factor,
                 group = interaction(x,y)), alpha = .15, size = .5)+# alpha = .15
  geom_point(data = testitall %>% filter(doy ==133), 
             aes(y = richness_Plants, x= red_3mean,
                 col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                   factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                 # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                 # size = ENS_Plants %>% cut(4) %>% factor,
                 group = interaction(x,y)), size = 2)+
  facet_grid(cols = vars(year))+theme_minimal()+
  labs(y = "Plant richness", col = "Plant richness",  x = "Red band reflectance at 3x resolution")+
  scale_color_viridis_d(direction = -1)+
  theme(panel.border = element_rect(color = "purple",
                                    fill = NA,
                                    linewidth =  1))

g2


g3 = ggplot()+
  geom_point(data = test2 %>% filter(doy ==133), 
             aes(y = richness_Plants, x= ARVI_3mean,
                 col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                   factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                 # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                 # size = ENS_Plants %>% cut(4) %>% factor,
                 group = interaction(x,y)), alpha = .15, size = .5)+# alpha = .15
  geom_point(data = testitall %>% filter(doy ==133), 
             aes(y = richness_Plants, x= ARVI_3mean,
                 col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                   factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                 # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                 # size = ENS_Plants %>% cut(4) %>% factor,
                 group = interaction(x,y)), size = 2,)+
  facet_grid(cols = vars(year))+theme_minimal()+
  labs(y = "Plant richness", col = "Plant richness", x = "ARVI at 3x resolution")+
  scale_color_viridis_d(direction = -1)+
  theme(panel.border = element_rect(color = "purple",
                                    fill = NA,
                                    linewidth =  1))

g3
#, rows = vars(Treatment_Category))
g4 = ggpubr::ggarrange(g, g2,g3, nrow=3, align = "hv")
ggsave(plot = g4, filename = "./figures/crossvalidation/paper_model_230726/important_features_133.jpeg", 
       width =9, height = 7)



#### partials ########

eight_feature_plant_richness_rf = h2o.loadModel("./processed_data/rf_models/paper_230703/parsimonious/richness_Plants_features_8")
h2o.varimp(testit)

testit <- h2o.explain(object = eight_feature_plant_richness_rf, 
                        newdata = as.h2o(test22),
                        top_n_features = 8)
  
figure = ggpubr::ggarrange(testit[[5]]$plots[[1]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[2]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[3]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[4]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[5]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[6]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[7]]+ rremove("ylab")+labs(title = ""),
                     testit[[5]]$plots[[8]]+ rremove("ylab")+labs(title = ""),
                     common.legend = T, 
                     legend = "bottom",
                     nrow = 2, ncol = 4, 
                     align = "hv", labels="auto")
  
  figure_out <- annotate_figure(figure, 
                                left = textGrob("Mean Response\n", rot = 90, 
                                                vjust = 1, gp = gpar(cex = 1)))
  
  ggsave(plot = figure_out, filename = paste0("./figures/crossvalidation/paper_model_230726/partials_8_feature_plant_richness.jpeg"), 
         height = 4, width = 8)

  
  two_feature_plant_richness_rf = h2o.loadModel("./processed_data/rf_models/paper_230703/parsimonious/richness_Plants_features_2")
  
  testit <- h2o.explain(object = two_feature_plant_richness_rf, 
                        newdata = as.h2o(test22),
                        top_n_features = 8)
  
  figure = ggpubr::ggarrange(testit[[5]]$plots[[1]]+ rremove("ylab")+labs(title = ""),
                             testit[[5]]$plots[[2]]+ rremove("ylab")+labs(title = ""),
                            
                             common.legend = T, 
                             legend = "bottom",
                             nrow = 1, ncol = 2, 
                             align = "hv", labels="auto")
  
  figure_out <- annotate_figure(figure, 
                                left = textGrob("Mean Response\n", rot = 90, 
                                                vjust = 1, gp = gpar(cex = 1)))
  
  ggsave(plot = figure_out, filename = paste0("./figures/crossvalidation/paper_model_230726/partials_2_feature_plant_richness.jpeg"), 
         height = 4, width = 8)
  
  # scale up ####
  
  h2o.init(nthreads = -1)
  ## load_data ## 
  pathplants = "./processed_data/planet/SR/buf_site_9/clean_plants_scale/plot_data/plots_datawide_2023_08_09"
  plant_path_out = map(paste0("_",1:9), ~ pathplants %>% gsub(pattern = "_9", replacement = .x))
  plants_data_wide = map(plant_path_out, read_rds) %>% do.call(rbind,.)
  
  pathpollinators = "./processed_data/planet/SR/buf_site_9/clean_pollinator_scale/plot_data/plots_datawide_2023_08_09"
  Poll_path_out = map(paste0("_",1:9), ~ pathpollinators %>% gsub(pattern = "_9", replacement = .x))
  pollinator_data_wide = map(Poll_path_out, read_rds) %>% do.call(rbind,.)
  
  ## load models ##
  Plant_richness_model = h2o.loadModel("./processed_data/rf_models/paper_230810/richness_Plants")#"./processed_data/rf_models/paper_230726/richness_Plants")
  Plant_ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230810/ENS_Plants")
  Plant_forb_model = h2o.loadModel("./processed_data/rf_models/paper_230810/percent_Forb")
  
  
  ##Plant parsimonios models 
  Plant_richness_model_pars = h2o.loadModel("./processed_data/rf_models/paper_230726/parsimonious/richness_Plants_Spring")#"./processed_data/rf_models/paper_230726/richness_Plant  s")
  Plant_richness_model_par_fall = h2o.loadModel("./processed_data/rf_models/paper_230726/parsimonious/richness_Plants_Fall")
  Plant_richness_model_pars_GNDVI = h2o.loadModel("./processed_data/rf_models/paper_230726/parsimonious/richness_Plants_1")#"./processed_data/rf_models/paper_230726/richness_Plants")
  Plant_richness_model_pars_GNDVI = h2o.loadModel("./processed_data/rf_models/paper_230726/parsimonious/")#"./processed_data/rf_models/paper_230726/richness_Plants")
  
  Plant_ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230810/ENS_Plants")
  
  #list plant models 
  list_plant_models = list(Plant_richness_model, Plant_ENS_model, Plant_richness_model_pars, Plant_richness_model_par_fall,Plant_richness_model_pars_GNDVI)
  
  ### 
  Bumblebees_richness_model = h2o.loadModel("./processed_data/rf_models/paper_230810/richness_Bumblebees")#"./processed_data/rf_models/paper_230726/richness_Plants")
  Bumblebees_ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230810/ENS_Bumblebees")
  
  Butterflies_richness_model = h2o.loadModel("./processed_data/rf_models/paper_230810/richness_Butterflies")#"./processed_data/rf_models/paper_230726/richness_Plants")
  Butterflies_ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230810/ENS_Butterflies")
  
  Bee_Groups_richness_model = h2o.loadModel("./processed_data/rf_models/paper_230726/richness_Bee_Groups")#"./processed_data/rf_models/paper_230726/richness_Plants")
  Bee_Groups_ENS_model = h2o.loadModel("./processed_data/rf_models/paper_230726/ENS_Bee_Groups")
  
  
  list_pollinator_models = list(Bumblebees_richness_model, Bumblebees_ENS_model,Butterflies_richness_model, 
                                Butterflies_ENS_model, Bee_Groups_richness_model, Bee_Groups_ENS_model)

plants_data_wide2 =plants_data_wide



plant_scale_predictllist = map(list_plant_models, 
                               ~ get_scale_prediction(.x, data_in = plants_data_wide2))

names(plant_scale_predictllist)= map(list_plant_models, function(H) H@model_id) %>% 
  unlist()

pollinator_scale_predictllist = map(list_pollinator_models, 
                               ~ get_scale_prediction(.x, data_in = pollinator_data_wide))

names(pollinator_scale_predictllist)= map(list_pollinator_models, function(H) H@model_id) %>% 
  unlist()

pollinator_scale_predictllist



plant_scale_stats = map2(plant_scale_predictllist, 
                         names(plant_scale_predictllist), 
                         response_stats_scale) %>% 
  rbindlist()

pollinator_scale_stats = map2(pollinator_scale_predictllist,
                              names(pollinator_scale_predictllist), 
                              response_stats_scale) %>% 
  rbindlist()
lux_arbor = response_stats_scale(model_in_list = plant_scale_predictllist$richness_Plants%>% 
                             filter(location %>% grepl(pattern="ux")), 
                     response_name = "richness_Plants")

site2_stats = response_stats_scale(model_in_list = plant_scale_predictllist$richness_Plants%>% 
                                  filter(image_site==2), 
                                 response_name = "richness_Plants")

#figure XX plant scale predition #### 
gps1 = RF_pairs_plot(df = plant_scale_predictllist$richness_Plants%>% filter(!is.na(prediction)),
                    lab_in = "plant richness",stats_in = plant_scale_stats %>% 
                      filter(Response=="richness_Plants")
)
# +
#   geom_point(data = plant_scale_predictllist$richness_Plants %>% 
#                filter(location %>% grepl(pattern="ux")), 
#              aes(x= response, y = prediction, col = 'Scale-up sites'))+
#   geom_point(data = plant_scale_predictllist$richness_Plants %>% 
#                filter(image_site==2), 
#              aes(x= response, y = prediction, col = 'shared imagery'))



gps2 = RF_pairs_plot(df = plant_scale_predictllist$ENS_Plants %>% filter(!is.na(prediction)),
                    lab_in = "plant Hill's number",stats_in = plant_scale_stats %>%
                      filter(Response=="ENS_Plants")
)
gpsfall = RF_pairs_plot(df = plant_scale_predictllist$richness_Plants_1 %>% filter(!is.na(prediction)),
                     lab_in = "plant richness, GNDVI",stats_in = plant_scale_stats %>%
                       filter(Response=="richness_Plants_1")
)

gps_out = ggpubr::ggarrange(gps1, gps2, ncol = 2,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/plant_scale_model_predictions.jpeg", 
       plot = gps_out, width = 8.5, height = 4, dpi = 700)


##figure SXX pollinator scale prediction #### 
gsp1 = RF_pairs_plot(df = pollinator_scale_predictllist$richness_Bee_Groups%>% filter(!is.na(prediction)),
                     lab_in = "bee groups richness",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="richness_Bee_Groups"), xshift = .73, yshit = 0
)
gsp2 = RF_pairs_plot(df = pollinator_scale_predictllist$richness_Bumblebees%>% filter(!is.na(prediction)),
                     lab_in = "bumblebees richness",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="richness_Bumblebees"), xshift = .73, yshit = 0
)
gsp3 = RF_pairs_plot(df = pollinator_scale_predictllist$richness_Butterflies%>% filter(!is.na(prediction)),
                     lab_in = "butterflies richness",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="richness_Butterflies"), xshift = .73, yshit = 0
)

gsph1 = RF_pairs_plot(df = pollinator_scale_predictllist$ENS_Bee_Groups%>% filter(!is.na(prediction)),
                     lab_in = "bee groups Hill's numaber",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="ENS_Bee_Groups"), xshift = .73, yshit = 0
)
gsph2 = RF_pairs_plot(df = pollinator_scale_predictllist$ENS_Bumblebees%>% filter(!is.na(prediction)),
                     lab_in = "bumblebees Hill's numaber",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="ENS_Bumblebees"), xshift = .73, yshit = 0
)
gsph3 = RF_pairs_plot(df = pollinator_scale_predictllist$ENS_Butterflies%>% filter(!is.na(prediction)),
                     lab_in = "butterflies Hill's numaber",stats_in = pollinator_scale_stats %>% 
                       filter(Response=="ENS_Butterflies"), xshift = .73, yshit = 0
)

gsp_out = ggpubr::ggarrange(gsp1, gsph1,gsp2, gsph2, gsp3, gsph3,
                            nrow = 3, ncol = 2,  align = 'hv')
ggsave(filename = "./figures/crossvalidation/paper_model_230726/pollinator_scale_model_predictions.jpeg", 
       plot = gsp_out, width = 8.5, height = 13, dpi = 700)


### scale up model 

data_frame = plants_data_wide2
predictors = predictors
response_listscale = names(plants_data_wide2)[c(6:7)]
# model_out_scale = map_h2o_rf_CV(data_frame = plants_data_wide2 %>% mutate(block = location %>%  as.numeric()),
#                           predictors = predictors,
#                           response_list = response_listscale, scope_label = "field_eval")

  data_in <- as.h2o(x = plants_data_wide2%>% 
                           # filter(richness_Plants>10) %>%
                           mutate(block = location %>%  factor()) )
  data_in_all <- as.h2o(x = rbind(test22 ,(plants_data_wide2 %>% mutate(block = 5+(location %>%  as.numeric())))))
  
  
  scaledata_drf_summary <-  h2o.randomForest(x = predictors,#,"Treatment_Category"
                                           y = "richness_Plants",
                                           ntrees = 900,
                                           min_rows = 100,
                                           sample_rate = .95,
                                           stopping_rounds = 20,
                                           seed = 10000,
                                           score_each_iteration = T,
                                           balance_classes = T,
                                           max_depth = 4,
                                           mtries = 16,
                                           model_id = "plant_richness_evaluation",
                                           fold_column = "block",
                                           keep_cross_validation_predictions= TRUE,
                                           training_frame = data_in %>% as.h2o)
  
Alldata_drf_summary <-  h2o.randomForest(x = predictors,#,"Treatment_Category"
                                             y = "richness_Plants",
                                             ntrees = 900,
                                             min_rows = 100,
                                             sample_rate = .95,
                                             stopping_rounds = 20,
                                             seed = 10000,
                                             score_each_iteration = T,
                                             balance_classes = T,
                                             max_depth = 4,
                                             mtries = 16,
                                             model_id = "plant_richness_evaluation",
                                             fold_column = "block",
                                             keep_cross_validation_predictions= TRUE,
                                             training_frame = data_in_all %>% as.h2o)
  
  path <- paste0(getwd(),"./processed_data/rf_models/paper_230810")
  mojo_destination <- h2o.saveModel(kbsdata_drf_summary, path = path,force = T,
                                    export_cross_validation_predictions = T)
  
  cvpreds <- h2o.getFrame(scaledata_drf_summary@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  # full_prediction = h2o.predict(scaledata_drf_summary, data_in)
  Nfold_data  = data_in %>% #data_in_split[[1]] %>% 
    as.data.frame() %>% mutate(response = richness_Plants)%>% 
    cbind(as.data.frame(cvpreds)) %>% dplyr::rename(prediction = predict) %>% 
    mutate(field = "Field, field model")
  # Nfold_data$full_prediction = as.vector(full_prediction)
  # Nfold_data$val = "Cross validation \nprediction"
  # # pout$val = "Full model"
  
  direct_prediction = plant_scale_predictllist$richness_Plants%>% filter(!is.na(prediction)) %>% 
    mutate(field = "Field, plot model")
  Nfold_data_select = rbind(Nfold_data %>% dplyr::select(field, prediction, response),
                            direct_prediction %>% as.data.frame()%>% 
                              dplyr::select(field, prediction, response), fill = T)
  
  Nfold_data_select$relerror = abs(Nfold_data_select$response-Nfold_data_select$prediction)/Nfold_data_select$response
  
  stats_field_plot = response_stats_scale(model_in_list =Nfold_data_select%>% filter(field =="Field, plot model"), 
                                   response_name = "Plant richness, fields - plot model")
   stats_field_field = response_stats_scale(model_in_list =Nfold_data_select%>% filter(field =="Field, field model"), 
                                   response_name = "Plant richness, fields - field model")
  
  field_plots_error = Nfold_data_select %>% filter(field =="Field, plot model")
  sn_all_field_plots = mean((field_plots_error$relerror), na.rm =T)/
    sd((field_plots_error$relerror), na.rm =T)
  
  
  Fields_field_error = Nfold_data_select %>% filter(field =="Field, field model")
  sn_all_Fields_field = mean((Fields_field_error$relerror), na.rm =T)/
    sd((Fields_field_error$relerror), na.rm =T)
  
  g_partmodel = ggplot(Nfold_data_select %>% filter(field!=T))+
    geom_point(aes(x= response, y= prediction, col = relerror*100))+
    geom_abline()+
    scale_color_viridis_c(
    
  )+scale_color_viridis_c(trans = 'log10')+theme_bw()+
    labs(col = "Relative \nerror(%)", 
         x = "Observed plant richness", 
         y = "predicted plant richness")+
    facet_wrap(~field)
  
  allcvpreds <- h2o.getFrame(Alldata_drf_summary@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  # all_full_prediction = h2o.predict(Alldata_drf_summary, data_in_all)
  allNfold_data  = data_in_all %>% #data_in_split[[1]] %>% 
    as.data.frame() %>% mutate(response = richness_Plants)%>% 
    cbind(as.data.frame(cvpreds)) %>% dplyr::rename(prediction = predict)
  # allNfold_data$all_full_prediction = as.vector(all_full_prediction)

  
  allNfold_data$relerror = abs(allNfold_data$response-allNfold_data$prediction)/allNfold_data$response
  allNfold_data$field = is.na(allNfold_data$image_site) %>% factor(labels = c("Plots","Fields"),
                                                     levels = c(T, F))
  stats_all_plots = response_stats_scale(model_in_list =allNfold_data %>% filter(field =="Plots"), 
                       response_name = "Plant richness, plots - full model")
  stats_all_fields = response_stats_scale(model_in_list =allNfold_data %>% filter(field =="Fields"), 
                                         response_name = "Plant richness, feilds - full model")
  stats_all = response_stats_scale(model_in_list =allNfold_data, 
                                          response_name = "Plant richness, all - full model")
  stats_table_fieds = rbind(stats_field_plot, stats_field_field, stats_all_fields, stats_all_plots, stats_all)
  plots_error = allNfold_data %>% filter(field =="Plots")
  sn_all_plots = mean((plots_error$relerror), na.rm =T)/
    sd((plots_error$relerror), na.rm =T)
  
  
  Fields_error = allNfold_data %>% filter(field =="Fields")
  sn_all_Fields = mean((Fields_error$relerror), na.rm =T)/
    sd((Fields_error$relerror), na.rm =T)
  
 g_fullmodel=  ggplot(rbind(Nfold_data_select, 
                            allNfold_data %>% 
                              dplyr::select(field, relerror, 
                                            prediction, response)) %>% filter(field!=T))+
    geom_point(aes(x= response, y= prediction, col = relerror*100))+
    geom_abline()+facet_grid(cols = vars(field))+
    scale_color_viridis_c(trans = 'log10')+theme_bw()+
    labs(col = "Relative \nerror(%)", 
         x = "Observed plant richness", 
         y = "Predicted plant richness")
  
 ggsave(filename = "./figures/crossvalidation/paper_model_230726/eval_fields_compare.jpeg",
        plot = g_fullmodel, width = 7, height = 3, dpi =300)
 
  scale_up_from_all = allNfold_data %>% filter(block>5) 
  
  
    stats_scale = response_stats_scale(model_in_list = Nfold_data, response_name = "richness_Plants")
    stats_all = response_stats_scale(model_in_list = allNfold_data, response_name = "richness_Plants")
    stats_all_PLOTS = response_stats_scale(model_in_list = allNfold_data %>% filter(block <6), response_name = "richness_Plants")
    stats_scale_up_all = response_stats_scale(model_in_list = scale_up_from_all, response_name = "richness_Plants")
  
    group_data_in = data_in_all %>% as.data.frame() %>% 
      group_by(Treatment_Category, trtmt_code, block, doy, year) %>% 
      summarise_all(~mean(.x, na.rm=T)) 

 
    g_all_out = ggarrange(gfullmodel, )
    testitall = test2_test %>% 
      group_by(Treatment_Category, trtmt_code, block, bin, year) %>% 
      summarise_all(~mean(.x, na.rm=T)) 
    plants_data_long = read_rds("./processed_data/planet/SR/buf_site_9/clean_plants_scale/plot_data/plots_datalong_2023_08_09")
    

    # Understanding the field evalution sites########
    plants_data_long2 = plants_data_long %>% 
      mutate(block = (location %>% as.numeric())+5)
    all_data_long = rbind(test2_test, plants_data_long2)
    plants_data_long_group = plants_data_long2 %>% 
      group_by(block, bin, year) %>% 
      summarise_all(~mean(.x, na.rm=T)) 
    
    
    all_long = rbind(plants_data_long_group,testitall)
    g_all = ggplot(data = all_long %>% filter(!is.na(richness_Plants)))+
      geom_line(aes(x = doy, y= GNDVI,
                    col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32, 64)) %>% 
                      factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                    # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                    # size = ENS_Plants %>% cut(4) %>% factor,
                    group = interaction(x,y)))+# ,alpha = .5
      facet_grid(cols = vars(year))+theme_minimal()+
      labs(x = "Day of the year", col = "Plant richness")+
      scale_color_viridis_d(direction = -1) 
    
    g_all2 = ggplot(data = all_long %>% filter(doy ==133))+
      # geom_point( 
      #   aes(y = richness_Plants, x= red_3mean,
      #       col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
      #         factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
      #       # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
      #       # size = ENS_Plants %>% cut(4) %>% factor,
      #       group = interaction(x,y)), alpha = .15, size = .5)+# alpha = .15
      geom_point(data = all_long %>% filter(doy ==133), 
                 aes(y = richness_Plants, x= red_3mean,
                     col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                       factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                     # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                     # size = ENS_Plants %>% cut(4) %>% factor,
                     shape= (block<6),
                     group = interaction(x,y)), size = 2)+
      facet_grid(cols = vars(year))+theme_minimal()+
      labs(y = "Plant richness", col = "Plant richness",  x = "Red band reflectance at 3x resolution")+
      scale_color_viridis_d(direction = -1)+
      theme(panel.border = element_rect(color = "purple",
                                        fill = NA,
                                        linewidth =  1))
    
    g_all2
    
    
    g3 = ggplot()+
      geom_point(data = all_data_long %>% filter(doy ==133), 
                 aes(y = richness_Plants, x= ARVI_3mean,
                     col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
                       factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                     # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                     # size = ENS_Plants %>% cut(4) %>% factor,
                     group = interaction(x,y)), alpha = .15, size = .5)+# alpha = .15
      geom_point(data = all_long %>% filter(doy ==133),
                 aes(y = richness_Plants, x= ARVI_3mean,
                     col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>%
                       factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                     # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                     # size = ENS_Plants %>% cut(4) %>% factor,
                     shape= (block<6),
                     group = interaction(x,y)), size = 2,)+
      facet_grid(cols = vars(year))+theme_minimal()+
      labs(y = "Plant richness", col = "Plant richness", x = "ARVI at 3x resolution")+
      scale_color_viridis_d(direction = -1)+
      theme(panel.border = element_rect(color = "purple",
                                        fill = NA,
                                        linewidth =  1))
    
    g3
    