
# 
# all_rast <- pmap(list(goodfiles, site_paths, buffsites),   function(H, G, J){
#   gc()
#   clean_planet(good_list = H, hard_path = G, buffsite = J)}) 
# 
# site_paths_2 = site_paths[[2]]
# 
# buffsites_2 = buffsites[[2]]
# buffsites_2_shape = st_read(buffsites_2)
# 
# plot(buffsites_2_shape$geometry)


# 1) set directory and paths  -----------------------
# 

tar_load(shape_plots)

shape_plots_kbs  = shape_plots %>% st_as_sf %>% filter(!is.na(block))

shape_plots_kbs$Treatment_Category = ifelse(shape_plots_kbs$trtmt_code == "G1" |
        shape_plots_kbs$trtmt_code  == "G2" |
        shape_plots_kbs$trtmt_code  == "G3", 
       "Annual Monoculture",
       ifelse( shape_plots_kbs$trtmt_code  == "G4" |
               shape_plots_kbs$trtmt_code  == "G5" |
               shape_plots_kbs$trtmt_code  == "G6" |
               shape_plots_kbs$trtmt_code  == "G7",
              "Low-diversity Perennial Grass",
              "High-diversity Perennial Polyculture"))


shape_plots_kbs_ext = shape_plots_kbs %>% st_bbox() %>% st_as_sfc() %>% vect

add_trtmt_code <- function(shape, trtmt_code){
  shape@data$trtmt_code <- trtmt_code
  return(shape)
}

## image locations 
buffsites2 <- list.files("./processed_data/extents/", pattern = "buf_site_..shp$", full.names = T)[2]

buffsites2

getvegindex = function(H){
  H$NDVI = (H[[4]]-H[[3]])/(H[[4]]+ H[[3]])
  H$SAVI = (H[[4]]-H[[3]])*(1+1)/(H[[4]]+ H[[3]]+1)
  H$ARVI = (H[[4]]-H[[3]]- (H[[1]]-H[[3]]))/(H[[4]]+ H[[3]]-(H[[1]]-H[[3]]))
  H$GNDVI = (H[[4]]-H[[2]])/(H[[4]]+H[[2]])
  H$CVI = (H[[4]])/(H[[2]])*(H[[3]])/(H[[2]])
  H$EVI = 2.5 * ((H[[4]] - H[[3]]) / (H[[4]] + 6 * H[[3]] - 7.5 * H[[1]] + 1))
  return(H)
  }

# site_vect <- lapply(buffsites, vect)

hard_path = "./processed_data/planet/SR/buf_site_1/01_RAW/"
subtit <- paste0("_", 1:8)
site_paths <- purrr::map(subtit, function(H) gsub(hard_path, pattern = "_1", replacement = H))

# 2) list all XM files -----------------------
H =site_paths[[1]]
raw_xml_files = purrr::map(site_paths, function(H){                            ### all here 
  temp_files = list.files(H, pattern =".xml$", full.names = T)
  
  raw_dates <- strsplit(temp_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  raw_xml_files = temp_files[lubridate::year(raw_dates)>2020]
})
all_files = NA
H = site_paths[[2]]
raw_xml_files = purrr::map(site_paths[2], function(H){                            ### all here 
  clean_dir = H %>% gsub(x = ., pattern = "01_RAW", "clean_kbs")
  if(!dir.exists(clean_dir)){dir.create(clean_dir)}
  
  test_files = list.files(H, pattern ="BGRN_SR_clip.tif$", full.names = T)
  temp_files = list.files(H, pattern ="_udm2_.*tif$", full.names = T)
  
  test_files_dates <- strsplit(test_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  temp_files_dates <- strsplit(temp_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  
  test_files = test_files[lubridate::year(test_files_dates)>2020&lubridate::year(test_files_dates)<2023]
  temp_files = temp_files[lubridate::year(temp_files_dates)>2020&lubridate::year(test_files_dates)<2023]
  
  library(tictoc)
  
  for(k in seq(test_files)){  
    tic()
    test = rast(test_files[k]) %>% terra::crop(shape_plots_kbs_ext) %>% terra::mask(shape_plots_kbs_ext)
    test_udm_2 = rast(temp_files[k])%>% terra::crop(shape_plots_kbs_ext)%>% terra::mask(shape_plots_kbs_ext)
  
    names(test_udm_2) = c("Clear map", "Snow map", "Shadow map", 
                        "Light haze map", "Heavy haze map", "Cloud map", 
                        "Confidence map","Unusable pixels")
  # plotRGB(test, stretch = "lin", b=1, g=2, r=3)
  
    percent_cells = (values(test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90)) %>% sum())/
      (terra::ncell(test_udm_2$`Clear map`))
    
  # plotRGB(test*test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90), stretch = "lin", b=1, g=2, r=3)
  # plot(test_udm_2)
  
  if(percent_cells>.5){
    all_files = c(all_files, test_files[k]) %>% na.omit
    temp_out = test*test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90)/10000
    
    ## build out all veg indixes. 
    temp_out = getvegindex(temp_out)
    writeRaster(x = temp_out, filename = file.path(clean_dir, basename(test_files[k])))
    
  }
  
  toc()
  
  return(clean_dir)
  }
})


## aggregate all #################

KBS_files = list.files(raw_xml_files[[1]], full.names = T)

getall = map(KBS_files, rast)
# getall_scale = map(getall, scale)   #decide what to do here 
getall_3 = map(getall, ~terra::aggregate(.x%>% mask(shape_plots_kbs %>% st_buffer(-3)), fact=3, fun=function(.x) mean(.x, na.rm=T))%>% 
                 terra::disagg(fact = 3))
getall_5 = map(getall, ~terra::aggregate(.x%>% mask(shape_plots_kbs %>% st_buffer(-3)), fact=5, fun=function(.x) mean(.x, na.rm=T)) %>% 
                 terra::disagg(fact = 5))

getall_3sd = map(getall, ~terra::aggregate(.x%>% mask(shape_plots_kbs %>% st_buffer(-3)), fact=3,fun=function(.x) sd(.x, na.rm=T)) %>% 
                 terra::disagg(fact = 3))
getall_5sd = map(getall, ~terra::aggregate(.x%>% mask(shape_plots_kbs %>% st_buffer(-3)), fact=5,fun=function(.x) sd(.x, na.rm=T)) %>% 
                 terra::disagg(fact = 5))



tar_load(dem_full)
dem_rast <- rast(dem_full)
# writeRaster(dem_rast$base_dem %>% 
#               terra::project( getall[[1]]$blue, method = "near") %>% 
#               terra::crop(., getall[[1]]$blue),
#             "./processed_data/DEM/DEM_9.tif", overwrite = T)
dem_rast_plots_WB = get_twi(dem_path = "./processed_data/DEM/DEM_9.tif", sitename = "kbs_9m")
dem_rast_plots = rast(dem_rast_plots_WB) %>% 
  terra::project( getall[[1]]$blue, method = "near") %>% 
  terra::crop(., getall[[1]]$blue)

dem_rast_plots_dt = dem_rast_plots %>% as.data.frame(.,xy=T) %>% 
  dplyr::select(-c("base_dem","DEM_breached", "DEM_smoothed"))
# test_all = terra::scale(getall[[1]])

getall_dt = map(getall, ~as.data.frame(.x,xy=T) %>% 
                  mutate(Date = .x@ptr@.xData$filenames() %>%  
                           strsplit(.,"_") %>% 
                           sapply(.,"[",8) %>% 
                           as.Date())) 
band_names = names(getall_dt[[1]])[3:12]

getall_3sd_dt = map2(getall_3sd, getall_3, function(H,G) as.data.frame(H/G,xy=T) %>% 
                  rename_at(.vars=(band_names),
                            function(H){paste0(H,"_3cv")}))

getall_3mean_dt = map(getall_3, ~as.data.frame(.x,xy=T) %>% 
                      rename_at(.vars=(band_names),
                                function(H){paste0(H,"_3mean")}))


getall_5sd_dt = map2(getall_5sd, getall_5, function(H,G) as.data.frame(H/G,xy=T) %>% 
                      rename_at(.vars=(band_names),
                                function(H){paste0(H,"_5cv")}))

getall_5mean_dt = map(getall_5, ~as.data.frame(.x,xy=T) %>% 
                        rename_at(.vars=(band_names),
                                  function(H){paste0(H,"_5mean")}))



feature_list = list(getall_dt, getall_3sd_dt, getall_3mean_dt,
                    getall_5sd_dt, getall_5mean_dt)
getall_out = pmap(feature_list,
     function(H1, H2, H3, H4, H5){merge(H1,H2) %>%
         merge(., H3) %>% merge(., H4,all.x=T) %>% 
         merge(., H5,all.x=T)}) %>% 
  data.table::rbindlist()

getall_out = getall_out %>% merge(dem_rast_plots_dt)


kbs_data_path =  site_paths %>% gsub(x = ., pattern = "01_RAW", "kbs_plots_data") 
if(!dir.exists(kbs_data_path)){dir.create(kbs_data_path)}
kbs_data_path = kbs_data_path[[2]] %>% paste0(., "plots_",Sys.time() %>%format("%Y_%m_%d"))

if(!file.exists(kbs_data_path)){
    ### Make the treatment and block rasters
  shape_plots_raster <- terra::rasterize(x = vect(shape_plots_kbs),
                                         y = getall[[1]],
                                         field = c("trtmt_code"))
  
  shape_plots_raster2 <- terra::rasterize(x = vect(shape_plots_kbs),
                                          y = getall[[1]],
                                          field = c("block"))
  
  shape_plots_raster3 <- terra::rasterize(x = vect(shape_plots_kbs),
                                          y = getall[[1]],
                                          field = c("Treatment_Category"))
  
  ### combine them 
  shape_plots_raster <- c(shape_plots_raster, shape_plots_raster2, shape_plots_raster3)
  
  ## convert them to a DF
  shape_plots_dt <- shape_plots_raster%>% 
    terra::as.data.frame(.,xy=T, na.rm=F) %>%
    filter(!is.na(trtmt_code)) %>% 
    as.data.table() 
  
  ### get richness data 
  rich_data <-  fread("./processed_data/total_plot_bioD.csv")
  
  ####### functional groups ############
  kbs_functional_data <-  readxl::read_excel	("./raw_data/BiodivPrediction_Plant_Composition_explore.xlsx",
                                    sheet ="kbs functional groups")
  names(kbs_functional_data) = c("species", "functional")
  kbs_functional_data %>% 
    pivot_wider(names_from =functional,
                values_from = species)
  
  plant_data <-  fread("./raw_data/BCSE_Plants.csv")
  
  plant_data$percent_grass = plant_data %>% 
    dplyr::select((kbs_functional_data %>% 
                     filter(functional=="Grass"))$species) %>% 
    apply(MARGIN = 1, FUN = sum)
  
  plant_data$percent_Legume= plant_data %>% 
    dplyr::select((kbs_functional_data %>% 
                     filter(functional=="Legume"))$species) %>% 
    apply(MARGIN = 1, FUN = sum)
  
  plant_data$percent_Woody= plant_data %>% 
    dplyr::select((kbs_functional_data %>% 
                     filter(functional=="Woody"))$species) %>% 
    apply(MARGIN = 1, FUN = sum)
  
  plant_data$percent_Forb = plant_data %>% 
    dplyr::select((kbs_functional_data %>% 
                     filter(functional=="Forb"))$species) %>% 
    apply(MARGIN = 1, FUN = sum)
  
  plant_data_function = plant_data %>% 
    dplyr::select(Treatment, Rep,Quadrat, 
                  percent_grass,percent_Legume, percent_Woody, percent_Forb) %>% 
    group_by(Treatment, Rep) %>% summarise_all(mean) %>% ungroup() %>% 
    dplyr::select( -Quadrat ) %>%
    rename (trtmt_code = Treatment,
            block = Rep) %>% 
    mutate(total = percent_grass + percent_Legume + percent_Woody +
             percent_Forb, 
           trtmt_code = trtmt_code %>% 
             strsplit(" ") %>% sapply(.,"[",1))

  #######
  
  richness_wide <- rich_data %>% pivot_wider(.,
                                             id_cols = c("Treatment", "Rep","Treatment_Category"),
                                             names_from = "Group",
                                             values_from = c("richness","ENS", "abundance", "richness_z"))
  richness_wide$Treatment <- richness_wide$Treatment %>% strsplit(.," ") %>% sapply(., "[", 1)
  
  richness_wide = richness_wide %>% 
    dplyr::select(names(richness_wide)[names(richness_wide) %>% 
                                         grepl(x = .,pattern="Treatment|Rep|Treatment_Category|Plants|Bee Groups|Bumblebees|Butterflies|Birds")]) 
  
  ### merge the shape file with richness data
  kbs_shapes <- merge(shape_plots_dt, richness_wide, 
                      by.x = c("trtmt_code", "block"),
                      all.x=T,  by.y =c("Treatment", "Rep")) %>% 
    mutate(Treatment_Category = ifelse(is.na(Treatment_Category.x),
                                       Treatment_Category.y %>% as.character(), 
                                       Treatment_Category.x%>% as.character()),
           Treatment_Category.y = NULL, 
           Treatment_Category.x = NULL)
  
  kbs_shapes = kbs_shapes %>% merge(plant_data_function)
  
  # kbs_shapes$Treatment_Category = ifelse()kbs_shapes$trtmt_code[is.na(kbs_shapes$Treatment_Category)]
  
  ### Merge the planet reflectance with the KBS plots   
  kbs_data <- merge(getall_out, kbs_shapes, fill = T) 
  names(kbs_data) <- names(kbs_data) %>% gsub(x = ., pattern= " ", replacement = "_")
  
  saveRDS(kbs_data, kbs_data_path)
}else{
  kbs_data = readRDS("./processed_data/planet/SR/buf_site_2/kbs_plots_data/plots_2023_07_26")
}

#### join all #####################3

test = kbs_data %>% left_join(rename(temp, Date = date)) 
# 
# test2 = test %>% 
#   group_by(Treatment_Category,trtmt_code, block, Date) %>% 
#   summarize_all(~mean(.x, na.rm=T)) 

dayrange =test$day %>% range
series_int =15
break_int = seq(dayrange[1]-1,dayrange[2]+1,series_int)
# write_rds(break_int,"./processed_data/scratch/KBS_break_int")

test2a = test %>% mutate(bin = cut(day, break_int)) 

test2 = test2a %>% 
  group_by(x,y,Treatment_Category, trtmt_code, block, bin, year) %>% 
  summarise_all(~mean(.x, na.rm=T))

int_table = as.data.frame(list(doy = (break_int+7)[1:19], bin = (test2$bin %>% unique())))

id_col_in = names(test2)[c(1:5,7,59:86)] #, 91:93
all_bands = names(test2)[c(8:17, 19:58)]
test2 = test2 %>% left_join(int_table) 

test2[,c("bin", "year")] %>% unique() %>% table

test2a = test2 %>% filter(year == 2021, bin == "(81,96]"|bin == "(111,126]")
test2c = test2a %>% group_by(x, y, trtmt_code, block, Treatment_Category) %>% 
  summarise_all(~mean(.x, na.rm=T)) %>% mutate(bin = "(96,111]",
                                               doy = int_table$doy[int_table$bin=="(96,111]"])

  ### make sure I have all the bins here ####
kbs_data_long_path = "./processed_data/planet/SR/buf_site_2/kbs_plots_data/plots_datalong_2023_07_26"
test2_test = rbind(test2, test2c) %>% mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
saveRDS(test2_test, kbs_data_long_path)

test21 = test2_test %>% filter(!is.na(doy)) %>% 
  pivot_wider(., 
              id_cols = id_col_in, 
              values_from = all_bands,
              names_from = doy)
kbs_data_wide_path = "./processed_data/planet/SR/buf_site_2/kbs_plots_data/plots_datawide_2023_07_26"
# saveRDS(test21, kbs_data_wide_path)
# test21 = read_rds(kbs_data_wide_path)
predictor_columns = (test21 %>% names)[c(7:9,35:984)]#982]

h2o.init(nthreads = -1)

test22 = test21 

test22$ENS_res_Bumblebees = as.vector(lm(ENS_Bumblebees~ENS_Plants, test21)$residuals)
test22$ENS_res_Butterflies = as.vector(lm(ENS_Butterflies~ENS_Plants, test21)$residuals)
test22$ENS_res_Bee_Groups = as.vector(lm(ENS_Bee_Groups~ENS_Plants, test21)$residuals)
test22$ENS_res_Birds = as.vector(lm(ENS_Birds~ENS_Plants, test21)$residuals)
test22$richness_res_Bumblebees = as.vector(lm(richness_Bumblebees~richness_Plants, test21)$residuals)
test22$richness_res_Butterflies = as.vector(lm(richness_Butterflies~richness_Plants, test21)$residuals)
test22$richness_res_Bee_Groups = as.vector(lm(richness_Bee_Groups~richness_Plants, test21)$residuals)
test22$richness_res_Birds = as.vector(lm(richness_Birds~richness_Plants, test21)$residuals)

test22$richness_res_Plants_treatment = as.vector(lm(richness_Plants~Treatment_Category, test21)$residuals)
test22$ENS_res_Plants_treatment = as.vector(lm(ENS_Plants~Treatment_Category, test21)$residuals)

lalal = test22  %>% mutate(Treatment_Category = Treatment_Category %>% factor)
plants_CV = cross_validate(data_in = lalal, #test22
               predictors = predictor_columns[predictor_columns %in% names(lalal)], 
               response = "Treatment_Category",
               label = "kbs_plants_bigset",scope_label = "ugh",
               fold_column = "block",
               max_depth = 5, 
               mtries = 20)

h2o.varimp(plants_CV[[2]])
test3 = test2 %>% 
  group_by(Treatment_Category, trtmt_code, block, bin, year) %>% 
  summarise_all(~mean(.x, na.rm=T))

test3meanall= test3 %>%
  group_by(bin, year) %>%
  dplyr::select(c(9:18)) %>%
  summarise_all(~mean(.x, na.rm=T))

test3sdall= test3 %>%
  group_by(bin, year)%>%
  dplyr::select(c(9:18)) %>%
  summarise_all(~sd(.x, na.rm=T))

names(test3meanall)[3:12]= paste0("MeanAll_", names(test3meanall)[3:12])

names(test3sdall)[3:12]= paste0("SDAll_", names(test3sdall)[3:12])



test4 = test3 %>% #left_join(test3meanall)%>% left_join(test3sdall) %>% 
  group_by(bin, year) %>% 
  mutate(NDVI_norm = (NDVI - MeanAll_NDVI)/SDAll_NDVI,
         GNDVI_norm = (GNDVI - MeanAll_GNDVI)/SDAll_GNDVI,
         EVI_norm = (EVI - MeanAll_EVI)/SDAll_EVI)



library(gganimate)

p =ggplot(testitall, 
       aes(y = richness_Plants, x= GNDVI_5cv,
           col = richness_Plants %>% cut(breaks=c(0,2,4,8,16,32,64)) %>% 
             factor(labels = c("0-2", "2-4", "4-8", "8-16", "16-32", "32-64")))) + 
  geom_point() + 
  # Here comes the gganimate code
  transition_states(
    doy,
    transition_length = 2,
    state_length = 1
  ) +
  enter_fade() + 
  exit_shrink() +
  ease_aes('sine-in-out')


ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  ease_aes('linear')




ggplot()+
  geom_line(data = test4, aes(x = day, y= NDVI_norm,
                              col = (percent_Woody/total) %>% cut(3),#ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor, 
                              # linetype = ENS_Plants %>% cut(breaks=c(0,2,4,8,16,32)) %>% factor,
                              # size = ENS_Plants %>% cut(4) %>% factor,
                              group = interaction(trtmt_code, block)))+
  facet_grid(cols = vars(year), rows = vars(Treatment_Category))

                                     
ggplot()+
  geom_line(data = test %>% sample_n(10000), 
             aes(x = TT, y= NDVI, col = Treatment_Category, group=paste0(x,y)), alpha =.2) 

ggplot()+geom_point(data = test, aes(x = TT, y= NDVI, col = lubridate::year(Date) %>% factor()))

ggplot()+geom_point(data = test, aes(x = lubridate::yday(Date),
                                     y= NDVI, col = lubridate::year(Date) %>% factor()))


getall_dt$Date = 
terra::time(getall) =  strsplit(KBS_files,"_") %>% sapply(.,"[",7) %>% as.Date()
temp_files_dates <- strsplit(all_files,"_") %>% sapply(.,"[",7) %>% as.Date()

plot(temp_files_dates)
plotRGB(rast(all_files[[40]]), stretch = "lin", b=1, g=2, r=3)
# 3) read XML files ------
# 
meta_data = map(raw_xml_files, function(H) (map(H, filter_xml))) #%>% do.call(rbind,.) %>% as.data.frame()

# Step 2: 
# 4) manual image selection. ------
# 
i=3
