

# 1) set directory and paths  -----------------------
# 
tar_load(shape_plots) ## add back to work flow 
tar_load(image_sites)
tar_load(dem_full)
dem_rast <- rast(dem_full)

## get images###
box_list_path = "./processed_data/scratch/box_list"
if(!file.exists(box_list_path)){
  point_path <-  "./raw_data/extent/scale_up_points.csv"
  plant_scale_path <-  "./raw_data/BiodivPrediction_Plant_Composition.xlsx"
  scale_points <-  read_csv(point_path)%>%
    st_as_sf(coords = c("Longitude","Latitude"), crs = st_crs(4326)) %>% 
    st_transform(.,shape_plots%>% crs)  ####################### hard coded
  
  
  library(readxl)
  plant_scale_data <-  read_xlsx(plant_scale_path)
  library(vegan)
  
  ## calcualte richness and ENS. 
  plant_scale_richness <-  plant_scale_data %>% group_by(wpt, location) %>% 
    summarize(richness = n(), abundance = sum(total_count), 
              ENS = exp(diversity(total_count)),
              shannon = diversity(total_count),
              simpson = diversity(total_count, index = "simpson")) %>% 
    ungroup()
  
  ### make spatial
  plant_scale_points <-  merge((scale_points  %>% 
                                  filter(Name %in% 
                                           (plant_scale_richness$wpt %>% unique()))),
                               plant_scale_richness, 
                               by.x = "Name", 
                               by.y = "wpt",
                               all.x = T) %>% 
    st_intersection(., image_sites %>% st_transform(.,crs(scale_points)))
  
  ### build the unit box 
  box <-  point_2_poly(x_vert = as.matrix(c(-.5,.5)),
                       y_vert =  as.matrix(c(-.5,.5)),
                       side_length =  10)#*3.28) ## keep in meters
  
  # write_sf(obj = st_as_sf(plant_scale_points), dsn="./processed_data/extents/plant_scale_point.shp")
  # write_sf(obj = st_as_sf(planet_scale_boxes_data), dsn="./processed_data/extents/plant_scale_boxes.shp")
  
  # 1.0) Build all the boxes and join with data. -------------------
  planet_scale_boxes_data = (box$x+plant_scale_points$geometry) %>%
    st_as_sf(crs = crs(shape_plots)) %>%
    st_join(.,plant_scale_points[,c("Name","location",
                                    "richness","ENS",
                                    "shannon", "simpson",
                                    "image_site")])
  
  # to vector in rast
  box <- vect(planet_scale_boxes_data) %>% 
    project(., crs(shape_plots)) #site_rast[[1]][[1]]))

  box_list <- split(box, box$image_site)
  names(box_list) <- paste0("buf_site_",box$image_site %>% unique())
  box_list_paths =  map2(box_list,names(box_list), function(H,G) {
    temp_path = file.path("./processed_data/scratch", paste0(G, "_boxes.shp"))
    terra::writeVector(H, temp_path, overwrite=TRUE)
    temp_path})
  
  box_data  = box %>% as.data.frame %>%  
    rename(.,"ENS_Plants" = "ENS", "richness_Plants" = "richness",
           "shannon_Plants" =  "shannon"  ,  "simpson_Plants" =  "simpson") %>% 
    mutate(., id = seq(1:n()))
  # split to lsit 
  write_rds(box_list_paths, "./processed_data/scratch/box_list")
  write_rds(box_data, "./processed_data/scratch/box_data")
  
  
  # saveRDS(scaleUP_data, scaleUP_data_path)
}else{
  box_list_paths = read_rds(box_list_path)
  box_list = map(box_list_paths, ~terra::vect(.x))
  box_data = read_rds("./processed_data/scratch/box_data")
}

box_data_out = map(box_list, function(H){
       a = st_as_sf(H) %>% st_transform(., st_crs(4236)) %>% st_centroid() 
       b = a %>% mutate(x = st_coordinates(a)[,1],
                        y = st_coordinates(a)[,2]) %>% 
         st_drop_geometry()
   }) %>% do.call(rbind,.) %>% dplyr::select(-c("f"))
write_csv(box_data_out, "./processed_data/scratch/Plant_div_location_summary.csv")

#11.1) BUild load and pollinator data #####====================

transect_data_all_path = "./processed_data/EXTRACTIONS/scale/ScaleUP_data_transects"
if(!file.exists(transect_data_all_path)){
  
  points <-  fread("./raw_data/extent/scale_up_points.csv")
  raw_data <-  fread("./raw_data/bd_prediction_pollinators.csv")
  
  
  raw_data <- raw_data %>% group_by(SiteWptNotes) %>% 
    summarise_all(., sum)
  raw_data <- raw_data%>% 
    mutate(richness_bumblebee = apply(raw_data[,3:8]>0, MARGIN = 1, FUN =sum),
           abundance_bumblebee = apply(raw_data[,3:8], MARGIN = 1, FUN =sum),
           ENS_bumblebee = exp(diversity(raw_data[,3:8])),
           richness_bee = apply(raw_data[,10:17]>0, MARGIN = 1, FUN =sum),
           abundance_bee = apply(raw_data[,10:17], MARGIN = 1, FUN =sum),
           ENS_bee = exp(diversity(raw_data[,10:17])),
           richness_butterfly = apply(raw_data[,c(18:38)]>0, MARGIN = 1, FUN =sum),
           abundance_butterfly = apply(raw_data[,c(18:38)], MARGIN = 1, FUN =sum),
           ENS_butterfly = exp(diversity(raw_data[,c(18:38)]))) %>% 
    as.data.table
  
  raw_data[,site := SiteWptNotes %>% strsplit(., ";") %>% sapply("[",1),]
  
  # plot(raw_data$richness, raw_data$ENS)
  ###11.1 get transects ################ 
  x = raw_data$SiteWptNotes %>% as.list
  x_numbers <- regmatches(x, gregexpr("[[:digit:]]+", x))
  coords <- regmatches(x, gregexpr("[[:digit:]]+.[[:digit:]]+", x))
  x_numbers[[1]][1]
  
  
  angle <- function(B,A){
    direction = atan(abs((B - A)[2] / (B - A)[1]))
    xsign = ifelse((B - A)[1]<0, -1,1)
    ysign = ifelse((B - A)[2]<0, -1,1)
    
    print(direction)
    if(xsign >0 &ysign>0){direction =direction
    }else if(xsign<0 & ysign>0){direction = pi -direction
    }else if(xsign<0 & ysign<0){direction = pi +direction
    }else if(xsign>0 & ysign<0){direction = 2*pi -direction
    }
    print(direction)
    return(direction)
  }
  
  # Set points 
  
  
  d = 100*3.28
  i = 3
  geometry <- list()
  for(i in 1:length(x_numbers)){
    print(i)
    point_length = length(x_numbers[[i]])
    
    if(point_length <=3){
      if(point_length == 3){
        x_numbers[[i]] <-  x_numbers[[i]][x_numbers[[i]]!="100"]
      }
      start =  x_numbers[[i]][1]
      end =  x_numbers[[i]][2]
      
      coord_start_pnt  = points[Name==start,c("Easting", "Northing")] 
      coord_end_pnt = points[Name==end,c("Easting", "Northing")]  #%>% as.numeric
      
      coord_start = coord_start_pnt %>% as.double()
      coord_end = coord_end_pnt %>% as.double()
      
      AN = angle(B = coord_end, A = coord_start)
      
      # if(AN>(pi/2)&AN<(2*pi)){
      coord_end2 = c(coord_start[1]+ (d)*cos(AN),
                     coord_start[2]+ (d)*sin(AN))
      # }else{
      #   coord_end2 = c(coord_start[1]+ (d)*sin(AN),
      #                  coord_start[2]+ (d)*cos(AN)) 
      # }
      
      temp_coords = as.data.frame(list("x" = c(coord_start[1], coord_end2[1]),
                                       "y" = c(coord_start[2], coord_end2[2])))
      
      temp <-  vect(temp_coords, geom = c("x", "y")) %>% as.lines
      crs(temp)  <- "epsg:5625"
      
      
      
    }else if(point_length==5){
      
      coord_start = coords[i] %>% unlist %>% as.numeric %>% as.list 
      names(coord_start) = c("y","x")
      coord_start$x = coord_start$x*-1
      temp_vect  = vect(coord_start %>% as.data.frame(), geom = c("x", "y"))
      crs(temp_vect) = "epsg:4326"
      temp_vect <- project(temp_vect, crs(geometry[[i-1]]))
      coord_start = crds(temp_vect)
      
      end =  x_numbers[[i]][5]
      coord_end = points[Name==end,c("Easting", "Northing")]  %>% as.numeric
      
      AN = angle(B = coord_end, A = coord_start)
      
      coord_end2 = c(coord_start[1]+ (d)*cos(AN),
                     coord_start[2]+ (d)*sin(AN))
      
      dist_check = sqrt((coord_start[1]-coord_end2[1])^2+(coord_start[2]-coord_end2[2])^2)
      
      temp_coords = as.data.frame(list("x" = c(coord_start[1], coord_end2[1]),
                                       "y" = c(coord_start[2], coord_end2[2])))
      
      temp <- vect(temp_coords, geom = c("x", "y")) %>% as.lines
      crs(temp)  <- "epsg:5625" ###
      
      
      # mapview(st_as_sf(temp))
    }else if(point_length==6){
      
      start =  x_numbers[[i]][1]
      end =  x_numbers[[i]][6]
      
      coord_start_pnt  = points[Name==start,c("Easting", "Northing")] 
      coord_end_pnt = points[Name==end,c("Easting", "Northing")]  #%>% as.numeric
      
      coord_start = coord_start_pnt %>% as.double()
      coord_end = coord_end_pnt %>% as.double()
      
      
      coord_mid = coords[i] %>% unlist %>% as.numeric %>% as.list 
      names(coord_mid) = c("y","x")
      coord_mid$x = coord_mid$x*-1
      temp_vect  = vect(coord_mid %>% as.data.frame(), geom = c("x", "y"))
      crs(temp_vect) = "epsg:4326"
      r = rast()
      crs(r) = "epsg:5625"
      temp_vect <- project(temp_vect, crs(r))
      coord_mid= crds(temp_vect)
      
      AN = angle(B = coord_end, A = coord_mid)
      
      d1 = sqrt((coord_end-coord_mid)[1]^2+(coord_end-coord_mid)[1]^2)
      
      # if(AN>(pi/2)&AN<(3*pi/2)){
      coord_end2 = c(coord_mid[1]+ (d-d1)*cos(AN),
                     coord_mid[2]+ (d-d1)*sin(AN))
      # }else{
      #   coord_end2 = c(coord_mid[1]+ (d-d1)*sin(AN),
      #                  coord_mid[2]+ (d-d1)*cos(AN)) 
      #              }
      # 
      temp_coords = as.data.frame(list("x" = c(coord_start[1],coord_mid[1], coord_end2[1]),
                                       "y" = c(coord_start[2], coord_mid[2], coord_end2[2])))
      temp<- vect(temp_coords, geom = c("x", "y")) %>% as.lines
      crs(temp)  <- "epsg:5625"
      # mapview(st_as_sf(vect(coord_end_pnt,geom = c("Easting", "Northing"), crs(temp)))) +
      #   mapview(st_as_sf(vect(coord_start_pnt,geom = c("Easting", "Northing"), crs(temp))))+
      #   mapview(st_as_sf(temp))
    }
    values(temp) = raw_data[i, c(41:50)]
    geometry = c(geometry, project(temp, crs(shape_plots)) %>% buffer(., 5))
    # rm(temp)
    
  }
  
  
  all_pollinator_transects <- terra::vect(geometry)
  all_pollinator_transects <- terra::intersect(all_pollinator_transects, 
                                               project(vect(image_sites), 
                                                       crs(all_pollinator_transects)))
  
  all_pollinator_transects_list <- split(all_pollinator_transects, all_pollinator_transects$image_site)
  
  
  ####
  names(all_pollinator_transects_list) <- paste0("buf_site_",1:9)
  all_pollinator_transects_list_paths =  map2(all_pollinator_transects_list,
                                              names(all_pollinator_transects_list), function(H,G) {
    temp_path = file.path("./processed_data/scratch", paste0(G, "_transects.shp"))
    terra::writeVector(H, temp_path, overwrite=TRUE)
    temp_path})
  ####
  write_rds(all_pollinator_transects_list_paths, "./processed_data/scratch/all_pollinator_transects_list")
  
  tryit <- map2(all_pollinator_transects_list, site_rast, function (H,G){
    # H = all_pollinator_transects_list[[5]]
    # G = site_rast[[5]]
    # J = G[[1]]
    # K = names(G)[[1]]
    temp = map2(G,names(G), function(J, K){
      temp_rast = crop(J,H) %>% mask(., H)
    })
    temp_dt <- rast_to_df(temp)
    shape_plots_dt <- terra::rasterize(x =  H,
                                       y = temp[[1]],
                                       field = c("site")) %>% 
      as.data.frame(.,xy=T)
    
    temp_dt= merge(temp_dt, shape_plots_dt)
    temp_dt= merge(temp_dt,H %>% as.data.frame(), by = "site")
    
  })
  transect_data_all <- tryit %>% 
    rbindlist(fill = T) %>% 
    dplyr::rename(.,
                  "richness_Bumblebees" =  "richness_bumblebee",
                  "abundance_Bumblebees" = "abundance_bumblebee" ,
                  "ENS_Bumblebees" = "ENS_bumblebee",
                  "richness_Bee_Groups" = "richness_bee",
                  "abundance_Bee_Groups" = "abundance_bee",
                  "ENS_Bee_Groups" = "ENS_bee",
                  "richness_Butterflies" = "richness_butterfly",
                  "abundance_Butterflies"= "abundance_butterfly" ,
                  "ENS_Butterflies" ="ENS_butterfly")
  

  saveRDS(transect_data_all, transect_data_all_path)
}else{ 
  all_pollinator_transects_list_paths = read_rds("./processed_data/scratch/all_pollinator_transects_list")
  all_pollinator_transects_list = map(all_pollinator_transects_list_paths, ~vect(.x))
  transect_data_all= read_rds(transect_data_all_path)
}

transects_data_out = map(all_pollinator_transects_list, function(H){
       a = st_as_sf(H) %>% st_transform(., st_crs(4236)) %>% st_centroid() 
       b = a %>% mutate(x = st_coordinates(a)[,1],
                        y = st_coordinates(a)[,2]) %>% 
         st_drop_geometry()
   }) %>% do.call(rbind,.) %>% dplyr::select(-c("f"))

names(transects_data_out)[1:9]= c("richness_Bumblebees" ,
                                  "abundance_Bumblebees" ,
                                  "ENS_Bumblebees" ,
                                  "richness_Bee_Groups",
                                  "abundance_Bee_Groups",
                                  "ENS_Bee_Groups" ,
                                  "richness_Butterflies",
                                  "abundance_Butterflies",
                                  "ENS_Butterflies")

write_csv(transects_data_out, "./processed_data/scratch/Pollinator_div_location_summary.csv")
## image locations 
getvegindex = function(H){
  H$NDVI = (H[[4]]-H[[3]])/(H[[4]]+ H[[3]])
  H$SAVI = (H[[4]]-H[[3]])*(1+1)/(H[[4]]+ H[[3]]+1)
  H$ARVI = (H[[4]]-H[[3]]- (H[[1]]-H[[3]]))/(H[[4]]+ H[[3]]-(H[[1]]-H[[3]]))
  H$GNDVI = (H[[4]]-H[[2]])/(H[[4]]+H[[2]])
  H$CVI = (H[[4]])/(H[[2]])*(H[[3]])/(H[[2]])
  H$EVI = 2.5 * ((H[[4]] - H[[3]]) / (H[[4]] + 6 * H[[3]] - 7.5 * H[[1]] + 1))
  return(H)
  }

hard_path = "./processed_data/planet/SR/buf_site_1/01_RAW/"
subtit <- paste0("_", 1:9)
site_paths <- purrr::map(subtit, function(H) gsub(hard_path, pattern = "_1", replacement = H))

# 2) list all XM files -----------------------
H =site_paths[[1]]

raw_xml_files = purrr::map(site_paths, function(H){                            ### all here 
  temp_files = list.files(H, pattern =".xml$", full.names = T)
  
  raw_dates <- strsplit(temp_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  raw_xml_files = temp_files[lubridate::year(raw_dates)>2020]
})
all_files = NA
all_files_pollinator  = NA
H = site_paths[[2]]
G = box_list[[2]]
I = all_pollinator_transects_list[[2]]
raw_xml_files = purrr::pmap(list(site_paths, box_list, all_pollinator_transects_list), function(H, G, I){                            ### all here 
  clean_dir = H %>% gsub(x = ., pattern = "01_RAW", "clean_plants_scale")
  clean_dir_pollinator = H %>% gsub(x = ., pattern = "01_RAW", "clean_pollinator_scale")
  if(!dir.exists(clean_dir)){dir.create(clean_dir)
  if(!dir.exists(clean_dir_pollinator)){dir.create(clean_dir_pollinator)}
  test_files = list.files(H, pattern ="BGRN_SR_clip.tif$", full.names = T)
  temp_files = list.files(H, pattern ="_udm2_.*tif$", full.names = T)
  
  test_files_dates <- strsplit(test_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  temp_files_dates <- strsplit(temp_files,"_") %>% sapply(.,"[",7) %>% as.Date()
  
  test_files = test_files[lubridate::year(test_files_dates)>2020&
                            lubridate::year(test_files_dates)<2023]
  temp_files = temp_files[lubridate::year(temp_files_dates)>
                            2020&lubridate::year(test_files_dates)<2023]
  
  library(tictoc)
  
  for(k in seq(test_files)){  
   
    test = rast(test_files[k]) %>% terra::crop(G)# %>% terra::mask(G)
    test_udm_2 = rast(temp_files[k])%>% terra::crop(G)#%>% terra::mask(G)
  
    names(test_udm_2) = c("Clear map", "Snow map", "Shadow map", 
                        "Light haze map", "Heavy haze map", "Cloud map", 
                        "Confidence map","Unusable pixels")
    # plotRGB(test, stretch = "lin", b=1, g=2, r=3)
  
    percent_cells = (values(test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90)) %>% sum())/
      (terra::ncell(test_udm_2$`Clear map`))
    
  # plotRGB(test*test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90), stretch = "lin", b=1, g=2, r=3)
  # plot(test_udm_2)
  
  if(percent_cells>.5){
    # all_files = c(all_files, test_files[k]) %>% na.omit
    temp_out = test*test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90)/10000
    
    ## build out all veg indixes. 
    temp_out = getvegindex(temp_out)
    outpath =  file.path(clean_dir, basename(test_files[k]))
    if(!file.exists(outpath)){
    writeRaster(x = temp_out, filename = file.path(clean_dir, basename(test_files[k])))
      }
    
  }
  
    test_pollinator = rast(test_files[k]) %>% terra::crop(I)# %>% terra::mask(G)
    test_udm_2_pollinator = rast(temp_files[k])%>% terra::crop(I)#%>% terra::mask(G)
    
    names(test_udm_2_pollinator) = c("Clear map", "Snow map", "Shadow map", 
                          "Light haze map", "Heavy haze map", "Cloud map", 
                          "Confidence map","Unusable pixels")
    # plotRGB(test, stretch = "lin", b=1, g=2, r=3)
    
    percent_cells = (values(test_udm_2_pollinator$`Clear map`*(test_udm_2_pollinator$`Confidence map`>90)) %>% sum())/
      (terra::ncell(test_udm_2_pollinator$`Clear map`))
    
    # plotRGB(test*test_udm_2$`Clear map`*(test_udm_2$`Confidence map`>90), stretch = "lin", b=1, g=2, r=3)
    # plot(test_udm_2)
    
    if(percent_cells>.5){
      # all_files_pollinator = c(all_files_pollinator, test_files[k]) %>% na.omit
      temp_out = test_pollinator*test_udm_2_pollinator$`Clear map`*(test_udm_2_pollinator$`Confidence map`>90)/10000
      
      ## build out all veg indixes. 
      temp_out = getvegindex(temp_out)
      outpath =  file.path(clean_dir_pollinator, basename(test_files[k]))
      if(!file.exists(outpath)){
        writeRaster(x = temp_out, filename = file.path(clean_dir_pollinator, basename(test_files[k])))
      }
      
    }  
    
 
  }
  }
  return(list(clean_dir, clean_dir_pollinator))
  
})


## aggregate all #################

plant_scale_files = sapply(raw_xml_files, "[", 1) %>% 
  map(.,~list.files(.x, full.names = T, pattern = "tif$")) 
pollinator_scale_files = sapply(raw_xml_files, "[", 2) %>% 
  map(.,~list.files(.x, full.names = T, pattern = "tif$")) 

file_list = pollinator_scale_files[[1]]
shape_in = all_pollinator_transects_list[[1]] 
label = shape_in$site %>%unique()
build_rast_set = function(file_list, shape_in, buffer=15, label){
  path_out =file.path(file_list %>% dirname() %>% unique(),"plot_data") 
  if(!dir.exists(path_out)){dir.create(path_out)}
  label= gsub(x = label, replacement = "_",pattern = " ")
  path_out = path_out %>% paste0(., "/data_",
                                 label,
                                 "_23_08_09")
  
  if(!file.exists(path_out)){
  getall = map(file_list, rast)
  # getall_scale = map(getall, scale)   #decide what to do here 
  base_r = getall[[1]][[1]]%>% crop(shape_in %>% buffer(buffer))
  tic()
  getall_3 = map(getall, ~terra::aggregate(.x %>% crop(base_r), 
                                           fact=3, fun="mean", na.rm=T)%>% 
                   terra::project(.,base_r, method = "near")) 
  toc()
  getall_5 = map(getall, ~terra::aggregate(.x %>% crop(base_r), 
                                           fact=5, fun="mean", na.rm=T)%>% 
                   terra::project(.,base_r, method = "near")) 
  
  getall_3sd = map(getall, ~terra::aggregate(.x%>% crop(base_r), 
                                             fact=3,fun="sd", na.rm=T)%>% 
                     terra::project(.,base_r, method = "near")) 
  getall_5sd = map(getall, ~terra::aggregate(.x%>% crop(base_r), 
                                             fact=5,fun="sd", na.rm=T)%>% 
                     terra::project(.,base_r, method = "near")) 
  
  
  dem_site_path = file.path("./processed_data/DEM",
                            paste0(file_list[[1]]  %>% strsplit("/") %>% sapply(.,"[",5),
                                   "-",label, ".tif"))
  ### to here...
  if(!file.exists(dem_site_path)){
    writeRaster(dem_rast$base_dem%>% 
                  terra::project(base_r), 
                dem_site_path, overwrite = T)
  }
  # writeRaster(dem_rast$base_dem, "./processed_data/DEM/DEM_9.tif", overwrite = T)
  dem_rast_plots_WB = get_twi(dem_path = dem_site_path, 
                              sitename = dem_site_path  %>% basename() %>% 
                                gsub(pattern = "[.]tif",x = .,replacement = ""))
  dem_rast_plots = rast(dem_rast_plots_WB) %>% 
    # terra::project( getall[[1]]$blue, method = "near") %>% 
    terra::crop(.,base_r)
  
  dem_rast_plots_dt = dem_rast_plots %>% as.data.frame(.,xy=T) %>% 
    dplyr::select(-c("base_dem","DEM_breached", "DEM_smoothed"))
  # test_all = terra::scale(getall[[1]])
  Date = map(file_list, ~.x %>% str_split("_") %>% sapply(.,"[",8) %>% as.Date)
  getall_dt = map2(getall, Date, ~as.data.frame(.x %>% 
                                           terra::project(.,base_r, method = "near"),
                                         xy=T) %>% 
                    mutate(Date = .y)) 
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
  
  shape_plots_dt <- map(names(shape_in), ~terra::rasterize(x = shape_in,
                                         y = base_r,
                                         field = .x)) %>% 
  rast %>% as.data.frame(xy = T)

  getall_out2 = shape_plots_dt %>% merge(getall_out)

  names(getall_out2) <- names(getall_out2) %>% gsub(x = ., pattern= " ", replacement = "_")
  
  saveRDS(getall_out2, path_out)
  
  getall_out2
}else{
  getall_out2 = readRDS(path_out)
}
}

### get plant df for scale up sites #####

plant_data = map2(plant_scale_files, box_list, function(H,G){ 
  if(length(G)>1){
    map(G %>% split(G$Name), function(K) build_rast_set(file_list =H, shape_in = K,
                                                        buffer= 30, label = K$location %>%unique()))%>% 
      do.call(rbind,.)
  }else{
    build_rast_set(file_list =H, shape_in = G,
                   buffer= 30, label = G$location %>%unique())
  } 
})%>% 
  do.call(rbind,.)

plant_data =plant_data %>% rename(richness_Plants =  richness,
                                            ENS_Plants= ENS)

### get pollinator df for scale up sites #####
pollinator_data = map2(pollinator_scale_files, all_pollinator_transects_list, function(H,G){ 
  if(length(G)>1){
    map(G %>% split(G$site), function(K) build_rast_set(file_list =H, shape_in = K,
                                                      buffer= 30, label = K$site %>%unique()))%>% 
      do.call(rbind,.)
  }else{
    build_rast_set(file_list =H, shape_in = G,
                   buffer= 30, label = G$site %>%unique())
  } 
})%>% 
  do.call(rbind,.)
pollinator_data =pollinator_data %>% rename(richness_Bumblebees =  richness_0,
                                            abundance_Bumblebees= abundance0,
                                            ENS_Bumblebees= ENS_bumbl0,
                                            richness_Bee_Groups= richness_1, 
                                            abundance_Bee_Groups = abundance1, 
                                            ENS_Bee_Groups= ENS_bee, 
                                            richness_Butterflies= richness_2, 
                                            abundance_Butterflies= abundance2, 
                                            ENS_Butterflies= ENS_butte0, 
                                            location = site)
names(plant_data$location)
#### join all #####################3
data_in = (plant_data %>% split(plant_data$image_site))[[2]]
data_in2 = (pollinator_data %>% split(pollinator_data$image_site))[[8]]
break_int = KBS_break_int = read_rds("./processed_data/scratch/KBS_break_int")
cli_data_series = temp
label = "clean_plants_scale"

plant_idcol =c("x","y","location","image_site","year",
               "richness_Plants","abundance_Plants" ,
               "ENS_Plants","DEM_flow_accum","DEM_slope",
               "DEM_TWI")
pollinator_idcol = c("x","y","location","image_site","year",
                     "richness_Bumblebees","abundance_Bumblebees" ,
                     "ENS_Bumblebees","richness_Bee_Groups",
                     "abundance_Bee_Groups","ENS_Bee_Groups",
                     "richness_Butterflies","abundance_Butterflies",
                     "ENS_Butterflies","DEM_flow_accum","DEM_slope",
                     "DEM_TWI")
allbands = names(plant_data)[c(11:20,22:61)]
build_periods = function(data_in, cli_data_series, break_int=NA, 
                         allbands= allbands,idcol = plant_idcol, 
                         label = "clean_plants_scale"){
  data_wide_path = file.path("./processed_data/planet/SR",
                             paste0("buf_site_",data_in$image_site %>% unique()),
                             label, "plot_data/plots_datawide_2023_08_09")
  data_long_path = file.path("./processed_data/planet/SR",
                             paste0("buf_site_",data_in$image_site %>% unique()),
                             label, "plot_data/plots_datalong_2023_08_09")
  
  
      if(T){#!file.exists(data_wide_path)
    test = data_in %>% left_join(rename(cli_data_series, Date = date)) 
    
    
    if(all(is.na(break_int))){
      dayrange =test$day %>% range
      series_int =15
      break_int = seq(dayrange[1]-1,dayrange[2]+1,series_int)
      
      }
    
    test2a = test %>% mutate(bin = cut(day, break_int)) 
    
    test2 = test2a %>% 
      group_by(x,y,location, image_site, bin, year) %>% 
      select_if(., is.numeric) %>% 
      summarise_all(~mean(.x, na.rm=T))
    
    int_table = as.data.frame(list(doy = (break_int+7)[1:19], 
                                   bin = levels(test2$bin)))# %>% unique()) %>% na.omit))
    
    id_col_in =  names(test2)[ names(test2)%in%idcol]# names(test2)[c(1:15,68:70)] #, 91:93
    all_bands = names(test2)[ names(test2)%in%all_bands]#names(test2)[17:67]
    test2 = test2 %>% left_join(int_table) 
    bin_list = test2$bin %>% unique() %>% na.omit()
    tableX = test2[,c("bin", "year")] %>% unique() %>% table
    missing_21 = names(tableX[,1])[tableX[,1]==0]
    missing_22 = names(tableX[,2])[tableX[,2]==0]
    
    missing_data = NA
    if(missing_21 %>% length>0 & missing_22 %>% length>0){
    missing_data = rbind(as.data.frame(list(year = 2021, bin = missing_21)),
          as.data.frame(list(year = 2022, bin= missing_22)))
    }else if(missing_21 %>% length>0){
      missing_data = as.data.frame(list(year = 2021, bin = missing_21))
      }else if(missing_22 %>% length>0){
        missing_data = as.data.frame(list(year = 2022, bin = missing_22))
        }
    
    
    if(any(!is.na(missing_data))){
    missing_data = full_join(test2 %>% ungroup %>% dplyr::select(id_col_in) %>% unique(), 
                             missing_data,relationship = "many-to-many") %>% left_join(int_table)
    
    test3_test = test2 %>% filter(!is.na(bin)) %>% 
      rbind(missing_data) %>% arrange(doy) %>% arrange(year)
    
  hmmm =   test3_test %>% group_by_at((id_col_in)) %>% 
      #fill(color, age, gender) %>% #default direction down
      fill(all_of(all_bands), .direction = "downup")
  tableX = hmmm[,c("bin", "year")] %>% unique() %>% table
    }else{
    hmmm = test2
  }
  ###################  
  # test2_test = test2 %>% filter(!is.na(bin))
    # 
    # if(missing_21 %>% length>0){
    #   for(i in seq(missing_21)){
    #     if(missing_21[[i]]== "(66,81]"){
    #       test2a = test2_test %>% filter(year == 2021,
    #                                      bin == bin_list[which(bin_list ==missing_21[[i]])+1])%>% 
    #         mutate(bin = missing_21[[i]], doy = int_table$doy[which(bin_list ==missing_21[[i]])])
    #       test2_test = rbind(test2_test, test2a) %>% mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #     }else if (missing_21[[i]]== "(336,351]"){
    #       test2a = test2_test %>% filter(year == 2021,
    #                                      bin == bin_list[which(bin_list ==missing_21[[i]])-1])%>% 
    #         mutate(bin = missing_21[[i]], doy = int_table$doy[which(bin_list ==missing_21[[i]])])
    #       test2_test = rbind(test2_test, test2a) %>% mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #     }else{
    #       test2a = test2_test %>% filter(year == 2021, 
    #                                      bin == bin_list[which(bin_list ==missing_21[[i]])-1]|
    #                                        bin == bin_list[which(bin_list ==missing_21[[i]])+1])
    #       test2c = test2a %>% group_by(x,y,location, image_site, year) %>% 
    #         summarise_all(~mean(.x, na.rm=T)) %>%
    #         mutate(bin = missing_21[[i]], doy = int_table$doy[which(bin_list ==missing_21[[i]])])
    #       test2_test = rbind(test2_test, test2c) %>% mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #     }
    #   }
    #   }
    # 
    # if(missing_22 %>% length>0){
    #   for(i in seq(missing_22)){
    #     if(missing_22[[i]]== "(66,81]"){
    #       test2a = test2_test %>% filter(year == 2022,
    #                                      bin == bin_list[which(bin_list ==missing_22[[i]])+1])%>% 
    #         mutate(bin = missing_22[[i]], doy = int_table$doy[which(bin_list ==missing_22[[i]])])
    #       test2_test = rbind(test2_test, test2a) %>% 
    #         mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #     }else if (missing_22[[i]]== "(336,351]"){
    #       test2a = test2_test %>% filter(year == 2022,
    #                                      bin == bin_list[which(bin_list ==missing_22[[i]])-1])%>% 
    #         mutate(bin = missing_22[[i]], doy = int_table$doy[which(bin_list ==missing_22[[i]])])
    #       test2_test = rbind(test2_test, test2a) %>% 
    #         mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #       }else{
    #     test2a = test2_test %>% filter(year == 2022, 
    #                               bin == bin_list[which(bin_list ==missing_22[[i]])-1]|
    #                                 bin == bin_list[which(bin_list ==missing_22[[i]])+1])
    #     test2c = test2a %>% group_by(x,y,location, image_site, year) %>% 
    #       summarise_all(~mean(.x, na.rm=T)) %>% 
    #       mutate(bin = missing_22[[i]], doy = int_table$doy[which(bin_list ==missing_22[[i]])])
    #     test2_test = rbind(test2_test, test2c) %>% mutate(Date = (as.Date(paste0(doy,"_",year),"%j_%Y")))
    #   }
    #   }
    # }
    # tableX = test2_test[,c("bin", "year")] %>% unique() %>% table
    ########################
   test21 = hmmm %>% filter(!is.na(bin)) %>% 
      pivot_wider(., 
                  id_cols = c(id_col_in), 
                  values_from = c(all_bands),
                  names_from = doy)
   saveRDS(hmmm, data_long_path)
   saveRDS(test21, data_wide_path)
      }else{
      test21 = readRDS(data_wide_path)
      }
}

data_in = (pollinator_data %>% split(pollinator_data$image_site))[[7]]
break_int = KBS_break_int = read_rds("./processed_data/scratch/KBS_break_int")
cli_data_series = temp
label = "clean_pollinator_scale"
pollinator_data_wide = map((pollinator_data %>% split(pollinator_data$image_site)), 
    ~build_periods(data_in = .x,
              cli_data_series = temp, 
              break_int=KBS_break_int, allbands = allbands,
              idcol = pollinator_idcol,
              label = "clean_pollinator_scale"))


pathpollinators = "./processed_data/planet/SR/buf_site_9/clean_pollinator_scale/plot_data/plots_datawide_2023_08_09"
Poll_path_out = map(paste0("_",1:9), ~ pathpollinators %>% gsub(pattern = "_9", replacement = .x))

pollinator_data_wide = map(Poll_path_out, read_rds) %>% do.call(rbind,.)

plants_data_wide = map((plant_data %>% split(plant_data$image_site)), 
                           ~build_periods(data_in = .x,
                                          cli_data_series = temp, 
                                          break_int=KBS_break_int,
                                          allbands = allbands,
                                          idcol = plant_idcol,
                                          label = "clean_plants_scale"))


pathplants = "./processed_data/planet/SR/buf_site_9/clean_plants_scale/plot_data/plots_datawide_2023_08_09"
plant_path_out = map(paste0("_",1:9), ~ pathplants %>% gsub(pattern = "_9", replacement = .x))

plants_data_wide = map(plant_path_out, ~read_rds(.x)) %>% do.call(rbind,.)

pathplants = "./processed_data/planet/SR/buf_site_9/clean_plants_scale/plot_data/plots_datalong_2023_08_09"
plant_path_out = map(paste0("_",1:9), ~ pathplants %>% gsub(pattern = "_9", replacement = .x))

plants_data_long = map(plant_path_out, read_rds) %>% do.call(rbind,.)
# ggplot(plants_data_long)+
#   geom_line(aes(x= doy, y = GNDVI_5mean, col = richness_Plants, group = location))+
#   theme_minimal()+scale_color_viridis_c()+facet_wrap(~year)
