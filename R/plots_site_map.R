# # Site map figure 
# kbs_shapes_WSG = shape_plots %>% st_as_sf %>% filter(!is.na(block)) %>% st_transform(st_crs(4326))
# # load july image 
# image_hard_path = "./processed_data/planet/SR/buf_site_2/01_RAW/5802311_1658620_2022-07-21_2477_BGRN_SR_clip.tif"
# site2_image = rast(image_hard_path)
# site2_image_wsg = site2_image %>% terra::project(., crs("epsg:4326"))
# df <- as.data.frame(site2_image_wsg, xy= TRUE)
# 
# 
# df <- df %>% filter(red != 0)   
# df <- df %>% rename(Red = red,   #Rename bands
#                     Green = green,
#                     Blue = blue)
# df2 = df %>% mutate(
#   Red = ifelse(Red>(mean(Red)+3*sd(Red)),(mean(Red)+3*sd(Red)),Red),
#   Green = ifelse(Green>(mean(Green)+3*sd(Green)),(mean(Green)+3*sd(Green)),Green),
#   Blue = ifelse(Blue>(mean(Blue)+3*sd(Blue)),(mean(Blue)+3*sd(Blue)),Blue),
# )
# 
# Map = ggplot()+                   #plot map
#   geom_raster(data = df2, aes(x = x, y =y, fill = rgb(r = Red,
#                          g = Green,
#                          b = Blue,
#                          maxColorValue = max((df2[,3:5])))), show.legend = FALSE) +
#   scale_fill_identity() + #coord_equal()+
#   theme_void()  
# g1 = geom_sf(data = kbs_shapes_WSG,
#              aes(), fill = "lightblue")
# Map+g1
# 
# ggplot()+g1
# ggsave(Map,                                            #save map
#        filename = paste0(here(), "/satellite_img.jpg"),
#        dpi = 200)
# 
# df <- as.data.frame(basemap, xy= TRUE)
# 
# plotRGB(site2_image, b=1, g=2, r=3, stretch ="lin")
# lines(shape_plots %>% st_as_sf %>% filter(!is.na(block)), col = "red")
# plot((box_list[[2]] %>% buffer(.,10)), col = "yellow", border = "yellow", add= T)
# plot((box_list[[2]] %>% buffer(.,10)), col = "yellow", border = "yellow", add= T)
# plot(all_pollinator_transects[all_pollinator_transects$image_site==2], col = "orange",border = "orange", add=T)
# 

################here############
airsct_image_file = "./raw_data/airscout/rgb/2021_05_13_1249001_geo.tif" 
# 
# airsct_image_file = "./raw_data/airscout/rgb/2021_07_02_1457241_geo.tif" 
airsct_image = rast(airsct_image_file)

plotRGB(airsct_image)

airsct_image_WSG = crop(airsct_image, kbs_shapes_WSG %>% st_buffer(50))
plotRGB(airsct_image_WSG)
lines(kbs_shapes_WSG)
airsct_image_WSG_dt = airsct_image_WSG %>% as.data.frame(xy=T)
names(airsct_image_WSG_dt)[3:5]= c("r","g","b")
AIR_MAP = ggplot()+ 
  geom_raster(data=airsct_image_WSG_dt, aes(x=x, y=y, fill=rgb(r,g,b, maxColorValue = 255))) + 
  scale_fill_identity()+
  geom_sf(data= kbs_shapes_WSG, 
                 aes(col = (kbs_shapes_WSG$trtmt_code  %>% factor(labels = c("Corn", "Prairie", "Sorghum", "Sorghum + covercrop",
                                                                            "Transitional Switchgrass", "Switchgrass",
                                                                            "Miscanthus", "Native Grasses", "Poplar" ,
                                                                            "Successional")))),
                 fill = NA)+labs(col= "Treatment") + theme_void()


ggsave(filename ="./figures/crossvalidation/paper_model_230726/Airscout_plot.jpeg", plot = AIR_MAP, 
       dpi = 700, width = 5, height = 3)

### states = 
states_sf =read_sf("D:/fowler53/BioD_KBS/raw_data/extent/cb_2018_us_state_500k/cb_2018_us_state_500k.shp")
# states_sf = states()
plot(states_sf["NAME"])

ggplot(states_sf %>% filter(STUSPS != "HI", STUSPS !="AK", STUSPS !="PR",
                            STUSPS !="VI",STUSPS !="MP",STUSPS !="AS",STUSPS !="GU"))+
  geom_sf(fill = NA)+theme_void()+
  geom_sf(data = box_list[[2]][1,] %>% st_as_sf(), size = 5, col = "blue") 
