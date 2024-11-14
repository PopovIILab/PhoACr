# Set working directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Install or call libraries
if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse, raster, sf, rnaturalearth, ggrepel, ggspatial)

# Read Natural Earth 2 map as raster object

natearth_map <- 
  raster::stack('NE2_HR_LC_SR_W_DR/NE2_HR_LC_SR_W_DR.tif') %>% # import tiff as rasterStack
  crop(extent(1, 85, 40, 65)) %>%                                   # subset to desired extent
  as.array %>%                                                      # convert to 3D array
  `/`(255) %>%                                                      # switch to proportions to meet rgb() requirements
  apply(c(1, 2), function(x) rgb(matrix(x, ncol = 3))) %>%          # collapse layers to RGB colors
  annotation_raster(1, 85, 40, 65)

# Create a map of world countries
world_countries <- ne_countries(scale = "small", returnclass = "sf")

# Set coordinates of cities
cities <- data.frame(city = c("Rostov-on-Don (2021-2023)", 
                              "Moscow (2015, 2021)",
                              "Vadum (2014)",
                              "Mustasaari (2016)*",
                              "Bad Segeberg (2007)",
                              "Novosibirsk (2021-2023)",
                              "Moensted (2016)**"),
                     lat = c(47.2357, 
                             55.6552, 
                             57.11778,
                             63.1167,
                             53.9370,
                             54.8324,
                             56.2646),
                     lng = c(39.7015, 
                             36.7436, 
                             9.857101,
                             21.7167,
                             10.3041,
                             83.3192,
                             9.1114))

cities <- st_as_sf(cities, coords = c("lng", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")

# Plot a map
ggplot(data = world_countries) +
  geom_sf()+
  natearth_map +
  annotate("text", x = 55, y = 54, 
           label = "Russia", color = "black", 
           alpha = 0.33, size = 6, fontface = "bold") +
  annotate("text", x = 23, y = 64, 
           label = "Finland", color = "black", 
           alpha = 0.33, size = 4, fontface = "bold") +
  geom_text_repel(data = cities,
                  aes(x = lng, y = lat, label = city),
                  size = 4,
                  fontface = "bold",
                  nudge_x = c(0, 0.15,-0.15, 0.3, 0.45, 0.6),
                  nudge_y = c(0.75, -0.75, 1, -1, 1.25, -1.25)) +
  geom_text_repel(data = world_countries %>% filter(!(admin %in% c("Kosovo"))),
                  stat = "sf_coordinates",
                  aes(label = name, geometry = geometry),
                  color = "black",
                  alpha = 0.33,
                  size = 4,
                  fontface = "bold") +
  geom_point(data = cities,
             aes(x = lng, y = lat),
             size = 3, 
             shape = 21,
             fill = "darkred") +
  coord_sf(xlim = c(1, 85), ylim = c(45, 65), expand = FALSE)+
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") +
  ylab("Latitude")

# Save it in 600 dpi
ggsave("images/geo-map.png", width = 16, height = 8, dpi=600)

# Special thanks to these guides:
# 1. https://stackoverflow.com/questions/69852503
# 2. https://r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
# 3. https://tsamsonov.github.io/r-geo-course/