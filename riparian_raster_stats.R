
library(rgdal)
library(raster)
library(sf)
library(tidyverse)
library(lubridate)
library(grid)
library(gridExtra)

phenology_tracker <- function(filename, figure_title)
{
  # *** Sentinel Vegetation Stats ***
  # Load Data
  raster_02_09 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/02_09_2016.tif")
  raster_02_19 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/02_19_2016.tif")
  raster_04_19 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/04_19_2016.tif")
  raster_07_28 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/07_28_2016.tif")
  raster_09_06 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/09_06_2016.tif")
  raster_10_06 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/10_06_2016.tif")
  raster_11_05 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/11_05_2016.tif")
  raster_11_25 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/11_25_2016.tif")
  raster_12_25 <- stack("D:/SERDP/Pendleton/Sentinel2/214b/stacked_10m/12_25_2016.tif")
  slope <- raster("D:/SERDP/Pendleton/LiDAR/products/terrain_slope_10m_utm.tif")
  chm <- raster("D:/SERDP/Pendleton/LiDAR/CHM/chm_10m_utm.tif")
  # Reproject Stream Segments
  shapefile <- readOGR(filename)
  shapefile_utm <- spTransform(shapefile, crs(raster_02_09[[1]]))
  # Crop by Extent
  raster_02_09 <- crop(raster_02_09, extent(shapefile_utm))
  raster_02_19 <- crop(raster_02_19, extent(shapefile_utm))
  raster_04_19 <- crop(raster_04_19, extent(shapefile_utm))
  raster_07_28 <- crop(raster_07_28, extent(shapefile_utm))
  raster_09_06 <- crop(raster_09_06, extent(shapefile_utm))
  raster_10_06 <- crop(raster_10_06, extent(shapefile_utm))
  raster_11_05 <- crop(raster_11_05, extent(shapefile_utm))
  raster_11_25 <- crop(raster_11_25, extent(shapefile_utm))
  raster_12_25 <- crop(raster_12_25, extent(shapefile_utm))
  slope <- crop(slope, extent(shapefile_utm))
  chm <- crop(chm, extent(shapefile_utm))
  # Mask by Riparian Zones
  raster_02_09 <- mask(raster_02_09, shapefile_utm)
  raster_02_19 <- mask(raster_02_19, shapefile_utm)
  raster_04_19 <- mask(raster_04_19, shapefile_utm)
  raster_07_28 <- mask(raster_07_28, shapefile_utm)
  raster_09_06 <- mask(raster_09_06, shapefile_utm)
  raster_10_06 <- mask(raster_10_06, shapefile_utm)
  raster_11_05 <- mask(raster_11_05, shapefile_utm)
  raster_11_25 <- mask(raster_11_25, shapefile_utm)
  raster_12_25 <- mask(raster_12_25, shapefile_utm)
  slope <- mask(slope, shapefile_utm)
  chm <- mask(chm, shapefile_utm)
  
  # ** Cell Stats **
  # Set up data frame
  veg_data <- data.frame(year = rep(2016,9),
                         month = c(2, 2, 4, 7, 9, 10, 11, 11, 12),
                         day = c(09, 19, 19, 28, 06, 06, 05, 25, 25))
  veg_data$day_of_year <- lubridate::yday(paste(veg_data$year, veg_data$month, veg_data$day, sep="-"))
  # Extract Stats
  veg_data$ndvi_05p <- c(cellStats(raster_02_09[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_02_19[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_04_19[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_07_28[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_09_06[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_10_06[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_05[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_25[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                         cellStats(raster_12_25[[10]], stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]])
  veg_data$ndvi_25p <- c(cellStats(raster_02_09[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_02_19[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_04_19[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_07_28[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_09_06[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_10_06[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_05[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_25[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                         cellStats(raster_12_25[[10]], stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]])
  veg_data$ndvi_50p <- c(cellStats(raster_02_09[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_02_19[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_04_19[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_07_28[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_09_06[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_10_06[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_05[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_25[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                         cellStats(raster_12_25[[10]], stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]])
  veg_data$ndvi_75p <- c(cellStats(raster_02_09[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_02_19[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_04_19[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_07_28[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_09_06[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_10_06[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_05[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_25[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                         cellStats(raster_12_25[[10]], stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]])
  veg_data$ndvi_95p <- c(cellStats(raster_02_09[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_02_19[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_04_19[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_07_28[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_09_06[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_10_06[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_05[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_11_25[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]],
                         cellStats(raster_12_25[[10]], stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]])
  
  phenology_plot <- ggplot(veg_data) + 
    geom_line(aes(x=day_of_year, y=ndvi_50p), col="black") + 
    geom_line(aes(x=day_of_year, y=ndvi_25p), col="green3") + 
    geom_line(aes(x=day_of_year, y=ndvi_75p), col="green3") + 
    geom_line(aes(x=day_of_year, y=ndvi_05p), col="red", linetype="dashed") + 
    geom_line(aes(x=day_of_year, y=ndvi_95p), col="red", linetype="dashed") + 
    scale_y_continuous(limits=c(0,0.8)) + 
    ggtitle(figure_title) + 
    xlab("Day of Year") + 
    ylab("NDVI") + 
    theme_bw()
  
  
  veg_data_fixed <- data.frame(quartile = c("05","25","50","75","95"))
  veg_data_fixed$chm <- c(cellStats(chm, stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                          cellStats(chm, stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                          cellStats(chm, stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                          cellStats(chm, stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                          cellStats(chm, stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]])
  veg_data_fixed$slope <- c(cellStats(slope, stat=function(s,...){quantile(s,probs=c(0.05),na.rm=TRUE)})[[1]],
                            cellStats(slope, stat=function(s,...){quantile(s,probs=c(0.25),na.rm=TRUE)})[[1]],
                            cellStats(slope, stat=function(s,...){quantile(s,probs=c(0.50),na.rm=TRUE)})[[1]],
                            cellStats(slope, stat=function(s,...){quantile(s,probs=c(0.75),na.rm=TRUE)})[[1]],
                            cellStats(slope, stat=function(s,...){quantile(s,probs=c(0.95),na.rm=TRUE)})[[1]])
  
  # Difference between Select Dates
  apr_feb_diff <- (raster_04_19[[10]]-raster_02_09[[10]])/(raster_04_19[[10]]+raster_02_09[[10]])
  apr_feb_diff_df <- as.data.frame(apr_feb_diff, xy=TRUE)
  
  
  feb_sep_diff <- (raster_02_09[[10]]-raster_11_05[[10]])
  feb_sep_diff_df <- as.data.frame(feb_sep_diff, xy=TRUE)
  
  chm_df <- as.data.frame(chm, xy=TRUE)
  slope_df <- as.data.frame(slope, xy=TRUE)

  output_data_steady <- apr_feb_diff_df
  names(output_data_steady) <- c("x","y","apr_feb")
  output_data_steady$feb_sep <- feb_sep_diff_df[,3]
  output_data_steady$chm <- chm_df[,3]
  output_data_steady$slope <- slope_df[,3]
  
  output_data_pheno <- as.data.frame(raster_02_09[[10]], xy=TRUE)
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_02_19[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_04_19[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_07_28[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_09_06[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_10_06[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_11_05[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_11_25[[10]]))
  output_data_pheno <- cbind(output_data_pheno, as.data.frame(raster_12_25[[10]]))
  
  return(list(phenology_plot,
              veg_data,
              veg_data_fixed,
              output_data_steady,
              output_data_pheno))
  
}


santa_margarita_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/santa_margarita.shp", "Santa Margarita Riparian Zone")
las_pulgas_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/las_pulgas.shp", "Las Pulgas Riparian Zone")
mountain_creek_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/mountain_creek.shp", "Mountain Creek Riparian Zone")
lowland_creek_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/lowland_creek.shp", "Lowland Creek Riparian Zones")
upland_shrubs_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/upland_shrubs.shp", "Upland Shrubs 1")
upland_shrubs_2_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/upland_shrubs_2.shp", "Upland Shrubs 2")
grassland_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/grassland.shp", "Grassland")
north_slope_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/north_slope.shp", "North Slope")
north_slope_2_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/north_slope_2.shp", "North Slope 2")
suburban_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/suburb.shp", "Suburban Zone")
suburban_2_out <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/suburban.shp", "Suburban Zone 2")


chaparral <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/chaparral.shp", "Chaparral")
chaparral_2 <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/chaparral_two.shp", "Chaparral (Second Segment)")
chaparral_3 <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/chaparral_3.shp", "Chaparral (Third Segment)")
upland_chaparral <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/upland_chaparral.shp", "Upland Chaparral")

coastal_floodplain_veg <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/coastal_floodplain_vegetation.shp", "Coastal Floodplain")
lowland_riparian_small <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/lowland_riparian_small.shp", "Lowland Creek Riparian Zone (2) (Second Segment)")
riparian_pine <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/riparian_pine.shp", "Riparian Pine")
seasonal_lake <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/seasonal_lake.shp", "Seasonal Lake")

coastal_disturbed_grassland <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/coastal_disturbed_grassland.shp", "Coastal Disturbed Grassland")
coastal_grassland <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/coastal_grassland.shp", "Coastal Grassland")
golf_lawn <- phenology_tracker("D:/SERDP/Pendleton/Riparian_Zone/manual/golf_lawn.shp", "Irrigated Turfgrass (Golf Course)")



riparian_plots <- grid.arrange(
             santa_margarita_out[[1]], las_pulgas_out[[1]], mountain_creek_out[[1]], lowland_creek_out[[1]], coastal_floodplain_veg[[1]], lowland_riparian_small[[1]], riparian_pine[[1]],
             top = textGrob(
               "Greenness Phenology in Pendleton Riparian Zones",
               gp = gpar(fontface = 2, fontsize = 15)),
             bottom = textGrob(
               "NDVI values within each habitat by month. Median values are in black, 25th and 75th percentiles in green, and 5th and 95th percentiles in dashed red.",
               gp = gpar(fontface = 1, fontsize = 9),
               hjust = 1,
               x = 1
             ))

upland <- grid.arrange(
             upland_shrubs_out[[1]], upland_shrubs_2_out[[1]], chaparral[[1]], chaparral_2[[1]], chaparral_3[[1]], upland_chaparral[[1]],
             top = textGrob(
               "Greenness Phenology in Pendleton Chaparral",
               gp = gpar(fontface = 2, fontsize = 15)),
             bottom = textGrob(
               "NDVI values within each habitat by month. Median values are in black, 25th and 75th percentiles in green, and 5th and 95th percentiles in dashed red.",
               gp = gpar(fontface = 1, fontsize = 9),
               hjust = 1,
               x = 1
             ))

upland <- grid.arrange(
  grassland_out[[1]], coastal_disturbed_grassland[[1]], coastal_grassland[[1]], golf_lawn[[1]],
  top = textGrob(
    "Greenness Phenology in Pendleton Grassland",
    gp = gpar(fontface = 2, fontsize = 15)),
  bottom = textGrob(
    "NDVI values within each habitat by month. Median values are in black, 25th and 75th percentiles in green, and 5th and 95th percentiles in dashed red.",
    gp = gpar(fontface = 1, fontsize = 9),
    hjust = 1,
    x = 1
  ))