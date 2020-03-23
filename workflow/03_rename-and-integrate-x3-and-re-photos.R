# Purpose: Rename the photo files within the RGB and RedEdge folders such that each file has a separate name

library(here)
library(tidyverse)
library(purrr)
library(exiftoolr)
library(sf)
library(viridis)

# A helper function to get arbitrary color palette for arbitrary range of color values in base R
scale_color_base <- function(value, colors=c("white", "black"), na.rm=FALSE, mapToRange=range(value, na.rm=na.rm), alpha=1, printRecast=FALSE)
{
  if (na.rm==FALSE & any(is.na(value)))
    stop("There are NAs in your vector. Using na.rm=TRUE will remove them for the mapping calculations, but add them back in in the final vector. Try that.")
  
  recast_value <- (value[!is.na(value)] - mapToRange[1]) / (diff(mapToRange))
  recast_value[recast_value < 0] <- 0
  recast_value[recast_value > 1] <- 1
  
  color_fnc <- colorRamp(colors=colors)
  plot_colors <- rep("NA", length(value))
  
  plot_colors[!is.na(value)] <- rgb(color_fnc(recast_value)/255, alpha=alpha)
  
  if (printRecast) {
    print(recast_value)
  }
  
  return (plot_colors)
}


all_sites <- list.files(here::here("data/data_drone/L0/flight-logs"))

# I manually iterated through the sites in order to monitor progress, since it involved moving around so much data

### IMPORTANT NOTE ###
# In order to not overload any one disk drive, I put each set of original photos onto the solid state drive one at a time
# That is, I took all the photos in their original folder structure configuration and copied them to the data/data_working/**current_site**/_x3_photos (or _re_photos)
# folder in order to run this script.

(current_site <- all_sites[32])

x3_dir <- paste0("data/data_working/", current_site, "_x3_photos")
re_dir <- paste0("data/data_working/", current_site, "_re_photos")

target_dir <- paste0("data/data_working/", current_site, "/", current_site, "_photos")

mission_footprint <- st_read(paste0("data/data_drone/L0/mission-footprint/site-bounds/", current_site, "_site-bounds.geojson"))
dem <- raster::raster("data/data_raw/srtm_30m.tif")

# Create new directory to store all the photos
dir.create(here::here(paste0("data/data_working/", current_site)))
dir.create(here::here(target_dir))

# Combine the RGB photos, and append an x3 to the filename
x3_photo_filepaths <- list.files(here::here(x3_dir), recursive = TRUE, full.names = TRUE, pattern = ".JPG")
x3_photo_old_names <- list.files(here::here(x3_dir), recursive = TRUE, pattern = ".JPG")
x3_photo_new_filepaths <- paste0(target_dir, "/x3_", str_replace_all(string = x3_photo_old_names, pattern = "/", replacement = "_"))

file.copy(from = x3_photo_filepaths, to = here::here(x3_photo_new_filepaths), overwrite = TRUE)

# The RedEdge photos need some extra curation help because the camera was just on a timer and so lots of extra pictures
# were taken (of the takeoff point, while the sUAS was taking off/landing, etc.)

# Before proceeding, add the photos of the calibration panel to a directory below the re_dir directory called "calibration/"
# The directory with the calibration photos (you've already manually added these photos to "here(re_dir)/calibration/" right?)
re_calibration_dir <- list.files(here::here(re_dir), full.names = TRUE, pattern = "calibration")
re_photo_dirs <- list.files(here::here(re_dir), full.names = TRUE)
re_nonCalibration_dirs <- re_photo_dirs[re_photo_dirs != re_calibration_dir]

# The list of all the RedEdge photos' filenames that aren't calibration panel photos
re_photo_filepaths <- list.files(re_nonCalibration_dirs, recursive = TRUE, full.names = TRUE, pattern = ".tif")

# Establish the takeoff point for the rededge mission as the location of the very first photo
# Note this assumes the first photo is taken of the takeoff point, not of the calibration panels
# (which were almost always taken some distance from the takeoff point). Ensure that the first
# photo is of the takeoff point (the plywood, in our case) before proceeding
takeoff_point <-
  re_photo_filepaths[1] %>% 
  exif_read(tags = c("GPSLongitude", "GPSLatitude", "GPSAltitude", "DateTimeOriginal", "Yaw", "Pitch", "Roll", "PressureAlt")) %>% 
  type_convert() %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326)

# Get the DEM elevation for the takeoff point location
takeoff_elev <- raster::extract(dem, takeoff_point, method = "bilinear")

# Calculate the offset between what the RedEdge altitude thinks it is and what the DEM thinks it is
takeoff_elev_offset <- takeoff_elev - takeoff_point$PressureAlt

# Use the exiftoolsr package to read the metadata of all of the rededge photos
# First filter any files that are of size 0 bytes (these are just bad triggers of the RedEdge and exiftools won't play nicely with them)
# Filter out photos that are missing any of the lat/lon or altitude

exif_data <-
  data_frame(SourceFile = list.files(here::here(re_dir), recursive = TRUE, full.names = TRUE, pattern = ".tif")) %>% 
  dplyr::mutate(re_photo_old_names = list.files(here::here(re_dir), recursive = TRUE, pattern = ".tif")) %>% 
  dplyr::mutate(DestFile = paste0(here::here(target_dir), "/re_", str_replace_all(string = re_photo_old_names, pattern = "/", replacement = "_"))) %>% 
  dplyr::filter(file.size(SourceFile) > 0) %>% 
  dplyr::filter(!str_detect(string = DestFile, pattern = "re_calibration")) %>% 
  left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude", "PressureAlt")), by = "SourceFile") %>% 
  type_convert() %>%
  dplyr::filter(!is.na(GPSLongitude), !is.na(GPSLatitude), !is.na(PressureAlt)) %>% 
  dplyr::select(SourceFile, DestFile, GPSLongitude, GPSLatitude, PressureAlt) 

# Convert the metadata to a spatial object, calculate the agl (above ground level altitude)
# as the difference between the RedEdge recorded altitude, the DEM altitude at that point
# and the offset between the DEM and the RedEdge altimeter

# Perform 2 checks: 1st, is the drone high enough? Give some wiggle room here and say any photo taken where 
# we calculate the drone was 90m or higher (agl) was within the mission parameters
# 2nd, check to see which photos are within the mission footprint (plus a 5 meter buffer) as defined by the
# flight logs from the RGB camera

# Then, "mission photos" are all photos that meet both of the above criteria

exif_data <-
  exif_data %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) 

exif_data <-
  exif_data %>% 
  mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>% 
  mutate(band_name = case_when(band_num == 1 ~ "blue",
                               band_num == 2 ~ "green",
                               band_num == 3 ~ "red",
                               band_num == 4 ~ "nir",
                               band_num == 5 ~ "re")) %>% 
  mutate(agl = PressureAlt - raster::extract(x = dem, y = ., method = "bilinear") + takeoff_elev_offset) %>% 
  mutate(suas_at_altitude = ifelse(agl > 80, yes = TRUE, no = FALSE)) %>% 
  mutate(suas_within_footprint = st_intersects(x = st_transform(., 3310), y = st_buffer(st_transform(mission_footprint, 3310), 10), sparse = FALSE)[ , 1]) %>% 
  mutate(mission_photo = suas_at_altitude & suas_within_footprint)

exif_data %>% 
  group_by(band_name) %>% 
  summarize(n = n())


dev.off()
par(mfrow = c(1, 2), mar = rep(0, 4))
plot_colors <- scale_color_base(value = exif_data$agl, colors = viridis(10))
plot(exif_data$geometry, type = "n")
plot(exif_data$geometry[exif_data$mission_photo], pch = 19, add = TRUE, col = plot_colors[exif_data$mission_photo])
plot(mission_footprint$geometry, add = TRUE)

plot(exif_data$geometry, type = "n")
plot(exif_data$geometry[!exif_data$mission_photo], pch = 19, add = TRUE, col = plot_colors[!exif_data$mission_photo])
plot(mission_footprint$geometry, add = TRUE)

photos_to_copy <-
  exif_data %>% 
  filter(mission_photo)

photos_to_copy %>% 
  group_by(band_name) %>% 
  summarize(n = n())


# Copy the photos from the working folder to the final photos folder using the new names
file.copy(from = photos_to_copy$SourceFile, to = photos_to_copy$DestFile, overwrite = TRUE)


# Add the calibration panel photos to the final photos directory
calibration_SourceDir <- paste0(re_dir, "/calibration")
calibration_SourceFile <- list.files(here::here(calibration_SourceDir), full.names = TRUE)
calibration_SourceNames <- list.files(here::here(calibration_SourceDir))
calibration_DestFile <- here::here(paste0(target_dir, "/re_calibration_", calibration_SourceNames))

file.copy(from = calibration_SourceFile, to = calibration_DestFile)

# Get the metadata for the photos and save to disk

# After manually curating the photos, get the list of all photos used for analysis
final_photo_list <- list.files(here::here(target_dir))

metadata <- data.frame(site = current_site, processed_photos = final_photo_list)

# Add the file size to the metadata
metadata <-
  metadata %>% 
  mutate(file_size_MB = file.size(here::here(paste0(target_dir, "/", metadata$processed_photos))) / 1e6)

write_csv(metadata, here::here(paste0("data/data_drone/L0/photos-metadata/", current_site, "_photos-metadata.csv")))

# These lines create the folders that will be used to house the intermediate Pix4D projects before they get merged
dir.create(here::here(paste0("data/data_working/", current_site, "/", current_site, "_x3")))
dir.create(here::here(paste0("data/data_working/", current_site, "/", current_site, "_re")))

