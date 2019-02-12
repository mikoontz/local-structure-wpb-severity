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


all_sites <- list.files(here::here("data/data_output/site_data"))

(current_site <- all_sites[36])

# Sites that require extra attention
# eldo_4k_2: two different take off points to maintain visual line of site and radio contact; need to calculate elevation for RedEdge photos separately
# sequ_4k_1: RedEdge GPS didn't work properly for eastern ~40% of survey area; need to process as separate subprojects with ground control points then merge (or throw away ~40% of survey area)
# sier_4k_1: RedEdge GPS didn't work properly for one flight in the middle of survey area: need to process as separate subprojects with ground control points then merge



# eldo_3k_1 (when current_site == all_sites[1]) uploaded to Pix4D Cloud on 2018-12-23 ## Needs more work; processed in two blocks
# eldo_3k_2 (when current_site == all_sites[2]) uploaded to Pix4D Cloud on 2018-12-26 ## Success!
# eldo_3k_3 (when current_site == all_sites[3]) uploaded to Pix4D Cloud on 2018-12-29 ## Currently cloud processing
# eldo_4k_1 (when current_site == all_sites[4]) uploaded to Pix4D Cloud on 2018- ## Standing by on upload until I can sort out multiple blocks from cloud processing 
# eldo_4k_2 (when current_site == all_sites[5]) uploaded to Pix4D Cloud on 2018- ## Still trying to get RE project to a single block
# eldo_4k_3 (when current_site == all_sites[6]) uploaded to Pix4D Cloud on 2018-
# eldo_5k_1 (when current_site == all_sites[7]) uploaded to Pix4D Cloud on 2018-
# eldo_5k_2 (when current_site == all_sites[8]) uploaded to Pix4D Cloud on 2018-
# eldo_5k_3 (when current_site == all_sites[9]) uploaded to Pix4D Cloud on 2018-
# sequ_4k_1 (when current_site == all_sites[10]) uploaded to Pix4D Cloud on 2018-

x3_dir <- paste0("data/data_working/", current_site, "_x3_photos")
re_dir <- paste0("data/data_working/", current_site, "_re_photos")

target_dir <- paste0("data/data_working/", current_site, "/", current_site, "_photos")

mission_footprint <- st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geojson"))
dem <- raster::raster("data/features/srtm_30m.tif")
plot_locations <- st_read(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations")))

# Create new directory to store all the photos
dir.create(here::here(paste0("data/data_working/", current_site)))
dir.create(here::here(target_dir))

# Combine the RGB photos, and append an x3 to the filename
x3_photo_filepaths <- list.files(here::here(x3_dir), recursive = TRUE, full.names = TRUE, pattern = ".JPG")
x3_photo_old_names <- list.files(here::here(x3_dir), recursive = TRUE, pattern = ".JPG")
x3_photo_new_filepaths <- paste0(target_dir, "/x3_", str_replace_all(string = x3_photo_old_names, pattern = "/", replacement = "_"))

# Start stan_3k_2 photo manipulation --------------------------------------
# Some photos have folks' houses in them (though I never flew over those houses). For privacy, I'm going to remove those photographs.

# mission_footprint_privacy <- st_read(here::here("data/data_output/site_data/stan_3k_2/stan_3k_2_mission-footprint/stan_3k_2_site-bounds_privacy.geoJSON"))
# 
# x3_exif_data <-
#   data_frame(SourceFile = list.files(here::here(x3_dir), recursive = TRUE, full.names = TRUE, pattern = "JPG")) %>% 
#   dplyr::mutate(x3_photo_old_names = list.files(here::here(x3_dir), recursive = TRUE, pattern = ".JPG")) %>% 
#   dplyr::mutate(DestFile = paste0(here::here(target_dir), "/x3_", str_replace_all(string = x3_photo_old_names, pattern = "/", replacement = "_"))) %>% 
#   dplyr::filter(file.size(SourceFile) > 0) %>% 
#   left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude", "AbsoluteAltitude", "GPSAltitude")), by = "SourceFile") %>% 
#   type_convert() %>%
#   dplyr::filter(!is.na(GPSLongitude), !is.na(GPSLatitude), !is.na(AbsoluteAltitude), !is.na(GPSAltitude)) %>% 
#   dplyr::select(SourceFile, DestFile, GPSLongitude, GPSLatitude, AbsoluteAltitude, GPSAltitude) 
# 
# x3_exif_data <-
#   x3_exif_data %>% 
#   st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) %>% 
#   mutate(agl = GPSAltitude - raster::extract(x = dem, y = ., method = "bilinear") + takeoff_elev_offset) %>% 
#   mutate(mission_photo = st_intersects(x = st_transform(., 3310), y = st_buffer(st_transform(mission_footprint_privacy, 3310), 5), sparse = FALSE)[ , 1])
# 
# dev.off()
# par(mfrow = c(1, 2), mar = rep(0, 4))
# plot_colors <- scale_color_base(value = x3_exif_data$agl, colors = viridis(10))
# plot(x3_exif_data$geometry, type = "n")
# plot(x3_exif_data$geometry[x3_exif_data$mission_photo], pch = 19, add = TRUE, col = plot_colors[x3_exif_data$mission_photo])
# plot(st_transform(plot_locations$geometry, 4326), add = TRUE, col = "red", pch = 19)
# plot(mission_footprint_privacy$geometry, add = TRUE)
# 
# plot(x3_exif_data$geometry, type = "n")
# plot(x3_exif_data$geometry[!x3_exif_data$mission_photo], pch = 19, add = TRUE, col = plot_colors[!x3_exif_data$mission_photo])
# plot(st_transform(plot_locations$geometry, 4326), add = TRUE, col = "red", pch = 19)
# plot(mission_footprint_privacy$geometry, add = TRUE)
# 
# x3_photos_to_copy <-
#   x3_exif_data %>% 
#   filter(mission_photo)
# 
# file.copy(from = x3_photos_to_copy$SourceFile, to = x3_photos_to_copy$DestFile, overwrite = TRUE)


# end stan_3k_2 -----------------------------------------------------------



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
plot(st_transform(plot_locations$geometry, 4326), add = TRUE, col = "red", pch = 19)
plot(mission_footprint$geometry, add = TRUE)

plot(exif_data$geometry, type = "n")
plot(exif_data$geometry[!exif_data$mission_photo], pch = 19, add = TRUE, col = plot_colors[!exif_data$mission_photo])
plot(st_transform(plot_locations$geometry, 4326), add = TRUE, col = "red", pch = 19)
plot(mission_footprint$geometry, add = TRUE)

photos_to_copy <-
  exif_data %>% 
  filter(mission_photo)

photos_to_copy %>% 
  group_by(band_name) %>% 
  summarize(n = n())


# Site eldo_4k_2 photo manipulation ---------------------------------------
# This site was the only one with 2 takeoff points (to keep the sUAS within visual line of sight and radio contact)
# Thus, there is a second takeoff point that should be accounted for
# Hmm, thinking further, there is really a "new" takeoff point for the first photo in each XXXXSET folder, where
# XXXX is an incrementing number left-padded with zeros. If the weather changes, the "pressure altitude" might
# be different for the first photo in each XXXXSET folder, because a new XXXXSET folder is created at the
# start of each new flight. A consideration for further development of this script (and a cleaner way to work
# with imagery that is flexible in its accounting for different takeoff points for each flight) is to grab the exif data
# from all the images within each XXXXET and base the elevation of each of those photos off the elevation
# of the pressure altitude recorded from the first photo within that set.

# re_photo_filepaths <- list.files(re_nonCalibration_dirs, recursive = TRUE, full.names = TRUE, pattern = ".tif")
# idx <- str_detect(string = re_photo_filepaths, pattern = "0010SET")
# tenSET_imgs <- re_photo_filepaths[idx]
# 
# takeoff_point <-
#   tenSET_imgs[1] %>%
#   exif_read(tags = c("GPSLongitude", "GPSLatitude", "GPSAltitude", "DateTimeOriginal", "Yaw", "Pitch", "Roll", "PressureAlt")) %>%
#   type_convert() %>%
#   st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326)
# 
# takeoff_elev <- raster::extract(dem, takeoff_point, method = "bilinear")
# 
# takeoff_elev_offset <- takeoff_elev - takeoff_point$PressureAlt
# 
# exif_data_tenSET <-
#   data_frame(SourceFile = list.files(here::here(re_dir), recursive = TRUE, full.names = TRUE, pattern = ".tif")) %>% 
#   dplyr::mutate(re_photo_old_names = list.files(here::here(re_dir), recursive = TRUE, pattern = ".tif")) %>% 
#   dplyr::filter(str_detect(string = SourceFile, pattern = "0010SET")) %>%
#   dplyr::mutate(DestFile = paste0(here::here(target_dir), "/re_", str_replace_all(string = re_photo_old_names, pattern = "/", replacement = "_"))) %>%
#   dplyr::filter(file.size(SourceFile) > 0) %>%
#   dplyr::filter(!str_detect(string = DestFile, pattern = "re_calibration")) %>%
#   left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude", "PressureAlt")), by = "SourceFile") %>%
#   type_convert() %>%
#   dplyr::filter(!is.na(GPSLongitude), !is.na(GPSLatitude), !is.na(PressureAlt)) %>%
#   dplyr::select(SourceFile, DestFile, GPSLongitude, GPSLatitude, PressureAlt)
# 
# exif_data_tenSET <-
#   exif_data_tenSET %>%
#   st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) %>%
#   mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>%
#   mutate(band_name = case_when(band_num == 1 ~ "blue",
#                                band_num == 2 ~ "green",
#                                band_num == 3 ~ "red",
#                                band_num == 4 ~ "nir",
#                                band_num == 5 ~ "re")) %>%
#   mutate(agl = PressureAlt - raster::extract(x = dem, y = ., method = "bilinear") + takeoff_elev_offset) %>%
#   mutate(suas_at_altitude = ifelse(agl > 90, yes = TRUE, no = FALSE)) %>%
#   mutate(suas_within_footprint = st_intersects(x = st_transform(., 3310), y = st_buffer(st_transform(mission_footprint, 3310), 5), sparse = FALSE)[ , 1]) %>%
#   mutate(mission_photo = suas_at_altitude & suas_within_footprint)
# 
# exif_data_tenSET %>%
#   group_by(band_name) %>%
#   summarize(n = n())
# 
# photos_to_copy_tenSET <-
#   exif_data_tenSET %>%
#   filter(mission_photo)
# 
# photos_to_copy_tenSET %>%
#   group_by(band_name) %>%
#   summarize(n = n())
# 
# photos_to_copy <-
#   photos_to_copy %>%
#   rbind(photos_to_copy_tenSET)


# End of eldo_4k_2 processing ---------------------------------------------


# start of sequ_4k_2 processing -------------------------------------------
# I revisited this site to capture imagery in the middle of the survey area because the RedEdge failed to fire the first.

# Note that all the photos are copied over to the same folder, but there should be TWO separate RE projects that get merged in order to use the reflectance panel images for the same day that the rest of the flight occurred.
# re_calibration_IMG_0000_1 through _5 are for the 1st RE subproject
# All re_card-1 photos are for subproject 1
# re_calibration_IMG_0002_1 through _5 are for the 2nd RE subproject
# All re_card-2 photos are for subproject 2

# end of sequ_4k_2 processing ---------------------------------------------


# sequ_4k_1 photo manipulation --------------------------------------------
# This site had a problem with the RedEdge GPS over the east ~40% of the survey area such that the time, latitude, longitude, and
# altitude weren't recorded for 40% of the photos. I marked ground control points on the already-created orthomosaic from processing
# just the RGB imagery and saved those 14 points to a shapefile. I will take the coordinates for those points and treat them as
# the true values for the conspicuous objects in the images. I'll mark the location of all of the 14 points in the RGB and in the
# multispectral photos following instructions here: https://support.pix4d.com/hc/en-us/articles/360000276046-How-to-import-and-mark-ground-control-points-GCPs-
# First I'll run step 1 on the RGB images, then I'll mark the 14 GCP's using the rayCloud option
# Next, I'll create a multispectral project and mark the 14 GCP's *BEFORE* running step 1, then I'll run step 1 and hopefully the 
# multispectral project will be able to use all of the images across the whole survey area.

# We're going to want to copy RedEdge photos that don't have the right EXIF metadata, so we'll have to modify our workflow a bit:
# Add a column for whether the photo was georeferenced properly or not. We'll include all the photos that we know we can say are mission photos
# and then all of the photos that weren't georeferenced. That way, we'll have done some of the subsetting of non-mission photos to 
# reduce the load of having to do so when there's no geolocation information
# 
# exif_data <-
#   exif_data %>% 
#   st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326, remove = FALSE) %>% 
#   mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>% 
#   mutate(band_name = case_when(band_num == 1 ~ "blue",
#                                band_num == 2 ~ "green",
#                                band_num == 3 ~ "red",
#                                band_num == 4 ~ "nir",
#                                band_num == 5 ~ "re")) %>% 
#   mutate(agl = PressureAlt - raster::extract(x = dem, y = ., method = "bilinear") + takeoff_elev_offset) %>% 
#   mutate(suas_at_altitude = ifelse(agl > 90, yes = TRUE, no = FALSE)) %>% 
#   mutate(suas_within_footprint = st_intersects(x = st_transform(., 3310), y = st_buffer(st_transform(mission_footprint, 3310), 5), sparse = FALSE)[ , 1]) %>% 
#   mutate(mission_photo = suas_at_altitude & suas_within_footprint) %>% 
#   mutate(georeferenced = ifelse(GPSLongitude == 0, yes = FALSE, no = TRUE))
# 
# exif_data %>% 
#   group_by(band_name) %>% 
#   summarize(n = n())
# 
# # ggplot(exif_data) +
# #   geom_sf(aes(color = agl)) +
# #   facet_grid(~ mission_photo) +
# #   scale_color_viridis_c()
# 
# photos_to_copy1 <-
#   exif_data %>% 
#   filter(mission_photo)
# 
# photos_to_copy2 <-
#   exif_data %>% 
#   filter(!georeferenced)
# 
# photos_to_copy1 %>% 
#   group_by(band_name) %>% 
#   summarize(n = n())
# 
# photos_to_copy2 %>% 
#   group_by(band_name) %>% 
#   summarize(n = n())
# 
# photos_to_copy <- rbind(photos_to_copy1, photos_to_copy2)
# 
# photos_to_copy %>% 
#   group_by(band_name) %>% 
#   summarize(n = n())

# We also can take the Ground Control Points shapefile that we marked in QGIS and turn it into a format that Pix4D can read
# (comma delimited text file with the point label, then the x coordinate, then the y coordinate)

# sequ_4k_1_gcp <- 
#   st_read(here::here("data/data_working/sequ_4k_1/sequ_4k_1_gcp/sequ_4k_1_gcp.shp")) %>% 
#   dplyr::mutate(id = paste0("gcp_", id))
# 
# Write the file as a CSV file and put the geometry column into two separate columns (one for X
# and one for Y); This also adds a column at the end-- not sure why
# st_write(sequ_4k_1_gcp, here::here("data/data_working/sequ_4k_1/sequ_4k_1_gcp/sequ_4k_1_gcp.csv"), layer_options = "GEOMETRY=AS_XY", delete_dsn = TRUE)
# 
# To remove the last column that was weirdly added by st_write(), we read the file back in
# as a data frame and just select the columns that we want.
# sequ_4k_1_gcp <- 
#   read.csv(here::here("data/data_working/sequ_4k_1/sequ_4k_1_gcp/sequ_4k_1_gcp.csv")) %>% 
#   dplyr::select(id, X, Y)
# 
# Then write it out as a file with no headers
# write.table(x = sequ_4k_1_gcp, file = here::here("data/data_working/sequ_4k_1/sequ_4k_1_gcp/sequ_4k_1_gcp.csv"), sep = ",", row.names = FALSE, col.names = FALSE)

# End of sequ_4k_1 photo manipulation --------------------------------------


# sier_4k_1 photo manipulation --------------------------------------------
# This is the Goat Mountain site
# Similar to sequ_4k_1 site, the RedEdge didn't record GPS for some of the photos. We will manually add the GPS locations
# in the same way that we did for sequ_4k_1

# exif_data <-
#   exif_data %>%
#   st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326, remove = FALSE) %>%
#   mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>%
#   mutate(band_name = case_when(band_num == 1 ~ "blue",
#                                band_num == 2 ~ "green",
#                                band_num == 3 ~ "red",
#                                band_num == 4 ~ "nir",
#                                band_num == 5 ~ "re")) %>%
#   mutate(agl = PressureAlt - raster::extract(x = dem, y = ., method = "bilinear") + takeoff_elev_offset) %>%
#   mutate(suas_at_altitude = ifelse(agl > 90, yes = TRUE, no = FALSE)) %>%
#   mutate(suas_within_footprint = st_intersects(x = st_transform(., 3310), y = st_buffer(st_transform(mission_footprint, 3310), 5), sparse = FALSE)[ , 1]) %>%
#   mutate(mission_photo = suas_at_altitude & suas_within_footprint) %>%
#   mutate(georeferenced = ifelse(GPSLongitude == 0, yes = FALSE, no = TRUE))
# 
# exif_data %>%
#   group_by(band_name) %>%
#   summarize(n = n())
# 
# # ggplot(exif_data) +
# #   geom_sf(aes(color = agl)) +
# #   facet_grid(~ mission_photo) +
# #   scale_color_viridis_c()
# 
# photos_to_copy1 <-
#   exif_data %>%
#   filter(mission_photo)
# 
# photos_to_copy2 <-
#   exif_data %>%
#   filter(!georeferenced)
# 
# photos_to_copy1 %>%
#   group_by(band_name) %>%
#   summarize(n = n())
# 
# photos_to_copy2 %>%
#   group_by(band_name) %>%
#   summarize(n = n())
# 
# photos_to_copy <- rbind(photos_to_copy1, photos_to_copy2)
# 
# photos_to_copy %>%
#   group_by(band_name) %>%
#   summarize(n = n())

# end sier_4k_1 photo manipulation ----------------------------------------





# Copy the photos from the working folder to the final photos folder using the new names
file.copy(from = photos_to_copy$SourceFile, to = photos_to_copy$DestFile, overwrite = TRUE)


# Add the calibration panel photos to the final photos directory
calibration_SourceDir <- paste0(re_dir, "/calibration")
calibration_SourceFile <- list.files(here::here(calibration_SourceDir), full.names = TRUE)
calibration_SourceNames <- list.files(here::here(calibration_SourceDir))
calibration_DestFile <- here::here(paste0(target_dir, "/re_calibration_", calibration_SourceNames))

file.copy(from = calibration_SourceFile, to = calibration_DestFile)

# After manually curating the photos, get the list of all photos used for analysis
final_photo_list <- list.files(here::here(target_dir))

metadata <- data.frame(site = current_site, processed_photos = final_photo_list)

# Add the file size to the metadata
metadata <-
  metadata %>% 
  mutate(file_size_MB = file.size(here::here(paste0(target_dir, "/", metadata$processed_photos))) / 1e6)

write_csv(metadata, here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_processed-photos.csv")))

# These lines create the folders that will be used to house the intermediate Pix4D projects before they get merged
dir.create(here::here(paste0("data/data_working/", current_site, "/", current_site, "_x3")))
dir.create(here::here(paste0("data/data_working/", current_site, "/", current_site, "_re")))

