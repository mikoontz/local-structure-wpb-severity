# Purpose: Write EXIF data to RedEdge photos that have no geolocation. Use information from the RGB flights
# and make a best guess as to the lat/lon of the photos with missing data. Then run it through Pix4D and
# use manual tie points/ground control points to clean everything up.

library(tidyverse)
library(sf)
library(exiftoolr)

# sequ_4k_1 ---------------------------------------------------------------

x3_photo_names <- list.files("data/data_working/sequ_4k_1/sequ_4k_1_photos/", pattern = "x3")
x3_photo_files <- list.files("data/data_working/sequ_4k_1/sequ_4k_1_photos/", pattern = "x3", full.names = TRUE)
x3_idx <- x3_photo_names >= "x3_100MEDIA_DJI_0801.JPG"

helper_photos <- 
  data_frame(x3_photo_names = x3_photo_names[x3_idx],
             SourceFile = x3_photo_files[x3_idx]) %>% 
  left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude")), 
            by = "SourceFile") %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) %>% 
  st_transform(32611)

helper_photos <-
  helper_photos %>% 
  bind_cols(as.data.frame(st_coordinates(.)) %>% setNames(c("lon", "lat")))

x3_north_endpoints <- c("x3_100MEDIA_DJI_0801.JPG",
                        "x3_100MEDIA_DJI_0868.JPG",
                        "x3_100MEDIA_DJI_0870.JPG",
                        "x3_100MEDIA_DJI_0940.JPG",
                        "x3_100MEDIA_DJI_0942.JPG",
                        "x3_101MEDIA_DJI_0013.JPG",
                        "x3_101MEDIA_DJI_0015.JPG",
                        "x3_101MEDIA_DJI_0087.JPG",
                        "x3_101MEDIA_DJI_0088.JPG",
                        "x3_101MEDIA_DJI_0159.JPG",
                        "x3_101MEDIA_DJI_0161.JPG",
                        "x3_101MEDIA_DJI_0231.JPG",
                        "x3_101MEDIA_DJI_0233.JPG",
                        "x3_101MEDIA_DJI_0302.JPG")

x3_south_endpoints <- c("x3_100MEDIA_DJI_0832.JPG",
                        "x3_100MEDIA_DJI_0834.JPG",
                        "x3_100MEDIA_DJI_0904.JPG",
                        "x3_100MEDIA_DJI_0906.JPG",
                        "x3_100MEDIA_DJI_0976.JPG",
                        "x3_100MEDIA_DJI_0978.JPG",
                        "x3_101MEDIA_DJI_0049.JPG",
                        "x3_101MEDIA_DJI_0051.JPG",
                        "x3_101MEDIA_DJI_0123.JPG",
                        "x3_101MEDIA_DJI_0125.JPG",
                        "x3_101MEDIA_DJI_0195.JPG",
                        "x3_101MEDIA_DJI_0197.JPG",
                        "x3_101MEDIA_DJI_0268.JPG",
                        "x3_101MEDIA_DJI_0268.JPG")

flight_details <-
  data_frame(transect = 23:36, 
             north = x3_north_endpoints, 
             south = x3_south_endpoints) %>% 
  tidyr::gather(key = "boundary", value = "x3_photo_names", -transect) %>% 
  dplyr::left_join(helper_photos)
  
re_photo_names <- list.files("data/data_working/sequ_4k_1/sequ_4k_1_photos/", pattern = "re")
re_photo_files <- list.files("data/data_working/sequ_4k_1/sequ_4k_1_photos/", pattern = "re", full.names = TRUE)

re_photos <-
  data_frame(re_photo_name = re_photo_names,
             SourceFile = re_photo_files) %>% 
  dplyr::filter(file.size(SourceFile) > 0) %>% 
  left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude", "PressureAlt")), by = "SourceFile") %>% 
  type_convert() %>%
  dplyr::filter(!is.na(GPSLongitude), !is.na(GPSLatitude), !is.na(PressureAlt)) %>% 
  dplyr::select(re_photo_name, SourceFile, GPSLongitude, GPSLatitude, PressureAlt) %>% 
  dplyr::mutate(core_name = substr(re_photo_name, start = 1, stop = nchar(re_photo_name) - 6)) %>% 
  dplyr::mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>% 
  mutate(band_name = case_when(band_num == 1 ~ "blue",
                               band_num == 2 ~ "green",
                               band_num == 3 ~ "red",
                               band_num == 4 ~ "nir",
                               band_num == 5 ~ "re"))

re_photos <-
  re_photos %>% 
  dplyr::mutate(georeferenced = ifelse(GPSLongitude == 0, yes = FALSE, no = TRUE)) %>% 
  dplyr::mutate(transect = NA)


re_north_endpoints <- c("re_card-2_0001SET_000_IMG_0036",
                        "re_card-2_0001SET_000_IMG_0138",
                        "re_card-2_0001SET_000_IMG_0144",
                        "re_card-2_0001SET_001_IMG_0046",
                        "re_card-2_0001SET_001_IMG_0052",
                        "re_card-2_0001SET_001_IMG_0154",
                        "re_card-2_0001SET_001_IMG_0159",
                        "re_card-2_0002SET_000_IMG_0050",
                        "re_card-2_0002SET_000_IMG_0056",
                        "re_card-2_0002SET_000_IMG_0158",
                        "re_card-2_0002SET_000_IMG_0163",
                        "re_card-2_0002SET_001_IMG_0067",
                        "re_card-2_0002SET_001_IMG_0073",
                        "re_card-2_0002SET_001_IMG_0173")

re_south_endpoints <- c("re_card-2_0001SET_000_IMG_0084",
                        "re_card-2_0001SET_000_IMG_0090",
                        "re_card-2_0001SET_000_IMG_0192",
                        "re_card-2_0001SET_000_IMG_0197",
                        "re_card-2_0001SET_001_IMG_0100",
                        "re_card-2_0001SET_001_IMG_0106",
                        "re_card-2_0001SET_002_IMG_0008",
                        "re_card-2_0001SET_002_IMG_0014",
                        "re_card-2_0002SET_000_IMG_0104",
                        "re_card-2_0002SET_000_IMG_0109",
                        "re_card-2_0002SET_001_IMG_0012",
                        "re_card-2_0002SET_001_IMG_0019",
                        "re_card-2_0002SET_001_IMG_0121",
                        "re_card-2_0002SET_001_IMG_0126")


# sier_4k_1 ---------------------------------------------------------------


x3_photo_names <- list.files("data/data_working/sier_4k_1/sier_4k_1_photos/", pattern = "x3")
x3_photo_files <- list.files("data/data_working/sier_4k_1/sier_4k_1_photos/", pattern = "x3", full.names = TRUE)

x3_idx <- (x3_photo_names >= "x3_100MEDIA_DJI_0487.JPG" & x3_photo_names <= "x3_100MEDIA_DJI_0702.JPG")

helper_photos <- 
  data_frame(x3_photo_names = x3_photo_names[x3_idx],
             SourceFile = x3_photo_files[x3_idx]) %>% 
  left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude")), 
            by = "SourceFile") %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) %>% 
  st_transform(32611)

helper_photos <-
  helper_photos %>% 
  bind_cols(as.data.frame(st_coordinates(.)) %>% setNames(c("lon", "lat")))

x3_north_endpoints <- c("x3_100MEDIA_DJI_0503.JPG",
                        "x3_100MEDIA_DJI_0505.JPG",
                        "x3_100MEDIA_DJI_0583.JPG",
                        "x3_100MEDIA_DJI_0584.JPG",
                        "x3_100MEDIA_DJI_0662.JPG",
                        "x3_100MEDIA_DJI_0664.JPG")

x3_south_endpoints <- c("x3_100MEDIA_DJI_0487.JPG",
                        "x3_100MEDIA_DJI_0543.JPG",
                        "x3_100MEDIA_DJI_0545.JPG",
                        "x3_100MEDIA_DJI_0623.JPG",
                        "x3_100MEDIA_DJI_0624.JPG",
                        "x3_100MEDIA_DJI_0702.JPG")

flight_details <-
  data_frame(transect = 6:11, 
             north = x3_north_endpoints, 
             south = x3_south_endpoints) %>% 
  tidyr::gather(key = "boundary", value = "x3_photo_names", -transect) %>% 
  dplyr::left_join(helper_photos)

re_photo_names <- list.files("data/data_working/sier_4k_1/sier_4k_1_photos/", pattern = "re")
re_photo_files <- list.files("data/data_working/sier_4k_1/sier_4k_1_photos/", pattern = "re", full.names = TRUE)

re_photos <-
  data_frame(re_photo_name = re_photo_names,
             SourceFile = re_photo_files) %>% 
  dplyr::filter(file.size(SourceFile) > 0) %>% 
  left_join(exif_read(.$SourceFile, tags = c("GPSLongitude", "GPSLatitude", "PressureAlt")), by = "SourceFile") %>% 
  type_convert() %>%
  dplyr::filter(!is.na(GPSLongitude), !is.na(GPSLatitude), !is.na(PressureAlt)) %>% 
  dplyr::select(re_photo_name, SourceFile, GPSLongitude, GPSLatitude, PressureAlt) %>% 
  dplyr::mutate(core_name = substr(re_photo_name, start = 1, stop = nchar(re_photo_name) - 6)) %>% 
  dplyr::mutate(band_num = str_sub(SourceFile, start = (nchar(SourceFile) - 4), end = (nchar(SourceFile) - 4))) %>% 
  mutate(band_name = case_when(band_num == 1 ~ "blue",
                               band_num == 2 ~ "green",
                               band_num == 3 ~ "red",
                               band_num == 4 ~ "nir",
                               band_num == 5 ~ "re"))

re_photos <-
  re_photos %>% 
  dplyr::mutate(georeferenced = ifelse(GPSLongitude == 0, yes = FALSE, no = TRUE)) %>% 
  dplyr::mutate(transect = NA)


re_north_endpoints <- c("re_card-1_0006SET_000_IMG_0063",
                        "re_card-1_0006SET_000_IMG_0070",
                        "re_card-1_0006SET_000_IMG_0183",
                        "re_card-1_0006SET_000_IMG_0191",
                        "re_card-1_0006SET_001_IMG_0102",
                        "re_card-1_0006SET_001_IMG_0108")

re_south_endpoints <- c("re_card-1_0006SET_000_IMG_0042",
                        "re_card-1_0006SET_000_IMG_0123",
                        "re_card-1_0006SET_000_IMG_0130",
                        "re_card-1_0006SET_001_IMG_0042",
                        "re_card-1_0006SET_001_IMG_0047",
                        "re_card-1_0006SET_001_IMG_0160")

transect_bounds <- 
  data_frame(transect = 6:11,
             north = re_north_endpoints,
             south = re_south_endpoints)

# Processing invariant of site --------------------------------------------


# re_photos_copy <- re_photos  
# re_photos <- re_photos_copy

# # Filter out the flight change photos
# re_photos <-
#   re_photos %>% 
#   dplyr::filter(!(core_name > "re_card-2_0001SET_002_IMG_0055" & core_name < "re_card-2_0002SET_000_IMG_0050"))

for (i in 1:nrow(transect_bounds)) {
  # For odd numbered transects (drone flew north to south)
  if (transect_bounds$transect[i] %% 2 == 1) {
    current_transect_idx <- re_photos$core_name >= transect_bounds$north[i] & re_photos$core_name <= transect_bounds$south[i]
    re_photos$transect[current_transect_idx] <- transect_bounds$transect[i]
    
    photos_in_this_transect <- sort(unique(re_photos$core_name[current_transect_idx]))
    num_of_photos_in_this_transect <- length(photos_in_this_transect)
    
    frac_seq <- 0:(num_of_photos_in_this_transect - 1) / (num_of_photos_in_this_transect - 1)
    
    x_north <- flight_details$lon[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "north"]
    x_south <- flight_details$lon[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "south"]
    
    y_north <- flight_details$lat[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "north"]
    y_south <- flight_details$lat[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "south"]
    
    x <- x_north + (x_south - x_north) * frac_seq
    y <- y_north + (y_south - y_north) * frac_seq
    
    new_coords <-
      data_frame(core_name = photos_in_this_transect, x = x, y = y)
    
    # Remember that all 5 of the bands from each shutter actuation need to be tagged with the
    # new geolocation in the EXIF metadata
    for (j in 1:nrow(new_coords)) {
      re_photos$GPSLongitude[which(re_photos$core_name %in% new_coords$core_name[j])] <- new_coords$x[j]
      re_photos$GPSLatitude[which(re_photos$core_name %in% new_coords$core_name[j])] <- new_coords$y[j]
    }
    
  } 
  # For even numbered transects (drone flew south to north)
  if (transect_bounds$transect[i] %% 2 == 0) {
    current_transect_idx <- re_photos$core_name >= transect_bounds$south[i] & re_photos$core_name <= transect_bounds$north[i]
    re_photos$transect[current_transect_idx] <- transect_bounds$transect[i]
    
    photos_in_this_transect <- sort(unique(re_photos$core_name[current_transect_idx]))
    num_of_photos_in_this_transect <- length(photos_in_this_transect)
    
    frac_seq <- 0:(num_of_photos_in_this_transect - 1) / (num_of_photos_in_this_transect - 1)
    
    x_north <- flight_details$lon[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "north"]
    x_south <- flight_details$lon[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "south"]
    
    y_north <- flight_details$lat[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "north"]
    y_south <- flight_details$lat[flight_details$transect == transect_bounds$transect[i] & flight_details$boundary == "south"]
    
    # The even transects fly south to north, so we reverse the order of the vector of new x and y
    # coordinates
    x <- rev(x_north + (x_south - x_north) * frac_seq)
    y <- rev(y_north + (y_south - y_north) * frac_seq)
    
    new_coords <-
      data_frame(core_name = photos_in_this_transect, x = x, y = y)
    
    # Remember that all 5 of the bands from each shutter actuation need to be tagged with the
    # new geolocation in the EXIF metadata
    for (j in 1:nrow(new_coords)) {
      re_photos$GPSLongitude[which(re_photos$core_name %in% new_coords$core_name[j])] <- new_coords$x[j]
      re_photos$GPSLatitude[which(re_photos$core_name %in% new_coords$core_name[j])] <- new_coords$y[j]
    }
    
  }
  print(paste0("Transect ", transect_bounds$transect[i], " complete."))
  
}

# re_photos_copy %>% 
#   dplyr::filter(GPSLongitude != 0)

re_photos_proper_geolocation <- 
  re_photos %>% 
  dplyr::filter(georeferenced) %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 4326) %>% 
  cbind(st_coordinates(.))


re_photos_hacked_geolocation <-
  re_photos %>% 
  dplyr::filter(GPSLongitude != 0) %>% 
  dplyr::filter(!georeferenced) %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 32611) %>% 
  st_transform(4326) %>% 
  cbind(st_coordinates(.))

re_photos_sf <-
  rbind(re_photos_proper_geolocation, re_photos_hacked_geolocation)

plot(st_geometry(re_photos_sf))

re_photos_copy_sf <- 
  re_photos_copy %>%
  dplyr::filter(GPSLongitude != 0) %>% 
  st_as_sf(coords = c("GPSLongitude", "GPSLatitude"), crs = 32611, remove = FALSE)


re_photos_df <-
  re_photos_sf %>% 
  as.data.frame() %>% 
  dplyr::select(re_photo_name, X, Y, PressureAlt)

re_photos_df

ggplot(re_photos_df, aes(x = X, y = Y, col = PressureAlt)) +
  geom_point() +
  scale_color_viridis_c() +
  coord_equal()

write.table(x = re_photos_df, file = here::here("data/data_working/sier_4k_1/sier_4k_1_re_geolocations.txt"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)



# From this stackOverflow answer: https://stackoverflow.com/questions/20923884/find-equidistant-points-between-two-coordinates
# frac is the sequence of fractional distances along the line between the north endpoint and the south endpoint. A sequence from (1 to the total number of photos between those points) - 1 divided by (total number of photos between those points - 1)

# x <- x_north + (x_south - x_north) * frac
# y <- y_north + (y_south - y_north) * frac