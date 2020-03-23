# Relocate the Pix4D processing reports, the DSMs, and the dense point clouds
# to their proper L1 home from deep within the default Pix4D file structure.
# Essentially, we are relocating all of the Level 1 products except the orthomosaics

# Recall that the processing happened in a separate folder on a solid-state hard drive with the
# following file structure (see also: workflow/07_integrate-x3-and-re-photos-for-one-project.txt)

# H:/ (a solid-state local drive)
# |---- my_project/
#       |---- readme.md
#       |---- my_project.Rproj
#       |---- data/
#             |---- data_drone/
#                   |---- L0/
#                   |---- L1/
#                   |---- L2/
#                   |---- L3a/
#                   |---- L3b/
#                   |---- L4/
#             |---- data_working/
#                   |---- "current_site".p4d
#                   |---- "current_site"_re_photos/
#                         |---- calibration/
#                   |---- "current_site"_x3_photos/
#                   |____ "current_site"/
#                         |---- "current_site"_photos/
#                         |---- "current_site"_re/
#                               |---- 1_initial/
#                               |---- 2_densification/
#                               |---- 3_dsm_ortho/
#                               |____ 4_index/
#                         |---- "current_site"_x3/
#                               |---- 1_initial/
#                               |---- 2_densification/
#                               |____ 3_dsm_ortho/
#                         |---- "current_site"_re.p4d  
#                         |---- "current_site"_x3.p4d
#                         |---- 1_initial/
#                         |---- 2_densification/
#                         |---- 3_dsm_ortho/
#                         |____ 4_index/      


source("workflow/01_make-processing-checklist.R")
all_sites <- sites_checklist$site

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# rehome pix4d processing reports -----------------------------------------

report_files <-
  tibble(site = all_sites,
         oldfilename = ifelse(site %in% merged_sites, 
                              yes = paste0(site, "_report.pdf"),
                              no = paste0(site, "_re_report.pdf")),
         oldpath = ifelse(site %in% merged_sites, 
                          yes = paste0("data/data_working/", site, "/1_initial/report/", oldfilename),
                          no = paste0("data/data_working/", site, "/", site, "_re/1_initial/report/", oldfilename)),
         newpath = paste0("data/data_drone/L1/pix4d-reports/", site, "_report.pdf"))

if(!dir.exists("data/data_drone/L1/pix4d-reports")) {
  dir.create("data/data_drone/L1/pix4d-reports", recursive = TRUE)
}

file.copy(from = report_files$oldpath, to = report_files$newpath)

# rehome dsm --------------------------------------------------------------

dsm_files <-
  tibble(site = all_sites,
         oldpath = ifelse(site %in% merged_sites, 
                          yes = paste0("data/data_working/", site, "/3_dsm_ortho/1_dsm/", site, "_dsm.tif"),
                          no = paste0("data/data_working/", site, "/", site, "_re/3_dsm_ortho/1_dsm/", site, "_re_dsm.tif")),
         newpath = paste0("data/data_drone/L1/dsm/", site, "_dsm.tif"))

if(!dir.exists("data/data_drone/L1/dsm")) {
  dir.create("data/data_drone/L1/dsm", recursive = TRUE)
}

for (i in 1:nrow(dsm_files)) {
  r <- raster::raster(dsm_files$oldpath[i])
  raster::writeRaster(x = r, filename = dsm_files$newpath[i])
  print(paste(dsm_files$site[i], "complete..."))
}


# dense point cloud files -------------------------------------------------

dense_point_cloud_files <-
  tibble(site = all_sites,
         oldpath = ifelse(site %in% merged_sites, 
                          yes = paste0("data/data_working/", site, "/2_densification/point_cloud/", site, "_densified_point_cloud.las"),
                          no = paste0("data/data_working/", site, "/", site, "_re/2_densification/point_cloud", site, "_re_Green_densified_point_cloud.las")),
         newpath = paste0("data/data_drone/L1/dense-point-cloud/", site, "_dense-point-cloud.las"))

if(!dir.exists("data/data_drone/L1/dense-point-cloud")) {
  dir.create("data/data_drone/L1/dense-point-cloud", recursive = TRUE)
}

file.copy(from = dense_point_cloud_files$oldpath, to = dense_point_cloud_files$newpath)