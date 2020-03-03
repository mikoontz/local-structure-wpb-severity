Analysis of forest structure data in yellow pine/mixed conifer forests of the Sierra Nevada Mountain Range and severity of western pine beetle in those forests.

Forest sites were chosen to have >40% ponderosa pine tree by basal area and >10% mortality of ponderosa pine tree by basal area.

The general workflow is:

1) Fly a grid pattern over a forested site using a DJI drone and the Map Pilot for iOS iPad app.
2) Take images every 2 seconds during the flight with lots of overlap between images.
3) Use photogrammetry (a.k.a. "structure from motion") software to turn the overlapping 2D images into a 3D point cloud, a 2D surface model, and a 2D orthomosaic. We use the Pix4D software and adjust processing parameters to best account for dense vegetation.
4) Use a point cloud manipulation software to classify the 3D point cloud into "ground" and "non-ground" points using a "cloth simulator filter". This filter also has the effect of using the classified "ground" points to generate a 2D continuous digital terrain model underneath the vegetation. We use CloudCompare for this.
5) Generate a "canopy height model" by subtracting away the CloudCompare generated terrain model from the Pix4D generated 2D surface model to leave just the surface representing the vegetation.
6) We use a variable window local maximum filter to identify the tree tops in the canopy height model. (In `R` using the `ForestTools` package)
7) We use a marker controlled watershed segmentation algorithm to identify the pixels in the orthomosaic that are associated with the individual trees. (Also in `R` and using the `ForestTools` package)
8) We extract the red, blue, and green values from the pixels of the orthomosaic representing each tree. We calculate the red-green index (RGI) as R / G and the red-green vegetation index as (R - G) / (R + G). We take the average and the 75th percentile of these 5 metrics across all pixels belonging to each tree.
9) Using QGIS, we overlay the polygons of the segmented crowns on top of the the RGB orthomosaic for the eldo_3k_1, eldo_5k_1, sequ_4k_1, and sequ_6k_1 sites. We made a copy of the segmented crowns polygons and added an attribute field called "live". We then panned around the orthomosaic and manually changed the "live" field for ~100 trees per site to 0 (for a dead tree) or 1 (for a live tree).
10) Using the ~400 manually classified crown segments and the orthomosaic for each site, we extract the R, G, B, RGI (RGI = R/G per pixel), and GBI (GBI = G/B per pixel) data from each pixel and take the mean value of each spectral signature across all the pixels within each crown polygon. Using the mean R, G, B, RGI, and GBI values per crown polygon as covariates, We fit a generalized linear model using a binomial family and a logit link to predict the probability that the tree represented by the crown segment was alive or dead.
11) Next, we extracted the R, G, B, RGI, and GBI means for all the crown polygons from each site's orthomosaic. We used the `velox` package in `R` to do this massive pixel value extraction from the raster files.
12) We used the logistic regression model that we fit above to predict the probability that each tree in the study (across all sites) was alive or dead. 

Some photos of the processing:

This is one of the data products from the Pix4D processing of the 2D aerial imagery. Here, the point cloud is visualized using the CloudCompare software.

![Original point cloud of a forested site, generated using Pix4D](figures/L1_eldo_3k_3_point-cloud_rgb-cloudcompare.png)

This is a portion of one of the other primary data products of the Pix4D processing-- an orthomosaic showing the top-down view of all objects in the scene with the spatial relationships amongst them preserved. This image is suitable for making measurments between the trees.

![An orthomosaic (a top-down view of all objects) of a forested site, generated using Pix4D](figures/L1_eldo_3k_3_ortho_rgb.png)

This is a view of just the points from the point cloud classified as "ground" after using the cloth simulator function in CloudCompare. Essentially, the trees from the previous image have been "deleted" off the landscape.

![Points from the point cloud classified as ground](figures/eldo_4k_2_no-trees.png)

We can use the same "cloth" that helped classify ground versus non ground to interpolate the terrain underneath the trees, even where the trees may have obscured the ground.
In 2 dimensions, the digital terrain model looks like this:

![2 dimensional digital terrain model](figures/L2_eldo_3k_3_dtm.png)

The Digital Surface Model is an output from Pix4D that represents the height of the surface, which includes the elevation of the ground plus the height of the vegetation.

![Digital Surface Model (DSM) of a forested site](figures/L1_eldo_3k_3_dsm.png)

By subtracting the Digital Terrain Model (DTM) away from the Digital Surface Model (DSM), we get a representation of the heights of all the vegetation-- a Canopy Height Model (CHM).

![A Canopy Height Model representing the height above the ground for each tree in the scene.](figures/L2_eldo_3k_3_chm.png)

The "variable window filter" algorithm detects tree tops as local maxima using the Canopy Height Model.

![Tree tops detected using the variable window filter algorithm on the Canopy Height Model.](figures/L3a_eldo_3k_3_ttops_cropped.png)

The tree top locations and the Canopy Height Model are then used to determine the spatial extent of each tree's crown using a marker controlled watershed segmentation algorithm.

![Crown spatial extent, or "segments", detected using marker controlled watershed segmentation algorithm on the Canopy Height Model and using the tree top locations.](figures/L3a_eldo_3k_3_crowns_cropped.png)


The classification of all trees in the study (whether they were alive or dead) based on the extracted pixel values from the orthomosaic for each crown segment generates a "forest stem map" for each site.

![A stem map showing the prediction of which trees are alive and which ones are dead based on extracted pixel values for each tree crown and the boosted logistic regression model.](figures/L3b_eldo_3k_3_live-dead.png)

We can use a different model to classify the live trees to species.

![A stem map showing the species of each tree based on extracted pixel values for each tree crown and a regularized discriminant analysis model. Ponderosa pine trees are the host to the western pine beetle, and all other species are non-hosts.](figures/L3b_eldo_3k_3_host-nonhost.png)