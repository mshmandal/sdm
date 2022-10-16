#############################################################################
# Author: Mohammad Shamim Hasan Mandal
# Date  : July 13 2022
# Name  : script.R
# Obj.  : Demo project: setup R project and make Maxent model
# Update: October 16 2022
#############################################################################
# SETUP R PROJECT

#start fresh
gc()            # clear memory
rm(list=ls())   # clear R session

# set work directory
getwd()
# setwd("D:/github/sdm") 

# create new directory to save data and results
dir.create("./data")
dir.create("./result")

#-------------------------------------------------------------------------------
# STEP 1: Install required libraries
#-------------------------------------------------------------------------------
# install.packages("raster",dependencies = T)
# install.packages("dismo",dependencies = T)
# install.packages("SDMtune",dependencies = T)

library(raster)
library(SDMtune)
library(dismo)

#-------------------------------------------------------------------------------
# Step 2 : Prepare data for model
#-------------------------------------------------------------------------------


## 2.1 GET BD SHAPEFILE and WORLCLIM data
ext = c(91,93,19,23)  # c(xmin,xmax,ymin,ymax)
bd  = getData('GADM',country='Bangladesh',level=2,path = "./data",download = F)

# to view
plot(bd,axes=T)

# crop the country shape file to study area extent
bd  = crop(bd,extent(ext))




## 2.2 Raster/Image data: GET WORLD CLIM DATA 

# browseURL("https://www.worldclim.org/data/bioclim.html") 
bio19 = getData('worldclim','bio',res=0.5,lon=90, lat=22,path = "./data",download=F)
bio19 = crop(bio19,extent(ext)) 

# HOW TO PLOT RASTER AND SHAPEFILE
plot(bio19[[1]],main="Annual mean temperature")
plot(bd,add=T)

# first four layers
plot(bio19[[1:4]])

# saving the results
pdf("./result/bio1-4.pdf",width = 6,height = 4)
plot(bio19[[1:4]])
dev.off()

# SAVING AND RE-READING
# writeRaster(bio19,filename = "./data/bio19_bd.grd",overwrite=T)
# bio19 = stack("./data/bio19_bd.grd")
 
#-------------------------------------------------------------------------------
# STEP 4: Prepare data for Maxent and run Maxent
#-------------------------------------------------------------------------------

# GET RANDOM POINTS
occ = dismo::randomPoints(n=20,bio19[[1]])
bg = dismo::randomPoints(n=1000,bio19[[1]])


# Prepare data
data <- SDMtune::prepareSWD(
  species = "Scientific Name", 
  p = occ, 
  a = bg,
  env = bio19, 
  )


# check the data object, what is inside
data@species
data@data[1:3,]
data@pa[1:6]
data@coords[1:6,]

# TRAIN A DEFAULT MAXENT MODEL
mx1 = train(method = "Maxent",fc="lqh",verbose = T,data = data)
names(bio19)

# Model results
plotResponse(mx1, var = c("bio1_29"), type = "cloglog")
plotROC(mx1)
auc(mx1)

# predict to 
map = predict(mx1, data=bio19, type="cloglog")

map # check

# plot the results
plotPred(map, 
 lt = "Habitat\nsuitability",
 colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
)


# Check SDM PACKAGE HELP
browseURL(url = "https://cran.r-project.org/web/packages/SDMtune/vignettes/basic-use.html")

#############################################################################
sessionInfo()
# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United Kingdom.utf8 
# [2] LC_CTYPE=English_United Kingdom.utf8   
# [3] LC_MONETARY=English_United Kingdom.utf8
# [4] LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.utf8    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dismo_1.3-5   SDMtune_1.1.5 raster_3.5-15 sp_1.4-7     
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.8.3        plyr_1.8.7          RColorBrewer_1.1-3 
#  [4] pillar_1.7.0        compiler_4.2.0      plotROC_2.3.0      
#  [7] rasterVis_0.51.2    tools_4.2.0         digest_0.6.29      
# [10] viridisLite_0.4.0   lifecycle_1.0.1     tibble_3.1.7       
# [13] gtable_0.3.0        lattice_0.20-45     png_0.1-7          
# [16] pkgconfig_2.0.3     rlang_1.0.2         cli_3.3.0          
# [19] DBI_1.1.2           parallel_4.2.0      rgdal_1.5-31       
# [22] hexbin_1.28.2       rJava_1.0-6         terra_1.5-21       
# [25] dplyr_1.0.9         stringr_1.4.0       generics_0.1.2     
# [28] vctrs_0.4.1         rgeos_0.5-9         grid_4.2.0         
# [31] tidyselect_1.1.2    glue_1.6.2          R6_2.5.1           
# [34] jpeg_0.1-9          fansi_1.0.3         latticeExtra_0.6-29
# [37] farver_2.1.0        ggplot2_3.3.6       purrr_0.3.4        
# [40] magrittr_2.0.3      scales_1.2.0        codetools_0.2-18   
# [43] ellipsis_0.3.2      assertthat_0.2.1    colorspace_2.0-3   
# [46] labeling_0.4.2      utf8_1.2.2          stringi_1.7.6      
# [49] munsell_0.5.0       crayon_1.5.1        zoo_1.8-10         
