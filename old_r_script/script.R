#############################################################################
# Author: Mohammad Shamim Hasan Mandal
# Date  : July 13 2022
# Name  : script.R
# Obj.  : Demo project: setup R project and make Maxent model
# Update: October 16 2022
#############################################################################
# SETUP R PROJECT

# Start fresh
# When we run a R project over and over again. Sometimes it is better to start
# from beginning. Clearing R memory and session variables

gc()            # clear memory
rm(list=ls())   # clear R session

# Set work directory
# This is an important step because we want to make sure our data and outputs
# are saved clearly inside our project folder. So, we first make sure that we
# are working inside our project directory. To check the current project
# directory use getwd() function. And to change the current working directory
# to a new one, use the setwd() function
getwd()
# setwd("D:/github/sdm") 

# Now we create two separate folder inside our project directory i.e. inside 
# the "D:/github/sdm" folder to save data and results
dir.create("./data")
dir.create("./result")


#-------------------------------------------------------------------------------
# STEP 1: Install required libraries
#-------------------------------------------------------------------------------
# For our project first we need to install few R packages. Make sure to run them
# for first time. Note: the dismo package is for running the maxent. 
# The maxent function requires rJava package. So, in your desktop/laptop Java
# should be installed, if not you will get error.

install.packages("raster",dependencies = T)
install.packages("dismo",dependencies = T)
install.packages("SDMtune",dependencies = T)

# Now load the packages
library(raster)
library(SDMtune)
library(dismo)

#-------------------------------------------------------------------------------
# Step 2 : Prepare data for model
#-------------------------------------------------------------------------------

# There are different ways of implementing a SDM.
# The following workflow is a practical way that I use. I do not state that
# this is the best way of most efficient way. However, I find that this works
# good for my projects.

# The most important and perhaps the hardest of all in making a SDM to decide
# study area, input predictors and then prepare for modelling. Here we, first
# decide an study area extent bounded by bounding box 
# 91,93,19,23 (min longitude, max longitude, min latitude, max latitude)
# Then we download a level2 administrative division shapefile from GADM website

## 2.1 GET BD SHAPEFILE and WORLCLIM data
ext = c(91,93,19,23)  # c(xmin,xmax,ymin,ymax)
# -30.506739, -30.530748, 143.158580, 143.128823
# aus_ext = c(-30.506739, -30.530748, 143.158580, 143.128823)
plot(ext(ext))

help(getData,raster)
bd  = getData('GADM',country='Bangladesh',level=2,path = "./data",download = F)
# aus=getData('GADM',country='Australia',level=2,path = "./data",download = T)
# plot(aus,axes=T)
bd_main = bd # saving in a new variable to use later

# crop the country shape file to study area extent
bd  = crop(bd,extent(ext))

# Make plot to visualize the data layers
par(mfrow=c(1,2))
plot(bd_main,axes=T, main="Bangladesh admin Level2")
plot(bd,axes=T, main="Cropped shapefile")
dev.off() # close the plot


## 2.2 Raster/Image data: GET WORLD CLIM DATA 

# The worldclim website provides the documentation of the bioclimatic variables
# and other important details about the data. There are other variables for 
# example elevation, windspeed etc. 
# Note several future climatic projection datasets are also avaialbe on
# the website. I recommend to check the website and read about the data
# browseURL("https://www.worldclim.org/data/bioclim.html") 


# Here we, download the 19 bioclimatic variables.
# We coose the resolution "0.5" degress which is roughly 1km2 grid size
# This is the spatial resolution in case of images. Because this is a big 
# image we can't download the global image, rather we need to give a longitude
# latidue coordinates as getData() function argument. Then an data will
# be downloaded sorrunding the coordinates (this is not the best description)
# but the function documentation is good enouhg for understanding.
# So, please read help(getData, raster) package documentation
# We also, provide a directory where to save the data as the path argument "path"
# The download=F argument specifies that don't download, if we alreay previously
# download

help(getData,raster)
bio19 = getData(
  name='worldclim',
  var='bio',
  res=0.5,
  lon=90, 
  lat=22,
  path = "./data",
  download=F
  )

# the number of layers
nlayers(bio19)

bio19 = crop(bio19,extent(ext))  # crop the data using the extent

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

# For maxent SDM
# We need occurrenct points which is sometimes called observed occurrences
# Then we also need where the species are not present, also called absense
# Logistic model, # Presense = 1, Absense  = 0
# But Maxent is popular because we do not necessarily need or have Absense data
# We can randomly genereate these points and using Maxent model try to 
# approximate the species distribution, please refer to 
# Phillips, S. J. (2005). A brief tutorial on Maxent. AT&T Research, 190(4), 231-259.

# Presence points
#occ = dismo::randomPoints(n=30,bio19[[1]])
# If you have your species occurrence points in a csv file where first two
# columns are "longitude" and "latitude" (order is important), then you can
# read them 
occ = read.csv("./data/occ.csv")


# Absense points
bg = dismo::randomPoints(n=1000,bio19[[1]])


# Prepare data SWD object which actually only requires for SDMtune model
# If we want to use other package we do not need this format. But, SDMtune
# package is very useful, later we will discover.
data <- SDMtune::prepareSWD(
  species = "Pathenium hys", 
  p = occ, 
  a = bg,
  env = bio19, 
  )


# check the data object, what is inside
data@species
data@data[1:10,1:3]
data@pa[1:6]
data@coords[1:6,]


#-------------------------------------------------------------------------------
# TRAIN A DEFAULT MAXENT MODEL
#-------------------------------------------------------------------------------
# the train() function of the SDMtune package is used to train a model
# The type of model to train is defined by the method argument. In our case
# it is "Maxent". Read the function documentation
help(train, SDMtune)

# For a maxent model we can specify some more argument, for example which 
# feature classes we want to use. The argument for that is "fc". The below code
# we only specify "lqh" which stands for linear, quardatic and hinge feature 
# classes. The regularization value in the original maxent implementation is
# 0.5. But increase of SDMtune is 1, so we can specify which one to use
# If we want we can tune this value i.e. find a value based on some algorithm
# search

mx1 = train(
  method = "Maxent",  # the name of the model
  fc="lqh",           # feature classes: linear, quadratic, hinge
  reg=0.5,            # regularization parameters
  verbose = T,        # show message during training
  data = data,        # the data variable
  iter = 1000         # number of iterations the logarithm will run
  )    


# Model results
# SDMtune has several functions to show the model results
# First we plot the response curves, for the first predictor variable
# we can do that for multiple variables also
plotResponse(mx1, var = c("bio1_29"), type = "logistic")

# Then we plot the ROC curve for our model which shows the traing AUC
auc(mx1)        # just to get the model AUC
plotROC(mx1)    # to plot the AUC curve


#-------------------------------------------------------------------------------
# APPLY MODEL PREDICTION
#-------------------------------------------------------------------------------
# Finally we predict our train model to the data to get probability of 
# occurrence in each pixel/gird in our study area
map = predict(mx1, data=bio19, type="logistic")

# The output is a raster layer, which is a image data class for raster package
# we can plot the map using plotPred() function. We could use the default plot()
# function or even ggplot2 package. But for now we will just use the plotPred 
# function.

# plot the results
dev.off() # clear the plotting window at first
plotPred(map, 
 lt = "Habitat\nsuitability",
 colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
)

# Here is another way of plotting the map
# crop the prediction to our study area
bd_pred = mask(map,bd)
plot(
  bd_pred,
  breaks=seq(0,1,0.2),
  col= rainbow(length(seq(0,1,0.2)))
 )
plot(bd,add=T) # overlay the shapefile on top of it

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

