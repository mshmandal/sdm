#############################################################################
# Author: Mohammad Shamim Hasan Mandal
# Date  : July 13 2022
# name  : script.R
# aim   : Demo project: setup R project and make Maxent model
#############################################################################
# SETUP R PROJECT

#start fresh
gc()
rm(list=ls())

# set work directory
setwd("C:/Git/sdm") 
library(raster)

dir.create("./data")
dir.create("./result")


#############################################################################
# GET BD SHAPEFILE and WORLCLIM data
ext = c(91,93,19,23)  # c(xmin,xmax,ymin,ymax)
bd  = getData('GADM',country='Bangladesh',level=2,path = "./data",download = F)
bd  = crop(bd,extent(ext))

# GET WORLD CLIM DATA # browseURL("https://www.worldclim.org/data/bioclim.html") 
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
 
#############################################################################
# MAXENT
library(SDMtune)
library(dismo)

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
