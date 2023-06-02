#################################################################################
###Retrieve data from MERRAclim and CHELSA for sub-Antractic/Antarctic islands###
#################################################################################

#extract BIO01, BIO04, BIO12, BIO15

#set workig directory
setwd("C:/Users/u21846708/Desktop/R_ACE_project")

library(raster)

#1. MERRAclim

#for reference, here is the link to the MERRAclim files: https://datadryad.org/stash/dataset/doi:10.5061/dryad.s2v81

#import raster and create RasterStack
rasStack_bio1 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db2/doi_10.5061_dryad.s2v81__v2/2_5m_mean_00s/2_5m_mean_00s_bio1.tif")
rasStack_bio4 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db2/doi_10.5061_dryad.s2v81__v2/2_5m_mean_00s/2_5m_mean_00s_bio4.tif")
rasStack_bio12 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db2/doi_10.5061_dryad.s2v81__v2/2_5m_mean_00s/2_5m_mean_00s_bio12.tif")
rasStack_bio15 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db2/doi_10.5061_dryad.s2v81__v2/2_5m_mean_00s/2_5m_mean_00s_bio15.tif")

#import geographical points
pointCoordinates <- data.frame(read.csv("C:/Users/u21846708/Desktop/R_ACE_project/geographical_coordinates.csv"))
#coordinates(pointCoordinates)= ~ longitude+ latitude

#extract vales from raster for the geographical points
rasValue_bio1 = extract(rasStack_bio1, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio4 = extract(rasStack_bio4, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio12 = extract(rasStack_bio12, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio15 = extract(rasStack_bio15, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude

#save data in csv file
combinePointValue=cbind(pointCoordinates,rasValue_bio1,rasValue_bio4,rasValue_bio12,rasValue_bio15)
combinePointValue
write.csv(combinePointValue,"extract_bio_ACE.csv", row.names = FALSE)

#2. CHELSA

#for reference, here is the link to the CHELSA files: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2F

#import raster and create RasterStack
rasStack_bio1 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db5/CHELSA_bio1_1981-2010_V.2.1.tif")
rasStack_bio4 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db5/CHELSA_bio4_1981-2010_V.2.1.tif")
rasStack_bio12 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db5/CHELSA_bio12_1981-2010_V.2.1.tif")
rasStack_bio15 <- raster("C:/Users/u21846708/Desktop/R_biogeography/extract_geochemical_data/db5/CHELSA_bio15_1981-2010_V.2.1.tif")

#import geographical points
pointCoordinates <- data.frame(read.csv("C:/Users/u21846708/Desktop/R_ACE_project/geographical_coordinates.csv"))
#coordinates(pointCoordinates)= ~ longitude+ latitude

#extract vales from raster for the geographical points
rasValue_bio1 = extract(rasStack_bio1, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio4 = extract(rasStack_bio4, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio12 = extract(rasStack_bio12, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude
rasValue_bio15 = extract(rasStack_bio15, pointCoordinates[,c(3,2)]) #column 3 reports longitude and column 2 reports latitude

#save data in csv file
combinePointValue=cbind(pointCoordinates,rasValue_bio1,rasValue_bio4,rasValue_bio12,rasValue_bio15)
combinePointValue
write.csv(combinePointValue,"extract_bio_ACE_CHELSA.csv", row.names = FALSE)
