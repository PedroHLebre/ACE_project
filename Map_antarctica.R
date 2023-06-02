pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop(x, " :Package not found")
  }
}

list.of.packages <- c("maptools","rgdal","sp", "oce","ocedata", "mapdata", "ggplot2", "RColorBrewer", "ggspatial", "ncdf4", "lubridate","scales")

# create list of installed packages
pkges = installed.packages()[,"Package"]
for (pk in list.of.packages) {
  pkgTest(pk)
}

# Load dataset with Lat Long coordinates of samples

indata=read.table (file.choose(), header = TRUE, row.names = 1, sep = ",")

latitude <- indata$Latitude
longitude <- indata$Longitude

# Longitude for this dataset goes from 0-360 but the transform library uses -180/180
# Convert from 360 to -180/180
lon180 <- ((longitude + 180) %% 360) - 180
indata$V2 <- lon180

dfcoords = cbind(lon180,latitude)      # coords in lon, lat order
sppoints = SpatialPoints(coords = indata)
spdf     = SpatialPointsDataFrame(coords = dfcoords, indata)

# Verify initial coordinates of spatial dataframe
coordsinit <- spdf@coords

# Define coordinate reference systems
crslonglat       = CRSargs(CRS("+init=epsg:4326")) # order is longitude latitude in R CRS
crsseaicepolster3412 = CRSargs(CRS("+init=epsg:3412"))

# Set initial CRS of spatial dataframe
proj4string(spdf) = CRS(crslonglat)
lat_lon_bbox       = spdf@bbox
print(lat_lon_bbox)

# Check the CRS 
crs_set = proj4string(spdf)

# Converts from existing lat-lon crs to polar stereo (3412)
spdfProjected = spTransform(spdf, CRS(crsseaicepolster3412))  
crs_projected = proj4string(spdfProjected)

coordsproj = spdfProjected@coords
bbox       = spdfProjected@bbox
print( bbox )


# define map extents
ylim <- c(-100,-50)
xlim <- c(-180,71)
data("coastlineWorldMedium") # included in ocedata
class(coastlineWorldMedium)

# set plot margins (bottom, left, top, right)
par(mar=c(2, 6, 2, 6))

## make a base
mapPlot(coastlineWorldMedium, 
        projection=crsseaicepolster3412,
        col="lightgray", 
        longitudelim=xlim, 
        latitudelim=ylim,
        main="Argo Float locations, lat/lon coordinates with oce map in EPSG:3412"
    )

## add points
mapPoints(longitude, latitude, col = "red", pch = 20, cex = 2)
