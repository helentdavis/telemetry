## telemetry-script.R
## R 3.5.1
## Date Created: 01 October 2018
## Intended Use: Written to analyze home range of raptors with telemetry units. 
## Home range calculations are based on Calabrese, J.M. and C.H. Fleming. 2016. ctmm: an R package for analyzing animal relocation dataas a continuous-time stochastic process. Methods in Ecology and Evolution (7) 1124-1132.
## Inputs: Input data must be in .csv format. Timestamp must be sorted oldest to newest in .csv prior to import in R.
## Variables: bird_name: unique name/code for individual bird
##            long_x: longitude
##            lat_y: latitude
##            timestamp: mm/dd/yyyy hh:mm:ss
##            sensor_type: gps, doppler, etc.
## Outputs: 95% homerange shapefile (.shp) and raster (.tif). Variogram plots can be saved manually using export in R studio.
## Instructions: Run script line-by-line from R Studio. Be aware of interactive variogram.fit() function, where GUESS input must be selected. Be aware as.telemetry() will default two-point equidistant projection unless otherwise specified.

# Install Required Packages 
install.packages(c("maptools","sp", "raster","move","ctmm"))

# Load Required Packages
packagelist <- c("maptools","sp", "raster", "move","ctmm")  
lapply(packagelist, require, character.only = TRUE)

#Create User Input Prompts
changewd <- function ()
{
  n <- readline(prompt = "Enter path to working directory:")
  return(n)
}
filename <- function ()
{
  n <- readline(prompt = "Enter dataset filename:")
  return(n)
}
shapename <- function ()
{
  n <- readline(prompt = "Enter output shape name:")
  return(n)
}

## Set working directory, Read in data
setwd(print(changewd()))
#verify successfully changed
print(getwd())
movedata <- read.csv(print(filename()))

#verify successfully changed
head(movedata)

#Reformat timestamp, remove NAs
movedata$time = as.POSIXct(movedata$timestamp, format="%m/%d/%Y %H:%M:%S", tz="GMT")
movedata <- movedata[!is.na(movedata$time),]

# Convert individual eagle data to move object. Only run function for eagle you want to analyze.
# Ignore warning unless dataset contains multiple eagles

movedata_move <- move(x=movedata$long_x, y=movedata$lat_y,
                      time=movedata$time,
                      proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                      data=movedata, animal=movedata$bird_name, 
                      sensor=movedata$sensor_type)

# Convert to telemetry object, identify needed projection to minimize distortion. Function will default two-point equidistant.
#Ignore warning if using individual eagle data
autokde <- as.telemetry(movedata_move) 

plot(autokde)

# Create variogram to assess autocorrelation and movement patterns. Reference Calabrese et al. 2016 to identify potential periodicity. 
gevario <- variogram(autokde)

plot(gevario)

# The variogram.fit function will fit a model to the variogram of our data.  Sliders in the plot window 
# will help select initial parameters that fit the data best. Either accept default parameters by selecting "GUESS"
# or move sliders until satisfied, then selected "GUESS"
# Guidance on determining proper slider settings can be found in Calabrese et al. 2016

gs <- variogram.fit(gevario, interactive=T)

# Use ctmm.select to select highest ranked model out of competing model set.
# Model parameters selected in variogram.fit will be considered

fitmods <- ctmm.select(autokde, CTMM=GUESS, verbose = T)
summary(fitmods)

# Models <2 delta AIC have similar support
# Select top model if no other models compete. 
# If additional models have similar support determine if differences are biologically meaningful or can be eliminated due to parsimony

top <- fitmods[[1]]

# Estimate and plot the UD.
# This process may take a few minutes.
akde.ou <- akde(autokde, CTMM=top)
plot(autokde, UD=akde.ou)

# Check home range size
summary(akde.ou)

# Define output filename
name <- print(shapename())
# Write raster with 95% homerange estimate
# Write shapefile with 95% small, mid, and large homerange estimates
writeShapefile(akde.ou, getwd(), file=name)
writeRaster(akde.ou, format="GTiff", filename=name)





