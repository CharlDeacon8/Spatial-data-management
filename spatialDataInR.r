#############
#Semi-automatic landcover classification of satellite images
#############

library(terra)

#Load usable bands individually

# Blue
b2 <- rast('band2.tif')
# Green
b3 <- rast('band3.tif')
# Red
b4 <- rast('band4.tif')
# Near Infrared (NIR)
b5 <- rast('band5.tif')

#Check bands if you want

b2
b3
b4
b5

#compare geometry of bands to see if they are similar(up to 3 bands simultaneously)

compareGeom(b2,b3)

#create a stack of the RGB layers

s <- c(b4, b3, b2)
s

#create a multiband image using all layers

filenames <- paste0('band', 1:11, ".tif") #Notice the number of the bands being left out
filenames

landsat <- rast(filenames)
landsat

writeRaster(landsat, filename="Multiband_landsat.tif", overwrite=TRUE)

#create true color image

landsatRGB <- c(b4, b3, b2)
plotRGB(landsatRGB, axes = TRUE, stretch = "lin", main = "Landsat True Color Composite")

writeRaster(landsatRGB, filename="RGBlandsat.tif", overwrite=TRUE)

#create false color image

landsatFCC <- c(b5, b4, b3)
plotRGB(landsatFCC, axes=TRUE, stretch="lin", main="Landsat False Color Composite")

writeRaster(landsatFCC, filename="FCClandsat.tif", overwrite=TRUE)

########################################




### LANDSAT CART Classification

#########################################

#load 6 bands from landsat image

library(sp)
library(rgdal)
library(terra)

### These two dataframes HAVE TO BE in the same CRS

raslist <- paste0('band', 2:7, ".tif")
landsat <- rast(raslist)
names(landsat) <- c('blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2')


samp = readOGR("trainingData32626.shp")
class(samp)

#samp_geom = geometry(samp)
#class(samp_geom)


plot(samp)
text(samp, samp$class)



#Generate random points within each polygon and convert to SpatVector object
set.seed(1)
# generate point samples from the polygons
ptsamp <- spsample(samp, 1000, type='random')
# add the land cover class to the points
ptsamp$class <- over(ptsamp, samp)$class
# We convert `ptsamp` to `SpatVector`
ptsamp <- vect(ptsamp)


# Use the x-y coordinates to extract the spectral values for the locations

xy <- as.matrix(geom(ptsamp)[,c('x','y')])
df <- extract(landsat, xy)

head(df) #Check if there are NA values

sampdata <- data.frame(class = ptsamp$class, df)

#Train classifier using 'class' from the created ptsamp object

library(rpart)
#Train the model
cartmodel <- rpart(as.factor(class)~., data = sampdata, method = 'class', minsplit = 7)

#Inspect and plot model structure

print(cartmodel)

plot(cartmodel, uniform=TRUE, main="Classification Tree")
text(cartmodel, cex = 1)

#Make predictions of the classified model across the extent of the entire image

classified <- predict(landsat, cartmodel, na.rm = TRUE)
classified

plot(classified)

#Plot final classification result, assigning distinct colours to each class

class <- c("plantation","old field","grassland","thicket", "water", "lake water", "road")
mycolor <- c('darkgreen', 'coral', 'chartreuse2', 'darkolivegreen4', 'blue', 'darkturquoise', 'burlywood4')
classdf <- data.frame(classvalue = c(1,2,3,4,5,6,7),
                      classnames = class,
                      color = mycolor
                      )
lulc <- app(classified, fun = which.max)
lulcc <- as.factor(lulc)

plot(lulcc)

#save as a .tiff file

writeRaster(lulcc, filename="classifiedLandsat.tif", overwrite=TRUE)



#Model evaluation

#fold cross evaluation
set.seed(99)

k <- 5
j <- sample(rep(1:k, each = round(nrow(sampdata))/k))

table(j)

#train and test the model 5 times

x <- list()
for (k in 1:5) {
    train <- sampdata[j!= k, ]
    test <- sampdata[j == k, ]
    cart <- rpart(as.factor(class)~., data=train, method = 'class',
                  minsplit = 7)
    pclass <- predict(cart, test, na.rm = TRUE)
    # assign class to maximum probablity
    pclass <- apply(pclass, 1, which.max)
    # create a data.frame using the reference and prediction
    x[[k]] <- cbind(test$class, as.integer(pclass))
}

#combine 5 generated list elements into a single dataframe, and compute a confusion matrix

y <- do.call(rbind, x)
y <- data.frame(y)
colnames(y) <- c('observed', 'predicted')

# confusion matrix
conmat <- table(y)

# change the name of the classes
colnames(conmat) <- classdf$classnames
rownames(conmat) <- classdf$classnames
print(conmat)

#Compute overall accuracy and the kappa metric
#Overall
# number of total cases/samples
n <- sum(conmat)
n

# number of correctly classified cases per class
diag <- diag(conmat)

OA <- sum(diag) / n
OA

#Kappa
# observed (true) cases per class
rowsums <- apply(conmat, 1, sum)
p <- rowsums / n

# predicted cases per class
colsums <- apply(conmat, 2, sum)
q <- colsums / n
expAccuracy <- sum(p*q)

kappa <- (OA - expAccuracy) / (1 - expAccuracy)
kappa

#Producer and user accuracy

# Producer accuracy
PA <- diag / colsums

# User accuracy
UA <- diag / rowsums
outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
outAcc

### SENTINAL2 CART Classification

#########################################

#load 6 bands from landsat image

library(sp)
library(rgdal)
library(terra)

### These two dataframes HAVE TO BE in the same CRS

raslist <- paste0('band', 1:7, ".tif")
landsat <- rast(raslist)
names(landsat) <- c('blue', 'green', 'red', 'NIR', 'WVP', 'SWIR1', 'SWIR2')


samp = readOGR("trainingData32726.shp")
class(samp)

#samp_geom = geometry(samp)
#class(samp_geom)


plot(samp)
text(samp, samp$class)



#Generate random points within each polygon and convert to SpatVector object
set.seed(1)
# generate point samples from the polygons
ptsamp <- spsample(samp, 1000, type='random')
# add the land cover class to the points
ptsamp$class <- over(ptsamp, samp)$class
# We convert `ptsamp` to `SpatVector`
ptsamp <- vect(ptsamp)


# Use the x-y coordinates to extract the spectral values for the locations

xy <- as.matrix(geom(ptsamp)[,c('x','y')])
df <- extract(landsat, xy)

head(df) #Check if there are NA values

sampdata <- data.frame(class = ptsamp$class, df)

#Train classifier using 'class' from the created ptsamp object

library(rpart)
#Train the model
cartmodel <- rpart(as.factor(class)~., data = sampdata, method = 'class', minsplit = 6)

#Inspect and plot model structure

print(cartmodel)

plot(cartmodel, uniform=TRUE, main="Classification Tree")
text(cartmodel, cex = 1)

#Make predictions of the classified model across the extent of the entire image

classified <- predict(landsat, cartmodel, na.rm = TRUE)
classified

plot(classified)

#Plot final classification result, assigning distinct colours to each class

class <- c("plantation","old field","grassland","thicket", "water", "lake water")
mycolor <- c('darkgreen', 'coral', 'chartreuse2', 'darkolivegreen4', 'blue', 'darkturquoise')
classdf <- data.frame(classvalue = c(1,2,3,4,5,6),
                      classnames = class,
                      color = mycolor
                      )
lulc <- app(classified, fun = which.max)
lulcc <- as.factor(lulc)

plot(lulcc)

#save as a .tiff file

writeRaster(lulcc, filename="classifiedLandsat.tif", overwrite=TRUE)



#Model evaluation

#fold cross evaluation
set.seed(99)

k <- 5
j <- sample(rep(1:k, each = round(nrow(sampdata))/k))

table(j)

#train and test the model 5 times

x <- list()
for (k in 1:5) {
    train <- sampdata[j!= k, ]
    test <- sampdata[j == k, ]
    cart <- rpart(as.factor(class)~., data=train, method = 'class',
                  minsplit = 7)
    pclass <- predict(cart, test, na.rm = TRUE)
    # assign class to maximum probablity
    pclass <- apply(pclass, 1, which.max)
    # create a data.frame using the reference and prediction
    x[[k]] <- cbind(test$class, as.integer(pclass))
}

#combine 5 generated list elements into a single dataframe, and compute a confusion matrix

y <- do.call(rbind, x)
y <- data.frame(y)
colnames(y) <- c('observed', 'predicted')

# confusion matrix
conmat <- table(y)

# change the name of the classes
colnames(conmat) <- classdf$classnames
rownames(conmat) <- classdf$classnames
print(conmat)

#Compute overall accuracy and the kappa metric
#Overall
# number of total cases/samples
n <- sum(conmat)
n

# number of correctly classified cases per class
diag <- diag(conmat)

OA <- sum(diag) / n
OA


#Kappa
# observed (true) cases per class
rowsums <- apply(conmat, 1, sum)
p <- rowsums / n

# predicted cases per class
colsums <- apply(conmat, 2, sum)
q <- colsums / n
expAccuracy <- sum(p*q)

kappa <- (OA - expAccuracy) / (1 - expAccuracy)
kappa

#Producer and user accuracy

# Producer accuracy
PA <- diag / colsums

# User accuracy
UA <- diag / rowsums
outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
outAcc


#############
#Alternative random forest classification
#############

ibrary(rgdal)
library(raster)
library(caret)


img <- brick("BioPatternGBR.tif")			###.TIF as a raster brick object
names(img) <- paste0("B", c(1,2,3))			###Change band names to B1, B2 and B3
plotRGB(img, r = 3, g = 1, b = 2)			###Check if the image loaded fine


trainData <- shapefile("Training.shp")			###Load training data into R

responseCol <- "Class"					###Tell R what column it should use




### Extract training pixel values

dfAll = data.frame(matrix(vector(), nrow = 0, ncol = length(names(img)) + 1))
for (i in 1:length(unique(trainData[[responseCol]]))){
category <- unique(trainData[[responseCol]])[i]
categorymap <- trainData[trainData[[responseCol]] == category,]
dataSet <- extract(img, categorymap)
if(is(trainData, "SpatialPointsDataFrame")){
dataSet <- cbind(dataSet, class = as.numeric(rep(category, nrow(dataSet))))
dfAll <- rbind(dfAll, dataSet[complete.cases(dataSet),])
}
if(is(trainData, "SpatialPolygonsDataFrame")){
dataSet <- dataSet[!unlist(lapply(dataSet, is.null))]
dataSet <- lapply(dataSet, function(x){cbind(x, class = as.numeric(rep(category, nrow(x))))})
df <- do.call("rbind", dataSet)
dfAll <- rbind(dfAll, df)
}
}



dfAll									###Check how many rows there are



### Model fitting and image classification

modFit_rf <- train(as.factor(class) ~ B1 + B2 + B3, method = "rf", data = dfAll)

###modFit_rf <- train(as.factor(class) ~ B1 + B2 + B3, method = "rf", data = dfAll, tuneGrid = data.frame(c(B1, B2, B3)))


beginCluster()

preds_rf <- clusterR(img, raster::predict, args = list(model = modFit_rf))

endCluster()


plot(preds_rf)

RF <- writeRaster(preds_rf, filename = "regions.tif", format = "GTiff", overwrite = TRUE)

cm_rf <- confusionMatrix(data = predict(modFit_rf, newdata = trainData),
                         trainData$Class)



#############
#Habitat suitability modelling from drone imagery (RGB bands)
#############

library(sp)
library(rgdal)			### used to export raster files from R
library(raster)			### used to read raster files into R
library(biomod2)		### used for habitat suitability models

# species occurrences
HSDtest<- read.delim("HSD.txt")

head(HSDtest)

# the name of studied species
myRespName <- 'Tf'

# the presence/absences data for our species
myResp <- as.numeric(HSDtest[,myRespName])

# the XY coordinates of species data
myRespXY <- HSDtest[,c("X","Y")]

myExpl<- stack(raster("red.tif", package="biomod2"), raster( "green.tif", package="biomod2"), raster( "blue.tif", package="biomod2"))

myExpl				### shows raster stack attributes

#This creates an environment for the data used in the model

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)

print(myBiomodData)		###Shows summary of each environmental variable (in this case RGB values)

myBiomodData				###Check data to see if it was imported correctly

plot(myBiomodData)			###Creates a plot of environmental data, but mainly grey (tiny pixels)

myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = c('SRE','CTA','RF','MARS','FDA'), NbRunEval=2, DataSplit=100, VarImport=3, models.eval.meth = c('TSS','ROC'), do.full.models=FALSE,  modeling.id="test")
				###This seems to be running fine

myBiomodModelOut

myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)

myBiomodModelEval["TSS","Testing.data","RF",,]			### 1
myBiomodModelEval["TSS","Testing.data","SRE",,]			### 0.25

myBiomodModelEval["ROC","Testing.data",,,]			###ROC scores of all variables

myBiomodModelEval["TSS","Testing.data",,,]				###TSS scores of all variables

# print variable importances
   get_variables_importance(myBiomodModelOut)			###Print the importance of variables

###ENSEMBLE SPECIES DISTRIBUTION MODELS USING TRITHEMIS FURVA AS AN EXAMPLE

###1 - Tf_AllData_Full_RF

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut, chosen.models = 'Tf_AllData_Full_RF', em.by='all', eval.metric = c('TSS'), eval.metric.quality.threshold = c(0.4), prob.mean = T, prob.cv = T, prob.ci = T, prob.ci.alpha = 0.05, prob.median = T, committee.averaging = T, prob.mean.weight = T, prob.mean.weight.decay = 'proportional')

 # print summary
myBiomodEM

# get evaluation scores
get_evaluations(myBiomodEM)

myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut, new.env = myExpl, proj.name = 'current', selected.models = 'Tf_AllData_Full_RF', binary.meth = 'TSS', compress = 'xz', clamping.mask = F, output.format = '.img') 		### output as .img (recommended)

myBiomodProj

 plot(myBiomodProj, str.grep = 'RF')  myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj

###2 - Tf_AllData_Full_SRE

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut, chosen.models = 'Tf_AllData_Full_SRE', em.by='all', eval.metric = c('TSS'), eval.metric.quality.threshold = c(0.4), prob.mean = T, prob.cv = T, prob.ci = T, prob.ci.alpha = 0.05, prob.median = T, committee.averaging = T, prob.mean.weight = T, prob.mean.weight.decay = 'proportional')

 # print summary
myBiomodEM

# get evaluation scores
get_evaluations(myBiomodEM)

myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut, new.env = myExpl, proj.name = 'current', selected.models = 'Tf_AllData_Full_SRE', binary.meth = 'TSS', compress = 'xz', clamping.mask = F, output.format = '.img')

myBiomodProj

 plot(myBiomodProj, str.grep = 'SRE')  myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj  
plot(myCurrentProj, str.grep = 'SRE')

#############
#Habitat suitability modelling from databased records and spatial datasets
#############

library(sp)
library(rgdal)			### used to export raster files from R
library(raster)			### used to read raster files into R
library(biomod2)		### used for habitat suitability models

#STEP 1 - Read all the required data into R

HSDtest<- read.delim("Ch1Test.txt")    ### This reads your text file into R, change this to the name of your text file (in orange) with presence/absence data

print(HSDtest)   ###shows a print of the text file with presence/absence data

# STEP 2 - Summon all the columns you want to use in the analysis

myRespName <- 'Al'   ###This is where you call the species you want to work with. The species name (i.e. the column name) in pink. When doing the next species, change the pink to the new species name, here, and in STEP 9.

myResp <- as.numeric(HSDtest[,myRespName])     ###This is where you assign values (0 or 1) to the species you want to work with

# the XY coordinates of species data
myRespXY <- HSDtest[,c("X","Y")]            ###This is where you assign the coordinates to the species you want to work with

# STEP 3 - Stack all the decomposed bands (.tif) into R. I am leaving the example of the ORTHOMOSAICS ONLY

myExpl<- stack(raster( "ACOCKS.tif", package="biomod2"), raster( "ERODE2.tif", package="biomod2"), raster( "MORPH.tif", package="biomod2"), raster( "RFL_SEA.tif", package="biomod2"), raster( "TempRegion.tif", package="biomod2"))

myExpl				### shows raster stack attributes

# STEP 4 - create an environment for the analysis to run

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)

plot(myBiomodData)			###Creates a plot of environmental data, but mainly grey (tiny pixels)

# STEP 5 - Run the model with SRE, CTA, RF, MARS, FDA ((although we are only really interested in RF))

myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = c('SRE','CTA','RF','MARS','FDA'), NbRunEval=2, DataSplit=100, VarImport=3, models.eval.meth = c('TSS','ROC'), do.full.models=FALSE,  modeling.id="test")

myBiomodModelOut         ### Please make sure that RF ran from this command

# STEP 6 - OPTIONAL BUT ADVISED: Shows the evaluation of RF

myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)

myBiomodModelEval["TSS","Testing.data","RF",,]

myBiomodModelEval["ROC","Testing.data", "RF",,]

# print variable importances
   get_variables_importance(myBiomodModelOut)			###Print the importance of variables

### STEP 7 - Ensemble species distribution models: this is where you create the actual data to take back to QGIS.

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut, chosen.models = 'Al_AllData_Full_RF', em.by='all', eval.metric = c('TSS'), eval.metric.quality.threshold = c(0.4), prob.mean = T, prob.cv = T, prob.ci = T, prob.ci.alpha = 0.05, prob.median = T, committee.averaging = T, prob.mean.weight = T, prob.mean.weight.decay = 'proportional')           ### This creates the Random Forests model for Trithemis furva (from the calculations in number 6)

myBiomodEM         ### Shows a summary of the model output

myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut, new.env = myExpl, proj.name = 'current', selected.models = 'Al_AllData_Full_RF', binary.meth = 'TSS', compress = 'xz', clamping.mask = F, output.format = '.img') 		### This creates a .img file of the model output, in a folder where all of your data is located. The file name should be (proj.current.SPECIESCODE), but please rename to something to identify it by later. From there, you can put it back into QGIS, and change the symbology to whatever you want it to be.

#############
#Generalised dissimilarity modelling (spatial plots as output)
#############

library(gdm)
library(raster)
library(rgdal)
library(vegan)
library(eulerr)


#STEP 1 - DATA SETUP

#Read in data
#Species data
sppTab<- read.delim("DFL.txt")

#Environmental data (rasters)(All these datasets should have the same extents, and the pixels should overlap precisely)
envRast<- stack(raster( "HUJan.tif", package="gdm"), raster( "HUFeb.tif", package="gdm"), raster( "HUMar.tif", package="gdm"), raster( "HUApr.tif", package="gdm"), raster( "HUMay.tif", package="gdm"), raster( "HUJun.tif", package="gdm"), raster( "HUJul.tif", package="gdm"), raster( "HUAug.tif", package="gdm"), raster( "HUSept.tif", package="gdm"), raster( "HUOct.tif", package="gdm"), raster( "HUNov.tif", package="gdm"), raster( "HUDec.tif", package="gdm"), raster( "HumJan.tif", package="gdm"), raster( "HumFeb.tif", package="gdm"), raster( "HumMar.tif", package="gdm"), raster( "HumApr.tif", package="gdm"), raster( "HumMay.tif", package="gdm"), raster( "HumJun.tif", package="gdm"), raster( "HumJul.tif", package="gdm"), raster( "HumAug.tif", package="gdm"), raster( "HumSep.tif", package="gdm"), raster( "HumOct.tif", package="gdm"), raster( "HumNov.tif", package="gdm"), raster( "HumDec.tif", package="gdm"), raster( "RadJan.tif", package="gdm"), raster( "RadFeb.tif", package="gdm"), raster( "RadMar.tif", package="gdm"), raster( "RadApr.tif", package="gdm"), raster( "RadMay.tif", package="gdm"), raster( "RadJun.tif", package="gdm"), raster( "RadJul.tif", package="gdm"), raster( "RadAug.tif", package="gdm"), raster( "RadSep.tif", package="gdm"), raster( "RadOct.tif", package="gdm"), raster( "RadNov.tif", package="gdm"), raster( "RadDec.tif", package="gdm"), raster( "RainJan.tif", package="gdm"), raster( "RainFeb.tif", package="gdm"), raster( "RainMar.tif", package="gdm"), raster( "RainApr.tif", package="gdm"), raster( "RainMay.tif", package="gdm"), raster( "RainJun.tif", package="gdm"), raster( "RainJul.tif", package="gdm"), raster( "RainAug.tif", package="gdm"), raster( "RainSep.tif", package="gdm"), raster( "RainOct.tif", package="gdm"), raster( "RainNov.tif", package="gdm"), raster( "RainDec.tif", package="gdm"), raster( "RCJan.tif", package="gdm"), raster( "RCFeb.tif", package="gdm"), raster( "RCMar.tif", package="gdm"), raster( "RCApr.tif", package="gdm"), raster( "RCMay.tif", package="gdm"), raster( "RCJun.tif", package="gdm"), raster( "RCJul.tif", package="gdm"), raster( "RCAug.tif", package="gdm"), raster( "RCSep.tif", package="gdm"), raster( "RCOct.tif", package="gdm"), raster( "RCNov.tif", package="gdm"), raster( "RCDec.tif", package="gdm"), raster( "TAJan.tif", package="gdm"), raster( "TAFeb.tif", package="gdm"), raster( "TAMar.tif", package="gdm"), raster( "TAApr.tif", package="gdm"), raster( "TAMay.tif", package="gdm"), raster( "TAJun.tif", package="gdm"), raster( "TAJul.tif", package="gdm"), raster( "TAAug.tif", package="gdm"), raster( "TASep.tif", package="gdm"), raster( "TAOct.tif", package="gdm"), raster( "TANov.tif", package="gdm"), raster( "TADec.tif", package="gdm"), raster( "TDRJan.tif", package="gdm"), raster( "TDRFeb.tif", package="gdm"), raster( "TDRMar.tif", package="gdm"), raster( "TDRApr.tif", package="gdm"), raster( "TDRMay.tif", package="gdm"), raster( "TDRJun.tif", package="gdm"), raster( "TDRJul.tif", package="gdm"), raster( "TDRAug.tif", package="gdm"), raster( "TDRSep.tif", package="gdm"), raster( "TDROct.tif", package="gdm"), raster( "TDRNov.tif", package="gdm"), raster( "TDRDec.tif", package="gdm"), raster( "TMAXJan.tif", package="gdm"), raster( "TMAXFeb.tif", package="gdm"), raster( "TMAXMar.tif", package="gdm"), raster( "TMAXApr.tif", package="gdm"), raster( "TMAXMay.tif", package="gdm"), raster( "TMAXJun.tif", package="gdm"), raster( "TMAXJul.tif", package="gdm"), raster( "TMAXAug.tif", package="gdm"), raster( "TMAXSep.tif", package="gdm"), raster( "TMAXOct.tif", package="gdm"), raster( "TMAXNov.tif", package="gdm"), raster( "TMAXDec.tif", package="gdm"), raster( "TMINJan.tif", package="gdm"), raster( "TMINFeb.tif", package="gdm"), raster( "TMINMar.tif", package="gdm"), raster( "TMINApr.tif", package="gdm"), raster( "TMINMay.tif", package="gdm"), raster( "TMINJun.tif", package="gdm"), raster( "TMINJul.tif", package="gdm"), raster( "TMINAug.tif", package="gdm"), raster( "TMINSep.tif", package="gdm"), raster( "TMINOct.tif", package="gdm"), raster( "TMINNov.tif", package="gdm"), raster( "TMINDec.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "ArtWet.tif", package="gdm"), raster( "LC.tif", package="gdm"), raster( "Frost.tif", package="gdm"), raster( "NatWet.tif", package="gdm"), raster( "PermNatRiv.tif", package="gdm"), raster( "PermTransRiv.tif", package="gdm"), raster( "SDR.tif", package="gdm"), raster( "SeasNatRiv.tif", package="gdm"), raster( "SeasTransRiv.tif", package="gdm"), raster( "Veg.tif", package="gdm"))

#Some environmental variables removed (for more processing efficiency)
envRast<- stack(raster( "RadJan.tif", package="gdm"), raster( "RadFeb.tif", package="gdm"), raster( "RadMar.tif", package="gdm"), raster( "RadApr.tif", package="gdm"), raster( "RadMay.tif", package="gdm"), raster( "RadJun.tif", package="gdm"), raster( "RadJul.tif", package="gdm"), raster( "RadAug.tif", package="gdm"), raster( "RadSep.tif", package="gdm"), raster( "RadOct.tif", package="gdm"), raster( "RadNov.tif", package="gdm"), raster( "RadDec.tif", package="gdm"), raster( "RainJan.tif", package="gdm"), raster( "RainFeb.tif", package="gdm"), raster( "RainMar.tif", package="gdm"), raster( "RainApr.tif", package="gdm"), raster( "RainMay.tif", package="gdm"), raster( "RainJun.tif", package="gdm"), raster( "RainJul.tif", package="gdm"), raster( "RainAug.tif", package="gdm"), raster( "RainSep.tif", package="gdm"), raster( "RainOct.tif", package="gdm"), raster( "RainNov.tif", package="gdm"), raster( "RainDec.tif", package="gdm"), raster( "TMAXJan.tif", package="gdm"), raster( "TMAXFeb.tif", package="gdm"), raster( "TMAXMar.tif", package="gdm"), raster( "TMAXApr.tif", package="gdm"), raster( "TMAXMay.tif", package="gdm"), raster( "TMAXJun.tif", package="gdm"), raster( "TMAXJul.tif", package="gdm"), raster( "TMAXAug.tif", package="gdm"), raster( "TMAXSep.tif", package="gdm"), raster( "TMAXOct.tif", package="gdm"), raster( "TMAXNov.tif", package="gdm"), raster( "TMAXDec.tif", package="gdm"), raster( "TMINJan.tif", package="gdm"), raster( "TMINFeb.tif", package="gdm"), raster( "TMINMar.tif", package="gdm"), raster( "TMINApr.tif", package="gdm"), raster( "TMINMay.tif", package="gdm"), raster( "TMINJun.tif", package="gdm"), raster( "TMINJul.tif", package="gdm"), raster( "TMINAug.tif", package="gdm"), raster( "TMINSep.tif", package="gdm"), raster( "TMINOct.tif", package="gdm"), raster( "TMINNov.tif", package="gdm"), raster( "TMINDec.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "LC.tif", package="gdm"), raster( "Frost.tif", package="gdm"), raster( "NatWet.tif", package="gdm"), raster( "PermNatRiv.tif", package="gdm"), raster( "SDR.tif", package="gdm"), raster( "SeasNatRiv.tif", package="gdm"), raster( "Veg.tif", package="gdm"))

#Some more environmental variables removed (for more processing efficiency)
envRast<- stack(raster( "RadJan.tif", package="gdm"), raster( "RadFeb.tif", package="gdm"), raster( "RadMar.tif", package="gdm"), raster( "RadApr.tif", package="gdm"), raster( "RadMay.tif", package="gdm"), raster( "RadJun.tif", package="gdm"), raster( "RadJul.tif", package="gdm"), raster( "RadAug.tif", package="gdm"), raster( "RadSep.tif", package="gdm"), raster( "RadOct.tif", package="gdm"), raster( "RadNov.tif", package="gdm"), raster( "RadDec.tif", package="gdm"), raster( "RainJan.tif", package="gdm"), raster( "RainFeb.tif", package="gdm"), raster( "RainMar.tif", package="gdm"), raster( "RainApr.tif", package="gdm"), raster( "RainMay.tif", package="gdm"), raster( "RainJun.tif", package="gdm"), raster( "RainJul.tif", package="gdm"), raster( "RainAug.tif", package="gdm"), raster( "RainSep.tif", package="gdm"), raster( "RainOct.tif", package="gdm"), raster( "RainNov.tif", package="gdm"), raster( "RainDec.tif", package="gdm"), raster( "TAJan.tif", package="gdm"), raster( "TAFeb.tif", package="gdm"), raster( "TAMar.tif", package="gdm"), raster( "TAApr.tif", package="gdm"), raster( "TAMay.tif", package="gdm"), raster( "TAJun.tif", package="gdm"), raster( "TAJul.tif", package="gdm"), raster( "TAAug.tif", package="gdm"), raster( "TASep.tif", package="gdm"), raster( "TAOct.tif", package="gdm"), raster( "TANov.tif", package="gdm"), raster( "TADec.tif", package="gdm"), raster( "TDRJan.tif", package="gdm"), raster( "TDRFeb.tif", package="gdm"), raster( "TDRMar.tif", package="gdm"), raster( "TDRApr.tif", package="gdm"), raster( "TDRMay.tif", package="gdm"), raster( "TDRJun.tif", package="gdm"), raster( "TDRJul.tif", package="gdm"), raster( "TDRAug.tif", package="gdm"), raster( "TDRSep.tif", package="gdm"), raster( "TDROct.tif", package="gdm"), raster( "TDRNov.tif", package="gdm"), raster( "TDRDec.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "LC.tif", package="gdm"), raster( "Frost.tif", package="gdm"), raster( "NatWet.tif", package="gdm"), raster( "PermNatRiv.tif", package="gdm"), raster( "SDR.tif", package="gdm"), raster( "SeasNatRiv.tif", package="gdm"), raster( "Veg.tif", package="gdm"))

#Some environmental variables removed (after covariation)
envRast<- stack(raster( "RadJan.tif", package="gdm"), raster( "RainJan.tif", package="gdm"), raster( "RainMay.tif", package="gdm"), raster( "TDRJun.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "LC.tif", package="gdm"), raster( "NatWet.tif", package="gdm"), raster( "PermNatRiv.tif", package="gdm"), raster( "SDR.tif", package="gdm"), raster( "SeasNatRiv.tif", package="gdm"), raster( "Veg.tif", package="gdm"))

#Significant environmental variables after backwards selection (species level)
envRast<- stack(raster( "RadJan.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "SDR.tif", package="gdm"))

#write.csv(envRast, file = "envRast.csv")

#STEP 2 - DEAL WITH THE BIASES REGARDING PRESENCE-ONLY DATA

#(CHOOSE A WEIGHTING METHOD BELOW)(Ferrier et al. 2007)

#Remove sites with less than x species records
gdmTab.rw <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envRast, sppFilter = 2)  ###JACCARD

gdmTab.rw <- na.omit(gdmTab.rw)		#Remove NA values

sum(is.na(gdmTab.rw))		#Make sure there are no NA values

#NOT USE
#Weight by site richness
gdmTab.rw <- formatsitepair(sppTab, bioFormat = 2, XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = "Record", predData = envRast, weightType = "richness")      ###Bray-curtis/Sørenson similarity (good for detecting underlying ecological gradients)

gdmTab.rw <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = "Record", predData = envRast, weightType = "richness")   ###Jaccard similarity (good for detecting underlying ecological gradients)

gdmTab.rw <- formatsitepair(sppTab, bioFormat = 2, dist = "raup", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = "Record", predData = envRast, weightType = "richness")  ###Raup-Crick similarity (good for dealing with unknown sample sizes)

gdmTab.rw$weights[1:5]

sum(is.na(gdmTab.rw))		#Make sure there are no NA values

gdmTab.rw <- na.omit(gdmTab.rw)		###Remove NA values

#write.csv(gdmTab.rw, file = "GDMTABLE.csv")

#Covariation (Check input data again)

COR <- cor(gdmTab.rw)

write.csv(COR, file = "covariation.csv")

#Model selection - significance determined by comparing the deviance explained by GDM fit to full model, to the deviance explained values of GDM fit to the permutated tables. Variable importance is quantified as percentage change in deviance between full model and permutated variables. Backwards elimination is used: At each step, the least important variable is dropped and the process continues until all non-significant predictors are removed.

#Set the time to Berlin time to remove time zone errors

Sys.setenv(TZ="Europe/Berlin")
Sys.getenv("TZ")

#Run backwards selection on ALL VARIABLES AFTER COVARIATION CHECKS to see which ones are important

modTest <- gdm.varImp(gdmTab.rw, geo = T, nPerm = 500, parallel = T, cores = 2, sampleSites = 0.5)
modTest
barplot(sort(modTest[[2]][,1], decreasing = T))		###From here, the significant variables are: Geographic distance, RadJan, Alt, SDR

#Save selection as an .RData file to restore later

save(modTest, file = "50%500p")
load("50%500p")

#Save entire workspace

save.image(file = "50%500p_workspace.RData")
load("50%500p_workspace.RData")

#STEP 3 - GDM ANALYSIS

gdm.1 <- gdm(gdmTab.rw, geo=T)

summary(gdm.1)		#Overview of the model - coefficients = 0 have no relationship with the biological pattern

str(gdm.1)				#Shorter summary

#Plots

length(gdm.1$predictors)			#To see how many panels there are


plot(gdm.1, plot.layout = c(4,3))

#Customization of I-spline plots - extract plotted values for each spline plot

gdm.1.splineDat <- isplineExtract(gdm.1)

str(gdm.1.splineDat)

summary(gdm.1.splineDat)

#Customised spline plots

plot(gdm.1.splineDat$x[, "Geographic"], gdm.1.splineDat$y[, "Geographic"], lwd = 3, type = "l", xlab = "Geographic distance", ylab = "Partial ecological distance")


#VISUALISTATION IN MULTIDIMENSIONAL SPACE

gdm.rast <- gdm(gdmTab.rw, geo=T)

gdm.rast <- na.omit(gdm.rast)

rastTrans <- gdm.transform(gdm.rast, envRast)				#Transform ispline plots to raster format

plot(rastTrans)


#Multidimensional shit

rastDat <- na.omit(getValues(rastTrans))

pcaSamp <- prcomp(rastDat)

pcaRast <- predict(rastTrans, pcaSamp, index = 1:3)

pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255

pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255

pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255


#Plot images (See which order looks the best for you...)
windows()
par(mfrow = c(2,3))


#RGB
plotRGB(pcaRast, r=1, g=2, b=3)

#Write to .tif file

library(rgdal)
RGB <- writeRaster(pcaRast, filename = "BioPatternRGB.tif", format = "GTiff", overwrite = TRUE)

#RBG
plotRGB(pcaRast, r=1, g=3, b=2)

###Write to .tif file

library(rgdal)
RBG <- writeRaster(pcaRast, filename = "BioPatternRBG.tif", format = "GTiff", overwrite = TRUE)

#GRB
plotRGB(pcaRast, r=2, g=1, b=3)

###Write to .tif file

library(rgdal)
GRB <- writeRaster(pcaRast, filename = "BioPatternGRB.tif", format = "GTiff", overwrite = TRUE)

#GBR
plotRGB(pcaRast, r=3, g=1, b=2)

###Write to .tif file

library(rgdal)
GBR <- writeRaster(pcaRast, filename = "BioPatternGBR.tif", format = "GTiff", overwrite = TRUE)

#BRG
plotRGB(pcaRast, r=2, g=3, b=1)

###Write to .tif file

library(rgdal)
BRG <- writeRaster(pcaRast, filename = "BioPatternBRG.tif", format = "GTiff", overwrite = TRUE)


#BGR
plotRGB(pcaRast, r=3, g=2, b=1)

###Write to .tif file

library(rgdal)
BGR <- writeRaster(pcaRast, filename = "BioPatternBGR.tif", format = "GTiff", overwrite = TRUE)

###Multidimensional plots (more spectral bands)

rastDat <- na.omit(getValues(rastTrans))

pcaSamp <- prcomp(rastDat)

pcaRast <- predict(rastTrans, pcaSamp, index = 1:5)


pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255

pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255

pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255

pcaRast[[4]] <- (pcaRast[[4]]-pcaRast[[4]]@data@min) /
  (pcaRast[[4]]@data@max-pcaRast[[4]]@data@min)*255

pcaRast[[5]] <- (pcaRast[[5]]-pcaRast[[5]]@data@min) /
  (pcaRast[[5]]@data@max-pcaRast[[5]]@data@min)*255



###Plot image

plotRGB(pcaRast, r=1, g=2, b=3)

###Write to .tif file

library(rgdal)
RF2 <- writeRaster(pcaRast, filename = "BioPattern(5 bands).tif", format = "GTiff", overwrite = TRUE)

###STEP 4 - VARIANCE PARTITIONING (SPECIES LEVEL)

###HERE, ALL COMBINATIONS OF ENVIRONMENTAL VARIABLES NEED TO BE FITTED TO CALCULATE THE EXPLAINED DEVIANCE

sppTab<- read.delim("DFL.txt")

########## 1 & 2

####Significant environmental variables after backwards selection
envALL<- stack(raster( "RadJan.tif", package="gdm"), raster( "Alt.tif", package="gdm"), raster( "SDR.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.ALL <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envALL, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.ALL <- na.omit(gdmTab.ALL)		###Remove NA values

sum(is.na(gdmTab.ALL))		###Make sure there are no NA values



gdm.1 <- gdm(gdmTab.ALL, geo=T)
explX1 <- gdm.1$explained
if(is.null( explX1 ) == TRUE){ explX1 <- 0}



gdm.2 <- gdm(gdmTab.ALL)
explX2 <- gdm.2$explained
if(is.null( explX2) == TRUE){ explX2 <- 0}




########## 3 & 4

####Significant environmental variables after backwards selection
envRA<- stack(raster( "RadJan.tif", package="gdm"), raster( "Alt.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.RA <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envRA, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.RA <- na.omit(gdmTab.RA)		###Remove NA values

sum(is.na(gdmTab.RA))		###Make sure there are no NA values


gdm.3 <- gdm(gdmTab.RA, geo=T)
explX3 <- gdm.3$explained
if(is.null( explX3) == TRUE){ explX3 <- 0}


gdm.4 <- gdm(gdmTab.RA)
explX4 <- gdm.4$explained
if(is.null( explX4 ) == TRUE){ explX4 <- 0}



########## 5 & 6


####Significant environmental variables after backwards selection
envRD<- stack(raster( "RadJan.tif", package="gdm"), raster( "SDR.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.RD <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envRD, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.RD <- na.omit(gdmTab.RD)		###Remove NA values

sum(is.na(gdmTab.RD))		###Make sure there are no NA values


gdm.5 <- gdm(gdmTab.RD, geo=T)
explX5 <- gdm.5$explained
if(is.null( explX5) == TRUE){ explX5 <- 0}


gdm.6 <- gdm(gdmTab.RD)
explX6 <- gdm.6$explained
if(is.null( explX6) == TRUE){ explX6 <- 0}



########## 7 & 8


####Significant environmental variables after backwards selection
envR<- stack(raster( "RadJan.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.R <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envR, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.R <- na.omit(gdmTab.R)		###Remove NA values

sum(is.na(gdmTab.R))		###Make sure there are no NA values


gdm.7 <- gdm(gdmTab.R, geo=T)
explX7 <- gdm.7$explained
if(is.null( explX7) == TRUE){ explX7 <- 0}


gdm.8 <- gdm(gdmTab.R)
explX8 <- gdm.8$explained
if(is.null( explX8) == TRUE){ explX8 <- 0}



##########  9 & 10

####Significant environmental variables after backwards selection
envAD<- stack(raster( "Alt.tif", package="gdm"), raster( "SDR.tif", package="gdm"))



###Remove sites with less than 2 species records
gdmTab.AD <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envAD, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.AD <- na.omit(gdmTab.AD)		###Remove NA values

sum(is.na(gdmTab.AD))		###Make sure there are no NA values


gdm.9 <- gdm(gdmTab.AD, geo=T)
explX9 <- gdm.9$explained
if(is.null( explX9) == TRUE){ explX9 <- 0}


gdm.10 <- gdm(gdmTab.AD)
explX10 <- gdm.10$explained
if(is.null( explX10) == TRUE){ explX10 <- 0}



############ 11 & 12

####Significant environmental variables after backwards selection
envA<- stack(raster( "Alt.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.A <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envA, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.A <- na.omit(gdmTab.A)		###Remove NA values

sum(is.na(gdmTab.A))		###Make sure there are no NA values


gdm.11 <- gdm(gdmTab.A, geo=T)
explX11 <- gdm.11$explained
if(is.null( explX11) == TRUE){ explX11 <- 0}


gdm.12 <- gdm(gdmTab.A)
explX12 <- gdm.12$explained
if(is.null( explX12) == TRUE){ explX12 <- 0}




############ 13 & 14

####Significant environmental variables after backwards selection
envD<- stack(raster( "SDR.tif", package="gdm"))


###Remove sites with less than 2 species records
gdmTab.D <- formatsitepair(sppTab, bioFormat = 2, dist = "jaccard", XColumn = "Long", YColumn = "Lat", sppColumn = "sppCol", siteColumn = NULL, predData = envD, sppFilter = 10, weightType="richness")  ###JACCARD

gdmTab.D <- na.omit(gdmTab.D)		###Remove NA values

sum(is.na(gdmTab.D))		###Make sure there are no NA values


gdm.13 <- gdm(gdmTab.D, geo=T)
explX13 <- gdm.13$explained
if(is.null( explX13) == TRUE){ explX13 <- 0}



gdm.14 <- gdm(gdmTab.D)
explX14 <- gdm.14$explained
if(is.null( explX14) == TRUE){ explX14 <- 0}


############# 15


##Load in a dataset with dummy variables to get geographical distance alone

gdmTabGeo<- read.delim("DUMMY2.txt")



gdm.15 <- gdm(gdmTabGeo, geo=T)
explX15 <- gdm.15$explained
if(is.null( explX15) == TRUE){ explX15 <- 0}








###Unique to each variable

R <- explX1 - explX9
vR <- R
if( R < 0){ R <- 0}


A <- explX1 - explX5
vA <- A
if( A < 0){ A <- 0}


S <- explX1 - explX3
vS <- S
if( S < 0){ S <- 0}


G <- explX1 - explX2
vG <- G
if( G < 0){ G <- 0}


###2 groups combined

RA <- explX1 - (R + A) - explX13
vRA <- RA
if( RA < 0){ RA <- 0}


RS <- explX1 - (R + S) - explX11
vRS <- RS
if( RS < 0){ RS <- 0}


RG <- explX1 - (R + G) - explX10
vRG <- RG
if( RG < 0){ RG <- 0}


AS <- explX1 - (A + S) - explX7
vAS <- AS
if( AS < 0){ AS <- 0}


AG <- explX1 - (A + G) - explX6
vAS <- AG
if( AG < 0){ AG <- 0}


SG <- explX1 - (S + G) - explX4
vAS <- SG
if( SG < 0){ SG <- 0}


#### groups combined

RAS <- explX1 - (R + A + S) - (RA + RS + AS) - explX15		###
vRAS <- RAS
if( RAS < 0){ RAS <- 0}


RSG <- explX1 - (R + S + G) - (RS + RG + SG) - explX12
vRSG <- RSG
if( RSG < 0){ RSG <- 0}


RAG <- explX1 - (R + A + G) - (RA + RG + AG) - explX14
vRAG <- RAG
if( RAG < 0){ RAG <- 0}


ASG <- explX1 - (A + S + G) - (AS + AG + SG) - explX8
vASG <- ASG
if( ASG < 0){ ASG <- 0}


RASG <- explX1 - (A + S + G) - (RA + RS +RG + AS + AG + SG) - (RAS + RSG + RAG + ASG)
vRASG <- RASG
if( RASG < 0){ RASG <- 0}





allParts <- (data.frame(R,A,S,G,RA,RS,RG,AS,AG,SG,RAS,RSG,RAG,ASG,RASG))

#round to 2 decimal spaces
allParts <- (round(allParts, digits=2))
allParts

#transpose data frame
allParts <- as.data.frame(t(allParts))

allocate<- (c("R","A","S","G","R&A", "R&S", "R&G", "A&S", "A&G", "S&G", "R&A&S", "R&S&G", "R&A&G","A&S&G","R&A&S&G"))

allocate

text <- paste(allocate,  "=", values$V1)
values <- data.frame(allocate, allParts)
values

onlyValues <- values[apply(values[2],1,function(z) !any(z==0)),]
onlyValues

VALUES <- onlyValues[2]
row.names(VALUES) <- onlyValues$allocate
data <- (t(VALUES))
data

borg <- c(data[,(1:(ncol(data)))])
borg


plot((euler(borg, shape = "ellipse")))

plot((euler(borg, shape = "ellipse")), quantities = TRUE)

###plot((euler(borg, shape = "ellipse")), labels = c("JanRad", "Alt", "SDR", "GeogDist"), quantities = TRUE)
