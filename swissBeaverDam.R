
####################Script for processing beaver dam location data, river characteristics and dam probability model

library(sf)

library(terra)

#load river network
rivJoinedData<-st_read("rivDataBafu.shp")

head(rivJoinedData)

#load presence and background data
snap2022<-st_read("presAbsSnapped20.shp")


#get coordinates of all data points for model evaluation later
coord22<-st_coordinates(snap2022)

snap2022$X<-coord22[,1]
snap2022$Y<-coord22[,2]

#small buffer to make sure dam intersects nearest stream
damBufferRiv<-st_buffer(snap2022,dist=1,endCapStyle = "ROUND")

#need this to track IDs later
damBufferRiv$damIDs<-1:nrow(damBufferRiv)

#set projections to be equal
damBufferRiv<-st_transform(damBufferRiv,st_crs(rivJoinedData))

#get intersection of dams and river network
interPointPred<-st_intersection(x=damBufferRiv,y=rivJoinedData)

#separate dams and background points
interPointPres<-interPointPred[interPointPred$pres==1,]
interPointBack<-interPointPred[interPointPred$pres==0,]

#drop geometry to do some data cleaning
dfPoint<-st_drop_geometry(interPointPred)


#drop ancillary fields
dfPoint<-dfPoint[,c(1,2,5:17)]  


#make sure all categorical data are numeric
dfPoint$pres<-as.numeric(dfPoint$pres)

dfPoint$A<-as.numeric(dfPoint$A)

#for data handling
library(dplyr)

#set up data frame to hold results - this section needed in case a dam intersects more than one stream

swissData<-data.frame(A=(rep(NA,length(unique(dfPoint$damIDs)))),pres=NA,X=NA,Y=NA,damIDs=NA,
                           discharge=NA,streamPower=NA,width=NA,terrainSlope=NA,gradient=NA,
                           urban100=NA,arable100=NA,pasture100=NA,wood100=NA,minSlope=NA)



#iterate over all variables and get mean for any point where intersecting streams are > 1
for(i in 1:ncol(swissData)){
  inter.i<-dfPoint[,i]
  
  interApp<-tapply(inter.i,dfPoint$damIDs,mean)
  print(i)
  swissData[,i]<-interApp
  
}

#set all data to numeric form
for(i in 1:ncol(swissData)){
  df.i<-swissData[,i]
  dfVector<-as.numeric(df.i)
  swissData[,i]<-dfVector
}


#dams
dataPres<-swissData[swissData$pres==1,]

#pseudo-absence
dataBack<-swissData[swissData$pres==0,]

#remove outliers (these were identified individually, the limits below remove them)
dataPres<-dataPres[dataPres$discharge<=0.7,]
dataPres<-dataPres[dataPres$width<20,]
dataPres<-dataPres[dataPres$gradient<0.3,]

# final data set
swissData<-rbind(dataPres,dataBack)

#anything with width == 0 is an error
swissData<-swissData[swissData$width>0,]


#object for number of presence and absence points to use as weights in regression
nPres<-nrow(swissData[swissData$pres==1,])
nBack<-nrow(swissData[swissData$pres==0,])


#add IDs for tracking rows
swissData$ID=1:nrow(swissData)

#data for training
dataTrain=sample_frac(swissData,0.8)

#data for testing 
dataTest=swissData[!swissData$ID%in%dataTrain$ID,]

#glm for the win (use quasibinomial to handle dispersion)
glmGeo<-glm(pres~streamPower+gradient+terrainSlope+minSlope+poly(width,2)+poly(pasture100,2)+
              poly(urban100,2)+poly(wood100,4),
            family = "quasibinomial",weights=(nPres/nBack)^(1-pres),data=dataTrain)

#check model results
summary(glmGeo)


#for model evaluation
library(precrec)

#make prediction to test data
pred <- predict(glmGeo,as.data.frame(dataTest),type="response")

#use precrec library's evalmod() function to obtain area under curve 
precrec_proc <- evalmod(scores = pred,labels = dataTest$pres,mode = "basic")

#extract proc thresholds, specificity and sensitivity outputs 

proc_df=data.frame(cbind(precrec_proc$score[[1]]$y,precrec_proc$sp[[1]]$x,precrec_proc$sn[[1]]$y))

colnames(proc_df)<-c("Tr","Specificity","Sensitivity")

#plot
plot(proc_df$Specificity,proc_df$Sensitivity)

#isolate threshold as that which maximises sensitivity and minimises specificity
TrMax<-proc_df[which.max(proc_df$Sensitivity-proc_df$Specificity),]

TrMax=TrMax$Tr


# set up new data frame for building response plots - below example is for response to stream width

#sequence of width values in range of damming data
width<-seq(0.1,15,
           length=1000)


# data frame with all other variables set to their median
newDataWidth<-data.frame(width=width,
                         gradient=median(swissData$gradient),
                         wood100=median(swissData$wood100),
                         pasture100=median(swissData$pasture100),
                         minSlope=median(swissData$minSlope),
                         terrainSlope=median(swissData$terrainSlope),
                         streamPower=median(swissData$streamPower),  
                         urban100=median(swissData$urban100))




#make a prediction
predWidth<-stats::predict(glmGeo,newdata = newDataWidth,type="response",se=T)


ucl<-predWidth$fit+1.96*predWidth$se.fit
lcl<-predWidth$fit-1.96*predWidth$se.fit

# final data frame for plotting

glmDataNew<-data.frame(newDataWidth,pred=(predWidth$fit),lcl=(lcl),
                       ucl=(ucl))




plot(glmDataNew$width,glmDataNew$pred,ylim=c(-0.002,0.15),xlim=c(0,15),cex=0.2,
     xlab="Width",ylab="Dam Prob")

polygon(c(glmDataNew$width,rev(glmDataNew$width)),c(glmDataNew$lcl,rev(glmDataNew$ucl)),
        col = "grey85", border = FALSE)


lines(glmDataNew$width,glmDataNew$pred,col="black",lwd=2)
lines(glmDataNew$width,glmDataNew$lcl,col="red",lwd=2,lty=2)
lines(glmDataNew$width,glmDataNew$ucl,col="red",lwd=2,lty=2)


#######################################################################
################Spatial Cross-validation##################


#build spatial layers from swissData 

presPoints<-swissData[swissData$pres==1,]
backPoints<-swissData[swissData$pres==0,]

#create sf object from coordinates
pres<-st_as_sf(presPoints,coords=c("X","Y"))

#background/pseudo-absence points
backPoints<-st_as_sf(backPoints,coords=c("X","Y"))

#set the crs
st_crs(pres)<-"epsg:2056"
st_crs(backPoints)<-"epsg:2056"


# make a grid for spatial blocking
area_fishnet_grid = st_make_grid(pres, c(100000, 100000), what = "polygons", square = T)

# plot
plot(area_fishnet_grid)

#check all points in correct place
plot(pres$geometry,add=T)


#get data into correct shape for analysis
fishnet_grid_sf = st_sf(area_fishnet_grid) %>%

#create id field to use as folds
mutate(grid_id = 1:length(lengths(area_fishnet_grid)))

#make sure only use data inside grid
fishInt<-st_contains(fishnet_grid_sf,pres$geometry)

#select points within
sel_logical = lengths(fishInt) > 0

#subset
fishnet_grid_Fin<-fishnet_grid_sf[sel_logical,]


# set folds
folds=fishnet_grid_Fin$grid_id

#combine pres/abs data the split data according to folds using the kfold() function.

dfST<-rbind(pres,backPoints)

#creat sf object and make sure projections all same
dfST<-st_as_sf(dfST,coords=c("X","Y"))

st_crs(dfST)<-st_crs(fishnet_grid_sf)

# now do the cross-validation
library(ranger) #for fast random forest
library(pROC) # for AUC plotting
rfK<-list() #list to hold AUC results

#set up formula for random forest model
rangerForm<-pres~streamPower+width+gradient+terrainSlope+minSlope+pasture100+urban100+wood100


for (i in folds) {
  
  #select fold
  fishNetTrain<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id!=i)
  
  #spatially subset the training data
  train <- st_intersection(fishNetTrain, dfST)#for presence values, select all folds which are not 'i' to train the model
  
  #spatially subset the test data
  fishNetTest<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id==i)
  test <- st_intersection(fishNetTest, dfST)  
  nPresTrain<-nrow(train[train$pres==1,])
  nBackTrain<-nrow(train[train$pres==0,])
  
  #rangerForm
  rfGeo<-ranger(formula = rangerForm,case.weights=(nPresTrain/nBackTrain)^(1-train$pres),data = as.data.frame(train), num.trees = 500, mtry = 2,min.node.size = 10)
  
  pred <- predict(rfGeo, as.data.frame(test),type="response")
  
  precrec_proc <- evalmod(scores = pred$predictions,labels = test$pres,mode = "prcroc")
  
  plot(precrec_proc)
 
  modauc <- precrec::auc(precrec::evalmod(scores = pred$predictions, 
                                          labels = test$pres))
  print(i)
  rfK[[i]] <- modauc$aucs[1]
  
  
}

mean(unlist(rfK))  




##############################spatial X-validation GLM#########################################################################################

glmK<-list()# list for glm results
par(mfrow=c(2,3))


for (i in folds) {
  
  
  fishNetTrain<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id!=i)
  
  train <- st_intersection(fishNetTrain, dfST)#for presence values, select all folds which are not 'i' to train the model
  
  
  fishNetTest<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id==i)
  
  test <- st_intersection(fishNetTest, dfST) 
  
  nPres<-nrow(train[train$pres==1,])
  
  nBack<-nrow(train[train$pres==0,])
  
  
  glmGeo<-glm(pres~streamPower+poly(width,2)+gradient+terrainSlope+minSlope
              +poly(pasture100,2)+poly(urban100,2)+poly(wood100,4),family = "quasibinomial",weights=(nPres/nBack)^(1-pres),
              data=train)
  
  
  pred <- predict(glmGeo, as.data.frame(test),type="response")
  
  precrec_proc <- evalmod(scores = pred,labels = test$pres,mode = "prcroc")
  
  plot(precrec_proc)
  
  predValsCurve <- roc(test$pres, pred)
 
  
  modauc <- precrec::auc(precrec::evalmod(scores = pred, 
                                          labels = test$pres))
  
  glmK[[i]] <- modauc$aucs[1]
  
  print(i)
  
}

mean(unlist(glmK),rm.na=T)


######==================================threshold estimation===========
glmTr<-list()# list for glm results
par(mfrow=c(2,3))


for (i in folds) {
  
  
  fishNetTrain<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id!=i)
  
  train <- st_intersection(fishNetTrain, dfST)#for presence values, select all folds which are not 'i' to train the model
  
  
  fishNetTest<-subset(fishnet_grid_Fin,fishnet_grid_Fin$grid_id==i)
  
  test <- st_intersection(fishNetTest, dfST) 
  
  nPres<-nrow(train[train$pres==1,])
  
  nBack<-nrow(train[train$pres==0,])
  
  
  glmGeo<-glm(pres~streamPower+poly(width,2)+gradient+terrainSlope+minSlope
              +poly(pasture100,2)+poly(urban100,2)+poly(wood100,4),family = "quasibinomial",weights=(nPres/nBack)^(1-pres),
              data=train)
  
  
  pred <- predict(glmGeo, as.data.frame(test),type="response")
  
  #use precrec library's evalmod() function to obtain area under curve 
  precrec_proc <- evalmod(scores = pred,labels = test$pres,mode = "basic")
  
  #extract proc thresholds, specificity and sensitivity outputs 
  
  proc_df=data.frame(cbind(precrec_proc$score[[1]]$y,precrec_proc$sp[[1]]$x,precrec_proc$sn[[1]]$y))
  
  colnames(proc_df)<-c("Tr","Specificity","Sensitivity")
  
  #plot
  plot(proc_df$Specificity,proc_df$Sensitivity)
  
  #isolate threshold as that which maximises sensitivity and minimises specificity
  TrMax<-proc_df[which.max(proc_df$Sensitivity-proc_df$Specificity),]
  
  TrMax=TrMax$Tr
  
  glmTr[[i]] <- TrMax
  
  
  
}

#mean optimal threshold


mean(unlist(glmTr))



