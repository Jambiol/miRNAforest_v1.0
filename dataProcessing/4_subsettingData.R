#Es fixa el directori de treball el directori on es troba l'arxiu .R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "3_oneHot_encode.RData")

set.seed(1980)
negData <- data[data$target == 0,]
posData <- data[data$target == 1,]

negIndex <- sample(1:dim(negData)[1], dim(negData)[1]/2)
posIndex <- sample(1:dim(posData)[1], dim(posData)[1]/2)

negData_features <- negData[negIndex,]
posData_features <- posData[posIndex,]
featuresData <- rbind(negData_features, posData_features)
write.csv(featuresData, file = "data/featuresDataNorm.csv", quote = FALSE, row.names = FALSE)

trainNegIndex <- sample(1:dim(negData_features)[1], round(dim(negData_features)[1]*.75))
trainPosIndex <- sample(1:dim(posData_features)[1], round(dim(posData_features)[1]*.75))

trainFeatures <- rbind(negData_features[trainNegIndex,], posData_features[trainPosIndex,])
testFeatures <- rbind(negData_features[-trainNegIndex,], posData_features[-trainPosIndex,])

write.csv(trainFeatures[,12:231], file = "data/trainFeaturesNorm.csv", quote = FALSE, row.names = FALSE)
write.csv(testFeatures[,12:231], file = "data/testFeaturesNorm.csv", quote = FALSE, row.names = FALSE)

negData_predict <- negData[-negIndex,]
posData_predict <- posData[-posIndex,]
predictData <- rbind(negData_predict, posData_predict)
write.csv(predictData, file = "data/predictDataNorm.csv", quote = FALSE, row.names = FALSE)

trainNegIndex <- sample(1:dim(negData_predict)[1], round(dim(negData_predict)[1]*.75))
trainPosIndex <- sample(1:dim(posData_predict)[1], round(dim(posData_predict)[1]*.75))

trainPredict <- rbind(negData_predict[trainNegIndex,], posData_predict[trainPosIndex,])
testPredict <- rbind(negData_predict[-trainNegIndex,], posData_predict[-trainPosIndex,])

write.csv(trainPredict, file = "data/trainPredictNorm.csv", quote = FALSE, row.names = FALSE)
write.csv(testPredict, file = "data/testPredictNorm.csv", quote = FALSE, row.names = FALSE)

save(data, featuresData, predictData, file = "4_subsettingData.RData")
rm(list=ls())
load(file = "4_subsettingData.RData")
