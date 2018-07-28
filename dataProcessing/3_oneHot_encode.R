#Es fixa el directori de treball el directori on es troba l'arxiu .R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "2_target.RData")

#Es defineix la funció oneHot, que implementarà la codificació one-hot encoding
oneHot <- function(x, length){
        
        #En primer lloc, la funció incorpora una llista que ens servirà de diccionari. En aquest llista, 
        #cada nucleòtid portarà associat un one-hot encoding. El caràcter 'o' ens servirà per omplir els
        #nucleòtids buit fins a 22 nt
        ntDict <- vector(mode="list", length=5)
        names(ntDict) <- c("A", "U", "C", "G", "N")
        
        #Es generen els 5 one-hot encoding
        ntbin <- rep(0,4)
        ntDict[[5]] <- ntbin
        for(i in 1:4){
                nt <- ntbin
                nt[i] <- 1
                ntDict[[i]] <- nt
        }
        
        #La seqüència d'entrada es divideix en els diferents caràcters, i afegim el caràcter o que ens servirà
        #per omplir els nucleòtids absents
        x <- toupper(x)
        ntseq <- strsplit(x, split="")[[1]]
        
        if(length(ntseq) < length){
                ntseq[(length(ntseq)+1):length] <- "N"
        }
        
        #La funció crea una nova llista, on es guarda el one-hot encoding del miRNA
        encodelist <- vector(mode="list", length=length(ntseq))
        for(i in 1:length(ntseq)){
                encodelist[[i]]<- ntDict[[ntseq[i]]]
        }
        
        #Finalment, el miRNA queda codificat en un vector de 22x4 variables
        as.vector(unlist(encodelist))
}

#miRNA one hot
mirna_seq <- data$miRNA_seq

SeqsEncode <- vector(mode="list", length=length(mirna_seq))
names(SeqsEncode) <- mirna_seq

len <- max(nchar(mirna_seq))
for(i in 1:length(mirna_seq)){
        SeqsEncode[[i]] <- oneHot(mirna_seq[i], len)
}

miRNAseqMatrix <- t(as.data.frame(SeqsEncode))
miRNA_oneHot <- as.data.frame(miRNAseqMatrix)
names(miRNA_oneHot) <- paste0("miR_OH", 1:length(names(miRNA_oneHot)))
rownames(miRNA_oneHot) <- NULL

#mRNA one hot
totalRegion_seq <- data$totalRegion

SeqsEncode <- vector(mode="list", length=length(totalRegion_seq))
names(SeqsEncode) <- totalRegion_seq

len <- max(nchar(totalRegion_seq))
for(i in 1:length(totalRegion_seq)){
        SeqsEncode[[i]] <- oneHot(totalRegion_seq[i], len)
}

mRNAseqMatrix <- t(as.data.frame(SeqsEncode))
mRNA_oneHot <- as.data.frame(mRNAseqMatrix)
rownames(mRNA_oneHot) <- NULL
names(mRNA_oneHot) <- paste0("mRNA_OH", 1:length(names(mRNA_oneHot)))

data <- cbind(data, miRNA_oneHot, mRNA_oneHot)
write.csv(data, file = "data/dataOneHot.csv", quote = FALSE, row.names = FALSE)

save(data, file = "dataOneHot.RData")
rm(list=ls())
load(file = "dataOneHot.RData")

set.seed(1980)
negData <- data[data$target == 0,]
posData <- data[data$target == 1,]

negIndex <- sample(1:dim(negData)[1], dim(negData)[1]/2)
posIndex <- sample(1:dim(posData)[1], dim(posData)[1]/2)

negData_features <- negData[negIndex,]
posData_features <- posData[posIndex,]
featuresData <- rbind(negData_features, posData_features)

trainNegIndex <- sample(1:dim(negData_features)[1], round(dim(negData_features)[1]*.75))
trainPosIndex <- sample(1:dim(posData_features)[1], round(dim(posData_features)[1]*.75))

trainFeatures <- rbind(negData_features[trainNegIndex,], posData_features[trainPosIndex,])
testFeatures <- rbind(negData_features[-trainNegIndex,], posData_features[-trainPosIndex,])

trainFeatures[1,12:191]
testFeatures[,12:191]

train_miR <- trainFeatures[,12:111]
test_miR <- testFeatures[,12:111]

train_mRNA <- trainFeatures[,112:191]
test_mRNA <- testFeatures[,112:191]

write.csv(train_mRNA, file = "data/train_mRNA.csv", quote = FALSE, row.names = FALSE)
write.csv(test_mRNA, file = "data/test_mRNA.csv", quote = FALSE, row.names = FALSE)

write.csv(trainFeatures[,12:191], file = "data/trainFeatures.csv", quote = FALSE, row.names = FALSE)
write.csv(testFeatures[,12:191], file = "data/testFeatures.csv", quote = FALSE, row.names = FALSE)

write.csv(train_miR, file = "data/train_miR.csv", quote = FALSE, row.names = FALSE)
write.csv(test_miR, file = "data/test_miR.csv", quote = FALSE, row.names = FALSE)
        
negData_predict <- negData[-negIndex,]
posData_predict <- posData[-posIndex,]
predictData <- rbind(negData_predict, posData_predict)

save(data, featuresData, predictData, file = "3_oneHot_NormWobble.RData")
rm(list=ls())
load(file = "3_oneHot_NormWobble.RData")

