#Es fixa el directori de treball el directori on es troba l'arxiu .R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "1_dataset.RData")

library(Biostrings)

totalRegion <- rep(NA, dim(dataset)[1])
match_8mer <- rep(0, dim(dataset)[1])
match_7mer_m8 <- rep(0, dim(dataset)[1])
match_7mer_A1 <- rep(0, dim(dataset)[1])
match_6mer <- rep(0, dim(dataset)[1])
GU_wobble <- rep(0, dim(dataset)[1])

miRNA_seq <- 'AAAAGCUGGGUUGAGAGGGCGA'
`3utr` <- 'GACCCCTGCAAATGGGACAGCCCCCCAGTTCATGAGGCCTGGCTATTGCAATATTTACTAGTAGAGGAACTCTATAGCAAGATGAAGAGGAAAAACAAACAAACAAACAAAAAAAAACACAAAAAAAGAAAGAATACTTTTTTATACCTCACTATGTTCTTTGAATATGTATTTTTCCTTTAAATTTCTGCCTTTAATTCTTTTGTTCCAAA'

dataset <- as.data.frame(cbind(miRNA_seq, `3utr`))
i <- 1
for(i in 1:dim(dataset)[1]){
        
        target <- dataset$`3utr`[i]
        target <- gsub('T', 'U', target)
        
        miRNAseq <- dataset$miRNA_seq[i]
        
        mer6 <- RNAString(substr(miRNAseq, 2, 7))
        mer6_rev <- as.character(reverseComplement(mer6))
        
        mer7_A1_rev <- paste0("A", mer6_rev)
        
        mer7_m8 <- RNAString(substr(miRNAseq, 2, 8))
        mer7_m8_rev <- as.character(reverseComplement(mer7_m8))
        
        mer8_rev <- paste0("A", mer7_m8_rev)
        
        if(grepl(mer8_rev, target)){
                pos <- gregexpr(mer8_rev, target)[[1]]
                pos <- pos[pos > 12]
                pos <- pos[pos+17 < nchar(target)]
                if(length(pos)>0){
                        match_8mer[i] <- 1
                        pos <- pos[1]
                        totalRegion[i] <- substr(target, pos-12, pos+17)
                }
        }else{
                if(grepl(mer7_m8_rev, target)){
                        pos <- gregexpr(mer7_m8_rev, target)[[1]]
                        pos <- pos[pos > 13]
                        pos <- pos[pos+16 < nchar(target)]
                        if(length(pos)>0){
                                match_7mer_m8[i] <- 1
                                pos <- pos[1]
                                totalRegion[i] <- substr(target, pos-13, pos+16)
                        }
                        
                }else{
                        if(grepl(mer7_A1_rev, target)){
                                pos <- gregexpr(mer7_A1_rev, target)[[1]]
                                pos <- pos[pos > 13]
                                pos <- pos[pos+16 < nchar(target)]
                                if(length(pos)>0){
                                        match_7mer_A1[i] <- 1
                                        pos <- pos[1]
                                        totalRegion[i] <- substr(target, pos-13, pos+16)
                                }
                        }else{
                                if(grepl(mer6_rev, target)){
                                        pos <- gregexpr(mer6_rev, target)[[1]]
                                        pos <- pos[pos > 14]
                                        pos <- pos[pos+15 < nchar(target)]
                                        if(length(pos)>0){
                                                match_6mer[i] <- 1
                                                pos <- pos[1]
                                                totalRegion[i] <- substr(target, pos-14, pos+15)
                                        }
                                }else{                                        
                                        totalRegion[i] <- NA
                                }
                        }
                }
        }
        if(is.na(totalRegion[i])){
                mer6_rev_GU <- gsub('A', '[AG]', mer6_rev)
                mer6_rev_GU <- gsub('C', '[CU]', mer6_rev_GU)
                
                mer7_A1_rev_GU <- paste0(mer6_rev_GU, "A")
                
                mer7_m8_rev_GU <- gsub('A', '[AG]', mer7_m8_rev)
                mer7_m8_rev_GU <- gsub('C', '[CU]', mer7_m8_rev_GU)
                
                mer8_rev_GU <- paste0(mer7_m8_rev_GU, "A")
                
                GU_dist <- c()
                coordinate <- c()
                theoricalRegion <- c()
                if(grepl(mer8_rev_GU, target)){
                        pos <- gregexpr(mer8_rev_GU, target)[[1]]
                        pos <- pos[pos > 12]
                        pos <- pos[pos+17 < nchar(target)]
                        if(length(pos)>0){
                                dist <- c()
                                for(j in 1:length(pos)){
                                        theorical_target <- substr(target, pos[j], pos[j]+7)
                                        dist[j] <- adist(mer8_rev, theorical_target)[1]
                                }
                                GU_dist[1] <- dist[which.min(dist)]
                                coordinate[1] <- pos[which.min(dist)]
                                theoricalRegion[1] <- substr(target, coordinate[1]-12, coordinate[1]+17)
                        }
                }
                if(grepl(mer7_m8_rev_GU, target)){
                        pos <- gregexpr(mer7_m8_rev_GU, target)[[1]]
                        pos <- pos[pos > 13]
                        pos <- pos[pos+16 < nchar(target)]
                        if(length(pos)>0){
                                dist <- c()
                                for(j in 1:length(pos)){
                                        theorical_target <- substr(target, pos[j], pos[j]+6)
                                        dist[j] <- adist(mer7_m8_rev, theorical_target)[1]
                                }
                                GU_dist[2] <- dist[which.min(dist)]
                                coordinate[2] <- pos[which.min(dist)]
                                theoricalRegion[2] <- substr(target, coordinate[2]-13, coordinate[2]+16)
                        }
                }
                if(grepl(mer7_A1_rev_GU, target)){
                        pos <- gregexpr(mer7_A1_rev_GU, target)[[1]]
                        pos <- pos[pos > 13]
                        pos <- pos[pos+16 < nchar(target)]
                        if(length(pos)>0){
                                dist <- c()
                                for(j in 1:length(pos)){
                                        theorical_target <- substr(target, pos[j], pos[j]+6)
                                        dist[j] <- adist(mer7_A1_rev, theorical_target)[1]
                                }
                                GU_dist[3] <- dist[which.min(dist)]
                                coordinate[3] <- pos[which.min(dist)]
                                theoricalRegion[3] <- substr(target, coordinate[3]-13, coordinate[3]+16)
                        }
                }
                if(grepl(mer6_rev_GU, target)){
                        pos <- gregexpr(mer6_rev_GU, target)[[1]]
                        pos <- pos[pos > 14]
                        pos <- pos[pos+15 < nchar(target)]
                        if(length(pos)>0){
                                dist <- c()
                                for(j in 1:length(pos)){
                                        theorical_target <- substr(target, pos[j], pos[j]+5)
                                        dist[j] <- adist(mer6_rev, theorical_target)[1]
                                }
                                GU_dist[4] <- dist[which.min(dist)]
                                coordinate[4] <- pos[which.min(dist)]
                                theoricalRegion[4] <- substr(target, coordinate[4]-14, coordinate[4]+15)
                        }
                }
                if(!is.null(coordinate)){
                        totalRegion[i] <- theoricalRegion[which.min(GU_dist)]
                        GU_wobble[i] <- as.numeric(GU_dist[which.min(GU_dist)])
                        
                        if(which.min(GU_dist) == 1){match_8mer[i] <- 1}
                        if(which.min(GU_dist) == 2){match_7mer_m8[i] <- 1}
                        if(which.min(GU_dist) == 3){match_7mer_A1[i] <- 1}
                        if(which.min(GU_dist) == 4){match_6mer[i] <- 1}
                }
        }
}

totalRegionData <- as.data.frame(cbind(match_8mer, match_7mer_m8, match_7mer_A1, match_6mer, GU_wobble, totalRegion))
head(totalRegionData)

data <- dataset[, -6]
data <- cbind(data, totalRegionData[, c(6,1:5)])
head(data)

data$totalRegion <- as.character(data$totalRegion)

sum(is.na(data[data$target == 0,]))
dim(data[data$target == 0,])

sum(is.na(data[data$target == 1,]))
dim(data[data$target == 1,])

data <- na.omit(data)

dim(data[data$target == 0,])
dim(data[data$target == 1,])

negData <- data[data$target == 0,]
posData <- data[data$target == 1,]

save(data, negData, posData, file = "2_target.RData")
rm(list=ls())
load(file = "2_target.RData")

