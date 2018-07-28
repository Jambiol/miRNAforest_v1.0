import rpy2.robjects as robjects

robjects.r('''
targetSeq <- function(miRNAseq, utr){

        library(Biostrings)

        totalRegion <- NA
        match_8mer <- 0
        match_7mer_m8 <- 0
        match_7mer_A1 <- 0
        match_6mer <- 0
        GU_wobble <- 0

        target <- gsub('T', 'U', utr)

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
                        match_8mer <- 1
                        pos <- pos[1]
                        totalRegion <- substr(target, pos-12, pos+17)
                }
        }else{
                if(grepl(mer7_m8_rev, target)){
                        pos <- gregexpr(mer7_m8_rev, target)[[1]]
                        pos <- pos[pos > 13]
                        pos <- pos[pos+16 < nchar(target)]
                        if(length(pos)>0){
                                match_7mer_m8 <- 1
                                pos <- pos[1]
                                totalRegion <- substr(target, pos-13, pos+16)
                        }

                }else{
                        if(grepl(mer7_A1_rev, target)){
                                pos <- gregexpr(mer7_A1_rev, target)[[1]]
                                pos <- pos[pos > 13]
                                pos <- pos[pos+16 < nchar(target)]
                                if(length(pos)>0){
                                        match_7mer_A1 <- 1
                                        pos <- pos[1]
                                        totalRegion <- substr(target, pos-13, pos+16)
                                }
                        }else{
                                if(grepl(mer6_rev, target)){
                                        pos <- gregexpr(mer6_rev, target)[[1]]
                                        pos <- pos[pos > 14]
                                        pos <- pos[pos+15 < nchar(target)]
                                        if(length(pos)>0){
                                                match_6mer <- 1
                                                pos <- pos[1]
                                                totalRegion <- substr(target, pos-14, pos+15)
                                        }
                                }else{                                        
                                        totalRegion <- NA
                                }
                        }
                }
        }
        if(is.na(totalRegion)){
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
                        totalRegion <- theoricalRegion[which.min(GU_dist)] 
                        GU_wobble <- GU_dist[which.min(GU_dist)]

                        if(which.min(GU_dist) == 1){match_8mer <- 1}
                        if(which.min(GU_dist) == 2){match_7mer_m8 <- 1}
                        if(which.min(GU_dist) == 3){match_7mer_A1 <- 1}
                        if(which.min(GU_dist) == 4){match_6mer <- 1}
                        }
        }
        
        GU_wobble <- sapply(GU_wobble, function(x) x/5)
        
        totalRegionData <- as.data.frame(cbind(match_8mer, match_7mer_m8, match_7mer_A1, 
                                               match_6mer, GU_wobble, totalRegion))
        return(totalRegionData)
}

oneHot <- function(miRNA, target){

        ntDict <- vector(mode="list", length=5)
        names(ntDict) <- c("A", "U", "C", "G", "N")
        
        ntbin <- rep(0,4)
        ntDict[[5]] <- ntbin
        for(i in 1:4){
                nt <- ntbin
                nt[i] <- 1
                ntDict[[i]] <- nt
        }
        
        ## miRNA
        miRNA <- as.character(miRNA)
        ntseq <- strsplit(miRNA, split="")[[1]]
        
        if(length(ntseq) < 25){
                ntseq[(length(ntseq)+1):25] <- "N"
        }
        
        encodelist <- vector(mode="list", length=length(ntseq))
        for(i in 1:length(ntseq)){
                encodelist[[i]]<- ntDict[[ntseq[i]]]
        }
        
        mirnaOneHot <- as.vector(unlist(encodelist))
        
        ## target
        totalRegion <- as.character(target[[6]])
        ntseq <- strsplit(totalRegion, split="")[[1]]
        
        encodelist <- vector(mode="list", length=length(ntseq))
        for(i in 1:length(ntseq)){
                encodelist[[i]]<- ntDict[[ntseq[i]]]
        }
        
        targetOneHot <- as.vector(unlist(encodelist))
        
        seed <- as.numeric(as.vector(unlist(target[1:5])))
        
        t(c(seed, mirnaOneHot, targetOneHot))
}
	''')

targetSeq = robjects.r['targetSeq']
oneHot = robjects.r['oneHot']