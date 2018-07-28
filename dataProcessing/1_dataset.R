#Es fixa el directori de treball el directori on es troba l'arxiu .R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##PROCESSAT DE LES DADES POSITIVES

#Es llegeix l'arxiu csv que conté les dades positives
positive <- read.csv("rawData/Positive_data.csv", header = TRUE)

require(biomaRt)
require(org.Hs.eg.db)

#S'obtenen les seq 3'UTR dels diferents transcripts
refseq <- unique(positive$mRNA_ID)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
utr3_pos <- getSequence(id=refseq, type="refseq_mrna", seqType="3utr", mart=ensembl)

#S'obté el nom del gen de cada transcript
gene_symbol <- getBM(attributes=c('hgnc_symbol', 'refseq_mrna'),filters = 'refseq_mrna', 
                     values = refseq, mart = ensembl)
positive <- merge(positive, gene_symbol, by.x = "mRNA_ID", by.y = "refseq_mrna")

#Es llegeix la seqüència dels miRNAs
library("seqinr")

mirna <- read.fasta("rawData/mature.fa")
hsa_mirna <- mirna[grepl("^hsa", names(mirna))]

miR_ID <- c()
miRNA_seq <- c()
for(i in 1:length(hsa_mirna)){
        miR_ID[i] <- attr(hsa_mirna[[i]], "name")
        seq <- c()
        for(j in 1:length(hsa_mirna[[i]])){
                seq[j] <- hsa_mirna[[i]][j]
        }
        miRNA_seq[i] <- paste(seq, sep="", collapse ="")
}
mirna_seq <- data.frame(cbind(miR_ID, miRNA_seq))

#Es creuen les dades, de manera que s'obté una taula amb cada combinació de miRNA i seq 3'UTR
positive.data <- merge(positive, utr3_pos, by.x = "mRNA_ID", by.y = "refseq_mrna")
positive.data <- merge(positive.data, mirna_seq, by = "miR_ID")
positive.data <- positive.data[,c(1, 5, 3, 2, 4)]
colnames(positive.data)[3] <- "gene_symbol"

positive.data <- positive.data[!(positive.data$`3utr` == "Sequence unavailable"),]
target <- rep(1, dim(positive.data)[1])
positive.data <- cbind(target, positive.data)

#Nombre de transcripts
length(unique(positive.data$mRNA_ID))

#Nombre de gens
length(unique(positive.data$gene_symbol))

#Nombre de miRNA
length(unique(positive.data$miR_ID))


##PROCESSAT DE LES DADES NEGATIVES

#Es llegeix l'arxiu csv que conté les dades negatives
negative <- read.csv("rawData/negativeTranscriptsTarbase.csv", sep = ";", header = TRUE)

require(biomaRt)
require(org.Hs.eg.db)

#S'obtenen les seq 3'UTR dels diferents transcripts
ensemblID <- unique(negative$EnsemblId)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
utr3_neg <- getSequence(id=ensemblID, type="ensembl_gene_id", seqType="3utr", mart=ensembl)

negative <- merge(negative, utr3_neg, by.x = "EnsemblId", by.y = "ensembl_gene_id")

#Es llegeix la seqüència dels miRNAs
library("seqinr")

mirna <- read.fasta("rawData/mature.fa")
hsa_mirna <- mirna[grepl("^hsa", names(mirna))]

miR_ID <- c()
miRNA_seq <- c()
for(i in 1:length(hsa_mirna)){
        miR_ID[i] <- attr(hsa_mirna[[i]], "name")
        seq <- c()
        for(j in 1:length(hsa_mirna[[i]])){
                seq[j] <- hsa_mirna[[i]][j]
        }
        miRNA_seq[i] <- paste(seq, sep="", collapse ="")
}
mirna_seq <- data.frame(cbind(miR_ID, miRNA_seq))

#Es creuen les dades, de manera que s'obté una taula amb cada combinació de miRNA i seq 3'UTR
negative.data <- merge(negative, mirna_seq, by.x = "miRNA", by.y = "miR_ID")
negative.data <- negative.data[,c(4, 1, 6, 3, 2, 5)]

#S'eliminen exemples que no tenen seq de 3'UTR o seq massa curtes
negative.data <- negative.data[!(negative.data$`3utr` == "Sequence unavailable"),]
negative.data <- negative.data[!(nchar(negative.data$`3utr`) < 22),]
dim(negative.data)

#Nombre de transcripts
length(unique(negative.data$`3utr`))

#Nombre de gens
length(unique(negative.data$gene_name))

#Nombre de miRNA
length(unique(negative.data$miRNA))

##Fusió de les dades positives i negatives
colnames(negative.data) <- colnames(positive.data)

dataset <- rbind(positive.data, negative.data)
dataset$target <- as.factor(dataset$target)
dataset$miR_ID <- as.character(dataset$miR_ID)
dataset$miRNA_seq <- as.character(dataset$miRNA_seq)
dataset$miRNA_seq <- toupper(dataset$miRNA_seq)

save(dataset, negative.data, positive.data, file = "1_dataset.RData")
rm(list=ls())
load(file = "1_dataset.RData")

write.csv(dataset, file = "data/1_dataset.csv", quote = FALSE, row.names = FALSE)
