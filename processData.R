## Library calls, load gene lists and cosmic signatures ##----------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(deconstructSigs)
library(data.table)
library(dplyr)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R")
load("C:/Users/robmc/Desktop/NUMB_files/data/signatures.genome.cosmic.v3.may2019.rda")
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"

## processData function ##----------------------
processData <- function(directory, filter = F, minMut = 30) {
  setwd(paste0(homeDir, "/data/", directory))
  
  sample <- fread("data_clinical_sample.txt", skip=4)
  patient <- fread("data_clinical_patient.txt", skip=4)
  clin <- merge(sample, patient, by = "PATIENT_ID")
  clin[clin == "[Not Available]"] <- NA
  clin <- Filter(function(x)!all(is.na(x)), clin)
  
  # filter format should be csv of sample names, with header as first line
  if (filter == T) {
    f <- read.csv("filter.csv")
    clin <- clin[clin$SAMPLE_ID %in% f[,1]]
  }
  
  if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
  } else { mutations <- fread("data_mutations.txt")}
  mutations <- mutations[mutations$Tumor_Sample_Barcode %in% clin$SAMPLE_ID, ]
  
  # Format chromosome column to match BSgenome format
  mutations$Chromosome <- sub("^", "chr", mutations$Chromosome)
  sigs <- mut.to.sigs.input(mut.ref = mutations, sample.id = "Tumor_Sample_Barcode", chr = "Chromosome",
                            pos = "Start_Position", ref = "Reference_Allele", 
                            alt = "Tumor_Seq_Allele2", bsg = BSgenome.Hsapiens.UCSC.hg19)
  finalMut <- sigs[rowSums(sigs)>minMut,]
  write.csv(finalMut, "tnm.csv")
  results <- vector("list", nrow(finalMut))
  names(results) <- row.names(finalMut)
  
  for(sID in row.names(finalMut)){
    results[[sID]] <- whichSignatures(finalMut,
                                      sample.id=sID,
                                      signatures.ref=signatures.genome.cosmic.v3.may2019,
                                      tri.counts.method="default",
                                      contexts.needed=TRUE)
  }
  expo <- do.call("rbind", sapply(results, "[", 1))
  Signature.unknown <- unlist(sapply(results, "[", 5))
  expo <- cbind(expo, Signature.unknown)
  # Create combined UV exposure signature column
  expo$UV_sig <- rowSums(expo[ , c("SBS7a", "SBS7b", "SBS7c", "SBS7d")])
  expo$SAMPLE_ID <- rownames(expo)
  expo <- expo[order(-expo$UV_sig),]
  expo$SAMPLE_ID <- gsub(".weights","",expo$SAMPLE_ID)
  write.csv(expo, "mutationalSignatures.csv", row.names = F)
  
  clinFinal <- clin[which(clin$SAMPLE_ID %in% expo$SAMPLE_ID), ]
  expo <- expo[order(expo$SAMPLE_ID),]
  clinFinal$UV_sig_value <- expo$UV_sig
  clinFinal$UV_sig[clinFinal$UV_sig_value >= 0.65] <- "High"
  clinFinal$UV_sig[clinFinal$UV_sig_value < 0.65 &
                     clinFinal$UV_sig_value > 0.2] <- "Moderate"
  clinFinal$UV_sig[clinFinal$UV_sig_value <= 0.2] <- "Low"
  clinFinal$Tumor_Sample_Barcode <- clinFinal$SAMPLE_ID
  
  mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B"
  mutations <- mutations[mutations$Hugo_Symbol %in% NERgenes, ]
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  
  clinMut <- as.data.frame(clinFinal$Tumor_Sample_Barcode)
  colnames(clinMut)[1] <- "Tumor_Sample_Barcode"
  for (i in NERgenes) {
    mut <- mutations[mutations$Hugo_Symbol == i]
    clinMut[[i]] <- 0
    clinMut[[i]][clinMut$Tumor_Sample_Barcode %in% mut$Tumor_Sample_Barcode] <- 1
  }
  write.csv(clinMut, "NER_mutSig.csv", row.names = F)
  
  cnv <- fread("data_cna.txt")
  cnv <- unique(cnv[!duplicated(cnv$Hugo_Symbol), ])
  cnv2 <- reshape2::melt(cnv, id.vars = "Hugo_Symbol")
  colnames(cnv2) <- c("Gene", "Sample_name", "CN")
  cnv2 <- cnv2[cnv2$Sample_name %in% clinFinal$Tumor_Sample_Barcode, ]
  cnv2$Gene[cnv2$Gene == "CUL4A" | cnv2$Gene == "CUL4B"] <- "CUL4A/B"
  
  clinCNV <- as.data.frame(clinFinal$Tumor_Sample_Barcode[clinFinal$Tumor_Sample_Barcode %in% cnv2$Sample_name])
  colnames(clinCNV)[1] <- "Tumor_Sample_Barcode"
  for (i in NERgenes) {
    cnv <- cnv2[cnv2$Gene == i, ]
    clinCNV[[i]] <- 0
    clinCNV[[i]][match(cnv$Sample_name, clinCNV$Tumor_Sample_Barcode)] <- cnv$CN
  }
  write.csv(clinCNV, "NER_cnvSig.csv", row.names = F)
  
  m <- data.frame(clinMut, row.names = 1)
  m <- m*0.1
  c <- data.frame(clinCNV, row.names = 1)
  
  both <- bind_rows(c %>% tibble::rownames_to_column(), 
                    m %>% tibble::rownames_to_column()) %>% 
    group_by(rowname) %>% 
    summarise_all(sum)
  
  both[both == 0 | both == 1 | both == 2] <- 0
  both[both == -0.9 | both == -1.9] <- 1
  both[both == 0.1] <- 1
  both[both == -1 | both == -2] <- 1
  both[both == 1.1 | both == 2.1] <- 1
  
  nerSig <- both[ , colnames(both) %in% c("CETN2", "GTF2H2", "ERCC8", "CDK7", "CCNH", "ERCC6", "RAD23B")]
  nerSig$NER_sig <- with(nerSig, rowSums(nerSig))
  nerSig$Tumor_Sample_Barcode <- both$rowname
  keep <- nerSig[,c("Tumor_Sample_Barcode", "NER_sig")]
  
  allDat <- merge(clinFinal, keep, by = "Tumor_Sample_Barcode")
  allDat <- allDat[order(-allDat$UV_sig_value),]
  write.csv(allDat, "clinical.csv", row.names = F)
}

## Execute code -------------------
processData("skcm_tcga", filter = T)