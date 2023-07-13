## Library calls, load gene lists and cosmic signatures ##----------------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(deconstructSigs)
library(data.table)
library(dplyr)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") # Load repair genes as list objects
load("C:/Users/robmc/Desktop/NUMB_files/data/signatures.genome.cosmic.v3.may2019.rda")
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files" # Create root data directory pointer 

## processData function ##----------------------
processData <- function(directory) {
  setwd(paste0(homeDir, "/data/", directory))
  setwd(homeDir)
  setwd(paste0(homeDir,"/data/skcm_tcga"))
  
  # Create initial clinical file by merging patient and sample data
  sample <- fread("data_clinical_sample.txt", skip=4)
  patient <- fread("data_clinical_patient.txt", skip=4)
  clin <- merge(sample, patient, by = "PATIENT_ID")
  # Remove clinical features which only contain NA values
  clin[clin == "[Not Available]"] <- NA
  clin <- Filter(function(x)!all(is.na(x)), clin)
  
  # Read in mutation data, accounting for both formats in cBioPortal
  if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
  } else { mutations <- fread("data_mutations.txt")}
  # Make sure all samples with mutations to be analyzed have associated clinical annotations
  mutations <- mutations[mutations$Tumor_Sample_Barcode %in% clin$SAMPLE_ID, ]
  
  # Format chromosome column to match BSgenome format
  mutations$Chromosome <- sub("^", "chr", mutations$Chromosome)
  # Prepare input object for deconstructSigs package
  sigs <- mut.to.sigs.input(mut.ref = mutations, sample.id = "Tumor_Sample_Barcode", chr = "Chromosome",
                            pos = "Start_Position", ref = "Reference_Allele", 
                            alt = "Tumor_Seq_Allele2", bsg = BSgenome.Hsapiens.UCSC.hg19)
  # Filter out samples with low mutation count as suggested by deconstructSigs package
  # This number can be changed, but further reading should be done to choose mutation cutoff
  finalMut <- sigs[rowSums(sigs)>50,]
  results <- vector("list", nrow(finalMut))
  names(results) <- row.names(finalMut)
  
  # Construct mutational profile for each tumor sample from already defined signatures (COSMIC in this case)
  for(sID in row.names(finalMut)){
    results[[sID]] <- whichSignatures(finalMut, # the matrix generated with mut.to.sigs.input
                                      sample.id=sID, # the current sample ID
                                      signatures.ref=signatures.genome.cosmic.v3.may2019, # the data.frame with the signatures that comes with deconstructSigs
                                      tri.counts.method="default", # which normalization method to use
                                      contexts.needed=TRUE) # set to TRUE if your input matrix contains counts instead of frequencies
  }
  # convert the exposures for each sample into a sample x signatures matrix
  expo <- do.call("rbind", sapply(results, "[", 1))
  # add the unknown value to the matrix such that the contributions add up to 1 per sample
  Signature.unknown <- unlist(sapply(results, "[", 5))
  expo <- cbind(expo, Signature.unknown)
  # Create combined UV exposure signature column
  expo$UV_sig <- rowSums(expo[ , c("SBS7a", "SBS7b", "SBS7c", "SBS7d")])
  expo$SAMPLE_ID <- rownames(expo)
  expo <- expo[order(-expo$UV_sig),]
  expo$SAMPLE_ID <- gsub(".weights","",expo$SAMPLE_ID)
  # Output imputed signatures with COSMIC signature names for downstream analysis / plot generation
  write.csv(expo, "mutationalSignatures.csv")
  
  # Create new UV signature annotation column based on % of mutational signature contribution by UV exposure
  clinFinal <- clin[which(clin$SAMPLE_ID %in% expo$SAMPLE_ID), ]
  clinFinal$UV_sig_value <- expo$UV_sig
  clinFinal$UV_sig[clinFinal$UV_sig_value >= 0.65] <- "High"
  clinFinal$UV_sig[clinFinal$UV_sig_value < 0.65 &
                     clinFinal$UV_sig_value > 0.2] <- "Moderate"
  clinFinal$UV_sig[clinFinal$UV_sig_value <= 0.2] <- "Low"
  clinFinal <- clinFinal[order(-clinFinal$UV_sig_value),]
  # Create tumor_sample_barcode column which is necessary for some maftools package functions
  clinFinal$Tumor_Sample_Barcode <- clinFinal$SAMPLE_ID
  
  mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B"  # Combine CUL4A and CUL4B for downstream analysis
  mutations <- mutations[mutations$Hugo_Symbol %in% NERgenes, ] # Isolate NER mutations
  mutations <- mutations[mutations$Variant_Classification != "Silent", ] # Filter out silent mutations
  
  # Parse through each sample with available UV signature data, and note mutation status of each NER gene (0 = WT, 1 = Mut)
  clinMut <- as.data.frame(clinFinal$Tumor_Sample_Barcode)
  colnames(clinMut)[1] <- "Tumor_Sample_Barcode"
  for (i in NERgenes) {
    mut <- mutations[mutations$Hugo_Symbol == i]
    clinMut[[i]] <- 0
    clinMut[[i]][clinMut$Tumor_Sample_Barcode %in% mut$Tumor_Sample_Barcode] <- 1
  }
  write.csv(clinMut, "NER_mutSig.csv", row.names = F)
  
  # Parse through each sample with available CNV data, and denote CNV status of each NER gene (-2,-1 = Del, 0 = WT, 1,2 = Amp)
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
  
  m <- clinMut[,2:26]
  rownames(m) <- clinMut$Tumor_Sample_Barcode
  m <- m*0.1
  
  c <- clinCNV[,2:26]
  rownames(c) <- clinCNV$Tumor_Sample_Barcode
  
  both <- bind_rows(c %>% tibble::rownames_to_column(), 
                    m %>% tibble::rownames_to_column()) %>% 
    # evaluate following calls for each value in the rowname column
    group_by(rowname) %>% 
    # add all non-grouping variables
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
  write.csv(allDat, "clinical.csv")
}