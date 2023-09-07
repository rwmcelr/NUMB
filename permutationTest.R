### Library calls/ scripts -------------------------------------
library(ggplot2)
library(ggrepel)
library(colorspace)
library(data.table)
library(reshape2)
library(dplyr)
library(qvalue)
homeDir<- "C:/Users/robmc/Desktop/NUMB_files"
'%notin%' <- Negate('%in%')

### Functions -------------------------------------
dataPrep <- function(dataDir) {
  dataDir <- "skcm_tcga"
  setwd(paste0(homeDir, "/data/", dataDir))
  
  cnv <- fread("data_cna.txt")
  cnv[cnv == 1 | cnv == 2] <- 0
  cnv[cnv == -1 | cnv == -2] <- 1
  cnv <- melt(cnv)
  cnvWide <- dcast(data = cnv,formula = variable~Hugo_Symbol, fun.aggregate = mean, value.var = "value")
  colnames(cnvWide)[1] <- "Tumor_Sample_Barcode"
  
  if(file.exists("data_mutations.txt")) { mutations <- fread("data_mutations.txt") } else {
    mutations <- fread("data_mutations_extended.txt")
  }
  
  mutations <- mutations[mutations$Variant_Classification != "Silent", ]
  mut <- distinct(mutations, End_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, .keep_all = T)
  mutWide <- dcast(data = mut,formula = Tumor_Sample_Barcode~Hugo_Symbol, fun.aggregate = length,
                   value.var = "Variant_Classification")
  
  # Keep only mutation and CNV data for which UV-signature and clinical data exist
  clin <- fread("clinical.csv")
  cnvFinal <- cnvWide[cnvWide$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
  mutFinal <- mutWide[mutWide$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
  
  # Some samples have mut data and no cnv data, splitting them off to make merging easier
  mutMerge <- mutFinal[mutFinal$Tumor_Sample_Barcode %in% cnvFinal$Tumor_Sample_Barcode, ]
  mutMerge <- mutMerge[!duplicated(mutMerge$Tumor_Sample_Barcode), ]
  mutNotMerge <- mutFinal[mutFinal$Tumor_Sample_Barcode %notin% cnvFinal$Tumor_Sample_Barcode, ]
  
  # Not all genes are shared between the datasets, splitting off unique CNV genes for easier merging
  rownames(cnvFinal) <- cnvFinal$Tumor_Sample_Barcode
  cnvInMut <- cnvFinal[ ,colnames(cnvFinal) %in% colnames(mutMerge)]
  cnvNotInMut <- cnvFinal[ ,colnames(cnvFinal) %notin% colnames(mutMerge)]
  # Splitting mutMerge into columns corresponding to CNV data and columns not corresponding
  rownames(mutMerge) <- mutMerge$Tumor_Sample_Barcode
  mutColCombine <- mutMerge[ ,colnames(mutMerge) %in% colnames(cnvInMut)]
  mutColNoCombine <- mutMerge[ ,colnames(mutMerge) %notin% colnames(cnvInMut)]
  
  # Summing alterations for genes shared by CNV data and mutation data
  tmp <- aggregate(. ~ Tumor_Sample_Barcode, rbind(mutColCombine,cnvInMut), sum)
  # Adding unique CNV data back to aggregated list
  cnvNotInMut$Tumor_Sample_Barcode <- rownames(cnvNotInMut)
  tmp2 <- merge(tmp, cnvNotInMut, by = "Tumor_Sample_Barcode")
  # Adding unique mutation data back to aggregated list
  mutColNoCombine$Tumor_Sample_Barcode <- rownames(mutColNoCombine)
  tmp3 <- merge(tmp2, mutColNoCombine, by = "Tumor_Sample_Barcode")
  
  everything <- plyr::rbind.fill(tmp3, mutNotMerge)
  rownames(everything) <- everything$Tumor_Sample_Barcode
  everything <- everything[,-1]
  everything[is.na(everything)] <- 0
  round <- round(nrow(everything)*.05)
  # Only keep genes with alterations in over 5% of samples (keeping everything creates an unreasonable computational burden)
  everythingSig <- everything[ ,colSums(everything) > round]
  everything$Tumor_Sample_Barcode <- rownames(everything)
  everythingSig$Tumor_Sample_Barcode <- rownames(everythingSig)
  UV <- clin[ , c("Tumor_Sample_Barcode","UV_sig")]
  UV <- UV[!duplicated(UV$Tumor_Sample_Barcode), ]
  sigUV <- merge(everythingSig, UV, by = "Tumor_Sample_Barcode")
  write.csv(sigUV, "MutsAndDels_over5percent.csv", row.names = F)

  finalDat <- sigUV[sigUV$UV_sig != "Moderate", ]
  write.csv(finalDat, "MutsAndDels_over5pct_UVfiltered.csv", row.names = F)
}

curve_ball <- function(m) {  # From https://doi.org/10.1038/ncomms5114
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}

permTest <- function(dataDir, num) {
  setwd(paste0(homeDir, "/data/", dataDir))
  permData <- fread("MutsAndDels_over5pct_UVfiltered.csv")
  row.names(permData) <- permData$Tumor_Sample_Barcode
  # Create matrix, and preserve original ratio of altered samples between UV high and UV low
  mat <- permData[ , -1]
  high <- mat[mat$UV_sig == "High", ]
  high <- high[, -"UV_sig"]
  highCols <- colSums(high)
  low <- mat[mat$UV_sig == "Low", ]
  low <- low[, -"UV_sig"]
  lowCols <- colSums(low)
  dif <- highCols - lowCols
  mat <- mat[, -"UV_sig"]
  h <- nrow(high)
  l <- nrow(low)
  
  x <- 1
  repeat {
    curved <- curve_ball(mat)
    curvedHigh <- curved[1:h, ]
    curvedHighCols <- colSums(curvedHigh)
    curvedLow <- curved[(h+1):(h+l), ]
    curvedLowCols <- colSums(curvedLow)
    curvedDif <- curvedHighCols - curvedLowCols
    
    if (x == 1) {
      results <- as.data.frame(curved[1,])
      results[,1] <- 0
      rownames(results) <- colnames(mat)
      colnames(results) <- "p"
    }
    # If randomized matrix has higher High-Low ratio for a gene, p value for that gene increases
    for (i in 1:length(dif)) {
      if (curvedDif[i] >= dif[i]) { results[i, ] <- results[i, ] + 1 }
    }

    # Stop permutation test when specified number of permutations are reached, and finalize p value statistic
    if (x == num){
      pVals <- results
      pVals$Gene <- rownames(pVals)
      pVals$p <- pVals$p / x
      write.csv(pVals, paste0(x,"_Permutations_Raw.csv"), row.names = F)
      break
    }
    x = x+1
    print(x)
    if (x%%100 == 0) { print(paste0(x, " permutations completed, ", num-x, " permutations left")) }
  }
}

plotResults <- function(dataDir) {
  setwd(paste0(homeDir, "/data/", dataDir))
  file.ls <- list.files(path=getwd(),pattern="Permutations_Raw.csv")
  results <- read.csv(file.ls)

  qobj <- qvalue(p = results$p)
  results$q <- qobj$qvalues
  results$qlog <- -log10(results$q + 0.00000001)
  results <- results[order(results$q),]

  # Highlight repair pathway specific genes, as well as melanoma control genes
  source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
  results$pathway <- ""
  results$pathway[results$Gene %in% NERgenes] <- "NER"
  results$pathway[results$Gene %in% melControl] <- "CTR"
  results$pathway[results$Gene %in% HRgenes] <- "HR"
  results$pathway[results$Gene %in% MMRgenes] <- "MMR"
  results$pathway[results$Gene %in% BERgenes] <- "BER"

  write.csv(results, "Permutation_Results_qValues.csv")
  
  results$rank <- seq.int(nrow(results))
  results <- results[order(-results$q),]
  results$rank[results$q == 0] <- max(results$rank[results$q == 0])
  results$rank <- factor(results$rank, levels = unique(results$rank))
  highlight <- results[results$pathway != "", ]

  resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
  if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
  setwd(resultsDir)

  # Plot -log transformed q values for each gene with previously described highlights
  permPlot <- ggplot(results, aes(x = rank, y = qlog)) +
    geom_point() +
    geom_point(data = highlight, aes(x = rank, y = qlog, color = factor(pathway)),
               size = 3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text=element_text(size = 8)) +
    geom_label_repel(aes(label=ifelse(pathway != "",as.character(Gene),'')),
                     box.padding = 0.5,
                     point.padding = 0.5,
                     max.overlaps = Inf,
                     segment.color = 'grey50')
  permPlot
  ggsave("Permutation_Test_Plot.pdf", width = 10, height = 5, units = "in")
}

### Data already prepared ###----------------------
dataPrep("skcm_tcga")
permTest("skcm_tcga", 100)
plotResults("skcm_tcga")
