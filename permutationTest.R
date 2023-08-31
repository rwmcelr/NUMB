### Library calls/ scripts -------------------------------------
library(maftools)
library(ggplot2)
library(colorspace)
library(data.table)
library(reshape2)
library(dplyr)
source("C:/Users/robmc/Desktop/NER/R/oncoPlease.R")
sourceDir <- "C:/Users/robmc/Desktop/NER"
skcmDir <- "C:/Users/robmc/Desktop/NER/data/skcm_perm/"
dfci19Dir <- "C:/Users/robmc/Desktop/NER/data/mel_dfci_2019/"
uclaDir <- "C:/Users/robmc/Desktop/NER/data/mel_ucla_2016/"
'%notin%' <- Negate('%in%')

### Data Preparation -------------------------------------
dataPrep <- function(dataDir) {
  setwd(dataDir)
  if(file.exists("scores.gistic")) {
    gisticLesions <- "all_lesions.conf_99.txt"
    gisticScores <- "scores.gistic"
    gisticAmp <- "amp_genes.conf_99.txt"
    gisticDel <- "del_genes.conf_99.txt"
    
    cnv <- readGistic(
      gisticAllLesionsFile = gisticLesions,
      #gisticAmpGenesFile = gisticAmp,
      gisticDelGenesFile = gisticDel,
      gisticScoresFile = gisticScores,
      cnLevel = "all",
      isTCGA = T,
      verbose = TRUE
    )
    summary <- cnv@data
    write.table(summary, "cnvSum.txt", row.names = F)
    cnvLong <- fread("cnvSum.txt")
    cnvWide <- dcast(data = cnvLong,formula = Tumor_Sample_Barcode~Hugo_Symbol, fun.aggregate = NULL,
                     value.var = "Variant_Classification")
  } else {
    cnv <- fread("data_cna.txt")
    cnv[cnv == 1 | cnv == 2] <- 0
    cnv[cnv == -1 | cnv == -2] <- 1
    cnv <- melt(cnv)
    cnvWide <- dcast(data = cnv,formula = variable~Hugo_Symbol, fun.aggregate = mean, value.var = "value")
    colnames(cnvWide)[1] <- "Tumor_Sample_Barcode"
  }
  
  if(file.exists("data_mutations.txt")) { mutations <- fread("data_mutations.txt") } else {
    mutations <- fread("data_mutations_extended.txt")
  }
  
  mutationsFiltered <- mutations[mutations$Variant_Classification != "Silent", ]
  mut <- distinct(mutationsFiltered, End_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, .keep_all = T)
  mutationsWide <- dcast(data = mut,formula = Tumor_Sample_Barcode~Hugo_Symbol, fun.aggregate = length,
                         value.var = "Variant_Classification")
  
  # Keep only mutation and CNV data for which UV-signature and clinical data exist
  clin <- fread("updatedClinical.csv")
  cnvFinal <- cnvWide[cnvWide$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
  mutationsWide$Tumor_Sample_Barcode <- substr(mutationsWide$Tumor_Sample_Barcode, 1, 12)
  mutFinal <- mutationsWide[mutationsWide$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
  
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
  everythingSig <- everything[ ,colSums(everything) > round]
  everything$Tumor_Sample_Barcode <- rownames(everything)
  everythingSig$Tumor_Sample_Barcode <- rownames(everythingSig)
  UV <- clin[ , c("Tumor_Sample_Barcode","UV_sig")]
  UV <- UV[!duplicated(UV$Tumor_Sample_Barcode), ]
  sigUV <- merge(everythingSig, UV, by = "Tumor_Sample_Barcode")
  sigUVhigh <- sigUV[sigUV$UV_sig == "High", ]
  sigUVlow <- sigUV[sigUV$UV_sig == "Low", ]
  write.csv(sigUV, "MutsAndDels_over5percent.csv", row.names = F)
  write.csv(sigUVhigh, "UVhigh_MutsAndDels_over5percent.csvgh_MutsAndDels_over5percent.csv", row.names = F)
  write.csv(sigUVlow, "UVlow_MutsAndDels_over5percent.csv", row.names = F)
  
  highlowsig <- rbind(sigUVhigh, sigUVlow)
  write.csv(highlowsig, "highLow_MutsAndDels_over5percent.csv", row.names = F)
  # highlow <- rbind(UVhigh, UVlow)
  # write.csv(highlow, "highLow_MutsAndDels_all.csv", row.names = F)
  setwd(sourceDir)
}

dataPrep(skcmDir)
dataPrep(dfci19Dir)
dataPrep(uclaDir)

### Data already prepared ###----------------------
curve_ball<-function(m){
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

perm <- function(dataDir, num) {
  setwd(dataDir)
  highlowsig <- fread("highLow_MutsAndDels_over5percent.csv")
  row.names(highlowsig) <- highlowsig$Tumor_Sample_Barcode
  m <- highlowsig[ , -1]
  high <- m[m$UV_sig == "High", ]
  high <- high[, -"UV_sig"]
  highCols <- colSums(high)
  low <- m[m$UV_sig == "Low", ]
  low <- low[, -"UV_sig"]
  lowCols <- colSums(low)
  originalDif <- highCols - lowCols
  m <- m[, -"UV_sig"]
  h <- nrow(high)
  l <- nrow(low)
  
  x <- 1
  repeat {
    curved <- curve_ball(m)
    curvedHigh <- curved[1:h, ]
    curvedHighCols <- colSums(curvedHigh)
    curvedLow <- curved[(h+1):(h+l), ]
    curvedLowCols <- colSums(curvedLow)
    curvedDif <- curvedHighCols - curvedLowCols
    
    if (x == 1) { 
      results <- as.data.frame(curved[1,])
      results[,1] <- 0 
      rownames(results) <- colnames(m)
      colnames(results) <- "p"
    }
    for (i in 1:length(originalDif)) {
      if (curvedDif[i] >= originalDif[i]) { results[i, ] <- results[i, ] + 1 }
    }
    
    if (x == num){
      pVals <- results
      pVals$Gene <- rownames(pVals)
      write.csv(pVals, paste0("pVals_",x,"_permutations.csv"))
      pVals$p <- pVals$p / x
      break
    }
    x = x+1
    print(x)
  }
}

perm(dfci19Dir, 1259)

setwd(dfci19Dir)
p1 <- read.csv("pVals_1240_permutations.csv", row.names = 1)
p2 <- read.csv("pVals_1241_permutations.csv", row.names = 1)
p3 <- read.csv("pVals_1252_permutations.csv", row.names = 1)
p4 <- read.csv("pVals_1253_permutations.csv", row.names = 1)
p5 <- read.csv("pVals_1257_permutations.csv", row.names = 1)
p6 <- read.csv("pVals_1258_permutations.csv", row.names = 1)
p7 <- read.csv("pVals_1259_permutations.csv", row.names = 1)
p8 <- read.csv("pVals_1260_permutations.csv", row.names = 1)

final <- p1
final$p <- 0
final$p <- p1$p + p2$p + p3$p + p4$p + p5$p + p6$p + p7$p + p8$p
final$p <- final$p / 10020

library(qvalue)
qobj <- qvalue(p = final$p)
final$q <- qobj$qvalues
final$qlog <- -log10(final$q + 0.00000001)
final <- final[order(final$q),]

source("C:/Users/robmc/Desktop/NER/R/oncoPlease.R")
final$plot <- 0
final$plot[final$Gene %in% NER] <- 5
final$plot[final$Gene %in% controlGenes] <- 4
final$plot[final$Gene %in% HR] <- 3
final$plot[final$Gene %in% MMR] <- 2
final$plot[final$Gene %in% BER] <- 1

final$pathway[final$plot == 5] <- "NER"
final$pathway[final$plot == 4] <- "CTR"
final$pathway[final$plot == 3] <- "HR"
final$pathway[final$plot == 2] <- "MMR"
final$pathway[final$plot == 1] <- "BER"

write.csv(final, "final_withQ.csv")

figure <- read.csv("final_withQ.csv")
figure <- figure[order(figure$q),]
figure$rank <- seq.int(nrow(figure))
figure <- figure[order(-figure$q),]
figure$rank[figure$q == 0] <- 849
figure$rank <- as.numeric(figure$rank)
figure$rank <- factor(figure$rank, levels = unique(figure$rank))
highlight <- figure[figure$plot != 0, ]

library(ggplot2)
library(ggrepel)
p <- ggplot(figure, aes(x = rank, y = qlog)) +
  geom_point() +
  geom_point(data = highlight, aes(x = rank, y = qlog, color = factor(pathway)),
             size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size = 8)) +
  geom_label_repel(aes(label=ifelse(plot != 0,as.character(Gene),'')),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   segment.color = 'grey50')
p
ggsave("figure.pdf", width = 10, height = 5, units = "in")



setwd(uclaDir)
p1 <- read.csv("pVals_10000_permutations.csv", row.names = 1)

final <- p1
final$p <- final$p / 10000

library(qvalue)
qobj <- qvalue(p = final$p)
final$q <- qobj$qvalues
final$qlog <- -log10(final$q + 0.00000001)
final <- final[order(final$q),]

source("C:/Users/robmc/Desktop/NER/R/oncoPlease.R")
final$plot <- 0
final$plot[final$Gene %in% NER] <- 5
final$plot[final$Gene %in% controlGenes] <- 4
final$plot[final$Gene %in% HR] <- 3
final$plot[final$Gene %in% MMR] <- 2
final$plot[final$Gene %in% BER] <- 1

final$pathway[final$plot == 5] <- "NER"
final$pathway[final$plot == 4] <- "CTR"
final$pathway[final$plot == 3] <- "HR"
final$pathway[final$plot == 2] <- "MMR"
final$pathway[final$plot == 1] <- "BER"

write.csv(final, "final_withQ.csv")

figure <- read.csv("final_withQ.csv")
figure <- figure[order(figure$q),]
figure$rank <- seq.int(nrow(figure))
figure <- figure[order(-figure$q),]
figure$rank[figure$q == 0] <- 849
figure$rank <- as.numeric(figure$rank)
figure$rank <- factor(figure$rank, levels = unique(figure$rank))
highlight <- figure[figure$plot != 0, ]

library(ggplot2)
library(ggrepel)
p <- ggplot(figure, aes(x = rank, y = qlog)) +
  geom_point() +
  geom_point(data = highlight, aes(x = rank, y = qlog, color = factor(pathway)),
             size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size = 8)) +
  geom_label_repel(aes(label=ifelse(plot != 0,as.character(Gene),'')),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   segment.color = 'grey50')
p
ggsave("figure.pdf", width = 10, height = 5, units = "in")
