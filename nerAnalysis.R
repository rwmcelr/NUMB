## Library calls, load gene lists and cosmic signatures ##----------------------
library(ggplot2)
library(colorspace)
library(data.table)
library(maftools)
library(reshape2)
library(dplyr)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"

## Function area  --------------------------------------------------------
generateOncoplot <- function(mutations, clin, cnv) {
  maf <- read.maf(maf = mutations,
                  clinicalData = clin,
                  cnTable = cnv,
                  isTCGA = F)
  
  
  uvcolors=c("#FF0000", "#FFA500", "#4F53B7")
  names(uvcolors)=c("High","Moderate","Low")
  anno_cols = list(UV_sig = uvcolors)
  feat <- names(anno_cols)
  
  pdf("NER_Oncoplot.pdf")
  oncoplot(maf = maf,
           genes = NERgenes, cohortSize = length(clin$Tumor_Sample_Barcode[]),
           clinicalFeatures = feat, annotationColor = anno_cols,
           removeNonMutated = F, draw_titv = F, drawRowBar = F,
           sampleOrder = clin$Tumor_Sample_Barcode,
           fontSize = 0.6, annotationFontSize = 0.5, legendFontSize = 0.8, anno_height = 1.8)
  dev.off()
}

analyzeNER <- function(clin, mutSig, cnvSig) {
  tempClin <- clin[, c("Tumor_Sample_Barcode", "UV_sig", "UV_sig_value", "NER_sig")]
  tempClin$NER_sig[tempClin$NER_sig > 1] <- 1
  pathDat <- as.data.frame(prop.table(table(tempClin$UV_sig, tempClin$NER_sig)))
  pathDat <- pathDat[pathDat$Var1 != "Moderate" & pathDat$Var2 == 1, ]
  
  # Code for normalizing % NER pathway altered for each level rather than whole cohort
  # tempClin$NER_sig <- tempClin$NER_sig / 7
  # highSum <- sum(tempClin$NER_sig[tempClin$UV_sig == "High"])
  # lowSum <- sum(tempClin$NER_sig[tempClin$UV_sig == "Low"])
  # n <- as.data.frame(table(tempClin$UV_sig))
  # pathDat <- data.frame(UV_Level = c("High", "Low"), Altered = c(highSum/n[1,2], lowSum/n[2,2]))
  
  nerPathPlot <- ggplot(pathDat, aes(x=Var1, y=Freq, fill = Var1)) +
    geom_col() +
    scale_y_continuous(limits = c(0, 0.6), expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    labs(y = "% Altered (NER Pathway)", x = "UV Signature Level", fill = "UV Signature Level") +
    scale_fill_manual(values = c("#FF0000", "#4F53B7")) +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  print(nerPathPlot)
  ggsave("NER_pcts.pdf", width = 4, height = 4, units = "in")
  
  m <- data.frame(mutSig, row.names = 1)
  m <- m*0.1
  c <- data.frame(cnvSig, row.names = 1)
  
  both <- bind_rows(c %>% tibble::rownames_to_column(), 
                    m %>% tibble::rownames_to_column()) %>% 
    group_by(rowname) %>% 
    summarise_all(sum)
  
  both[both == 0 | both == 1 | both == 2] <- 0
  both[both == -0.9 | both == -1.9] <- 1
  both[both == 0.1] <- 1
  both[both == -1 | both == -2] <- 1
  both[both == 1.1 | both == 2.1] <- 1
  
  nerSig <- both[ , colnames(both) %in% c("rowname", "CETN2", "GTF2H2", "ERCC8", "CDK7", "CCNH", "ERCC6", "RAD23B")]
  colnames(nerSig)[1] <- "Tumor_Sample_Barcode"
  geneClin <- merge(tempClin, nerSig, by = "Tumor_Sample_Barcode")
  geneClin <- geneClin[ , -c(1,3,4)]
  geneDat <- data.frame()
  for (i in 2:ncol(geneClin)) {
    temp <- geneClin[, c(1,i)]
    toAdd <- data.frame(UV_Level = c("High", "Low"), Gene = rep(colnames(temp)[2], 2),
                        Altered = c(sum(temp[,2][temp$UV_sig == "High"])/length(temp$UV_sig),
                                    sum(temp[,2][temp$UV_sig == "Low"])/length(temp$UV_sig)))
    # Code for examining % of samples altered by UV level, rather than as a % of total samples in cohort
    # toAdd <- data.frame(UV_Level = c("High", "Low"), Gene = rep(colnames(temp)[2], 2),
    #                     Altered = c(sum(temp[,2][temp$UV_sig == "High"])/length(temp$UV_sig[temp$UV_sig == "High"]),
    #                                 sum(temp[,2][temp$UV_sig == "Low"])/length(temp$UV_sig[temp$UV_sig == "Low"])))
    geneDat <- rbind(geneDat, toAdd)
  }
  
  nerGenePlot <- ggplot(geneDat, aes(x=reorder(Gene,-Altered), y=Altered, fill = UV_Level)) +
    geom_col(position = "dodge") +
    scale_y_continuous(limits = c(0, 0.6), expand = c(0,0)) +
    labs(y = "% of Altered Samples", x = NULL, fill = "UV Signature Level") +
    scale_fill_manual(values = c("#FF0000", "#4F53B7")) +
    scale_x_discrete(expand = c(0.005, 0)) +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  print(nerGenePlot)
  ggsave("NER_gene_pcts.pdf", width = 8, height = 4, units = "in")
}

nerAnalysis <- function(dataDir) {
  setwd(paste0(homeDir, "/data/", dataDir))
  
  clin <- read.csv("clinical.csv")
  clin <- clin[order(-clin$UV_sig_value),]
  
  if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
  } else { mutations <- fread("data_mutations.txt") }
  mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B" 
  
  cnvSig <- read.csv("NER_cnvSig.csv")
  mutSig <- read.csv("NER_mutSig.csv")
  
  cnv <- fread("data_cna.txt")
  cnv <- unique(cnv[!duplicated(cnv$Hugo_Symbol), ])
  cnv <- cnv[, -(2:3)]
  cnv2 <- reshape2::melt(cnv, id.vars = "Hugo_Symbol")
  colnames(cnv2) <- c("Gene", "Sample_name", "CN")
  cnv2$Gene[cnv2$Gene == "CUL4A" | cnv2$Gene == "CUL4B"] <- "CUL4A/B"
  cnv2$CN[cnv2$CN == -2 | cnv2$CN == -1] <- "Del"
  cnv2$CN[cnv2$CN == 2] <- "Amp"
  cnv2 <- cnv2[cnv2$CN != 0, ]
  cnv2 <- cnv2[cnv2$CN != 1, ]
  cnv2 <- cnv2[cnv2$Sample_name %in% clin$Tumor_Sample_Barcode, ]
  
  resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
  if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
  setwd(resultsDir)
  
  generateOncoplot(mutations, clin, cnv2)
  analyzeNER(clin, mutSig, cnvSig)
}
## Execute code -------------------
nerAnalysis("skcm_tcga")
