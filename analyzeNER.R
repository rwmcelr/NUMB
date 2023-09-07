## Library calls, load gene lists and cosmic signatures ##----------------------
library(ggplot2)
library(colorspace)
library(data.table)
library(maftools)
library(reshape2)
library(dplyr)
library(survival)
library(survminer)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"

## Function area  --------------------------------------------------------
generateOncoplot <- function(mutations, clin, cnv) {
  # Leverage maftools package to generate oncoplot, maf object can take some time to load
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

analyzeNER <- function(clin, combinedSig) {
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

  # Plot % of samples with NER pathway alteration by UV level (excluding moderate)
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

  nerSig <- combinedSig[ , colnames(combinedSig) %in% c("rowname", "CETN2", "GTF2H2", "ERCC8", "CDK7", "CCNH", "ERCC6", "RAD23B")]
  nerSig$Tumor_Sample_Barcode <- rownames(nerSig)
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

  # Plot % of samples with alteration for each specific gene in NER signature by UV level (excluding moderate)
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

nerSurvival <- function(clin) {
  # Create multi-condition strata for each combination of UV level and NER proficiency status (excluding UV moderate)
  surv <- clin
  surv <- surv[surv$UV_sig != "Moderate", ]
  surv$Combination[surv$UV_sig == "High" & surv$NER_sig == 0] <- "UV High, NER Proficient"
  surv$Combination[surv$UV_sig == "High" & surv$NER_sig > 0] <- "UV High, NER Deficient"
  surv$Combination[surv$UV_sig == "Low" & surv$NER_sig == 0] <- "UV Low, NER Proficient"
  surv$Combination[surv$UV_sig == "Low" & surv$NER_sig > 0] <- "UV Low, NER Deficient"
  surv$OS_MONTHS <- as.numeric(surv$OS_MONTHS)
  surv$OS_STATUS <- as.numeric(substr(surv$OS_STATUS, 1, 1))

  # Plot survival for each possible combination status
  survival <- ggsurvplot(fit = surv_fit(Surv(OS_MONTHS, OS_STATUS) ~ Combination, data = surv),
             xlab = "Months",
             ylab = "Overall survival probability",
             pval = T,
             legend = "bottom")
  ggpar(survival, font.legend = list(size = 7))
  ggsave("Survival_NER.pdf", width = 8.62, height = 5.47, units = "in")
}

# Load data and execute plotting functions
nerAnalysis <- function(dataDir, subDir) {
  setwd(paste0(homeDir, "/data/", dataDir))
  
  clin <- read.csv("clinical.csv")
  clin <- clin[order(-clin$UV_sig_value),]
  
  if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
  } else { mutations <- fread("data_mutations.txt") }
  mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B" 
  
  combinedSig <- read.csv("NER_mutCNVsig.csv", row.names = 1)
  
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
  
  if (subDir != "") {
    resultsDir <- paste0(homeDir,"/results/",dataDir,"/",subDir,"/")
    if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
    setwd(resultsDir)
  } else {
    resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
    if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
    setwd(resultsDir)
  }
  
  generateOncoplot(mutations, clin, cnv2)
  analyzeNER(clin, combinedSig)
  nerSurvival(clin)
}
## Execute code -------------------
nerAnalysis("skcm_tcga")
