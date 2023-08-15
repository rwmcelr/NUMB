library(ggplot2)
library(colorspace)
library(reshape2)
library(dplyr)
library(data.table)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"
## Function area  --------------------------------------------------------
plotSigs <- function(sigs, clin) {
  sig <- sigs[,-c(7:10, 69)] # Remove UV signature for now, will be re-added later
  hr_defect <- c("SBS3")
  mmr_defect <- c("SBS6", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
  ber_defect <- c("SBS30", "SBS36")
  indirect_uv <- c("SBS38")
  chemotherapy <- c("SBS25", "SBS31", "SBS35", "SBS86", "SBS87")
  all_named <- c(hr_defect, mmr_defect, ber_defect, indirect_uv, chemotherapy)
  
  sig$Other <- rowSums(sig[colnames(sig) %notin% all_named])
  sig$HR_Defect <- rowSums(sig[colnames(sig) %in% hr_defect])
  sig$MMR_Defect <- rowSums(sig[colnames(sig) %in% mmr_defect])
  sig$BER_Defect <- rowSums(sig[colnames(sig) %in% ber_defect])
  sig$Indirect_UV <- rowSums(sig[colnames(sig) %in% indirect_uv])
  sig$Chemotherapy <- rowSums(sig[colnames(sig) %in% chemotherapy])
  sig$UV_Signature <- sigs$UV_sig
  
  # Set tumor sample barcode to factor to maintain UV signature decreasing order
  sig$Tumor_Sample_Barcode <- as.character(rownames(sig))
  sig$Tumor_Sample_Barcode <- factor(sig$Tumor_Sample_Barcode, levels = unique(sig$Tumor_Sample_Barcode))
  
  # Make plot data frame and color palette
  temp <- sig[,c("Tumor_Sample_Barcode","UV_Signature","HR_Defect","MMR_Defect","BER_Defect",
                 "Indirect_UV","Chemotherapy","Other")]
  sigPlotDat <- melt(temp, id.vars = "Tumor_Sample_Barcode")
  pal <- diverging_hcl(7, palette=" Berlin")
  pal[1] <- "#FFFF00"
  pal[7] <- "#808080"
  
  # Create the plot, using signatures to dictate fill of the bars
  cosmicPlot <- ggplot(sigPlotDat, aes(x=Tumor_Sample_Barcode, y=value, fill=variable)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    scale_fill_manual(values = pal) +
    theme(axis.text.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle("TCGA_SKCM COSMIC Signatures") +
    labs(fill="COSMIC Signature") +
    theme(legend.position = "bottom") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill=guide_legend(nrow=4,byrow=TRUE))
  print(cosmicPlot)
  ggsave(paste0(dataDir,"_COSMICsigs.pdf"), width = 8, height = 8, units = "in")
  
  clin$UV_sig <- factor(clin$UV_sig, levels = unique(clin$UV_sig))
  UVbyLplot <- ggplot(clin, aes(x=UV_sig, y=UV_sig_value, fill=UV_sig)) +
    geom_boxplot() + 
    labs(x="UV Signature Level", y="UV Signature (%)") + 
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(UVbyLplot)
  ggsave(paste0(dataDir,"_AvgUVsigByLevel.pdf"))
}
  
plotTNM <- function(clin, tnm) {
  tnm2 <- as.data.frame(tnm$X)
  tnm2$C.A <- rowSums(tnm[ , c(2:17)])
  tnm2$C.G <- rowSums(tnm[ , c(18:33)])
  tnm2$T.A <- rowSums(tnm[ , c(50:65)])
  tnm2$T.C <- rowSums(tnm[ , c(66:81)])
  tnm2$T.G <- rowSums(tnm[ , c(82:97)])
  tnm2$C.Tdipyr <- rowSums(tnm[ , c(35,37:41,43,45:49)])
  tnm2$C.Tother <- rowSums(tnm[ , c(34,36,42,44)])
  colnames(tnm2) <- c("Tumor_Sample_Barcode", "C>A", "C>G", "T>A", 
                      "T>C", "T>G", "C>T at Dipyrimidine Site", "C>T Other")
  
  UVvals <- clin[, c("Tumor_Sample_Barcode", "UV_sig_value", "UV_sig")]
  temp <- merge(tnm2, UVvals, by = "Tumor_Sample_Barcode")
  temp <- temp[order(-temp$UV_sig_value),]
  temp2 <- temp[,1:8]
  temp2$Tumor_Sample_Barcode <- as.character(temp2$Tumor_Sample_Barcode)
  temp2$Tumor_Sample_Barcode <- factor(temp2$Tumor_Sample_Barcode, levels = unique(temp2$Tumor_Sample_Barcode))
  tnmPlotDat <- reshape2::melt(temp2)
  pal2 <- diverging_hcl(8, palette="Purple-Green")
  
  tnmPlot <- ggplot(tnmPlotDat, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
    geom_bar(position="fill", stat="identity", width = 1) +
    scale_fill_manual(values = pal2) +
    theme(axis.text.x=element_blank(),
          axis.title.y=element_blank()) +
    labs(fill="Nucleotide Change") +
    theme(legend.position = "bottom") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  print(tnmPlot)
  ggsave(paste0(dataDir,"_TNM.pdf"), width = 8, height = 8, units = "in")
  
  CTpct <- temp %>% mutate(CT.Dipyr.pct = .[,7]/rowSums(select(., -c(1,9,10))))
  CTpct$UV_sig <- factor(CTpct$UV_sig, levels=unique(CTpct$UV_sig))
  CTbreakdownPlot <- ggplot(CTpct, aes(fill=UV_sig, y=CT.Dipyr.pct, x=UV_sig)) + 
    geom_boxplot() + 
    labs(x="UV Signature Level", y="C>T at dipyrimidine site (%)") + 
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(CTbreakdownPlot)
  ggsave(paste0(dataDir,"_CTdipyrPctByLevel.pdf"))
}

plotSNVs <- function(clin, mutations) {
  # Average SNV plot
  muts <- mutations[mutations$Variant_Classification != "Silent",]
  temp <- as.data.frame(table(muts$Tumor_Sample_Barcode))
  colnames(temp)[1] <- "Tumor_Sample_Barcode"
  tempClin <- merge(clin, temp, by = "Tumor_Sample_Barcode")
  tempClin <- tempClin[order(-tempClin$UV_sig_value),]
 
  tempClin$UV_sig <- factor(tempClin$UV_sig, levels = unique(tempClin$UV_sig))
  AvgSNVplot <- ggplot(tempClin, aes(x = UV_sig, y = Freq, fill = UV_sig)) +
    geom_boxplot() +
    scale_y_continuous(limits = quantile(tempClin$Freq, c(0.1, 0.9))) +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1)) +
    labs(x = "UV Signature Level", y = "Average # of SNVs") +
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(AvgSNVplot)
  ggsave(paste0(dataDir,"_AverageSNVbyLevel.pdf"))
  
  # Individual SNV plot
  tempClin$Tumor_Sample_Barcode <- factor(tempClin$Tumor_Sample_Barcode, levels = unique(tempClin$Tumor_Sample_Barcode))
  
  library(ggbreak)
  sampleSNVplot <- ggplot(tempClin, aes(fill = UV_sig, y=Freq, x=Tumor_Sample_Barcode)) +
    geom_bar(stat="identity", width = 1, position = "dodge") +
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7", "#808080")) +
    theme(axis.text.x=element_blank()) +
    labs(fill="UV Signature Level", x = "Sample", y = "Total SNV") +
    theme(legend.position = "none") +
    scale_y_break(c(3000, 15000 ) ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  print(sampleSNVplot)
  ggsave(paste0(dataDir,"_AverageSNVperSample.pdf"), width = 7, height = 2, units = "in")
  
  tempClin$index <- seq.int(nrow(tempClin)) 
  tempClin$TMBlog10 <- log10(tempClin$TMB_NONSYNONYMOUS)
  tempClin$TMBlog10[tempClin$TMBlog10 < 0] <- 0
  
  mean <- tempClin %>% 
    group_by(UV_sig) %>% 
    summarise(mean_val = mean(TMBlog10))
  
  TMBplot <- ggplot(tempClin, aes(x=-UV_sig_value, y = TMBlog10, color = UV_sig)) +
    geom_point() +
    scale_color_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0.05)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(x = "UV Signature (Decreasing)") +
    geom_hline(data= mean, aes(yintercept = mean_val,col=UV_sig)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Mean lines in this graph need to be adjusted in illustrator to only cover their respective UV level
  TMBplot
  ggsave("TMB.pdf", width = 12, height = 2, units = "in")
}

## Execute code -------------------
dataDir <- "skcm_tcga"
setwd(paste0(homeDir, "/data/", dataDir))

clin <- read.csv("clinical.csv")
sigs <- read.csv("mutationalSignatures.csv", row.names = "SAMPLE_ID")
tnm <- read.csv("tnm.csv")
if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
} else { mutations <- fread("data_mutations.txt")}
mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B" 

resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
setwd(resultsDir)

plotSigs(sigs)
plotTNM(clin, tnm)
