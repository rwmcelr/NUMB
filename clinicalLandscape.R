library(ggplot2)
library(colorspace)
library(reshape2)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") # Load repair genes as list objects
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files" # Create root data directory pointer
## SIGNATURE/TITV GRAPHS --------------------------------------------------------
clinScape <- function(dataDir, resultsDir) {
  # Load data, then change to results directory to avoid cluttering data directory
  dataDir <- "skcm_tcga"
  setwd(paste0(homeDir, "/data/", dataDir))
  clin <- read.csv("clinical.csv")
  sigs <- read.csv("mutationalSignatures.csv", row.names = "SAMPLE_ID")
  tnm <- read.csv("tnm.csv")
  # setwd(resultsDir)
  
  # Group all repair signature defects, Indirect UV, and Chemotherapy as their own levels, and set all other signatures to "other"
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
  ggsave("tcga_skcm_COSMICsigs.pdf", width = 6.5, height = 2, units = "in")
  
  # Condense trinucelotide matrix, and find C>T at dipyrimidine occurrence rate
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
  
  # Sort tnm data by decreasing UV, then prep for plotting
  UVvals <- clin[, c("Tumor_Sample_Barcode", "UV_sig_value")]
  temp <- merge(tnm2, UVvals, by = "Tumor_Sample_Barcode")
  temp <- temp[order(-temp$UV_sig_value),]
  temp2 <- temp[,1:8]
  temp2$Tumor_Sample_Barcode <- as.character(temp2$Tumor_Sample_Barcode)
  temp2$Tumor_Sample_Barcode <- factor(temp2$Tumor_Sample_Barcode, levels = unique(temp2$Tumor_Sample_Barcode))
  tnmPlotDat <- melt(temp2)
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
  ggsave("tcga_skcm_titv.pdf", width = 6.5, height = 2, units = "in")
}