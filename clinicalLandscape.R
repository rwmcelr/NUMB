library(ggplot2)
library(colorspace)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") # Load repair genes as list objects
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files" # Create root data directory pointer
## SIGNATURE/TITV GRAPHS --------------------------------------------------------
processData <- function(dataDir, resultsDir) {
  dataDir <- "skcm_tcga"
  setwd(paste0(homeDir, "/data/", dataDir))
  clin <- read.csv("clinical.csv") # Read in generated clinical file
  sigs <- read.csv("mutationalSignatures.csv", row.names = "SAMPLE_ID") # Read in imputed signature file

  # Create a plot showing mutation signature by decreasing UV-sig
  sigs$Tumor_Sample_Barcode <- as.character(sigs$Tumor_Sample_Barcode)
  sigs$Tumor_Sample_Barcode <- factor(sigs$Tumor_Sample_Barcode, levels = unique(sigs$Tumor_Sample_Barcode))
  
  # Transforming the expo dataframe to create a stacked barplot with ggplot2
  dat <- melt(sigs, id.vars = "Tumor_Sample_Barcode")
  dat2 <- transform(dat, Sample = factor(Tumor_Sample_Barcode, levels = unique(sigs$Tumor_Sample_Barcode)))
  
  # Creating a color palette, making Signature 7 distinct from the rest
  c2 <- diverging_hcl(12, palette="Berlin")
  c2[1] <- "#808080"
  c2[7] <- "#FFFF00"
  c2[2] <- "#358856" #HR Defect
  c2[5] <- "#06FCED" #BER Defect
  c2[6] <- "#FF8B3D" #Indirect UV
  c2[9] <- "#06FCB4" #MMR Defect
  c2[11] <- "#8E06FC" #Chemotherapy
  
  # Create the plot, using signatures to dictate fill of the bars
  cosmicPlot <- ggplot(dat2, aes(x=Sample, y=value, fill=variable)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    scale_fill_manual(values = c2) +
    theme(axis.text.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle("TCGA_SKCM COSMIC Signatures") +
    labs(fill="COSMIC Signature") +
    #theme(legend.position = "bottom") +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill=guide_legend(nrow=4,byrow=TRUE))
  print(cosmicPlot)
  ggsave("tcga_skcm_COSMICsigs.pdf", width = 6.5, height = 2, units = "in")
  
  # Generate trinucleotide matrix plot
  tnm <- trinucleotideMatrix(maf = maf, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19', prefix = "chr")
  tnm <- as.data.frame(tnm[[1]])
  tnm$Tumor_Sample_Barcode <- row.names(tnm)
  tnm2 <- as.data.frame(tnm$Tumor_Sample_Barcode)
  tnm2$C.A <- rowSums(tnm[ , c(1:16)])
  tnm2$C.G <- rowSums(tnm[ , c(17:32)])
  tnm2$C.T <- rowSums(tnm[ , c(33:48)])
  tnm2$T.A <- rowSums(tnm[ , c(49:64)])
  tnm2$T.C <- rowSums(tnm[ , c(65:80)])
  tnm2$T.G <- rowSums(tnm[ , c(81:96)])
  tnm2$C.Tdipyr <- rowSums(tnm[ , c(34,36:40,42,44:48)])
  tnm2$C.Tother <- rowSums(tnm[ , c(33,35,41,43)])
  colnames(tnm2) <- c("Tumor_Sample_Barcode", "C>A", "C>G", "C>T", "T>A", 
                      "T>C", "T>G", "C>T at Dipyrimidine Site", "C>T Other")
  
  UVvals <- clin[, c("Tumor_Sample_Barcode", "UV_sig_value")]
  tnm2 <- tnm2[which(tnm2$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode),]
  temp <- merge(tnm2, UVvals, by = "Tumor_Sample_Barcode")
  temp <- temp[order(-temp$UV_sig_value),]
  tnm <- temp[,1:9]
  tnm$Tumor_Sample_Barcode <- as.character(tnm$Tumor_Sample_Barcode)
  tnm$Tumor_Sample_Barcode <- factor(tnm$Tumor_Sample_Barcode, levels = unique(tnm$Tumor_Sample_Barcode))
  dat <- melt(tnm)
  dat <- dat[which(dat$variable != "C>T"),]
  c1 <- diverging_hcl(8, palette="Purple-Green")
  
  plot <- ggplot(dat, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
    geom_bar(position="fill", stat="identity", width = 1) +
    scale_fill_manual(values = c1) +
    theme(axis.text.x=element_blank(),
          axis.title.y=element_blank()) +
    labs(fill="Nucleotide Change") +
    #theme(legend.position = "bottom") +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  print(plot)
  ggsave("tcga_skcm_titv.pdf", width = 6.5, height = 2, units = "in")
}