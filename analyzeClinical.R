## Library calls, load gene lists and cosmic signatures ##----------------------
library(ggplot2)
library(colorspace)
library(reshape2)
library(dplyr)
library(data.table)
library(survival)
library(survminer)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"

## Function area  --------------------------------------------------------
plotSigs <- function(sigs, clin) {
  # Name signatures to be highlited during plotted, anything not named here will be shown as "Other"
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
  
  # Plot signature breakdown per sample using signatures to dictate fill of the bars, sorted by decreasing UV signature level
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
  ggsave("COSMICsigs.pdf", width = 8, height = 8, units = "in")

  # Plot average UV signature as % of total signature by level
  clin$UV_sig <- factor(clin$UV_sig, levels = unique(clin$UV_sig))
  UVbyLplot <- ggplot(clin, aes(x=UV_sig, y=UV_sig_value, fill=UV_sig)) +
    geom_boxplot() + 
    labs(x="UV Signature Level", y="UV Signature (%)") + 
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(UVbyLplot)
  ggsave("Avg_Sig7_By_Level.pdf")
  
  print("Signature graphs generated successfully!")
}
  
plotTNM <- function(clin, tnm) {
  # Condense trinucleotide contexts into standard SNV notation
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

  # Plot SNV class breakdown per sample with further detailing of C>T mutations with respect to occurrence at dipyrimidine sites
  # sorted by decreasing UV signature level
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
  ggsave("TNM.pdf", width = 8, height = 8, units = "in")

  # Plot C>T transversions at dipyrimidine occurrence (in context of all mutation classes) by UV level
  CTpct <- temp %>% mutate(CT.Dipyr.pct = .[,7]/rowSums(select(., -c(1,9,10))))
  CTpct$UV_sig <- factor(CTpct$UV_sig, levels=unique(CTpct$UV_sig))
  CTbreakdownPlot <- ggplot(CTpct, aes(fill=UV_sig, y=CT.Dipyr.pct, x=UV_sig)) + 
    geom_boxplot() + 
    labs(x="UV Signature Level", y="C>T at dipyrimidine site (%)") + 
    scale_fill_manual(values = c("#FF0000", "#FFA500", "#4F53B7")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(CTbreakdownPlot)
  ggsave("Dipyr_Pct_By_Level.pdf")
  
  print("TNM graphs generated successfully!")
}

plotSNVs <- function(clin, mutations) {
  # Plot average SNV by UV level
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
  ggsave("Average_SNV_By_Level.pdf")
  
  # Plot SNV per sample, sorted by descending UV signature level
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
  ggsave("SNV_Per_Sample.pdf", width = 7, height = 2, units = "in")
  
  tempClin$index <- seq.int(nrow(tempClin)) 
  tempClin$TMBlog10 <- log10(tempClin$TMB_NONSYNONYMOUS)
  tempClin$TMBlog10[tempClin$TMBlog10 < 0] <- 0
  
  mean <- tempClin %>% 
    group_by(UV_sig) %>% 
    summarise(mean_val = mean(TMBlog10))

  # Plot TMB per sample, sorted by decreasing UV signature level
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
  
  print("SNV graphs generated successfully!")
}

plotFeatures <- function(clin) {
  # Group ages, groups are semi-arbitrary
  clin$Age_Group[clin$AGE < 30] <- "< 30"
  clin$Age_Group[clin$AGE >= 30 &
                   clin$AGE <= 50] <- "30-50"
  clin$Age_Group[clin$AGE > 50 &
                   clin$AGE <= 80] <- "51-80"
  clin$Age_Group[clin$AGE > 80] <- "80+"
  agecolors <- sequential_hcl(4, palette="Viridis")

  # Plot prevalence of UV-signature high for each age group
  ageTable <- prop.table(table(clin$Age_Group, clin$UV_sig), margin = 1)
  age <- as.data.frame(ageTable[,1])
  colnames(age) <- "High"
  agePlot <- ggplot(age, aes(x = rownames(age), y = High))+
    geom_col(aes(fill = rownames(age)), width = 0.6, position = "stack") +
    labs(fill = "Age Group", y = "Prevalence of Sig 7 High (%)", x = "Age Group") +
    coord_flip() +
    scale_fill_manual(values = agecolors) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
  agePlot
  ggsave("HighByAge.pdf", width = 6, height = 3, units = "in")
  
  # Stage plot (SKCM specific, needs to be formatted for different stage naming conventions to be applicable to other cohorts)
  if ("AJCC_PATHOLOGIC_TUMOR_STAGE" %in% colnames(clin)) {
    clin$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub('[ABC]', '', clin$AJCC_PATHOLOGIC_TUMOR_STAGE)
    clin$Tumor_Stage <- clin$AJCC_PATHOLOGIC_TUMOR_STAGE
    clin$Tumor_Stage[clin$Tumor_Stage == "Stage III" |
                                     clin$Tumor_Stage == "Stage IV"] <- "High Stage"
    clin$Tumor_Stage[is.na(clin$Tumor_Stage) == F & 
                                     clin$Tumor_Stage != "High Stage"] <- "Low Stage"

    # Plot prevalence of UV-signature high for each stage group
    stageDat <- table(clin$Tumor_Stage, clin$UV_sig)
    stageTable <- stageDat[,-3]
    stage <- as.data.frame(prop.table(stageTable, margin = 2))
    stagePlot <- ggplot(stage, aes(x = Var2, y = Freq)) + 
      geom_col(aes(fill = Var1), width = 0.6, position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = c("#E3242B", "#0A1172")) +
      labs(fill = "Stage", x = "UV Level", y = "Stage Prevalence (%)") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
    stagePlot
    ggsave("HighByStage.pdf", width = 6, height = 3, units = "in")
  }
  
  # Create feature list from user selected inputs, only allowing inputs which exist within clinical file
  features <- list()
  input <- ""
  while (input != "done") {
    input <- readline(prompt = "Enter name of clinical feature for analysis exactly as it appears in clinical info. \nEnter 'help' for a list of clinical features, or enter 'done' to exit: ")
    if (input == "help") {
      print(colnames(clin))
    } else if (input == "done") {
      features <- append(features, "UV_sig")
      print("Exiting, rest of function will now complete.")
    } else if (input %notin% colnames(clin)) { 
      print("Feature not found, please try again!")
    } else {
      features <- append(features, input)
      print("Feature added successfully!")
    }
  }

  # Create random palettes for each user selected feature, matching Age_Group palette to previous plot if selected
  cols <- hcl_palettes("sequential")
  for (i in 1:length(features)) {
    tClin <- clin[,c(features[[i]], "UV_sig_value")] 
    tClin$index <- 1:nrow(tClin)
    if (ceiling(length(unique(tClin[,1]))/4) > 1) {
      hVal <- 1.2+ceiling(length(unique(tClin[,1]))/4)*.1
    } else { hVal <- 1.2 }
    
    if (colnames(tClin)[1] == "Age_Group") {
      tPal <- sequential_hcl(4, palette="Viridis")
      hVal <- 1.2
    } else if (colnames(tClin)[1] != "UV_sig") {
      tPal <- sequential_hcl(length(unique(tClin[,1])), palette=sample(rownames(cols),1, replace = F))
    } else {
      tPal <- c("#FF0000", "#4F53B7", "#FFA500")
    }

    # Plot feature status per sample, sorted by decreasing UV signature level
    featureGraph <- ggplot(tClin, aes(x = index)) +
      geom_bar(aes(fill = tClin[,1])) +
      labs(fill = colnames(tClin)[1], x = "Signature 7 (decreasing)", y = NULL) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(values = tPal) +
      theme(legend.position = "bottom") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8),
            axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
    featureGraph
    ggsave(paste0(colnames(tClin)[1], ".pdf"), width = 8, height = hVal, units = "in")
  }
  print("Feature graphs generated successfully!")
}

plotControls <- function(clin, mutations) {
  mutFilter <- mutations[mutations$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
  mutFilter <- mutFilter[mutFilter$Hugo_Symbol %in% melControl, ]
  mutFilter <- mutFilter[mutFilter$Variant_Classification != "Silent", ]

  # Plot mutation status by UV level for each melanoma control gene, with specific conditions for BRAF
  for (i in 1:length(melControl)) {
    gene <- melControl[i]
    geneMut <- mutFilter[mutFilter$Hugo_Symbol == gene, ]
    mutated <- unique(mutFilter$Tumor_Sample_Barcode[mutFilter$Hugo_Symbol == gene])
    if (gene == "BRAF") {
      brafV <- unique(geneMut$Tumor_Sample_Barcode[geneMut$HGVSp_Short == "p.V600E"])
    }
    
    tempClin <- clin[, c("Tumor_Sample_Barcode", "UV_sig")]
    tempClin$Status <- "WT"
    tempClin$Status[tempClin$Tumor_Sample_Barcode %in% mutated] <- "Mut"
    if (gene == "BRAF") {
      tempClin$Status[tempClin$Status == "Mut"] <- "Other Mut"
      tempClin$Status[tempClin$Tumor_Sample_Barcode %in% brafV] <- "V600E"
      pal <- c("#228B22", "#960019", "#777B7E")
    } else {
      pal <- c("#228B22", "#960019")
    }
    
    mutTable <- as.data.frame(prop.table(table(tempClin$UV_sig, tempClin$Status), margin = 1))
    mutTable$Var1 <- factor(mutTable$Var1, levels = c("High", "Moderate", "Low"))
    if (gene == "BRAF") {
      mutTable$Var2 <- factor(mutTable$Var2, levels = c("WT", "V600E", "Other Mut"))
    } else {
      mutTable$Var2 <- factor(mutTable$Var2, levels = c("WT", "Mut"))
    }
    
    genePlot <- ggplot(mutTable, aes(x=Var1, y=Freq, fill=Var2)) +
      geom_col(position="dodge", color="black") +
      scale_fill_manual(values=pal) +
      scale_x_discrete(expand = c(0.005, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    genePlot
    ggsave(paste0(gene,".pdf"), width = 6, height = 3, units = "in")
  }
  print("Control gene graphs generated successfully!")
  return(clin)
}

plotSurvival <- function(clin) {
  # Plot survival curve for cohort, stratifying by UV signature level
  clin$OS_MONTHS <- as.numeric(clin$OS_MONTHS)
  clin$OS_STATUS <- as.numeric(substr(clin$OS_STATUS, 1, 1))
  
  ggsurvplot(fit = surv_fit(Surv(OS_MONTHS, OS_STATUS) ~ UV_sig, data = clin),
             xlab = "Months",
             ylab = "Overall survival probability",
             pval = T,
             palette = c("#FF0000", "#4F53B7", "#FFA500"))
  ggsave("Survival.pdf")
  
  print("Survival graph generated successfully!")
}

plot10gene <- function(clin, mutations) {
  # Plot 10 gene mutational signature (previously reported in another paper) by UV signature level 
  tenGene <- c("LRP1B", "GPR98", "XIRP2", "PKHD1L1", "USH2A", "DNAH9", "PCDH15", "DNAH10", "TP53", "PCDHAC1")
  initial <- mutations[mutations$Hugo_Symbol %in% tenGene, ]
  dat <- initial[ , c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")]
  dat$Variant_Classification[dat$Variant_Classification != "Silent"] <- 1
  dat$Variant_Classification[dat$Variant_Classification == "Silent"] <- 0
  dat2 <- as.data.frame(reshape(dat, idvar = "Tumor_Sample_Barcode", timevar = "Hugo_Symbol", direction = "wide"))
  dat2[is.na(dat2)] <- 0
  dat2[,2:11] <- as.numeric(unlist(dat2[,2:11]))
  final <- data.frame(Tumor_Sample_Barcode = dat2$Tumor_Sample_Barcode, Sig = apply(dat2[,2:11], 1, sum))
  
  t <- merge(clin, final, all.x = T)
  t$Sig[is.na(t$Sig)] <- 0
  t$Sig <- t$Sig / 10
  
  plotDat <- t[, c("UV_sig", "UV_sig_value", "Sig")]
  plotDat <- plotDat[order(-plotDat$UV_sig_value),]
  plotDat$UV_sig <- factor(plotDat$UV_sig, levels = unique(plotDat$UV_sig))

  corr <- cor.test(x=plotDat$UV_sig_value, y=plotDat$Sig, method = 'spearman')
  tenGenePlot <- ggplot(plotDat, aes(x=-UV_sig_value, y=Sig)) +
    geom_point(aes(color=UV_sig)) +
    annotate("text", x=-0.2, y=0.3, label= paste0("Spearman p-value: ",round(corr$p.value, digits = 4))) +
    annotate("text", x=-0.2, y=0.25, label= paste0("Spearman rho: ",round(corr$estimate, digits = 4))) +
    scale_x_reverse(expand = c(0.003,0)) +
    scale_y_discrete(expand = c(0, 0.003)) +
    scale_color_manual(values=c("#FF0000", "#FFA500", "#4F53B7")) +
    theme_classic() +
    theme(legend.position = "none")
  print(tenGenePlot)
  ggsave("10_Gene_Plot.pdf", height = 3, width = 6, units = "in")
}

# Load required data (all generated by processData.R) and execute plot generation
generatePlots <- function(dataDir, subDir = "") {
  setwd(paste0(homeDir, "/data/", dataDir))
  clin <- read.csv("clinical.csv")
  clin <- clin[order(-clin$UV_sig_value),]
  sigs <- read.csv("mutationalSignatures.csv", row.names = "SAMPLE_ID")
  tnm <- read.csv("tnm.csv")
  if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
  } else { mutations <- fread("data_mutations.txt") }
  mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B" 
  
  if (subDir != "") {
    resultsDir <- paste0(homeDir,"/results/",dataDir,"/",subDir,"/")
    if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
    setwd(resultsDir)
  } else {
    resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
    if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
    setwd(resultsDir)
  }
  
  plotSigs(sigs, clin)
  plotTNM(clin, tnm)
  plotSNVs(clin, mutations)
  plotFeatures(clin)
  plotControls(clin, mutations)
  plotSurvival(clin)
  plot10gene(clin, mutations)
  
  print("All plots generated successfully!")
}

## Execute code -------------------
generatePlots("skcm_tcga")
