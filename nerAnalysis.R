## Library calls, load gene lists and cosmic signatures ##----------------------
library(ggplot2)
library(colorspace)
library(data.table)
library(maftools)
library(dplyr)
source("C:/Users/robmc/Desktop/NUMB/repairGenes.R") 
'%notin%' <- Negate('%in%')
homeDir <- "C:/Users/robmc/Desktop/NUMB_files"

## Function area  --------------------------------------------------------
dataDir <- "skcm_tcga"
setwd(paste0(homeDir, "/data/", dataDir))
clin <- read.csv("clinical.csv")
clin <- clin[order(-clin$UV_sig_value),]
if (!file.exists("data_mutations.txt")) { mutations <- fread("data_mutations_extended.txt") 
} else { mutations <- fread("data_mutations.txt") }
mutations$Hugo_Symbol[mutations$Hugo_Symbol == "CUL4A" | mutations$Hugo_Symbol == "CUL4B"] <- "CUL4A/B" 
cnv <- fread("data_cna.txt")
cnvSig <- read.csv("NER_cnvSig.csv")
mutSig <- read.csv("NER_mutSig.csv")

resultsDir <- paste0(homeDir,"/results/",dataDir,"/")
if (!dir.exists(resultsDir)) {  dir.create(resultsDir)  }
setwd(resultsDir)

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

maf <- read.maf(maf = mutations,
                clinicalData = clin,
                cnTable = cnv2,
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

tempClin <- clin[, c("Tumor_Sample_Barcode", "UV_sig", "UV_sig_value", "NER_sig")]
tempClin$NER_sig[tempClin$NER_sig > 1] <- 1
alteredTable <- prop.table(table(tempClin$UV_sig, tempClin$NER_sig), margin = 1)

p <- ggplot(temp, aes(x = temp[,2], y = temp[,1], fill = temp[,2])) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#FF0000", "#0000FF")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(y = "UV Signature", legend = "NER Signature") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size = 8))
p
ggsave(paste0(gene,".pdf"), height = 6, width = 6, units = "in")

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