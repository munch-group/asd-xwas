library(dplyr)
library(ggplot2)

setwd("/faststorage/jail/project/ChrXh2/shannon/data/")

data <- read.table("1000G/20130606_g1k.ped", header = T, sep = "\t")

# If the column data$Population contains:
# CHB JPT CHS CDX or KHV: code as EAS
# CEU TSI FIN GBR or IBS: code as EUR
# YRI LWK MAG MSL ESN ASW ACB or GWD: code as AFR
# MXL PUR CLM or PEL: code as AMR 
# GIH PJL BEB STU or ITU: code as SAS

data$Super_Population <- with(data, ifelse(Population %in% c("CHB", "JPT", "CHS", "CDX", "KHV"), "EAS",
                                           ifelse(Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS"), "EUR",
                                                  ifelse(Population %in% c("YRI", "LWK", "MAG", "MSL", "ESN", "ASW", "ACB", "GWD"), "AFR",
                                                         ifelse(Population %in% c("MXL", "PUR", "CLM", "PEL"), "AMR",
                                                                ifelse(Population %in% c("GIH", "PJL", "BEB", "STU", "ITU"), "SAS", NA))))))

# Merge with 1000G genetic data to ensure that the IDs match up
genetic <- read.table("1000G/1000G.fam", header = F)
merged <- merge(data, genetic, by.x = "Individual.ID", by.y= "V2")

# Save the output 
G1000 <- merged[, c("V1", "Individual.ID", "Gender", "Population", "Super_Population")]
colnames(G1000) <- c("FID", "IID", "Gender", "Population", "Super_Population")
write.table(G1000, file = "1000G/1000G_populations.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# Read in the eigenvectors
data <- read.table("PCA_iPSYCH2015_1000G_danes.eigenvec", header = F)
colnames(data) <- c("FID", "IID", 
                    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
# Read in iPSYCH IDs
iPSYCH <- read.table("/faststorage/jail/project/ChrXh2/data/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn.fam", header = F)
colnames(iPSYCH) [1:2] <- c("FID", "IID") 

# Add pop_group column based on Super_Population
G1000$pop_group <- with(G1000, ifelse(Super_Population == "EAS", "East Asia",
                                      ifelse(Super_Population == "EUR", "Europe",
                                             ifelse(Super_Population == "AFR", "Africa",
                                                    ifelse(Super_Population == "AMR", "Admixed America",
                                                           ifelse(Super_Population == "SAS", "South Asia", NA))))))

# Merge iPSYCH data for visualization
data <- merge(data, G1000[, c("IID", "pop_group")], by = "IID", all.x = TRUE)
data$pop_group[data$IID %in% iPSYCH$IID] <- "iPSYCH2015"

# Define group colours
group_colors <- c("iPSYCH2015" = "black", "Africa" = "#F5C710", "Admixed America" = "#2297E6", 
                  "East Asia" = "#CD0BBC", "South Asia" = "#61D04F", "Europe" = "#DF536B")


# Calculate mean and sd for EUR group
eur_mean_PC1 <- mean(data$PC1[data$pop_group == "Europe"], na.rm = TRUE)
eur_sd_PC1 <- sd(data$PC1[data$pop_group == "Europe"], na.rm = TRUE)
eur_mean_PC2 <- mean(data$PC2[data$pop_group == "Europe"], na.rm = TRUE)
eur_sd_PC2 <- sd(data$PC2[data$pop_group == "Europe"], na.rm = TRUE)

# Plot PC1 against PC2, colour according to ancestry grouping
# Plot the black points (iPSYCH2015) below the coloured points (1000 Genomes population groups)
# Add lines for the PC cutoffs
png("../results/iPSYCH_1000G_danes_PCA_v1.png", width = 7, height = 6, units = "in", res = 300)
ggplot(data) +
  geom_point(data = subset(data, pop_group == "iPSYCH2015"), aes(x = PC1, y = PC2, color = pop_group), size = 2, alpha = 0.3) +
  geom_point(data = subset(data, pop_group != "iPSYCH2015"), aes(x = PC1, y = PC2, color = pop_group), size = 2, alpha = 0.3) +
  scale_color_manual(name = ' ', values = group_colors,
                     breaks = c("East Asia", "Admixed America", "South Asia", "Europe", "Africa", "iPSYCH2015")) +
  geom_vline(xintercept = eur_mean_PC1 + 6 * eur_sd_PC1, linetype = "dashed", color = "gray62") +
  geom_vline(xintercept = eur_mean_PC1 - 6 * eur_sd_PC1, linetype = "dashed", color = "gray62") +
  geom_hline(yintercept = eur_mean_PC2 + 6 * eur_sd_PC2, linetype = "dashed", color = "gray62") +
  geom_hline(yintercept = eur_mean_PC2 - 6 * eur_sd_PC2, linetype = "dashed", color = "gray62") +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_family = "Arial")
dev.off()

# Plot the coloured points (1000 Genomes population groups) below the black points (iPSYCH2015)
png("../results/iPSYCH_1000G_danes_PCA_v2.png", width = 7, height = 6, units = "in", res = 300)
ggplot(data) +
  geom_point(data = subset(data, pop_group != "iPSYCH2015"), aes(x = PC1, y = PC2, color = pop_group), size = 2, alpha = 0.3) +
  geom_point(data = subset(data, pop_group == "iPSYCH2015"), aes(x = PC1, y = PC2, color = pop_group), size = 2, alpha = 0.3) +
  scale_color_manual(name = ' ', values = group_colors,
                     breaks = c("East Asia", "Admixed America", "South Asia", "Europe", "Africa", "iPSYCH2015")) +
  geom_vline(xintercept = eur_mean_PC1 + 6 * eur_sd_PC1, linetype = "dashed", color = "gray62") +
  geom_vline(xintercept = eur_mean_PC1 - 6 * eur_sd_PC1, linetype = "dashed", color = "gray62") +
  geom_hline(yintercept = eur_mean_PC2 + 6 * eur_sd_PC2, linetype = "dashed", color = "gray62") +
  geom_hline(yintercept = eur_mean_PC2 - 6 * eur_sd_PC2, linetype = "dashed", color = "gray62") +
  labs(x = "PC1", y = "PC2") +
  theme_classic(base_family = "Arial")
dev.off()

# Exclusion filtering based on PC1 and PC2 thresholds
PC1_upper <- eur_mean_PC1 + 6 * eur_sd_PC1
PC1_lower <- eur_mean_PC1 - 6 * eur_sd_PC1
PC2_upper <- eur_mean_PC2 + 6 * eur_sd_PC2
PC2_lower <- eur_mean_PC2 - 6 * eur_sd_PC2

# Count how many individuals would be excluded
exclude <- length(which(data$pop_group == "iPSYCH2015")) - length(which(data$PC1<PC1_upper & data$PC1>PC1_lower & data$PC2<PC2_upper & data$PC2>PC2_lower & data$pop_group == "iPSYCH2015"))

if(exclude > 0) {
  print("iPSYCH2015 individuals fall outside of 6 SD of the EUR mean, consider excluding")
  keep <- which(data$PC1 < PC1_upper & data$PC1 > PC1_lower & data$PC2 < PC2_upper & data$PC2 > PC2_lower & data$pop_group == "iPSYCH2015")
  write.table(data[keep, c(2,1) ], file = "iPSYCH_danes_EUR.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  ("No iPSYCH2015 individuals fall outside of 6 SD of the EUR mean for PC1 and PC2 ")
}
