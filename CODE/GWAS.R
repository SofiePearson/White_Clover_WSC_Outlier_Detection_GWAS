## GWAS for WSC, leaf area and other nutritional attributes:

## Extract samples with phenotypes out of vcf file:
# vcftools --remove-indels --keep WSC_pheno_samples.txt --recode --recode-INFO-all --vcf ../ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > WSC_pheno_samples.vcf
## Kept 605 individuals

library(rrBLUP)
library(vcfR)
library(plyr)
library(calibrate)
library(readr)
library(dplyr)
library(tibble)
library(calibrate)
library(adegenet)
library(poppr)
library(reshape2)
library(ggplot2)
library(scales)

########################
## Import data
########################
pheno_traits_geno <- read.vcfR("WSC_pheno_samples.vcf")

########################
## create matrix for genotype data
########################
## For rrBLUP, data needs to be in 1 = homo, 0 = hetero, -1 = homo
pheno_traits_geno_gt <- extract.gt(pheno_traits_geno, element = "GT")
pheno_traits_geno_gt <- revalue(pheno_traits_geno_gt, c("0/0"=1, "1/1"=-1, "0/1"=0, "1/0"=0))

## Missing data encoded as NA
pheno_traits_geno_gt[is.na(pheno_traits_geno_gt)] <- NA

pheno_traits_geno_gt_num=matrix(as.numeric(pheno_traits_geno_gt),nrow=dim(pheno_traits_geno_gt)[1])
colnames(pheno_traits_geno_gt_num)=colnames(pheno_traits_geno_gt)
rownames(pheno_traits_geno_gt_num)=rownames(pheno_traits_geno_gt)

## swap rows and columns around
pheno_traits_num_rrBLUP_format <- t(pheno_traits_geno_gt_num)

## Find order of samples
pheno_traits_pop <- gsub("-.$","",colnames(pheno_traits_geno@gt),fixed = FALSE)
pheno_traits_pop <- pheno_traits_pop[-1]
write.table(pheno_traits_pop, file = "Pheno_sample_order.txt")

## Find samples that have both phenotype and genotype information.
## Then save phenotype file with the following headers: gid, WSC, SSS, Ash, CP, NDF, ADF, Lipid and LA
## For samples with LA and no NIRS data - add "NA" and vice versa.
## Sort samples in the same order as VCf file.
## And save the phenotype file as "GWAS_rrBLUP_pheno_format.txt"
## Remove the gid column and save as "myrrBLUP_pheno.txt"

## Remove sample names and Chr 17 SNPs then save as "myrrBLUP_geno.txt"
str(pheno_traits_num_rrBLUP_format)
pheno_traits_num_rrBLUP_format_table <- as_tibble(pheno_traits_num_rrBLUP_format)
myrrBLUP_chr17 <- pheno_traits_num_rrBLUP_format_table %>% dplyr::select(-contains("S17"))
write.table(myrrBLUP_chr17, file = "myrrBLUP_geno.txt", sep = "\t")
dim(myrrBLUP_chr17)
# [1]   605 14626

########################
## Import required data files
########################
myrrBLUP <- as.matrix(read.table(file = "myrrBLUP_geno.txt", header = TRUE)) ## Genotype data
myrrBLUP_pheno <- as.matrix(read.table(file = "myrrBLUP_pheno.txt", header = TRUE)) ## Phenotype data
## Check that the files loaded in correctly
dim(myrrBLUP)
#[1]  605 14626
dim(myrrBLUP_pheno)
#[1] 605   8

########################
## Impute missing markers using A.mat command:
########################
imputeEM = A.mat(myrrBLUP, max.missing = 0.5, impute.method = "EM", return.imputed = T)
# [1] "A.mat converging:"
# [1] 0.0258
# [1] 0.00606
markers_imputeEM = imputeEM$imputed
dim(markers_imputeEM)
#[1]  605 5757
write.table(markers_imputeEM, file = "markers_imputeEM.txt", sep = "\t")
## 5757 markers left that were imputed

## swap rows and columns around
markers_imputeEM_swap <- t(markers_imputeEM)

## Add in SNP_ID
colnames(markers_imputeEM_swap) <- pheno_traits_pop

## Add in Chromosome and position information
SNP_ID <- rownames(markers_imputeEM_swap)
SNP_ID <- as.data.frame(SNP_ID)
CHR <- SNP_ID %>% tidyr::separate(SNP_ID, c("Chromosome","Position"))
CHR$Chromosome <- gsub("S","", as.character(CHR$Chromosome))
GWAS_rrBLUP_geno_format_imputed_EM <- cbind(SNP_ID,CHR,markers_imputeEM_swap)

## Save as "GWAS_rrBLUP_geno_format_imputed_EM.txt"
write.table(GWAS_rrBLUP_geno_format_imputed_EM, file = "GWAS_rrBLUP_geno_format_imputed_EM.txt", quote = F, sep = "\t", col.names= T, row.names = F)

########################
## PCA for GWAS samples
########################

## input 5757 SNP and 605 sample data from above

GWAS_PCA <- read.csv("PCA_samples_SNPs.csv")
GWAS_PCA[GWAS_PCA == 1] <- "AA"
GWAS_PCA[GWAS_PCA == 0] <- "AC"
GWAS_PCA[GWAS_PCA == -1] <- "CC"

head(GWAS_PCA[,c(1:10)])
str(GWAS_PCA)
ind <- as.character(GWAS_PCA$sample) ## individual labels
population <- as.character(GWAS_PCA$Pool) ## K=15 clustering results
locus <- GWAS_PCA[, -c(1,2)]

gen = df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
gen$tab[1:5, 1:10]
gen
#' /// GENIND OBJECT /////////
#'   
#'   // 605 individuals; 5,757 loci; 11,514 alleles; size: 29.7 Mb
#' 
#' // Basic content
#' @tab:  605 x 11514 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-2)
#' @loc.fac: locus factor for the 11514 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: df2genind(X = locus, sep = "", ind.names = ind, pop = population, 
#'                  ploidy = 2)
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 83-137)

x <- tab(gen, NA.method = "mean")
pca1 <- dudi.pca(x, scannf = FALSE, nf = 10)

percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,2),
        names.arg = round(percent, 1))

## Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

## Add a column containing individuals
ind_coords$Ind = indNames(gen)

## Add a column with the Pool IDs
ind_coords$Pool = gen$pop

## Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3, Axis4,Axis5,Axis6,Axis7,Axis8,Axis9,Axis10) ~ Pool, data = ind_coords, FUN = mean)

## Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Pool", suffix = c("",".cen"))

## Define colour palette
cols = c("#42d4f4","#e6194b","#ffe119","#4363d8","#3cb44b")


## Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

## Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

## Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Pool), show.legend = FALSE)+
  ## points
  geom_point(aes(fill = Pool), shape = 21, size = 3, show.legend = FALSE)+
  ## centroids
  geom_label(data = centroid, aes(label = Pool, fill = Pool), size = 4, show.legend = FALSE)+
  ## colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  ## custom labels
  labs(x = xlab, y = ylab)+
    ## custom theme
  ggtheme

## Save as "GWAS_samples_PCA_center_names.tiff" w = 956  H = 796 (Supplementary Figure 3, top)

########################
## Test the GWAS
########################
pheno <- as.data.frame(read.table("GWAS_rrBLUP_pheno_format.txt", header = TRUE))
# GWAS_rrBLUP_geno_format_imputed_EM <- read_delim("GWAS_rrBLUP_geno_format_imputed_EM.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

geno2 <- GWAS_rrBLUP_geno_format_imputed_EM ## OR use "GWAS_rrBLUP_geno_format_imputed_EM" (can skip the remove row names step on next line)
rownames(geno2) <- c() ## remove row names
geno2 <- as.data.frame(geno2)
geno2$Chromosome <- as.numeric(geno2$Chromosome)
geno2$Position <- as.numeric(geno2$Position)

## Run GWAS:
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
# [1] "GWAS for trait: WSC"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: SSS"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: Ash"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: CP"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: NDF"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: ADF"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: Lipid"
# [1] "Variance components estimated. Testing markers."
# [1] "GWAS for trait: LA"
# [1] "Variance components estimated. Testing markers."

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute.txt", sep = "\t")

########################
## Run GWAS seperately for each trait
########################
## Get data in correct format:
## Using the "GWAS_rrBLUP_pheno_format.txt" file, split so there are 8 text files. Each one has "gid" and "[trait]" headers
## e.g. 
##   gid	  WSC
## 40-06	20.05
## 40-05	19.47
## and save each trait as "myrrBLUP_WSC.txt" "myrrBLUP_SSS.txt" etc

########################
## Run GWAS for SSS
########################
pheno <- as.data.frame(read.table("myrrBLUP_SSS.txt", header = TRUE))
hist(pheno$SSS, xlab = (expression("SSS  (g kg"^"-1"*' DM)')),col = "white", main = NULL) ## (Supplementary Figure 3, bottom)

#GWAS_rrBLUP_geno_format_imputed_EM <- read_delim("GWAS_rrBLUP_geno_format_imputed_EM.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#geno2 <- GWAS_rrBLUP_geno_format_imputed_EM
#geno2 <- as.data.frame(geno2)

par(mfrow=c(1,1))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_SSS.png" W = 350 H = 400

## Save output from the GWAS. This will be used to create Manhattan plot
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_SSS.txt", sep = "\t")

## To create Manhattan input files:
## have to change column names to "SNP", CHR", "BP" and "P"
## and convert -log10(p) back into just p for each trait
## Use previous manhattan plot code to create the manhattan plots
## The cut off for 5% is 5.06 on the -log10 scale (0.05/5757 = 8.68508E-06, then -log10(8.68508E-06) = 5.061226225)

SSS <- rrBLUP_GWAS_EM_impute
colnames(SSS) <- c("SNP", "CHR", "BP", "P")
SSS$P <- 10^(-(SSS$P))
head(SSS$P)
write.table(SSS, "GWAS_manhattan_SSS.txt", sep = "\t")

########################
## Run GWAS for WSC
########################
pheno <- as.data.frame(read.table("myrrBLUP_WSC.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_WSC.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_WSC.txt", sep = "\t")

## Create Manhattan input file
WSC <- rrBLUP_GWAS_EM_impute
colnames(WSC) <- c("SNP", "CHR", "BP", "P")
WSC$P <- 10^(-(WSC$P))
write.table(WSC, "GWAS_manhattan_WSC.txt", sep = "\t")

########################
## Run GWAS for Ash
########################
pheno <- as.data.frame(read.table("myrrBLUP_Ash.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_Ash.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_Ash.txt", sep = "\t")

## Create Manhattan input file
Ash <- rrBLUP_GWAS_EM_impute
colnames(Ash) <- c("SNP", "CHR", "BP", "P")
Ash$P <- 10^(-(Ash$P))
write.table(Ash, "GWAS_manhattan_Ash.txt", sep = "\t")

########################
## Run GWAS for CP
########################
pheno <- as.data.frame(read.table("myrrBLUP_CP.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_CP.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_CP.txt", sep = "\t")

## Create Manhattan input file
CP <- rrBLUP_GWAS_EM_impute
colnames(CP) <- c("SNP", "CHR", "BP", "P")
CP$P <- 10^(-(CP$P))
write.table(CP, "GWAS_manhattan_CP.txt", sep = "\t")

########################
## Run GWAS for NDF
########################
pheno <- as.data.frame(read.table("myrrBLUP_NDF.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_NDF.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_NDF.txt", sep = "\t")

## Create Manhattan input file
NDF <- rrBLUP_GWAS_EM_impute
colnames(NDF) <- c("SNP", "CHR", "BP", "P")
NDF$P <- 10^(-(NDF$P))
write.table(NDF, "GWAS_manhattan_NDF.txt", sep = "\t")

########################
## Run GWAS for ADF
########################
pheno <- as.data.frame(read.table("myrrBLUP_ADF.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_ADF.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_ADF.txt", sep = "\t")

## Create Manhattan input file
ADF <- rrBLUP_GWAS_EM_impute
colnames(ADF) <- c("SNP", "CHR", "BP", "P")
ADF$P <- 10^(-(ADF$P))
write.table(ADF, "GWAS_manhattan_ADF.txt", sep = "\t")

########################
## Run GWAS for Lipid
########################
pheno <- as.data.frame(read.table("myrrBLUP_Lipid.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_Lipid.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_Lipid.txt", sep = "\t")

## Create Manhattan input file
Lipid <- rrBLUP_GWAS_EM_impute
colnames(Lipid) <- c("SNP", "CHR", "BP", "P")
Lipid$P <- 10^(-(Lipid$P))
write.table(Lipid, "GWAS_manhattan_Lipid.txt", sep = "\t")

########################
## Run GWAS for LA
########################
pheno <- as.data.frame(read.table("myrrBLUP_LA.txt", header = TRUE))
rrBLUP_GWAS_EM_impute <- GWAS(pheno, geno2, n.PC = 2, min.MAF = 0.03, plot = TRUE)
## save QQ plot as "rrBLUP_QQ_impute_EM_LA.png" W = 350 H = 400

## Save output from the GWAS
write.table(rrBLUP_GWAS_EM_impute, file = "GWAS_results_EM_impute_LA.txt", sep = "\t")

## Create Manhattan input file
LA <- rrBLUP_GWAS_EM_impute
colnames(LA) <- c("SNP", "CHR", "BP", "P")
LA$P <- 10^(-(LA$P))
write.table(LA, "GWAS_manhattan_LA.txt", sep = "\t")

## Copied the output from each "GWAS_results_EM_impute_[TRAIT].txt" file into spreadsheet called "p-values_from_GWAS" in "GWAS_outliers.xlsx" excel file
## saved the manhattan input files in "Manhattan_input_for_each_trait" spreadsheet in "GWAS_outliers.xlsx" excel file
## Keep the Q-Q plots produced above (combined to create Supplementary Figure 11). But use the following manhattan plot code to create the manhattan plots
## The cut off for 5% is 5.06 on the -log10 scale (0.05/5757 = 8.68508E-06, then -log10(8.68508E-06) = 5.061226225) (see "Manhattan_input_for_each_trait")
## Also look at SNPs above 3 on the -log10(p) axis
## Used output to create Table 2 and Supplementary Table 11

# WSC <- read.table("GWAS_manhattan_WSC.txt", header = TRUE)
# SSS <- read.table("GWAS_manhattan_SSS.txt", header = TRUE)
# NDF <- read.table("GWAS_manhattan_NDF.txt", header = TRUE)
# ADF <- read.table("GWAS_manhattan_ADF.txt", header = TRUE)
# CP <- read.table("GWAS_manhattan_CP.txt", header = TRUE)
# Lipid <- read.table("GWAS_manhattan_Lipid.txt", header = TRUE)
# Ash <- read.table("GWAS_manhattan_Ash.txt", header = TRUE)
# LA <- read.table("GWAS_manhattan_LA.txt", header = TRUE)

## Run the following manhattan function: 
## NOTE: Need to change the suggestive line and genomewide line values based on number of markers (SUGGESTIVE = -log10(p) > 3 AND GENOMEWIDE = 0.05)

manhattan_p_value <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                              col=c("#e6194b", "#fabebe", "#f58231","#ffe119", "#aaffc3", "#3cb44b", "#46f0f0", "#4363db"), chrlabs=NULL,
                              suggestiveline=-log10(0.001), genomewideline=-log10(8.68508E-06), 
                              highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = FALSE, ...) {
  
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=19, cex = 2,
                   xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))),
                   xlab=xlabel, ylab=expression(-log[10](italic(p)-value)))
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=col[1], ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=19, cex = 2, ...))
      icol=icol+1
    }
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
  
  # Highlight snps from a character vector
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col="#000000", pch=13, cex = 2,...)) 
  }
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), 
           textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.7), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.7, ...)
    }
  }  
  par(xpd = FALSE)
}

par(mfrow=c(4,2))
manhattan_p_value(WSC,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(A)~~"Water-soluble carbohydrate"), adj = 0, cex.main = 2)

manhattan_p_value(Ash,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(E)~~"Ash"), adj = 0, cex.main = 2)

manhattan_p_value(SSS,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(B)~~"Soluble sugars and starches"), adj = 0, cex.main = 2)

manhattan_p_value(NDF,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(F)~~"Neutral detergent fibre"), adj = 0, cex.main = 2)

manhattan_p_value(CP,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(C)~~"Crude protein"), adj = 0, cex.main = 2)

manhattan_p_value(ADF,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(G)~~"Acid detergent fibre"), adj = 0, cex.main = 2)

manhattan_p_value(Lipid,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(D)~~"Lipid"), adj = 0, cex.main = 2)

manhattan_p_value(LA,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(H)~~"Leaf area"), adj = 0, cex.main = 2)

## save as "GWAS_Manhattan_plots.TIFF" w = 1600 H = 900

########################
## Create final plots
########################

## Highlight the twelve SNPs ID greater than -log10(p) of 3 for WSC and SSS:
par(mfrow=c(2,1))
manhattan_p_value(WSC,annotatePval = FALSE, ylim =c(0,8), highlight = c("S1_2338028","S6_16865077","S11_7689073","S5_47903593","S1_13746102","S9_31736793","S3_51874121","S3_51874123","S4_30841647","S9_23070656","S1_102715","S8_27636527"))
title(expression(bold(A)~~"Water-soluble carbohydrate"), adj = 0, cex.main = 1.5)

manhattan_p_value(SSS,annotatePval = FALSE, ylim =c(0,8), highlight = c("S1_2338028","S6_16865077","S11_7689073","S5_47903593","S1_13746102","S9_31736793","S3_51874121","S3_51874123","S4_30841647","S9_23070656","S1_102715","S8_27636527"))
title(expression(bold(B)~~"Soluble sugars and starches"), adj = 0, cex.main = 1.5)

## save as "GWAS_Manhattan_plots_12_SNPs_highlighted.TIFF" w = 1600 H = 900 (Figure 6)
## save as "GWAS_Manhattan_plots_12_SNPs_highlighted.PDF" A4 landscape

par(mfrow=c(3,2))
manhattan_p_value(CP,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(A)~~"Crude protein"), adj = 0, cex.main = 2)

manhattan_p_value(NDF,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(D)~~"Neutral detergent fibre"), adj = 0, cex.main = 2)

manhattan_p_value(Lipid,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(B)~~"Lipid"), adj = 0, cex.main = 2)

manhattan_p_value(ADF,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(E)~~"Acid detergent fibre"), adj = 0, cex.main = 2)

manhattan_p_value(Ash,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(C)~~"Ash"), adj = 0, cex.main = 2)

manhattan_p_value(LA,annotatePval = FALSE, ylim =c(0,8))
title(expression(bold(F)~~"Leaf area"), adj = 0, cex.main = 2)

## save as "GWAS_Manhattan_plots_6.TIFF" w = 1600 H = 900 (Supplementary Figure 12)
## save as "GWAS_Manhattan_plots_6.PDF" A4 landscape