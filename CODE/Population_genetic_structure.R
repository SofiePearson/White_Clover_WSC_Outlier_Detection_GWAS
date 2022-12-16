library(readr)
library(ggplot2)
library(ggpubr)
library(vcfR)
library(adegenet)
library(hierfstat)
library(dplyr)
library(poppr)

############################################################################
## Check positive S9 samples across GBS libraries and duplicate samples using PCA
########################
## Group samples based on pool and have "control" samples seperate
data <- read.vcfR(file = "1113_samples_and_S9_no_multiallelic_dp5to150_mm0.2_maf0.03.vcf")
my_genlight <- vcfR2genlight(data)
print(my_genlight)
ploidy(my_genlight) <- 2
## View sample order
my_genlight$ind.names
write.table(my_genlight$ind.names, file = "genlight_names.txt", sep = "\t")
## make pop file:
## assigned population names and Pool names to the samples and saved a .txt file as:
## "1113_samples_and_S9_sample_pop_Pool.txt"
pop <- read.table(file = "1113_samples_and_S9_sample_pop_Pool.txt", header = TRUE) ## File containing list of sample names in one column and list of corresponding populations in the second column. The samples should be ordered the same as the vcf file.
## Check sample names are the same between the vcf file and pop file:
all(colnames(data@gt)[-1] == pop$sample)
# [1] TRUE
## Confirm the ID's are ordered the same in my_genlight and pop:
my_genlight$ind.names
pop$sample

## Assign pop info to genlight object
pop(my_genlight) <- pop$Pool
my_genlight

## Sort pools so the order is: WNZLL, WNZSL, WUSLL, FNZLL, FNZSL and then Control
levels(my_genlight$pop)
# [1] "Control-1" "Control-2" "Control-3" "FNZLL"     "FNZSL"     "WNZLL"     "WNZSL"     "WUSLL"  
my_genlight$pop <- factor(my_genlight$pop,
                          levels(my_genlight$pop)[c(6,7,8,4,5,1,2,3)])
levels(my_genlight$pop)
# [1] "WNZLL"     "WNZSL"     "WUSLL"     "FNZLL"     "FNZSL"     "Control-1" "Control-2" "Control-3"
## View example of SNP names
my_genlight$loc.names[1:10]
# [1] "S1_94724"  "S1_102715" "S1_102724" "S1_178752" "S1_182660" "S1_182664" "S1_182668" "S1_182672"
# [9] "S1_182673" "S1_182678"

## Compute the PCA on the SNP genotypes and plot it:
pca1 <- glPca(my_genlight, nf=4) # nf = number of PC axes to retain (here, 4) - this takes a while
pca1 # prints summary

barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$pop <- pop(my_genlight)

cols <- c("#e6194b", "#ffe119", "#3cb44b", "#4363d8", "#42D4F4", "black", "#f032e6", "#f58231") # Same colours as phenotype figures for each pool, used black for S9, magenta for LM control and orange for LE control)

set.seed(9)
p <- ggplot(pca1.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=3)
p <- p + scale_color_manual(values = cols) 
p <- p + theme_bw()
p

## Add % variance to plot:
sum(pca1$eig[1:4]) ## Only used first 4 PCs
# 27.44863
View(pca1) ## The first two eigenvalues correspond to PC1 and PC2
## PC1 = (8.0000058/27.44863) x 100
PC1var <- 29.15
## PC2 = (6.7159123/27.44863) x 100
PC2var <- 24.47

set.seed(9)
p <- ggplot(pca1.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=4,alpha = 0.6)
p <- p + scale_color_manual(values = cols, name = "") 
p <- p + theme_bw(base_size = 18)
p <- p + xlab(paste0("PC1 (",PC1var,"% variance)"))
p <- p + ylab(paste0("PC2 (",PC2var,"% variance)"))
p
## Save as "PCA.tiff" w = 956  H = 796

## Order the plotted data so the three controls are plotted last and come out on top of the other points
## Also re-label key so controls dont have the "-" inbetween
set.seed(9)
p <- ggplot(pca1.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=4,alpha = 0.6)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-1'),
                    size =4)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-2'),
                    size =4)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-3'),
                    size =4)
p <- p + scale_color_manual(values = cols, name = "", labels = c("WNZLL", "WNZSL", "WUSLL", "FNZLL", "FNZSL", "Control 1", "Control 2", "Control 3")) 
p <- p + theme_bw(base_size = 18)
p <- p + xlab(paste0("PC1 (",PC1var,"% variance)"))
p <- p + ylab(paste0("PC2 (",PC2var,"% variance)"))
p
## Save as "PCA_FINAL.tiff" w = 956  H = 796

## Change the shape depending on Generation and Selection

## Low-End  square 15
## Low-Mid  circle 19
## P        asterisk 8
## High-Mid diamond 18
## High-End triangle 17

## want Control-1 as a circle so in the text file change Control-1 to LM

## Assign pop info to genlight object
pop <- read.table(file = "1113_samples_and_S9_sample_pop_Pool_Gen.txt", header = TRUE) 
pop(my_genlight) <- pop$Gen
my_genlight

## sort Generation
levels(my_genlight$pop)
# [1] "HE"        "HM"        "LE"        "LM"        "P" 
my_genlight$pop <- factor(my_genlight$pop,
                          levels(my_genlight$pop)[c(3,4,5,2,1)])
levels(my_genlight$pop)
# [1] "LE" "LM" "P"  "HM" "HE"
pca1.scores$Gen <- pop(my_genlight)
head(pca1.scores)
#              PC1        PC2       PC3         PC4   pop Gen
# 40-02 -2.0679006 -0.5842705 0.1620086  0.13588941 FNZSL   P
# 38-44  0.5906719  0.6616968 1.7318241 -1.52611201 WNZSL  LE
# 38-45  1.1991103  0.9572461 0.8562993 -1.21483213 WNZSL  LE
# 38-46  1.1613051  1.0014304 1.3663491 -1.24796882 WNZSL  LE
# 38-47  0.9421142  0.9597866 1.0023154 -1.22106630 WNZSL  LE
# 40-06 -1.2366726 -0.9418451 0.2094897 -0.03502158 FNZSL   P

set.seed(9)
p <- ggplot(pca1.scores, aes(x=PC1, y=PC2, colour=pop, shape = Gen)) 
p <- p + geom_point(size=4,alpha = 0.6)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-1'),
                    size =4)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-2'),
                    size =4)
p <- p + geom_point(data = subset(pca1.scores, pop == 'Control-3'),
                    size =4)
p <- p + scale_color_manual(values = cols, name = "", labels = c("WNZLL", "WNZSL", "WUSLL", "FNZLL", "FNZSL", "Control 1", "Control 2", "Control 3")) 
p <- p + theme_bw(base_size = 18)
p <- p + xlab(paste0("PC1 (",PC1var,"% variance)"))
p <- p + ylab(paste0("PC2 (",PC2var,"% variance)"))
p <- p + scale_shape_manual(name = "", 
                            labels = c("Low-End","Low-Mid","Parent","High-Mid","High-End"), 
                            values = c(15,19,8,18,17))
p
## Save as "PCA_col_pool_shape_gen.tiff" w = 956  H = 796

############################################################################
## Pairwise FST matrix
########################
## Input data
data <- read.vcfR(file = "ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf") 
my_genind <- vcfR2genind(data)

########################
## Pairwise FST matrix - a priori (i.e. 24 pops)
########################
pop <- read.table(file = "24_SAMPLE_POP.txt", header = TRUE) ## File containing list of sample names in one column and list of corresponding populations in the second column. The samples should be ordered the same as the vcf file.
## Assign the pop column to genind object
pop(my_genind) <- pop$pop 

## Check genind object is correct
my_genind 
#' /// GENIND OBJECT /////////
#'   
#'   // 1,113 individuals; 14,743 loci; 29,486 alleles; size: 133.1 Mb
#' 
#' // Basic content
#' @tab:  1113 x 29486 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-2)
#' @loc.fac: locus factor for the 29486 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: adegenet::df2genind(X = t(x), sep = sep)
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 41-48)

## Convert genind to hierfstat format
my_hierfstat <- genind2hierfstat(my_genind) 

## Create pairwise FST matrix for 24 populations
genet.dist(my_hierfstat, method = "WC84") ## This took 4-5 hours!

## Used output to create Supplementary Table 4

########################
## Pairwise FST matrix - K-means determined genetic clusters
########################
all(colnames(data@gt)[-1] == pop$sample)
# [1] TRUE
my_genlight <- vcfR2genlight(data)
ploidy(my_genlight) <- 2 ## Although white clover is a tetraploid the subgenomes are distinct enough to be able to call it a diploid
pop(my_genlight) <- pop$pop
my_genlight
## Principal components analysis
vcf.pca <- glPca(my_genlight) ## This takes a while
## Will get asked how many PCs to retain
# choose number of axes to retain:
#   800 ## Can retain as many as you want at this step. 800 looked like 90%

## K-means clustering algorithm. I ran this 10 times and chose the clustering output with the lowest BIC score from the 10 runs. Each time retain 800 PCs.

grp1 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca)
# Choose the number of PCs to retain:
#   800 ## Kept about 90%
# Choose number of clusters:
#   12 ## Choose the lowest BIC value
write.table(grp1$grp, file = "grp1_clustering_K12-1.txt")
## Save a copy of the BIC graph (BIC_K12-1.png) H = 815 W = 815

grp2 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp2$grp, file = "grp2_clustering_K12-2.txt")
## Save a copy of the BIC graph (BIC_K12-2.png) H = 815 W = 815

grp3 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp3$grp, file = "grp3_clustering_K12-3.txt")
## BIC_K12-3 H = 815 W = 815

grp4 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp4$grp, file = "grp4_clustering_K11-1.txt")
## BIC_K11-1 H = 815 W = 815

grp5 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp5$grp, file = "grp5_clustering_K12-4.txt")
## BIC_K12-4 H = 815 W = 815

grp6 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp6$grp, file = "grp6_clustering_K11-2.txt")
## BIC_K11-2 H = 815 W = 815 (Supplementary Figure 7)

grp7 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp7$grp, file = "grp7_clustering_K12-5.txt")
## BIC_K12-5 H = 815 W = 815

grp8 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp8$grp, file = "grp8_clustering_K12-6.txt")
## BIC_K12-6 H = 815 W = 815

grp9 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp9$grp, file = "grp9_clustering_K12-7.txt")
## BIC_K12-7 H = 815 W = 815

grp10 <- find.clusters(my_genlight, max.n.clust = 40, glPca = vcf.pca, n.pca = 800)
write.table(grp10$grp, file = "grp10_clustering_K12-8.txt")
## BIC_K12-8 H = 815 W = 815

## Use grp6
grp6$stat ## K=11 6538.781 (Had lowest BIC value from all the runs)
# grp6 <- read.table("grp6_clustering_K11-2.txt") ## Can read in from saved text file
## To make it compatible with the code used below, change the column name:
# names(grp6)[names(grp6) == 'x'] <- 'grp'

table(pop(my_genlight), grp6$grp) ## See which individuals group into each genetic cluster

## Assign cluster numbers names to each sample:
grp6$grp[grp6$grp == "1"] <- "WNZSL-L"
grp6$grp[grp6$grp == "2"] <- "PARENT"
grp6$grp[grp6$grp == "3"] <- "FNZSL-L"
grp6$grp[grp6$grp == "4"] <- "FNZLL-H"
grp6$grp[grp6$grp == "5"] <- "WUSLL-L"
grp6$grp[grp6$grp == "6"] <- "WUSLL-H"
grp6$grp[grp6$grp == "7"] <- "WNZLL-H"
grp6$grp[grp6$grp == "8"] <- "FNZSL-H"
grp6$grp[grp6$grp == "9"] <- "FNZLL-L"
grp6$grp[grp6$grp == "10"] <- "WNZSL-H"
grp6$grp[grp6$grp == "11"] <- "WNZLL-L"
table(pop(my_genlight), grp6$grp) ## Check everything was assigned correctly

## Created file with both original population and cluster info from above:           
K11_group <- read.table(file = "11_SAMPLE_CLUSTER_POP.txt", header = T)

all(colnames(grp6$grp)[-1] == K11_group$grp)
# [1] TRUE
all(colnames(my_genlight@pop)[-1] == K11_group$pop)
# [1] TRUE

## Order the original populations after importing data
K11_group$pop <- as.factor(K11_group$pop)
levels(K11_group$pop)
# [1] "FNZLL-High-End" "FNZLL-High-Mid" "FNZLL-Low-End"  "FNZLL-Low-Mid"  "FNZLL-Parent"  
# [6] "FNZSL-High-End" "FNZSL-High-Mid" "FNZSL-Low-End"  "FNZSL-Low-Mid"  "FNZSL-Parent"  
#[11] "WNZLL-High-End" "WNZLL-High-Mid" "WNZLL-Low-End"  "WNZLL-Low-Mid"  "WNZLL-Parent"  
#[16] "WNZSL-High-End" "WNZSL-High-Mid" "WNZSL-Low-End"  "WNZSL-Low-Mid"  "WUSLL-High-End"
#[21] "WUSLL-High-Mid" "WUSLL-Low-End"  "WUSLL-Low-Mid"  "WUSLL-Parent"  

K11_group$pop <- factor(K11_group$pop, levels(K11_group$pop)[c(15,11,13,12,14,16,18,17,19,24,20,22,21,23,5,1,3,2,4,10,6,8,7,9)])
levels(K11_group$pop)
# [1] "WNZLL-Parent"   "WNZLL-High-End" "WNZLL-Low-End"  "WNZLL-High-Mid" "WNZLL-Low-Mid" 
# [6] "WNZSL-High-End" "WNZSL-Low-End"  "WNZSL-High-Mid" "WNZSL-Low-Mid"  "WUSLL-Parent"  
#[11] "WUSLL-High-End" "WUSLL-Low-End"  "WUSLL-High-Mid" "WUSLL-Low-Mid"  "FNZLL-Parent"  
#[16] "FNZLL-High-End" "FNZLL-Low-End"  "FNZLL-High-Mid" "FNZLL-Low-Mid"  "FNZSL-Parent"  
#[21] "FNZSL-High-End" "FNZSL-Low-End"  "FNZSL-High-Mid" "FNZSL-Low-Mid" 

## Order the clusters
K11_group$cluster <- as.factor(K11_group$cluster)
levels(K11_group$cluster)
# [1] "FNZLL-H" "FNZLL-L" "FNZSL-H" "FNZSL-L" "PARENT"  "WNZLL-H" "WNZLL-L" "WNZSL-H" "WNZSL-L" "WUSLL-H"
#[11] "WUSLL-L"
K11_group$cluster <- factor(tablevalues$cluster, levels(K11_group$cluster)[c(5,6,7,8,9,10,11,1,2,3,4)])
levels(K11_group$cluster)
# [1] "PARENT"  "WNZLL-H" "WNZLL-L" "WNZSL-H" "WNZSL-L" "WUSLL-H" "WUSLL-L" "FNZLL-H" "FNZLL-L" "FNZSL-H"
#[11] "FNZSL-L"

K11_group <- as.data.frame(K11_group)
head(K11_group)
#  sample cluster           pop
#1  40-02  PARENT  FNZSL-Parent
#2  38-44 WNZSL-L WNZSL-Low-End
#3  38-45 WNZSL-L WNZSL-Low-End
#4  38-46 WNZSL-L WNZSL-Low-End
#5  38-47 WNZSL-L WNZSL-Low-End
#6  40-06  PARENT  FNZSL-Parent
data <- table(K11_group$cluster, K11_group$pop) ## Check table design
print(data)

mylab <-c("PARENT","WNZLL-H","WNZLL-L","WNZSL-H","WNZSL-L","WUSLL-H","WUSLL-L","FNZLL-H","FNZLL-L","FNZSL-H","FNZSL-L")
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

table.value(table(K11_group$pop,K11_group$cluster), col.lab=make.italic(mylab), clabel.col = 0.8, clabel.row = 0.8)
## Save as "inf_11_vs_ori_24_ordered.tiff" (Figure 3)

## K-means/DAPC clustering used as populations for pairwise FST:
pop(my_genind) <- K11_group$cluster
my_genind
my_hierfstat2 <- genind2hierfstat(my_genind) 
genet.dist(my_hierfstat2, method = "WC84") ## This took a while (few hours)

## Used output to create Supplementary Table 5

############################################################################
## Run DAPC using K-means determined clusters
########################
set.seed(999) ## This should make it reproducable
system.time(dapc_cv_K11_overall <- xvalDapc(tab(my_genlight, NA.method = "mean"), K11_group$cluster, training.set = 0.9, result = "overall", n.pca = 1:50, n.rep = 100, parallel = "snow", ncpus = 4L))
# user  system elapsed 
# 197.96  559.50 3441.17 

## Save the DAPC cross-validation as "xvaldapc_K11" W = 815 and H = 815

dapc_cv_K11_overall[-1] ## To get the values for RMSE and MSA
# $`Median and Confidence Interval for Random Chance`
#       2.5%        50%      97.5% 
# 0.07314973 0.09132942 0.10772536 

# $`Mean Successful Assignment by Number of PCs of PCA`													
#	        1	        2	        3	        4	        5	        6	        7	        8     	  9	       10	       11	       12	       13
#	0.4670909	0.7308182	0.8493636	0.9265455	0.9350909	0.9928182	0.9928182	0.9969091	0.9964545	0.9949091	0.9971818	0.9958182	0.9950909
#	       14	       15	       16	       17	       18	       19	       20	       21	       22	       23	       24	       25	       26
#	0.9945455	0.9969091	0.9960909	0.9971818	0.9968182	0.9968182	0.9956364	0.9956364	0.9956364	0.9973636	0.9971818	0.9972727	0.9974545
#	   27	       28	       29	       30	       31	       32	       33	   34	       35	       36	       37	       38	       39
#	0.998	0.9981818	0.9985455	0.9977273	0.9986364	0.9986364	0.9977273	0.999	0.9986364	0.9984545	0.9986364	0.9981818	0.9979091
#	       40	       41	       42	       43	       44	       45	       46	   47	       48	       49	       50		
#	0.9980909	0.9985455	0.9976364	0.9988182	0.9977273	0.9980909	0.9988182	0.999	0.9987273	0.9988182	0.9984545		
# $`Number of PCs Achieving Highest Mean Success`
# [1] "34"

# $`Root Mean Squared Error by Number of PCs of PCA`											
#	          1	          2	          3	          4	          5	          6	          7	          8	          9	         10	         11
#	0.533855461	0.270748885	0.153119086	0.077406302	0.068779201	0.010718024	0.010640636	0.005749596	0.006363636	0.008907235	0.005529784
#	         12	         13	         14	         15	         16	         17	         18	         19	         20	         21	      22
#	0.007496556	0.008430562	0.008528029	0.005749596	0.006741999	0.005677271	0.005961308	0.005961308	0.007158189	0.007385489	0.007606
#	         23	         24	         25	         26	         27	         28	         29	         30	         31	         32	         33
#	0.006363636	0.005378254	0.006030227	0.006030227	0.004453618	0.004453618	0.003636364	0.006098367	0.003748278	0.003748278	0.006363636
#	         34	         35	         36	         37	         38	        39	         40	         41	         42	         43	         44
#	0.003277774	0.004165978	0.004165978	0.003962635	0.005454545	0.00522233	0.005061604	0.004810457	0.006298367	0.003277774	0.005821022
#	         45	         46	         47	         48	         49	         50
#	0.004359847	0.004545455	0.003520894	0.004264014	0.003520894	0.004165978					
# $`Number of PCs Achieving Lowest MSE`
# [1] "43"

## Used output to create Supplementary Figure 8

## Run DAPC using just 6 PCs                                                
set.seed(999)
dapc <- dapc(my_genlight, K11_group$cluster, n.pca = 6, glPca = vcf.pca) 
# Choose the number discriminant functions to retain (>=1): 
#   6

## Plot DAPC scatter                                                
scatter(dapc)
rainbow <- c("#aaffc3", "#000000","#f58231", "#e6beff", "#4363db", "#ffe119","#fabebe","#46F0F0", "#3cb44b","#e6194b", "#911eb4")
scatter(dapc, bg="white", pch=20, cell=0, cstar=0, col=rainbow, solid=.4, cex=3,clab=0, leg=TRUE)

########################
## FINAL PLOT
########################
## Remove the population key and instead have the PCA eigenvalues present in the DAPC scatter of x = 1 and y = 2
scatter(dapc, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=F, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,1], dapc$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,2], pch=4,
       cex=3, lwd=3, col=rainbow)
## Save as DAPC_scatter_K11_FINAL  W= 956 H = 796

## Use the following code for the DA Eigenvalues in the top right corner as cant get crosses in right spot once added in
scatter(dapc, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright", posi.pca = "bottomright",       ## Add in own population labels in powerpoint
        leg=F, posi.leg = "bottomright", inset.solid = 1)   ## Make inset solid
## Saved as DAPC_scatter_K11_FINAL_eigen w = 956  H = 796                                            

## Analyse how much percent of genetic variance is explained by each axis
percent = dapc$eig/sum(dapc$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,40),
        names.arg = round(percent, 1))
## DF1 = 24.6% and DF2 = 19.3%                                                

## Two figures combined to make Figure 4

########################
## Explore other DF combinations
########################  
## Save all dimensions as w = 956  H = 796.
## 1 and 3
scatter(dapc, xax = 1, yax = 3,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,1], dapc$grp.coord[,3], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,3], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_1-3
scatter(dapc, xax = 1, yax = 3,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright")
## K11_1-3-eig

## 1 and 4
scatter(dapc, xax = 1, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,1], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_1-4
scatter(dapc, xax = 1, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright")
## K11_1-4-eig

## 1 and 5
scatter(dapc, xax = 1, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,1], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_1-5
scatter(dapc, xax = 1, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright")
## K11_1-5-eig

## 1 and 6
scatter(dapc, xax = 1, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,1], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_1-6
scatter(dapc, xax = 1, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright")
## K11_1-6-eig

## 2 and 3
scatter(dapc, xax = 2, yax = 3,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,2], dapc$grp.coord[,3], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,2], dapc$grp.coord[,3], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_2-3
scatter(dapc, xax = 2, yax = 3,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_2-3-eig

## 2 and 4
scatter(dapc, xax = 2, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "topleft")
par(xpd=TRUE)
points(dapc$grp.coord[,2], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,2], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_2-4
scatter(dapc, xax = 2, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "topleft") 
## K11_2-4-eig

## 2 and 5
scatter(dapc, xax = 2, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,2], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,2], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_2-5
scatter(dapc, xax = 2, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_2-5-eig

## 2 and 6
scatter(dapc, xax = 2, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,2], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,2], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_2-6
scatter(dapc, xax = 2, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_2-6-eig

## 3 and 4
scatter(dapc, xax = 3, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,3], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,3], dapc$grp.coord[,4], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_3-4
scatter(dapc, xax = 3, yax = 4,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_3-4-eig

## 3 and 5
scatter(dapc, xax = 3, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,3], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,3], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_3-5
scatter(dapc, xax = 3, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_3-5-eig

## 3 and 6
scatter(dapc, xax = 3, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,3], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,3], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_3-6
scatter(dapc, xax = 3, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topleft",
        leg=TRUE, posi.leg = "bottomright") 
## K11_3-6-eig

## 4 and 5
scatter(dapc, xax = 4, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,4], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,4], dapc$grp.coord[,5], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_4-5
scatter(dapc, xax = 4, yax = 5,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topleft",
        leg=TRUE, posi.leg = "bottomright") 
## K11_4-5-eig

## 4 and 6
scatter(dapc, xax = 4, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,4], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,4], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_4-6
scatter(dapc, xax = 4, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_4-6-eig

## 5 and 6
scatter(dapc, xax = 5, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4, clab=0,
        mstree=F, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg = "bottomright")
par(xpd=TRUE)
points(dapc$grp.coord[,5], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=8, col="black")
points(dapc$grp.coord[,5], dapc$grp.coord[,6], pch=4,
       cex=3, lwd=3, col=rainbow)
## K11_5-6
scatter(dapc, xax = 5, yax = 6,
        ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=rainbow, solid=.5, cex=4,
        posi.da="topright",
        leg=TRUE, posi.leg = "bottomright") 
## K11_5-6-eig

############################################################################
## AMOVA in poppr
########################
## Organise strata information
amova_K11_pop <- read.table(file ="AMOVA_11_CLUSTER_POP.txt", header = TRUE)                                                  
head(amova_K11_pop)
pop(my_genind) <- amova_K11_pop$Cluster_Pop                                                
strata(my_genind) <- data.frame(pop(my_genind))                                            
splitStrata(my_genind) <- ~Cluster/Pop 
head(strata(my_genind))

## Remove loci with missing values greater than 5% from genind object

## poppr.amova (which can be used for significance testing) will not use loci and samples that have greater than 5% missing values
## Remove loci with more than 5% missing. As long as there are more than 30 markers, AMOVA will yield acceptable results https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3797491/
my_genind_remove <- missingno(my_genind, type = "loci", cutoff = 0.05, quiet = FALSE, freq = FALSE)   
## 14401 loci contained missing values greater than 5%, therefore only 342 loci will be used.                                              

########################
## Compare within and among population variation at K = 24
########################                                                               
genindamova_na_pop <- poppr.amova(my_genind_remove, ~Pop, within = F) ## 19.3% (between populations) 80.7% (within populations)                                             
genindamova_na_pop
# $call
# ade4::amova(samples = xtab, distances = xdist, structures = xstruct)
#
# $results
#                  Df   Sum Sq   Mean Sq
# Between samples   23 10623.66 461.89806
# Within samples  1089 41486.34  38.09581
# Total           1112 52110.00  46.86151
#
# $componentsofcovariance
#                                Sigma         %
# Variations  Between samples  9.139375  19.34866
# Variations  Within samples  38.095814  80.65134
# Total variations            47.235189 100.00000
#
# $statphi
#                         Phi
# Phi-samples-total 0.1934866

## More variation observed within populations than among populations  

## significance testing for variation among 24 populations:                                                 
set.seed(999)
genindamova_na_pop_signif <- randtest(genindamova_na_pop, nrepet = 9999)
plot(genindamova_na_pop_signif)                                             
genindamova_na_pop_signif                                               
# Monte-Carlo test
# Call: as.randtest(sim = res, obs = sigma[1])
#
# Observation: 9.139375 
#
# Based on 999 replicates
# Simulated p-value: 0.001 
# Alternative hypothesis: greater 
#
#      Std.Obs   Expectation      Variance 
# 3.213351e+02 -1.066866e-03  8.091294e-04  

########################
## Compare within and among population variation at: among K = 11 clusters, among populations (K = 24) within clusters, within clusters
######################## 
genindamova_na <- poppr.amova(my_genind_remove, ~Cluster/Pop, within = F)  ## 15.4% (between cluster), 5.4% (between samples within cluster), 79.2% (within clusters)
genindamova_na
# $call
# ade4::amova(samples = xtab, distances = xdist, structures = xstruct)
# 
# $results
#                                  Df    Sum Sq   Mean Sq
# Between Cluster                  10  8937.738 893.77385
# Between samples Within Cluster   17  2206.825 129.81324
# Within samples                 1085 40965.433  37.75616
# Total                          1112 52109.997  46.86151
# 
# $componentsofcovariance
#                                                Sigma          %
# Variations  Between Cluster                 7.341626  15.402428
# Variations  Between samples Within Cluster  2.567595   5.386708
# Variations  Within samples                 37.756159  79.210864
# Total variations                           47.665380 100.000000
# 
# $statphi
#                           Phi
# Phi-samples-total   0.2078914
# Phi-samples-Cluster 0.0636745
# Phi-Cluster-total   0.1540243

## Do significance testing:
set.seed(999)
genindamova_na_signif <- randtest(genindamova_na, nrepet = 9999)
plot(genindamova_na_signif)                                             
genindamova_na_signif                                             
# class: krandtest lightkrandtest 
# Monte-Carlo tests
# Call: randtest.amova(xtest = genindamova_na, nrepet = 9999)
# 
# Number of tests:   3 
# 
# Adjustment method for multiple comparisons:   none 
# Permutation number:   9999 
# Test       Obs   Std.Obs   Alter Pvalue
# 1  Variations within samples 37.756159 -276.3777    less  1e-04
# 2 Variations between samples  2.567595   40.0929 greater  1e-04
# 3 Variations between Cluster  7.341626   13.2417 greater  1e-04                                               

## Significant for all measures  

## Used AMOVA output above for Supplementary Table 6                                             

############################################################################
## AMOVA in pegas
########################

########################
## Compare within and among population variation at K = 24
########################                                             

## Import population information    
amova_24pop <- read.table(file = "AMOVA_24_POP.txt", header = TRUE) ## This is just "24_SAMPLE_POP.txt" but without the sample column   

## Create copy of my_genlight
my_genlight2 <- my_genlight  

## Assign strata to genlight                                               
strata(my_genlight2) <- amova_24pop
head(strata(my_genlight2))
table(strata(my_genlight2))

## Create euclidean distance matrix                                                
K24_dist <- dist(my_genlight2)
K24_stra <- strata(my_genlight2)

## Run amova using the 24 population structure                                                
K24amova <- pegas::amova(K24_dist ~Pop, data = K24_stra)
K24amova
# Analysis of Molecular Variance
#
# Call: pegas::amova(formula = K24_dist ~ Pop, data = K24_stra)
#
#            SSD       MSD   df
# Pop    981378.4 42668.626   23
# Error 3075147.2  2823.827 1089
# Total 4056525.6  3647.955 1112
#
# Variance components:
#        sigma2 P.value
# Pop    859.26       0
# Error 2823.83        
#
# Phi-statistics:
# Pop.in.GLOBAL 
#      0.233299 
#
# Variance coefficients:
#        a 
# 46.37103
## % is calculated by adding the two sigma values together (859.26+2823.83) (equals: 3683.09)
## then for among populations % : (859.26/3683.09)*100 equals 23.3
## and within populations % : (2823.83/3683.09)*100 equals 76.7
## More variation observed within populations than among populations                                                  

########################
## Compare within and among population variation at: among K = 11 clusters, among populations (K = 24) within clusters, within clusters
########################                                                
amova_K11_pop <- read.table(file ="AMOVA_11_CLUSTER_POP.txt", header = TRUE)
my_genlight3 <- my_genlight
strata(my_genlight3) <- amova_K11_pop
head(strata(my_genlight3))
table(strata(my_genlight3))
splitStrata(my_genlight3) <- ~Cluster/Pop
head(strata(my_genlight3))
table(strata(my_genlight3))
K11_dist <- dist(my_genlight3)
K11_stra <- strata(my_genlight3)
K11amova <- pegas::amova(K11_dist ~Cluster/Pop, data = K11_stra)
K11amova
# Analysis of Molecular Variance
# 
# Call: pegas::amova(formula = K11_dist ~ Cluster/Pop, data = K11_stra)
# 
#               SSD       MSD   df
# Cluster  797180.6 79718.063   10
# Pop      184197.8 14169.059   13
# Error   3075147.2  2823.827 1089
# Total   4056525.6  3647.955 1112
# 
# Variance components:
#          sigma2 P.value
# Cluster  669.39        
# Pop      245.32       0
# Error   2823.83        
# 
# Phi-statistics:
#   Cluster.in.GLOBAL (Phi_CT)     Pop.in.GLOBAL (Phi_ST)    Pop.in.Cluster (Phi_SC) 
# 0.17905051                 0.24466866                 0.07992959 
# 
# Variance coefficients:
#   a        b        c 
# 46.24756 46.53154 97.82013

## % is calculated by adding the three sigma values together (669.39+245.32+2823.83) (equals: 3738.54)
## then for among clusters % : (669.39/3738.54)*100 equals 17.9
## then for among populations % : (245.32/3738.54)*100 equals 6.6                                                       
## and within clusters % : (2823.83/3738.54)*100 equals 75.5

## Both AMOVA results are similar. pegas::amova can handle missing data but poppr.amova will not use loci and samples that have greater than 5% missing values