############################################################################
## PCAdapt outlier detection
########################

## Use vcf files that were created for DAPC and BayeScan and KGD analyses

## Convert into .PLINK format first:
# /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Pcadapt/
# conda activate plink
# plink --vcf /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WNZLL_samples_DAPC_grouping.vcf --maf 0.03 --make-bed --out WNZLL_samples_DAPC_grouping --indiv-sort a
# plink --vcf /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WNZSL_samples_DAPC_grouping.vcf --maf 0.03 --make-bed --out WNZSL_samples_DAPC_grouping --indiv-sort a
# plink --vcf /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WUSLL_samples_DAPC_grouping.vcf --maf 0.03 --make-bed --out WUSLL_samples_DAPC_grouping --indiv-sort a
# plink --vcf /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/FNZLL_samples_DAPC_grouping.vcf --maf 0.03 --make-bed --out FNZLL_samples_DAPC_grouping --indiv-sort a
# plink --vcf /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/FNZSL_samples_DAPC_grouping.vcf --maf 0.03 --make-bed --out FNZSL_samples_DAPC_grouping --indiv-sort a

## NUMBER OF SNPS LEFT AFTER MAF 0.03, PLINK CONVERSION AND REMOVING CHR 17 SNPs:
## WNZLL 11,061 AND 188 SAMPLES
## WNZSL 11,479 AND 186 SAMPLES
## WUSLL 11,171 AND 195 SAMPLES
## FNZLL 10,979 AND 182 SAMPLES
## FNZSL 10,976 AND 184 SAMPLES

library(pcadapt)
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Pcadapt")

## IMPORT PLINK FILES AND RUN PCADAPT:
WNZLL <- read.pcadapt("WNZLL_samples_DAPC_grouping.bed", type = "bed", type.out = c("bed", "matrix"))
WLL <- pcadapt(input = WNZLL, K = 20, min.maf = 0.03)

WNZSL <- read.pcadapt("WNZSL_samples_DAPC_grouping.bed", type = "bed", type.out = c("bed", "matrix"))
WSL <- pcadapt(input = WNZSL, K = 20, min.maf = 0.03)

WUSLL <- read.pcadapt("WUSLL_samples_DAPC_grouping.bed", type = "bed", type.out = c("bed", "matrix"))
WUS <- pcadapt(input = WUSLL, K = 20, min.maf = 0.03)

FNZLL <- read.pcadapt("FNZLL_samples_DAPC_grouping.bed", type = "bed", type.out = c("bed", "matrix"))
FLL <- pcadapt(input = FNZLL, K = 20, min.maf = 0.03)

FNZSL <- read.pcadapt("FNZSL_samples_DAPC_grouping.bed", type = "bed", type.out = c("bed", "matrix"))
FSL <- pcadapt(input = FNZSL, K = 20, min.maf = 0.03)

########################
## changed screeplot look by increasing point size, making background white and adding in titles
## see "scree_plotSP" functions below
########################

## dependents:
library(ggplot2)
library(plotly)
library(ggthemes)
library(ggpubr)
library(gridExtra)
#' Principal Components Analysis Scree Plot
#'
#' \code{scree_plot} plots the scee plot associated with the principal components analysis performed on the dataset.
#' NB : \code{pcadapt} has to be run on the dataset in order to get an output readable by \code{plot.screePlot}
#'
#' @param x an output from \code{pcadapt} containing the singular values.
#' @param K an integer specifying the number of components to take into account in the scree plot.
#'
#' @examples
#' ## see ?fastpcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot geom_line guides ggtitle
#' @importFrom plotly plot_ly layout
#' 
#' @export
#'
scree_plotSP = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         colour = I("#000000"), 
                         size = I(4),
                         xlab = "Principal Component", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line(size = 0.8, colour = I("#000000")) + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(expression(bold(A)~~"WNZLL"))+
      theme_bw() +
      theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15), plot.title = element_text(size = 20))
    print(p0)
  }
}
scree_plotSP(WLL, 20)
WLL_scree <- scree_plotSP(WLL, 20)

scree_plotSP = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         colour = I("#000000"), 
                         size = I(4),
                         xlab = "Principal Component", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line(size = 0.8, colour = I("#000000")) + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(expression(bold(B)~~"WNZSL"))+
      theme_bw() +
      theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15), plot.title = element_text(size = 20))
    print(p0)
  }
}
scree_plotSP(WSL, 20)
WSL_scree <- scree_plotSP(WSL, 20)

scree_plotSP = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         colour = I("#000000"), 
                         size = I(4),
                         xlab = "Principal Component", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line(size = 0.8, colour = I("#000000")) + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(expression(bold(C)~~"WUSLL"))+
      theme_bw() +
      theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15), plot.title = element_text(size = 20))
    print(p0)
  }
}
scree_plotSP(WUS, 20)
WUS_scree <- scree_plotSP(WUS, 20)

scree_plotSP = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         colour = I("#000000"), 
                         size = I(4),
                         xlab = "Principal Component", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line(size = 0.8, colour = I("#000000")) + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(expression(bold(D)~~"FNZLL"))+
      theme_bw() +
      theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15), plot.title = element_text(size = 20))
    print(p0)
  }
}
scree_plotSP(FLL, 20)
FLL_scree <- scree_plotSP(FLL, 20)

scree_plotSP = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         colour = I("#000000"), 
                         size = I(4),
                         xlab = "Principal Component", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line(size = 0.8, colour = I("#000000")) + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(expression(bold(E)~~"FNZSL"))+
      theme_bw() +
      theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15), plot.title = element_text(size = 20))
    print(p0)
  }
}
scree_plotSP(FSL, 20)
FSL_scree <- scree_plotSP(FSL, 20)

ggarrange(WLL_scree,FLL_scree,WSL_scree,FSL_scree,WUS_scree, ncol = 2, nrow = 3, font.label = list(family = "sans"))
## Save as PDF A4 portrait and 840 x 1185 TIFF (Supplementary Figure 2)

########################
## create popfile for EACH pool by using the first column of the .fam file "WNZLL_samples_dp5to150_mm0.2_no_multiallelic.fam" 
########################
## READ IN EACH POOL'S POPFILE:

WLLpop <- read.table("WLLpopfile2.txt", header = FALSE)
WLLpop <- as.character(WLLpop$V1)
WLLpop <- as.factor(WLLpop)
levels(WLLpop)
# [1] "WNZLL-High-End" "WNZLL-High-Mid" "WNZLL-Low-End"  "WNZLL-Low-Mid" 
WLLpop <- factor(WLLpop,
                 levels(WLLpop)[c(1,2,4,3)])
levels(WLLpop)
# [1] "WNZLL-High-End" "WNZLL-High-Mid" "WNZLL-Low-Mid"  "WNZLL-Low-End" 

WSLpop <- read.table("WSLpopfile2.txt", header = FALSE)
WSLpop <- as.character(WSLpop$V1)
WSLpop <- as.factor(WSLpop)
levels(WSLpop)
# [1] "WNZSL-High-End" "WNZSL-High-Mid" "WNZSL-Low-End"  "WNZSL-Low-Mid" 
WSLpop <- factor(WSLpop,
                 levels(WSLpop)[c(1,2,4,3)])
levels(WSLpop)
# [1] "WNZSL-High-End" "WNZSL-High-Mid" "WNZSL-Low-Mid"  "WNZSL-Low-End" 

WUSpop <- read.table("WUSpopfile2.txt", header = FALSE)
WUSpop <- as.character(WUSpop$V1)
WUSpop <- as.factor(WUSpop)
levels(WUSpop)
# [1] "WUSLL-High-End" "WUSLL-High-Mid" "WUSLL-Low-End"  "WUSLL-Low-Mid" 
WUSpop <- factor(WUSpop,
                 levels(WUSpop)[c(1,2,4,3)])
levels(WUSpop)
# [1] "WUSLL-High-End" "WUSLL-High-Mid" "WUSLL-Low-Mid"  "WUSLL-Low-End" 

FLLpop <- read.table("FLLpopfile2.txt", header = FALSE)
FLLpop <- as.character(FLLpop$V1)
FLLpop <- as.factor(FLLpop)
levels(FLLpop)
# [1] "FNZLL-High-End" "FNZLL-High-Mid" "FNZLL-Low-End"  "FNZLL-Low-Mid" 
FLLpop <- factor(FLLpop,
                 levels(FLLpop)[c(1,2,4,3)])
levels(FLLpop)
# [1] "FNZLL-High-End" "FNZLL-High-Mid" "FNZLL-Low-Mid"  "FNZLL-Low-End" 

FSLpop <- read.table("FSLpopfile2.txt", header = FALSE)
FSLpop <- as.character(FSLpop$V1)
FSLpop <- as.factor(FSLpop)
levels(FSLpop)
# [1] "FNZSL-High-End" "FNZSL-High-Mid" "FNZSL-Low-End"  "FNZSL-Low-Mid" 
FSLpop <- factor(FSLpop,
                 levels(FSLpop)[c(1,2,4,3)])
levels(FSLpop)
# [1] "FNZSL-High-End" "FNZSL-High-Mid" "FNZSL-Low-Mid"  "FNZSL-Low-End" 

########################
## CREATE CUSTOM SCORE PLOTS:
########################
## make five different colour palettes:
cbpalette1 <- c("#dfc27d", "#a6611a", "#018571", "#80cdc1")
cbpalette2 <- c("#f1b6da", "#d01c8b", "#008837", "#a6dba0")
cbpalette3 <- c("#fdb863", "#e66101", "#5e3c99", "#b2abd2")
cbpalette4 <- c("#f4a582", "#ca0020", "#0571b0", "#92c5de")
cbpalette5 <- c("#c2a5cf", "#7b3294", "#4dac26", "#b8e186")

## TO VIEW THE PC% VALUES:
head(WLL$singular.values)
# 0.3497738 0.1769604 = 34.98% for PC1 and 17.69% for PC2
head(WSL$singular.values)
# 0.3107311 0.1799233 = 31.07 and 17.99
head(WUS$singular.values)
# 0.3219666 0.2101511 = 32.20 and 21.02
head(FLL$singular.values)
# 0.3592923 0.1998444 = 35.93 and 19.98
head(FSL$singular.values)
# 0.3517550 0.2210601 = 35.18 and 22.11

## dependents:
library(ggplot2)
# install.packages("magrittr")
library(magrittr)
# install.packages("ggthemes")
library(ggthemes)

## need to manually write in PC% values for the PC axes (line 47)
## PC for each pool
## WLL:   PC 1 = 35.0% and PC 2 = 17.7% and PC 3 = 17.4%
## WSL:	PC 1 = 31.1% and PC 2 = 18.0% and PC 3 = 18.0%
## WUS:	PC 1 = 32.2% and PC 2 = 21.0% and PC 3 = 18.2%
## FLL:	PC 1 = 35.9% and PC 2 = 20.0% and PC 3 = 18.9%
## FSL:	PC 1 = 35.2% and PC 2 = 22.1% and PC 3 = 18.2%
## and change the Pool name:

score_plotSP = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j], .name_repair = "minimal")    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot(df, aes_string("PC_i", "PC_j")) + 
      geom_point() + 
      ggtitle(expression(bold(A)~~"WNZLL")) +
      labs(x = paste0("PC", i, " (35.0%)"), y = paste0("PC", j, " (17.7%)"))+
      theme_bw() +
      theme(legend.text = element_text(size = 10), axis.text=element_text(size = 10), 
            axis.title=element_text(size = 15), plot.title = element_text(size = 20)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))
    
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(colour = factor(df$Pop)), size = 3)
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            colors = col,
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}
WLL_score <- score_plotSP(WLL, i = 1, j = 2, pop = WLLpop, col = cbpalette1)

score_plotSP = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j], .name_repair = "minimal")    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot(df, aes_string("PC_i", "PC_j")) + 
      geom_point() + 
      ggtitle(expression(bold(B)~~"WNZSL")) +
      labs(x = paste0("PC", i, " (31.1%)"), y = paste0("PC", j, " (18.0%)"))+
      theme_bw() +
      theme(legend.text = element_text(size = 10), axis.text=element_text(size = 10), 
            axis.title=element_text(size = 15), plot.title = element_text(size = 20)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))
    
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(colour = factor(df$Pop)), size = 3)
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            colors = col,
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}
WSL_score <- score_plotSP(WSL, i = 1, j = 2, pop = WSLpop, col = cbpalette2)

score_plotSP = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j], .name_repair = "minimal")    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot(df, aes_string("PC_i", "PC_j")) + 
      geom_point() + 
      ggtitle(expression(bold(C)~~"WUSLL")) +
      labs(x = paste0("PC", i, " (32.2%)"), y = paste0("PC", j, " (21.0%)"))+
      theme_bw() +
      theme(legend.text = element_text(size = 10), axis.text=element_text(size = 10), 
            axis.title=element_text(size = 15), plot.title = element_text(size = 20)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))
    
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(colour = factor(df$Pop)), size = 3)
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            colors = col,
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}
WUS_score <- score_plotSP(WUS, i = 1, j = 2, pop = WUSpop, col = cbpalette3)

score_plotSP = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j], .name_repair = "minimal")    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot(df, aes_string("PC_i", "PC_j")) + 
      geom_point() + 
      ggtitle(expression(bold(D)~~"FNZLL")) +
      labs(x = paste0("PC", i, " (20.0%)"), y = paste0("PC", j, " (18.9%)"))+
      theme_bw() +
      theme(legend.text = element_text(size = 10), axis.text=element_text(size = 10), 
            axis.title=element_text(size = 15), plot.title = element_text(size = 20)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))
    
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(colour = factor(df$Pop)), size = 3)
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            colors = col,
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}
FLL_score <- score_plotSP(FLL, i = 1, j = 2, pop = FLLpop, col = cbpalette4)

score_plotSP = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j], .name_repair = "minimal")    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot(df, aes_string("PC_i", "PC_j")) + 
      geom_point() + 
      ggtitle(expression(bold(E)~~"FNZSL")) +
      labs(x = paste0("PC", i, " (35.2%)"), y = paste0("PC", j, " (22.1%)"))+
      theme_bw() +
      theme(legend.text = element_text(size = 10), axis.text=element_text(size = 10), 
            axis.title=element_text(size = 15), plot.title = element_text(size = 20)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))
    
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(colour = factor(df$Pop)), size = 3)
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            colors = col,
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}
FSL_score <- score_plotSP(FSL, i = 1, j = 2, pop = FSLpop, col = cbpalette5)

legend_a <- get_legend(WLL_score, position = NULL)
legend_b <- get_legend(WSL_score, position = NULL)
legend_c <- get_legend(WUS_score, position = NULL)
legend_d <- get_legend(FLL_score, position = NULL)
legend_e <- get_legend(FSL_score, position = NULL)

a <- as_ggplot(legend_a)
b <- as_ggplot(legend_b)
c <- as_ggplot(legend_c)
d <- as_ggplot(legend_d)
e <- as_ggplot(legend_e)

key <- ggarrange(a,d,b,e,c, ncol = 2, nrow = 3)

ggarrange(WLL_score, FLL_score, WSL_score, FSL_score,WUS_score,key, ncol = 2, nrow = 3, common.legend = F, legend = "none")
## Save as PDF A4 portrait and 840 x 1185 TIFF (Figure 5)

## Assess the score plots for PC 2 and 3 (change the PC % values in code above before running)
# WLL_2_3 <- score_plotSP(WLL, i = 2, j = 3, pop = WLLpop, col = cbpalette1)
# WSL_2_3 <- score_plotSP(WSL, i = 2, j = 3, pop = WSLpop, col = cbpalette2)
# WUS_2_3 <- score_plotSP(WUS, i = 2, j = 3, pop = WUSpop, col = cbpalette3)
# FLL_2_3 <- score_plotSP(FLL, i = 2, j = 3, pop = FLLpop, col = cbpalette4)
# FSL_2_3 <- score_plotSP(FSL, i = 2, j = 3, pop = FSLpop, col = cbpalette5)
# 
# legend_a <- get_legend(WLL_2_3, position = NULL)
# legend_b <- get_legend(WSL_2_3, position = NULL)
# legend_c <- get_legend(WUS_2_3, position = NULL)
# legend_d <- get_legend(FLL_2_3, position = NULL)
# legend_e <- get_legend(FSL_2_3, position = NULL)
# 
# a <- as_ggplot(legend_a)
# b <- as_ggplot(legend_b)
# c <- as_ggplot(legend_c)
# d <- as_ggplot(legend_d)
# e <- as_ggplot(legend_e)
# 
# key <- ggarrange(a,d,b,e,c, ncol = 2, nrow = 3)
# 
# ggarrange(WLL_2_3, FLL_2_3, WSL_2_3, FSL_2_3,WUS_2_3,key, ncol = 2, nrow = 3, common.legend = F, legend = "none")
## Save as PDF A4 portrait and 840 x 1100 TIFF (Supplementary Figure 9)

########################
## FIND SNPS THAT ARE DRIVING THE SEPARATION OF POPULATIONS BASED ON PC 1
########################
## Because the first axis explain the majority of the variation (scree plot K = 1) run PCAdapt again, this time using a K value of 1
WLLK1 <- pcadapt(input = WNZLL, K = 1, min.maf = 0.03)
WSLK1 <- pcadapt(input = WNZSL, K = 1, min.maf = 0.03)
WUSK1 <- pcadapt(input = WUSLL, K = 1, min.maf = 0.03)
FLLK1 <- pcadapt(input = FNZLL, K = 1, min.maf = 0.03)
FSLK1 <- pcadapt(input = FNZSL, K = 1, min.maf = 0.03)

write.table(WLLK1$pvalues, file = "WLL_p_values.txt",sep="\t")
write.table(WSLK1$pvalues, file = "WSL_p_values.txt",sep="\t")
write.table(WUSK1$pvalues, file = "WUS_p_values.txt",sep="\t")
write.table(FLLK1$pvalues, file = "FLL_p_values.txt",sep="\t")
write.table(FSLK1$pvalues, file = "FSL_p_values.txt",sep="\t")

## CREATE MANHATTAN INPUT FILES WITH FOLLOWING HEADINGS: "CHR" "SNP" "BP" and "P" + REMOVE ANY SNPS FROM CHR 17
## FILES CALLED "[POOL]_MAN_INPUT.txt"

## IMPORT MANHATTAN INPUT FILES:
WLLK1_MAN <- read.table("WLL_MAN_INPUT.txt", sep="\t",header=TRUE)
WSLK1_MAN <- read.table("WSL_MAN_INPUT.txt", sep="\t",header=TRUE)
WUSK1_MAN <- read.table("WUS_MAN_INPUT.txt", sep="\t",header=TRUE)
FLLK1_MAN <- read.table("FLL_MAN_INPUT.txt", sep="\t",header=TRUE)
FSLK1_MAN <- read.table("FSL_MAN_INPUT.txt", sep="\t",header=TRUE)

## IDENTIFY OUTLIERS using Bonferroni correction at an alpha of 5%:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("qvalue")
library(qvalue)

qvalWLLK1 <- qvalue(WLLK1$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qvalWLLK1 < alpha)
length(outliers)
# 963
## Bonferroni correction
padj <- p.adjust(WLLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 165
outliers_WLLK1 <- get.pc(WLLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WLLK1,file="outliers_WLLK1a-0.05.txt",sep="\t")

qvalWSLK1 <- qvalue(WSLK1$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qvalWSLK1 < alpha)
length(outliers)
# 537
## Bonferroni correction
padj <- p.adjust(WSLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 94
outliers_WSLK1 <- get.pc(WSLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WSLK1,file="outliers_WSLK1a-0.05.txt",sep="\t")

qvalWUSK1 <- qvalue(WUSK1$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qvalWUSK1 < alpha)
length(outliers)
# 588
## Bonferroni correction
padj <- p.adjust(WUSK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 119
outliers_WUSK1 <- get.pc(WUSK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WUSK1,file="outliers_WUSK1a-0.05.txt",sep="\t")

qvalFLLK1 <- qvalue(FLLK1$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qvalFLLK1 < alpha)
length(outliers)
# 955
## Bonferroni correction
padj <- p.adjust(FLLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 203
outliers_FLLK1 <- get.pc(FLLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_FLLK1,file="outliers_FLLK1a-0.05.txt",sep="\t")

qvalFSLK1 <- qvalue(FSLK1$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qvalFSLK1 < alpha)
length(outliers)
# 562
## Bonferroni correction
padj <- p.adjust(FSLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 99
outliers_FSLK1 <- get.pc(FSLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_FSLK1,file="outliers_FSLK1a-0.05.txt",sep="\t")

########################
## IDENTIFY OUTLIERS using Bonferroni correction at an alpha of 1%:
########################

qvalWLLK1 <- qvalue(WLLK1$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qvalWLLK1 < alpha)
length(outliers)
# 527
## Bonferroni correction
padj <- p.adjust(WLLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 117
outliers_WLLK1 <- get.pc(WLLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WLLK1,file="outliers_WLLK1a-0.01.txt",sep="\t")

qvalWSLK1 <- qvalue(WSLK1$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qvalWSLK1 < alpha)
length(outliers)
# 261
## Bonferroni correction
padj <- p.adjust(WSLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 69
outliers_WSLK1 <- get.pc(WSLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WSLK1,file="outliers_WSLK1a-0.01.txt",sep="\t")

qvalWUSK1 <- qvalue(WUSK1$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qvalWUSK1 < alpha)
length(outliers)
# 281
## Bonferroni correction
padj <- p.adjust(WUSK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 81
outliers_WUSK1 <- get.pc(WUSK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_WUSK1,file="outliers_WUSK1a-0.01.txt",sep="\t")

qvalFLLK1 <- qvalue(FLLK1$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qvalFLLK1 < alpha)
length(outliers)
# 553
## Bonferroni correction
padj <- p.adjust(FLLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 145
outliers_FLLK1 <- get.pc(FLLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_FLLK1,file="outliers_FLLK1a-0.01.txt",sep="\t")

qvalFSLK1 <- qvalue(FSLK1$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qvalFSLK1 < alpha)
length(outliers)
# 294
## Bonferroni correction
padj <- p.adjust(FSLK1$pvalues, method = "bonferroni")
outliersB <- which(padj < alpha)
length(outliersB)
# 59
outliers_FSLK1 <- get.pc(FSLK1, outliersB) #This tells you the SNP and the PC (match the SNPs by their order)
write.table(outliers_FSLK1,file="outliers_FSLK1a-0.01.txt",sep="\t")

########################
## Manhattan plots:
########################
library(qqman)

## The cut off mark proposed by qqman:
## Genomewideline is set to -log10(5e-8), which is about 7.30 on the -log10 scale. 
## The p=5e-8 number comes from the notion that there are ~1 million independent tests in the genome, 
## so by a Bonferroni correction, 0.05/1e6=5e-8, 
## the widely accepted holy grail p-value for "genome-wide significance."

## Find the Bonferroni cut off thresholds for each pool (P-values):
## FDR THRESHOLD [0.05 OR 0.01]/ NUMBER OF SNPS USED FOR EACH POOL

# 0.05

# WLL = 4.52039E-06
# FLL = 4.55415E-06
# WSL = 4.35578E-06
# FSL = 4.55539E-06
# WUS = 4.47588E-06

# 0.01

# WLL = 9.04077E-07
# FLL = 9.1083E-07
# WSL = 8.71156E-07
# FSL = 9.11079E-07
# WUS = 8.95175E-07

## IMPORT THE OUTLIER SNPS FROM EACH POOL (IN COMMON BETWEEN MULTIPLE POOLS AND P-VALUE ABOVE FDR = 0.05)

WLL_outliers <- c("1_34133860","1_50382136","2_47264284","3_28211144","4_39100380","5_29652780","7_43806856","12_33367740","13_4850703","14_125557","15_6515376","15_6515380","16_35905618")
WSL_outliers <- c("1_34133860","2_47264284","4_58491291","5_29652780","9_36683987","9_42641599","11_9345053","11_9345077","12_16032423","15_6515376","15_6515380","15_7508099","15_7508130","15_7508139")
WUS_outliers <- c("1_50382136","2_6673787","4_3910146","4_9733285","4_29835602","4_36148372","4_83057071","5_35660916","6_31429353","6_31429365","9_36683987","9_42641599","12_33367740","15_7508099","15_7508130","15_7508139","16_760795","16_35905618")
FLL_outliers <- c("2_6673787","3_28211144","4_3910146","4_9733285","4_29835602","4_83057071","6_31429353","6_31429365","7_5692751","7_43806856","11_9345053","11_9345077","12_3437942","12_16032423","12_17533776","14_125557","15_17315253","16_760795")
FSL_outliers <- c("1_50382136","4_36148372","4_39100380","4_58491291","5_35660916","7_5692751","12_3437942","12_17533776","13_4850703","15_17315253")

## Used output above to create Supplementary Table 7

## RUN THE FOLLOWING MANHATTAN FUNCTION: 
## NOTE: NEED TO CHANGE THE SUGGESTIVE LINE AND GENOMEWIDE LINE VALUES FOR EACH POOL (SUGGESTIVE = 0.05 AND GENOMEWIDE = 0.01)
library(calibrate)
manhattan_p_value <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                              col=c("#e6194b", "#fabebe", "#f58231","#ffe119", "#aaffc3", "#3cb44b", "#46f0f0", "#4363db"), chrlabs=NULL,
                              suggestiveline=-log10(4.47588E-06), genomewideline=-log10(8.95175E-07), 
                              highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  
  # Not sure why, but package check will warn without this.
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

par(mfrow=c(3,2))

manhattan_p_value(WLLK1_MAN[is.na(WLLK1_MAN$P)==FALSE,],annotatePval = FALSE, annotateTop = FALSE,highlight = WLL_outliers, ylim =c(0,40))
title(expression(bold(A)~~"WNZLL"), adj = 0, cex.main = 2)
## CHANGE LINE VALUES FOR EACH POOL

manhattan_p_value(FLLK1_MAN[is.na(FLLK1_MAN$P)==FALSE,],annotatePval = FALSE, annotateTop = FALSE,highlight = FLL_outliers, ylim =c(0,40))
title(expression(bold(D)~~"FNZLL"), adj = 0, cex.main = 2)

manhattan_p_value(WSLK1_MAN[is.na(WSLK1_MAN$P)==FALSE,],annotatePval = FALSE, annotateTop = FALSE,highlight = WSL_outliers, ylim =c(0,40))
title(expression(bold(B)~~"WNZSL"), adj = 0, cex.main = 2)

manhattan_p_value(FSLK1_MAN[is.na(FSLK1_MAN$P)==FALSE,],annotatePval = FALSE, annotateTop = FALSE,highlight = FSL_outliers, ylim =c(0,40))
title(expression(bold(E)~~"FNZSL"), adj = 0, cex.main = 2)

manhattan_p_value(WUSK1_MAN[is.na(WUSK1_MAN$P)==FALSE,],annotatePval = FALSE, annotateTop = FALSE,highlight = WUS_outliers, ylim =c(0,40))
title(expression(bold(C)~~"WUSLL"), adj = 0, cex.main = 2)
## save image width: 1600 and height:850
## P_VALUE_MANHATTAN_PLOTS_PCADAPT.tiff

########################
## LIST OF OUTLIER SNPS AT FDR = 0.05 ARE STORED IN "FDR_0.05_COMMON_UPSET" SHEET
## These SNPs are outliers ID at FDR = 0.05 threshold
## NAME OF WHICH POOL AND POINT OF SELECTION CYCLE AS HEADER
## Save data as csv AND IMPORT INTO R TO CREATE UPSET PLOT
########################
library(readr)
FDR_0.05_COMMON_UPSET <- read_csv("FDR_0.05_COMMON_UPSET.csv")
View(FDR_0.05_COMMON_UPSET)

ListInput <- list(WNZLL = FDR_0.05_COMMON_UPSET$WNZLL[1:165],
                  WNZSL = FDR_0.05_COMMON_UPSET$WNZSL[1:94],
                  WUSLL = FDR_0.05_COMMON_UPSET$WUSLL[1:119],
                  FNZLL = FDR_0.05_COMMON_UPSET$FNZLL[1:203],
                  FNZSL = FDR_0.05_COMMON_UPSET$FNZSL[1:99])
library(UpSetR)
upset(fromList(ListInput), order.by = "freq", nsets = 5, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per pool", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2, 2), 
      set_size.show = TRUE, set_size.scale_max = 250)
## save image width: 1000 and height: 700 
"Upset_WSC_pool_PCAdapt_outliers_FDR_5.png"

############################################################################
## BayeScan analysis
########################
## install bayescan on hpc:
# pwd
# dataset/GBS_WClover_WSC/archive/Sofie/Bayescan
## also in:
# home/pearsons/conda-envs/
## Download zip:
# wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
## Extract files:
# unzip BayeScan2.1.zip
## Find program:
# cd BayeScan2.1/binaries/ (dataset/GBS_WClover_WSC/archive/Sofie/Bayescan/BayeScan2.1/binaries)
# cd BayeScan2.1/binaries/ (home/pearsons/conda-envs/BayeScan2.1/binaries)
## Allow execution of programs:
# chmod u+x BayeScan2.1_linux64bits
## Run program:
# ./BayeScan2.1_linux64bits [name].bsc
## eg:
# ./BayeScan2.1_linux64bits WNZLL_C.bsc

## Run BayeScan analysis using K = 11 clusters as determined by DAPC analysis (exclude parental K = 10)
## Files are stored in hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan
# conda activate vcftools
# cd dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/
# vcftools --remove-indels --keep WNZLL_samples_DAPC_grouping.txt --recode --recode-INFO-all --vcf ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > WNZLL_samples_DAPC_grouping.vcf
# vcftools --remove-indels --keep WNZSL_samples_DAPC_grouping.txt --recode --recode-INFO-all --vcf ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > WNZSL_samples_DAPC_grouping.vcf
# vcftools --remove-indels --keep WUSLL_samples_DAPC_grouping.txt --recode --recode-INFO-all --vcf ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > WUSLL_samples_DAPC_grouping.vcf
# vcftools --remove-indels --keep FNZLL_samples_DAPC_grouping.txt --recode --recode-INFO-all --vcf ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > FNZLL_samples_DAPC_grouping.vcf
# vcftools --remove-indels --keep FNZSL_samples_DAPC_grouping.txt --recode --recode-INFO-all --vcf ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf --stdout > FNZSL_samples_DAPC_grouping.vcf

## Create Bayescan input files
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan")

# library(devtools)
# devtools::install_github("rstudio/httpuv")
library(httpuv)
library(adegenet)
library(hierfstat)
library(pegas)
library(vcfR)
library(ape)
library(RColorBrewer)
library(poppr)
library(plyr)

## vcf files with populations I want to compare:
WNZLL_vcf <- read.vcfR("WNZLL_samples_DAPC_grouping.vcf")
#188 SAMPLES, 14743 SNPs
WNZSL_vcf <- read.vcfR("WNZSL_samples_DAPC_grouping.vcf")
#186 SAMPLES, 14743 SNPs
WUSLL_vcf <- read.vcfR("WUSLL_samples_DAPC_grouping.vcf")
#195 SAMPLES, 14743 SNPs
FNZLL_vcf <- read.vcfR("FNZLL_samples_DAPC_grouping.vcf")
#182 SAMPLES, 14743 SNPs
FNZSL_vcf <- read.vcfR("FNZSL_samples_DAPC_grouping.vcf")
#184 SAMPLES, 14743 SNPs

########################
## Find order of samples in each vcf file and assign population info to create popfiles
########################
WNZLL_vcfsamples=gsub("-.$","",colnames(WNZLL_vcf@gt),fixed = FALSE)
WNZLL_vcfsamples=WNZLL_vcfsamples[-1]
write.table(WNZLL_vcfsamples, file = "WNZLL_samples.txt", sep = "\t")

WNZSL_vcfsamples=gsub("-.$","",colnames(WNZSL_vcf@gt),fixed = FALSE)
WNZSL_vcfsamples=WNZSL_vcfsamples[-1]
write.table(WNZSL_vcfsamples, file = "WNZSL_samples.txt", sep = "\t")

WUSLL_vcfsamples=gsub("-.$","",colnames(WUSLL_vcf@gt),fixed = FALSE)
WUSLL_vcfsamples=WUSLL_vcfsamples[-1]
write.table(WUSLL_vcfsamples, file = "WUSLL_samples.txt", sep = "\t")

FNZLL_vcfsamples=gsub("-.$","",colnames(FNZLL_vcf@gt),fixed = FALSE)
FNZLL_vcfsamples=FNZLL_vcfsamples[-1]
write.table(FNZLL_vcfsamples, file = "FNZLL_samples.txt", sep = "\t")

FNZSL_vcfsamples=gsub("-.$","",colnames(FNZSL_vcf@gt),fixed = FALSE)
FNZSL_vcfsamples=FNZSL_vcfsamples[-1]
write.table(FNZSL_vcfsamples, file = "FNZSL_samples.txt", sep = "\t")
## Assign population labels (L or H) to each sample in the above .txt files and save as [pool]_popfile.txt

########################
## READ POP FILES
########################
WNZLL_C_vcfpops <- read.table(file = "WNZLL_popfile.txt", header = TRUE) 
WNZSL_C_vcfpops <- read.table(file = "WNZSL_popfile.txt", header = TRUE)
WUSLL_C_vcfpops <- read.table(file = "WUSLL_popfile.txt", header = TRUE)
FNZLL_C_vcfpops <- read.table(file = "FNZLL_popfile.txt", header = TRUE)
FNZSL_C_vcfpops <- read.table(file = "FNZSL_popfile.txt", header = TRUE)

## Convert vcf file to genind to hierfstat to bayescan format
WNZLL_C_GI <- vcfR2genind(WNZLL_vcf)
ploidy(WNZLL_C_GI) <- 2
strata(WNZLL_C_GI) <- WNZLL_C_vcfpops
setPop(WNZLL_C_GI) <- ~Pop
popNames(WNZLL_C_GI)
WNZLL_C_GI
WNZLL_C_hierfstat <- genind2hierfstat(WNZLL_C_GI)
write.bayescan(dat = WNZLL_C_hierfstat, diploid = TRUE, fn = "WNZLL.bsc")

WNZSL_C_GI <- vcfR2genind(WNZSL_vcf)
ploidy(WNZSL_C_GI) <- 2
strata(WNZSL_C_GI) <- WNZSL_C_vcfpops
setPop(WNZSL_C_GI) <- ~Pop
popNames(WNZSL_C_GI)
WNZSL_C_GI
WNZSL_C_hierfstat <- genind2hierfstat(WNZSL_C_GI)
write.bayescan(dat = WNZSL_C_hierfstat, diploid = TRUE, fn = "WNZSL.bsc")

WUSLL_C_GI <- vcfR2genind(WUSLL_vcf)
ploidy(WUSLL_C_GI) <- 2
strata(WUSLL_C_GI) <- WUSLL_C_vcfpops
setPop(WUSLL_C_GI) <- ~Pop
popNames(WUSLL_C_GI)
WUSLL_C_GI
WUSLL_C_hierfstat <- genind2hierfstat(WUSLL_C_GI)
write.bayescan(dat = WUSLL_C_hierfstat, diploid = TRUE, fn = "WUSLL.bsc")

FNZLL_C_GI <- vcfR2genind(FNZLL_vcf)
ploidy(FNZLL_C_GI) <- 2
strata(FNZLL_C_GI) <- FNZLL_C_vcfpops
setPop(FNZLL_C_GI) <- ~Pop
FNZLL_C_GI
FNZLL_C_hierfstat <- genind2hierfstat(FNZLL_C_GI)
write.bayescan(dat = FNZLL_C_hierfstat, diploid = TRUE, fn = "FNZLL.bsc")

FNZSL_C_GI <- vcfR2genind(FNZSL_vcf)
ploidy(FNZSL_C_GI) <- 2
strata(FNZSL_C_GI) <- FNZSL_C_vcfpops
setPop(FNZSL_C_GI) <- ~Pop
FNZSL_C_GI
FNZSL_C_hierfstat <- genind2hierfstat(FNZSL_C_GI)
write.bayescan(dat = FNZSL_C_hierfstat, diploid = TRUE, fn = "FNZSL.bsc")

## BAYESCAN INPUT FILES STORED IN \\hpcsamba\dataset\GBS_WClover_WSC\archive\Sofie\Analyses\Bayescan
########################
## RUN BAYESCAN ON HPC:
########################
#cd home/pearsons/conda-envs/BayeScan2.1/binaries
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WNZLL.bsc
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WNZSL.bsc
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WUSLL.bsc
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/FNZLL.bsc
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/FNZSL.bsc

## COPY BAYESCAN RESULTS INTO: \\hpcsamba\dataset\GBS_WClover_WSC\archive\Sofie\Analyses\Bayescan\Results
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/Results")

########################
## FIND OUTLIERS
########################
## Run the following function before the ***plot_bayescan("WNZLL_fst.txt",0,FDR=0.01)*** command

plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}
# Typical usage:
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type
# > plot_bayescan("output_fst.txt",0,FDR=0.05)
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",0,FDR=0.05)
# results$outliers
# results$nb_outliers
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/Results")

plot_bayescan("WNZLL_fst.txt",0,FDR=0.01)
plot_bayescan("WNZSL_fst.txt",0,FDR=0.01)
plot_bayescan("WUSLL_fst.txt",0,FDR=0.01)
plot_bayescan("FNZLL_fst.txt",0,FDR=0.01)
plot_bayescan("FNZSL_fst.txt",0,FDR=0.01)

plot_bayescan("WNZLL_fst.txt",0,FDR=0.05)
plot_bayescan("WNZSL_fst.txt",0,FDR=0.05)
plot_bayescan("WUSLL_fst.txt",0,FDR=0.05)
plot_bayescan("FNZLL_fst.txt",0,FDR=0.05)
plot_bayescan("FNZSL_fst.txt",0,FDR=0.05)

## SNPs are outliers if the q-value from the [Pool]_fst.txt file is less than 0.01 (at a FDR alpha of 1%)
## Create manhattan plot input files in excel - need CHR, BP, P and SNP as headers
## remove "S" infront of chromosome number and remove SNPs from chromosome 17
## See below the manhattan_q_value function

## plotting posterior distribution
mydata=read.table("WNZLL.sel",colClasses="numeric")
parameter="Fst1"
plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
## you can plot population specific Fst coefficient by setting
## parameter="Fst1"
## you also have access to the likelihood with:
## parameter="logL"

# install.packages("boa")
## if you have the package "boa" installed, you can very easily obtain Highest Probability 
## Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
library(boa)
boa.hpd(mydata[[parameter]],0.05)
# Lower Bound Upper Bound 
#   0.1507167   0.1633356 

########################
## PLOT MANHATTAN PLOTS FOR Q-VALUE AND FST
########################
library(calibrate)
manhattan_q_value <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                              col=c("#e6194b", "#fabebe", "#f58231","#ffe119", "#aaffc3", "#3cb44b", "#46f0f0", "#4363db"), chrlabs=NULL,
                              suggestiveline=-log10(0.05), genomewideline=-log10(0.01), 
                              highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  
  # Not sure why, but package check will warn without this.
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
                   xlab=xlabel, ylab=expression(-log[10](italic(q)-value)))
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

########################
## Find chromosome and position info to create Manhattan plots:
########################
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan")

## chromosomal and bp position
head(WNZLL_C_GI$loc.n.all)
write.table(WNZLL_C_GI$loc.n.all, file = "WNZLL_SNP_order.txt", sep = "\t")
write.table(WNZSL_C_GI$loc.n.all, file = "WNZSL_SNP_order.txt", sep = "\t")
write.table(WUSLL_C_GI$loc.n.all, file = "WUSLL_SNP_order.txt", sep = "\t")
write.table(FNZLL_C_GI$loc.n.all, file = "FNZLL_SNP_order.txt", sep = "\t")
write.table(FNZSL_C_GI$loc.n.all, file = "FNZSL_SNP_order.txt", sep = "\t")

## Results from BayeScan are in [pool]_fst.txt files
## from these fst files, use the "qval" and "fst" columns for the manhattan plots
## Created manhattan input files and saved as "[pool]_MAN_P_INPUT.txt" and "[pool]_MAN_FST_INPUT.txt"
## Each one has either the qvalue or FST in the "P" column

setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/Results")
## Import chromosome, bp, SNP and q-values ("P") for each pool
WLL_q <- read.table(file = "WNZLL_MAN_P_INPUT.txt", header = T)
WSL_q <- read.table(file = "WNZSL_MAN_P_INPUT.txt", header = T)
WUS_q <- read.table(file = "WUSLL_MAN_P_INPUT.txt", header = T)
FLL_q <- read.table(file = "FNZLL_MAN_P_INPUT.txt", header = T)
FSL_q <- read.table(file = "FNZSL_MAN_P_INPUT.txt", header = T)

########################
## MAKE LIST OF OUTLIER SNPS YOU WANT TO HIGHLIGHT FOR EACH MANHATTAN PLOT
########################
## Excel spreadsheet saved as "Bayescan_FDR-0.01_results_Dec_2019" in \\hpcsamba\dataset\GBS_WClover_WSC\archive\Sofie\Analyses\Bayescan\Results

## Copy common outliers into a character vector for plotting:

## These SNPs are outliers in more than 2 pools (and have an FST of greater than 0.3)
# WNZLL_1_q_outliers <- c("2_23112313","9_1044750","13_4850703","5_28434929","11_46932936","11_63153863","4_13559491","16_32428574","8_40904996","8_40905003","8_40905002")
# WNZSL_1_q_outliers <- c("4_13559491","8_40904996","8_40905002","8_40905003","4_71509072")
# WUSLL_1_q_outliers <- c("4_13559491","9_1044750","2_40720445","11_4462899","5_28434929")
# FNZLL_1_q_outliers <- c("11_4462899","11_63153863","11_46932936","4_71509072")
# FNZSL_1_q_outliers <- c("13_4850703","16_32428574","2_23112313","2_40720445")

WNZLL_5_q_outliers <- c("2_23112313","9_1044750","13_4850703","14_7499208","5_28434929","11_46932936","11_63153863","4_13559491","16_32428574","13_23675501","8_40904996","8_40905003","8_40905002","1_3522737","12_19238640")
WNZSL_5_q_outliers <- c("4_13559491","8_40904996","8_40905002","8_40905003","4_71509072","11_17480539","12_19238640","12_16032423")
WUSLL_5_q_outliers <- c("4_13559491","9_1044750","11_21109408","11_21109404","2_40720445","11_4462899","5_28434929","2_32222083","2_14186624","2_14186629","13_23675463","5_39504991")
FNZLL_5_q_outliers <- c("11_4462899","11_63153863","11_46932936","4_71509072","12_16032423","1_3522737","2_14186629","2_14186624","14_7499208","11_17480539","13_23675501","13_23675463","2_32222083")
FNZSL_5_q_outliers <- c("13_4850703","16_32428574","2_23112313","5_39504991","2_40720445","11_21109404","11_21109408")

########################
## Create Manhattan plot:
########################
par(mfrow=c(3,2))
manhattan_q_value(WLL_q,annotatePval = FALSE, suggestiveline = T, highlight = WNZLL_5_q_outliers, ylim =c(0,5))
title(expression(bold(A)~~"WNZLL"), adj = 0, cex.main = 2)
manhattan_q_value(FLL_q,annotatePval = FALSE, suggestiveline = T, highlight = FNZLL_5_q_outliers, ylim =c(0,5))
title(expression(bold(D)~~"FNZLL"), adj = 0, cex.main = 2)
manhattan_q_value(WSL_q,annotatePval = FALSE, suggestiveline = T, highlight = WNZSL_5_q_outliers, ylim =c(0,5))
title(expression(bold(B)~~"WNZSL"), adj = 0, cex.main = 2)
manhattan_q_value(FSL_q,annotatePval = FALSE, suggestiveline = T, highlight = FNZSL_5_q_outliers, ylim =c(0,5))
title(expression(bold(E)~~"FNZSL"), adj = 0, cex.main = 2)
manhattan_q_value(WUS_q,annotatePval = FALSE, suggestiveline = T, highlight = WUSLL_5_q_outliers, ylim =c(0,5))
title(expression(bold(C)~~"WUSLL"), adj = 0, cex.main = 2)
# save image width: 1600 and height:850
# Q_VALUE_MANHATTAN_PLOTS.tiff

########################
## Find outliers in common with KGD and PCAdapt analyses
########################
## UPSET PLOT INFO BELOW:

############################################################################################
## IN EXCEL, OPEN UP THE "POOL_Cx_SNPs.txt" FILE WITH THE CORRESPONDING "RESULTS.txt" FILE
## save both in "Bayescan_outlier_SNP_ID_results.xlsx" with a separate tab for each grouping (both FDR in same tab)
## THE COLUMNS WITH THE YELLOW HIGHLIGHTED HEADING ARE THE ONES THAT WILL BE COMPARED TO FIND THE SNP NAME

############################################################################################
## LIST OF OUTLIER SNPS AT FDR = 0.01 ARE STORED IN "compiled-C2_C4_C3_C6-FDR_0.01" SHEET
## NAME OF WHICH POOL AND POINT OF SELECTION CYCLE AS HEADER
## Save data as csv AND IMPORT INTO R TO CREATE UPSET PLOT
library(readr)
COMPILED_C2_C4_C3_C6_FDR_0_01 <- read_csv("COMPILED-C2-C4-C3-C6-FDR-0.01.csv")
View(COMPILED_C2_C4_C3_C6_FDR_0_01)

ListInput <- list(WNZLL_C2 = COMPILED_C2_C4_C3_C6_FDR_0_01$WNZLL_C2[1:65],
                  WNZLL_C4 = COMPILED_C2_C4_C3_C6_FDR_0_01$WNZLL_C4[1:86],
                  WNZSL_C2 = COMPILED_C2_C4_C3_C6_FDR_0_01$WNZSL_C2,
                  WNZSL_C4 = COMPILED_C2_C4_C3_C6_FDR_0_01$WNZSL_C4[1:69],
                  WUSLL_C2 = COMPILED_C2_C4_C3_C6_FDR_0_01$WUSLL_C2[1:76],
                  WUSLL_C4 = COMPILED_C2_C4_C3_C6_FDR_0_01$WUSLL_C4[1:68],
                  FNZLL_C3 = COMPILED_C2_C4_C3_C6_FDR_0_01$FNZLL_C3[1:50],
                  FNZLL_C6 = COMPILED_C2_C4_C3_C6_FDR_0_01$FNZLL_C6[1:48],
                  FNZSL_C3 = COMPILED_C2_C4_C3_C6_FDR_0_01$FNZSL_C3[1:58],
                  FNZSL_C6 = COMPILED_C2_C4_C3_C6_FDR_0_01$FNZSL_C6[1:65]
)
library(UpSetR)
upset(fromList(ListInput), order.by = "freq", nsets = 10, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per population pair", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2, 2), 
      set_size.show = TRUE, set_size.scale_max = 100)
## save image width: 1771 and height: 796
"Upset_WSC_pool_selection_cycles_BayeScan_outliers.png"

## combine the selection cycle SNPs and have the outliers just for each pool
## SAVED CSV AS "COMPILED-POOL-FDR-0.01.CSV"
## IMPORT DATA INTO R
COMPILED_POOL_FDR_0_01 <- read_csv("COMPILED-POOL-FDR-0.01.csv")
View(COMPILED_POOL_FDR_0_01)

ListInputPOOL <- list(WNZLL = COMPILED_POOL_FDR_0_01$WNZLL[1:120],
                      WNZSL = COMPILED_POOL_FDR_0_01$WNZSL[1:127],
                      WUSLL = COMPILED_POOL_FDR_0_01$WUSLL[1:117],
                      FNZLL = COMPILED_POOL_FDR_0_01$FNZLL[1:87],
                      FNZSL = COMPILED_POOL_FDR_0_01$FNZSL[1:111])
upset(fromList(ListInputPOOL), order.by = "freq", nsets = 5, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per pool", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2.5, 2), 
      set_size.show = TRUE, set_size.scale_max = 150)
## save image width: 1092 and height: 620
"Upset_WSC_pool_BayeScan_outliers.png"

## 43 OUTLIERS IN COMMON BETWEEN MULTIPLE POOLS

## FOUND LIST OF OUTLIER SNPS FROM PCADAPT ANALYSIS (SAVED IN "PCADAPT_FINAL_SCRIPT_WSC.bash")
## made combination of pcadapt and bayescan outliers and import into R
PCADAPT_AND_BAYESCAN_OUTLIERS <- read_csv("PCADAPT_AND_BAYESCAN_OUTLIERS.csv", 
                                          col_types = cols(X11 = col_skip(), X12 = col_skip(),X13 = col_skip(), X14 = col_skip(), X15 = col_skip(), X16 = col_skip(), X17 = col_skip()))
## make upset plot for bayescan and pcadapt
ListInputPCADAPT_BAYESCAN <- list(
  WNZLL_PCADAPT = PCADAPT_AND_BAYESCAN_OUTLIERS$WNZLL_PCADAPT[1:77],
  WNZLL_BAYESCAN = PCADAPT_AND_BAYESCAN_OUTLIERS$WNZLL_BAYESCAN[1:120],
  WNZSL_PCADAPT = PCADAPT_AND_BAYESCAN_OUTLIERS$WNZSL_PCADAPT[1:44],
  WNZSL_BAYESCAN = PCADAPT_AND_BAYESCAN_OUTLIERS$WNZSL_BAYESCAN[1:127],
  WUSLL_PCADAPT = PCADAPT_AND_BAYESCAN_OUTLIERS$WUSLL_PCADAPT[1:57],
  WUSLL_BAYESCAN = PCADAPT_AND_BAYESCAN_OUTLIERS$WUSLL_BAYESCAN[1:117],
  FNZLL_PCADAPT = PCADAPT_AND_BAYESCAN_OUTLIERS$FNZLL_PCADAPT[1:83],
  FNZLL_BAYESCAN = PCADAPT_AND_BAYESCAN_OUTLIERS$FNZLL_BAYESCAN[1:87],
  FNZSL_PCADAPT = PCADAPT_AND_BAYESCAN_OUTLIERS$FNZSL_PCADAPT[1:29],
  FNZSL_BAYESCAN = PCADAPT_AND_BAYESCAN_OUTLIERS$FNZSL_BAYESCAN[1:111]
)
upset(fromList(ListInputPCADAPT_BAYESCAN), 
      sets = c("WNZLL_PCADAPT", "WNZLL_BAYESCAN", 
               "WNZSL_PCADAPT", "WNZSL_BAYESCAN",
               "WUSLL_PCADAPT", "WUSLL_BAYESCAN",
               "FNZLL_PCADAPT", "FNZLL_BAYESCAN",
               "FNZSL_PCADAPT","FNZSL_BAYESCAN"),
      group.by = "sets", order.by = "freq", 
      keep.order = TRUE, nsets = 10, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per pool", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2, 2), 
      set_size.show = TRUE, set_size.scale_max = 150)
## still playing with the above parameters

upset(fromList(ListInputPCADAPT_BAYESCAN), order.by = "freq", 
      nsets = 10, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per pool", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2, 2), 
      set_size.show = TRUE, set_size.scale_max = 150)
## this above one is okay for now. save image width: 1771 and height: 796
"Upset_WSC_pool_PCADAPT_and_BAYESCAN_outliers.png"

upset(fromList(ListInputPCADAPT_BAYESCAN), 
      nsets = 10, point.size = 3.5, 
      line.size = 2, mainbar.y.label = "Number of outlier SNPs", 
      sets.x.label = "Outlier SNPs per pool", 
      text.scale = c(2, 1.7, 1.5, 1.5, 2, 2), 
      set_size.show = TRUE, set_size.scale_max = 150)

###################################################################################
## Run 24 populations all together:
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses")
## vcf files with populations I want to compare:
ALL_vcf <- read.vcfR("ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf")
## 1113 SAMPLES, 14743 SNPs

## Find order of samples in each vcf file and assign population info to create popfiles
## Sample order found in "1113_sample_order.txt"
## Assign population labels to each sample in the above .txt files and save as 1113_popfile.txt

## READ POP FILES
ALL_pops <- read.table(file = "1113_popfile.txt", header = TRUE) 

## Convert vcf file to genind to hierfstat to bayescan format
ALL_vcf_GI <- vcfR2genind(ALL_vcf)
ploidy(ALL_vcf_GI) <- 2
strata(ALL_vcf_GI) <- ALL_pops
setPop(ALL_vcf_GI) <- ~pop
popNames(ALL_vcf_GI)
ALL_vcf_GI
ALL_hierfstat <- genind2hierfstat(ALL_vcf_GI)
setwd("C:/Sofie")
write.bayescan(dat = ALL_hierfstat, diploid = TRUE, fn = "1113.bsc")
## Copy bsc file to HPC (Sofie/Analyses/Bayescan/1113.bsc)
# pwd = /home/pearsons/conda-envs/BayeScan2.1/binaries
# set prior odds at 100
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/1113.bsc -pr_odds 100 -out_freq
## Check FST distribution

###################################################################################
## Try WNZLL pool with parental plants added in:

## Run 5 WNZLL populations together:
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan")
## vcf files with populations I want to compare:
WLL_vcf <- read.vcfR("WNZLL_POOL.vcf")
## 235 SAMPLES, 14743 SNPs

## Find order of samples in each vcf file and assign population info to create popfiles
## Sample order found in "1113_sample_order.txt"
## Assign population labels to each sample in the above .txt files and save as 1113_popfile.txt

## READ POP FILES
WLL_pops <- read.table(file = "WNZLL_popfile_5.txt", header = TRUE) 

## Convert vcf file to genind to hierfstat to bayescan format
WLL_vcf_GI <- vcfR2genind(WLL_vcf)
ploidy(WLL_vcf_GI) <- 2
strata(WLL_vcf_GI) <- WLL_pops
setPop(WLL_vcf_GI) <- ~pop
popNames(WLL_vcf_GI)
WLL_vcf_GI
WLL_hierfstat <- genind2hierfstat(WLL_vcf_GI)
setwd("C:/Sofie")
write.bayescan(dat = WLL_hierfstat, diploid = TRUE, fn = "WNZLL_POOL.bsc")
## Copy bsc file to HPC (Sofie/Analyses/Bayescan/WNZLL_POOL.bsc)
# pwd = /home/pearsons/conda-envs/BayeScan2.1/binaries
# set prior odds at 100
# ./BayeScan2.1_linux64bits /dataset/GBS_WClover_WSC/archive/Sofie/Analyses/Bayescan/WNZLL_POOL.bsc
## Check FST distribution

############################################################################
## KGD-FST analysis
########################

## setwd("/bifo/scratch/GBS_WClover_WSC/Sofie/fst_test")
## Final data and analysis stored in archive:
setwd("//hpcsamba/dataset/GBS_WClover_WSC/archive/Sofie/FST")

gform <- "Tassel"	
genofile <- "ALL_DATA_nofail_no_pos_or_neg_dp5to150_mm0.2_maf0.03_no_multiallelic.vcf.ra.tab"
source("//hpcsamba/dataset/Ryegrass_Bulks/active/bin/KGD/GBS-Chip-Gmatrix.R") #Run this script
library('plotly')
controls <- which(grepl("Blank", seqID, ignore.case=TRUE)==TRUE) ##identifying positive and negative controls
controls <- append(controls, which(grepl("S9-", seqID, ignore.case=TRUE)==TRUE)) ##identifying positive and negative controls
samp.remove(controls, keep=FALSE) ##removing controls, not needed if previously done 
Gfull      <- calcG() ##calculating G Matrix
#subset_SNPs <- which(HWdis > -0.05) ##Filtering low HWdis SNPs

info.file <- read.table("24_SAMPLE_POP.txt", header=TRUE) ##Reading in metadata file with sample and population information
populations <- rep(NA, length(seqID)) ##Generating vector of populations per sample
populations <- info.file$pop[match(seqID, info.file$sample)] ##Generating vector of populations per sample

#groups <- populations ##Needed for plotly PCA
#samp.info <- list(populations=populations, seqID=seqID) ##Needed for plotly PCA
#groups <- as.vector(populations) ##Needed for plotly PCA
#GHWdgm.05  <- calcG(subset_SNPs, "HWdgm.05_HWdgm.05", npc=3,withPlotly=T,plotly.group=groups, samp.info=samp.info) ##Needed for plotly PCA
GHWdgm.05  <- calcG(subset_SNPs, "HWdgm.05_HWdgm.05", npc=3) ##Re-calculating GRM with filtered data
writeG(Guse=GHWdgm.05, outname="clover_wsc_Nov18.vcf.depth_5_250", outtype=1:6) ##Writing GRM and principle components to files

source("//hpcsamba/dataset/Ryegrass_Bulks/active/bin/KGD/GBS-PopGen.R") ##Sourcing pop gen functions
het <- heterozygosity() ##calculating observed and expected heterozygosity
popmaf(populations=populations) ##Generating MAF distribution plots per population
fst.pairwise <- Fst.GBS.pairwise(populations=populations, sortlevels=TRUE) ##Calculating pairwise Fst estimates between each population
#up to here
pos_sub1 <- pos[] ##Vector of SNP position from filtered SNPs
chrom_sub1 <- chrom[] ##Vector of each chromosme SNP position from filtered SNPs 
manhatplot(fst.pairwise[1,2,], chrom_sub1, pos_sub1, "FNZLL-C3-H_FNZLL-C3-L", qdistn=qchisq, df=2) ##Generating manhattan plot

## Made manhattan plot order chromosomes numerically
manhatplot_SP <- function(value, chrom, pos, plotname, qdistn=qunif, keyrot=0, ...) {
  chromcol <- colourby(chrom)
  colkey(chromcol,"chrom",srt=keyrot)
  plotord <- order(as.numeric(chrom),as.numeric(pos))
  png(paste0(plotname,"-Manhat.png"),width=1200)
  plot(value[plotord], col=chromcol$sampcol[plotord],xlab="Position",ylab=substitute(value),cex=0.8, xaxt="n")
  dev.off()
  png(paste0(plotname,"-QQ.png"))
  qqplot(qdistn(ppoints(length(value)),...), y=value, xlab="Theoretical quantiles", ylab=paste0(substitute(value)," quantiles"), 
         sub="Line for mid 98% of values", col=chromcol$sampcol[order(value)])
  qqline(value,col=2, distribution = function(p) qdistn(p, ...),prob=c(0.01,0.99))
  dev.off()
}

levels(populations)
# [1] "FNZLL-C3-H" "FNZLL-C3-L" "FNZLL-C6-H" "FNZLL-C6-L" "FNZLL-P"    "FNZSL-C3-H" "FNZSL-C3-L" "FNZSL-C6-H" "FNZSL-C6-L"
# [10] "FNZSL-P"    "WNZLL-C2-H" "WNZLL-C2-L" "WNZLL-C4-H" "WNZLL-C4-L" "WNZLL-P"    "WNZSL-C2-H" "WNZSL-C2-L" "WNZSL-C4-H"
# [19] "WNZSL-C4-L" "WUSLL-C2-H" "WUSLL-C2-L" "WUSLL-C4-H" "WUSLL-C4-L" "WUSLL-P" 
manhatplot_SP(fst.pairwise[1,2,], chrom_sub1, pos_sub1, "FNZLL-C3-H_FNZLL-C3-L", qdistn=qchisq, df=2) ##Generating manhattan plot
manhatplot_SP(fst.pairwise[3,4,], chrom_sub1, pos_sub1, "FNZLL-C6-H_FNZLL-C6-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[6,7,], chrom_sub1, pos_sub1, "FNZSL-C3-H_FNZSL-C3-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[8,9,], chrom_sub1, pos_sub1, "FNZSL-C6-H_FNZSL-C6-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[11,12,], chrom_sub1, pos_sub1, "WNZLL-C2-H_WNZLL-C2-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[13,14,], chrom_sub1, pos_sub1, "WNZLL-C4-H_WNZLL-C4-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[16,17,], chrom_sub1, pos_sub1, "WNZSL-C2-H_WNZSL-C2-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[18,19,], chrom_sub1, pos_sub1, "WNZSL-C4-H_WNZSL-C4-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[20,21,], chrom_sub1, pos_sub1, "WUSLL-C2-H_WUSLL-C2-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise[22,23,], chrom_sub1, pos_sub1, "WUSLL-C4-H_WUSLL-C4-L", qdistn=qchisq, df=2)

## Get FST for each SNP in each pairwise comparison
dim(fst.pairwise)
# [1]    24    24 14743
library(reshape2)
df <- dcast(melt(fst.pairwise), Var1 + Var2~Var3)
write.table(df, file = "fst_pairwise_out.txt")

########################
## USE K = 11 POPULATIONS
########################
info.file1 <- read.table("11_SAMPLE_POP.txt", header=TRUE) ## Reading in metadata file with sample and population information
populations1 <- rep(NA, length(seqID)) ##Generating vector of populations per sample
populations1 <- info.file1$pop[match(seqID, info.file1$sample)] ## Generating vector of populations per sample
popmaf(populations=populations1) ##Generating MAF distribution plots per population
fst.pairwise1 <- Fst.GBS.pairwise(populations=populations1, sortlevels=TRUE) ## Calculating pairwise Fst estimates between each population

pos_sub1 <- pos[] ## Vector of SNP position from filtered SNPs
chrom_sub1 <- chrom[] ## Vector of each chromosme SNP position from filtered SNPs 

manhatplot_SP <- function(value, chrom, pos, plotname, qdistn=qunif, keyrot=0, ...) {
  chromcol <- colourby(chrom)
  colkey(chromcol,"chrom",srt=keyrot)
  plotord <- order(as.numeric(chrom),as.numeric(pos))
  png(paste0(plotname,"-Manhat.png"),width=1200)
  plot(value[plotord], col=chromcol$sampcol[plotord],xlab="Position",ylab=expression('Pairwise  F'[ST]),cex=0.8, xaxt="n")
  dev.off()
  png(paste0(plotname,"-QQ.png"))
  qqplot(qdistn(ppoints(length(value)),...), y=value, xlab="Theoretical quantiles", ylab=paste0(substitute(value)," quantiles"), 
         sub="Line for mid 98% of values", col=chromcol$sampcol[order(value)])
  qqline(value,col=2, distribution = function(p) qdistn(p, ...),prob=c(0.01,0.99))
  dev.off()
}

levels(populations1)
# "FNZLL-H" "FNZLL-L" "FNZSL-H" "FNZSL-L" "PARENT"  "WNZLL-H" "WNZLL-L" "WNZSL-H" "WNZSL-L" "WUSLL-H" "WUSLL-L"
manhatplot_SP(fst.pairwise1[1,2,], chrom_sub1, pos_sub1, "FNZLL-H_FNZLL-L", qdistn=qchisq, df=2) ##Generating manhattan plot
manhatplot_SP(fst.pairwise1[3,4,], chrom_sub1, pos_sub1, "FNZSL-H_FNZSL-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise1[6,7,], chrom_sub1, pos_sub1, "WNZLL-H_WNZLL-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise1[8,9,], chrom_sub1, pos_sub1, "WNZSL-H_WNZSL-L", qdistn=qchisq, df=2)
manhatplot_SP(fst.pairwise1[10,11,], chrom_sub1, pos_sub1, "WUSLL-H_WUSLL-L", qdistn=qchisq, df=2)

## Get FST for each SNP in each pairwise comparison
dim(fst.pairwise1)
# [1]    11    11 14743
library(reshape2)
df1 <- dcast(melt(fst.pairwise1), Var1 + Var2~Var3)
write.table(df1, file = "fst_pairwise_out_K11.txt")
## in excel fill in the var pop names and keep the comparisons you want
## find order of SNPs - chromosome and position and write table
write.table(chrom, file = "chrom_order.txt")
write.table(pos, file = "pos_order.txt")

## Import chromosome, bp, SNP and FST values for each pool
FLL_FST <- read.table("FNZLL_FST.txt",sep="\t",header=TRUE)
FSL_FST <- read.table("FNZSL_FST.txt",sep="\t",header=TRUE)
WLL_FST <- read.table("WNZLL_FST.txt",sep="\t",header=TRUE)
WSL_FST <- read.table("WNZSL_FST.txt",sep="\t",header=TRUE)
WUS_FST <- read.table("WUSLL_FST.txt",sep="\t",header=TRUE)

## MAKE LIST OF OUTLIER SNPS YOU WANT TO HIGHLIGHT FOR EACH MANHATTAN PLOT
## These SNPs are present in more than 2 pools and have an FST of greater than 0.3
WNZLL0.3 <- c("1_3522737","1_52318951","1_88536713","1_88536715","1_88536722","1_88536745","1_88536753","2_17285940","2_19713712","2_22339748","2_23112313","2_24478676","2_24478681","2_32471335","2_32471346","2_50131173","2_50131186","2_50131221","2_50131225","3_8924597","3_8924618","3_24629067","3_24629078","3_24629098","3_28211144","3_48617602","3_48617607","3_52167304","3_54016436","4_13559491","4_24997513","4_33506228","4_47136437","4_47136450","4_47136456","4_54982767","4_65854669","4_71509072","4_72861142","4_75238813","5_5810604","5_5810606","5_5810616","5_5810638","5_23376003","5_37306891","5_40927453","5_40927454","5_45087490","6_937289","6_9463432","6_12367754","6_15128286","6_17993297","6_17993317","6_19904807","6_26636343","6_28860092","6_28860101","7_26536796","7_46864277","8_20902710","8_21857020","8_22608025","8_23866166","8_40904996","8_40905002","8_40905003","8_49745906","8_49745930","9_1044750","9_11025810","9_20443858","9_20479288","9_20666911","9_28140350","9_28140353","9_29022193","9_30605540","9_50695342","9_50695348","9_50695373","9_58721639","9_58721645","9_58721646","9_60771287","9_60771332","10_302864","10_20927968","11_9270652","11_9345053","11_9345077","11_18271224","11_18271235","11_24485515","11_26653519","11_42909236","11_46136083","11_46932936","11_58320018","11_61697899","11_63153863","12_2164767","12_12639439","12_24136974","12_50027346","12_50027362","13_4850703","13_4850721","13_10076378","13_10076420","13_10676963","13_23675477","13_23675501","14_6628207","14_7499208","14_29883839","14_29964501","14_31922837","15_6515376","15_6515380","15_14594474","16_16208938","16_16208955","16_32428574","16_32428586","16_32428616","16_34749250")
WNZSL0.3	<- c("1_62358850","2_6663484","2_18189884","2_31152286","2_45181105","3_12212339","3_12212364","3_30076983","3_49659978","3_49659979","3_54016436","4_9512082","4_13559491","4_19134918","4_21710417","4_30841647","4_37385889","4_37385894","4_37385918","4_54982767","4_71509072","5_21556617","6_12367754","6_15128286","6_24241000","6_24241025","6_26636343","7_26536786","7_26536796","7_27080552","7_35302392","7_46864277","8_15066059","8_15066082","8_23866139","8_23866166","8_24770464","8_40904996","8_40905002","8_40905003","8_46683920","9_1044750","9_20443858","9_20479288","9_20666911","9_50695342","9_50695348","9_50695373","9_58748228","9_59794335","10_31677849","11_9345053","11_9345077","11_17480539","11_21109384","11_21109404","11_21109408","12_12639439","12_16032423","12_24136974","12_33878524","12_50027346","12_50027362","13_4850721","13_5205472","13_5205476","13_5205477","13_32324441","13_32324479","13_40514506","14_29964501","15_6515376","15_6515380","16_14595671","16_16208938","16_16208955")
WUSLL0.3	<- c("1_25533350","1_26165134","1_62358850","1_71060882","1_71060883","1_71060888","1_79119252","1_88536713","1_88536715","1_88536722","1_88536745","1_88536753","2_6673787","2_14186624","2_14186629","2_22339748","2_47562237","2_47562240","2_50131173","2_50131186","2_50131221","2_50131225","3_24325289","3_30076983","3_52167304","4_9733285","4_13559491","4_19134918","4_33506228","4_36494940","4_36494946","4_36494976","4_47136437","4_47136450","4_47136456","4_65854669","4_72861142","5_5810616","5_5810638","5_31890940","5_37306891","5_40927453","5_40927454","6_9463432","6_14364301","6_15128286","6_24241000","6_24241025","6_31429353","6_31429365","6_35669244","6_35669268","6_35669280","7_22745064","7_35221900","7_35302392","7_52124106","8_20902710","8_22608025","8_23866139","9_1044750","9_17663260","9_17663267","9_20479288","9_39773614","9_41449789","9_41449808","9_58748228","10_302864","10_28507355","10_31677849","11_4462899","11_21109384","11_21109404","11_21109408","11_26653519","12_2164767","12_33878524","13_4850703","13_10076378","13_10076420","13_23675463","13_23675477","14_31922837","15_11569524","16_18188475","16_32428575","16_32428586","16_32428616","16_34749250")
FNZLL0.3	<- c("1_3522737","1_25533350","1_52318951","1_62358850","1_68649268","1_71060882","1_71060883","1_71060888","1_88536713","1_88536715","1_88536722","1_88536745","1_88536753","2_6663484","2_6673787","2_14186624","2_14186629","2_17285940","2_18189884","2_19713712","2_31152286","2_32471335","2_32471346","2_45181105","2_47562237","2_47562240","2_50131225","3_8924597","3_8924618","3_12212339","3_12212364","3_24629067","3_24629078","3_24629098","3_28211144","3_49659978","3_49659979","3_61266795","4_9512082","4_9733285","4_13559491","4_21710417","4_24997513","4_36494940","4_36494946","4_36494976","4_47136437","4_71509072","5_21556617","5_40927453","5_40927454","6_11698503","6_17993297","6_17993345","6_24241000","6_24241025","6_31429353","6_31429365","6_35669244","6_35669268","6_35669280","7_27080552","7_35221900","7_35302392","7_52124106","8_20902710","8_21857020","8_22608025","8_24770464","8_46683920","8_49745906","8_49745930","8_49778173","9_1044750","9_11025820","9_11025821","9_17663260","9_17663267","9_30605540","9_36684034","9_39773614","9_58721639","9_58721645","9_58721646","10_6118067","10_6118068","10_20927961","10_20927968","10_28507355","11_4462899","11_9002661","11_9345053","11_9345077","11_14899937","11_17480539","11_18271224","11_18271235","11_46136083","11_46932936","11_58477409","11_61697899","11_63153863","12_3437942","12_16032423","12_33878524","12_50027346","12_50027362","13_10583141","13_10676963","13_14007722","13_15943291","13_23675463","13_23675501","13_32324441","13_32324479","13_40514506","14_7499208","14_29883839","15_6428131","15_11569524","15_17315253")
FNZSL0.3 <- c("1_26165134","1_62358850","1_68649268","1_79119252","2_19713712","2_23112313","2_24478676","2_24478681","3_12212339","3_12212364","3_24325289","3_30076983","3_48617602","3_48617607","3_49659978","3_49659979","3_61266795","4_30841647","4_37385889","4_37385894","4_37385918","4_65854669","4_75238813","5_5810604","5_5810606","5_5810616","5_5810638","5_23376003","5_31890940","5_45087490","6_937289","6_11698503","6_14364301","6_17993297","6_17993317","6_17993345","6_19904807","6_28860092","6_28860101","7_22745064","7_26536786","7_46864277","8_15066059","8_15066082","8_40904996","8_40905002","8_40905003","8_49778173","9_11025810","9_11025820","9_11025821","9_28140350","9_28140353","9_29022193","9_30605540","9_36684034","9_41449789","9_41449808","9_50695342","9_50695348","9_50695373","9_59794335","9_60771287","9_60771332","10_6118067","10_6118068","10_20927961","10_20927968","11_4462899","11_9002661","11_9270652","11_14899937","11_21109404","11_21109408","11_24485515","11_42909236","11_58320018","11_58477409","12_3437942","12_12639439","13_4850703","13_5205472","13_5205476","13_5205477","13_10583141","13_14007722","13_15943291","13_23675477","13_23675501","14_6628207","15_6428131","15_14594474","15_17315253","16_14595671","16_18188475","16_32428574","16_32428575")

## Plot manhattan plot (run function at bottom of script)
manhattanSP(FLL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = FNZLL0.3, ylim = c(0,1), logp=F, main = "Manhattan plot of pairwise FST (>0.3) for FNZLL")
manhattanSP(FSL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = FNZSL0.3, ylim = c(0,1), logp=F, main = "Manhattan plot of pairwise FST (>0.3) for FNZSL")
manhattanSP(WLL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WNZLL0.3, ylim = c(0,1), logp=F, main = "Manhattan plot of pairwise FST (>0.3) for WNZLL")
manhattanSP(WSL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WNZSL0.3, ylim = c(0,1), logp=F, main = "Manhattan plot of pairwise FST (>0.3) for WNZSL")
manhattanSP(WUS_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WUSLL0.3, ylim = c(0,1), logp=F, main = "Manhattan plot of pairwise FST (>0.3) for WUSLL")
## save image width: 1818 and height: 636
## FNZLL_MANHATTAN_FST.tiff ETC

## OTHER COLOUR COMBINATIONS
## c("dark blue", "cornflowerblue") #rrBLUP GWAS colours
## c("gray10", "gray60") #greyscale
## c("#e6194b", "#fabebe", "#f58231","#ffe119", "#aaffc3", "#3cb44b", "#46f0f0", "#4363db") #colourful

par(mfrow=c(3,2))
manhattanSP(WLL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WNZLL0.3, ylim = c(0,1.1), logp=F)
title(expression(bold(A)~~"WNZLL"), adj = 0, cex.main = 2)
manhattanSP(FLL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = FNZLL0.3, ylim = c(0,1.1), logp=F)
title(expression(bold(D)~~"FNZLL"), adj = 0, cex.main = 2)
manhattanSP(WSL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WNZSL0.3, ylim = c(0,1.1), logp=F)
title(expression(bold(B)~~"WNZSL"), adj = 0, cex.main = 2)
manhattanSP(FSL_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = FNZSL0.3, ylim = c(0,1.1), logp=F)
title(expression(bold(E)~~"FNZSL"), adj = 0, cex.main = 2)
manhattanSP(WUS_FST,annotatePval = FALSE, suggestiveline = FALSE, highlight = WUSLL0.3, ylim = c(0,1.1), logp=F)
title(expression(bold(C)~~"WUSLL"), adj = 0, cex.main = 2)
## save image width: 1600 and height:850
## FST_MANHATTAN_PLOTS.tiff