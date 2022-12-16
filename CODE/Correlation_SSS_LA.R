library(readr)
library(ggpubr)
library(ggplot2)
library(MASS)

########################
## Input LA and SSS datasheets (n = 447 and 600, respectively)
########################
## Correlation analysis includes 297 plants from 5 pools
## Transformed SSS and Total_WSC data to g/kg Dry Matter, i.e. multiplied raw data by 10
x <- read.table("Correlation_SSS_LA.txt", header = TRUE)
WLL <- read.table("Correlation_SSS_LA_WNZLL.txt", header = TRUE)
WSL <- read.table("Correlation_SSS_LA_WNZSL.txt", header = TRUE)
WUS <- read.table("Correlation_SSS_LA_WUSLL.txt", header = TRUE)
FLL <- read.table("Correlation_SSS_LA_FNZLL.txt", header = TRUE)
FSL <- read.table("Correlation_SSS_LA_FNZSL.txt", header = TRUE)

########################
## Test for data normality in each dataset
########################
## Shapiro.Wilk Normailty test
shapiro.test(x$NovLeafArea)
# W = 0.97248, p-value = 1.823e-05
shapiro.test(x$SSS)
# W = 0.99357, p-value = 0.2368 *
shapiro.test(WLL$NovLeafArea)
# W = 0.97635, p-value = 0.294 *
shapiro.test(WLL$SSS)
# W = 0.97745, p-value = 0.3303 *
shapiro.test(WSL$NovLeafArea)
# W = 0.91829, p-value = 0.000655
shapiro.test(WSL$SSS)
# W = 0.97881, p-value = 0.3804 *
shapiro.test(WUS$NovLeafArea)
# W = 0.98052, p-value = 0.486 *
shapiro.test(WUS$SSS)
# W = 0.97303, p-value = 0.2313 *
shapiro.test(FLL$NovLeafArea)
# W = 0.96888, p-value = 0.1285 *
shapiro.test(FLL$SSS)
# W = 0.97452, p-value = 0.2414 *
shapiro.test(FSL$NovLeafArea)
# W = 0.91797, p-value = 0.000636
shapiro.test(FSL$SSS)
# W = 0.97344, p-value = 0.2143 *
## * = fail to reject null hypothesis of normal distribution.
## All SSS datasets follow a normal distribution but 3 leaf area datasets do not:
## All, WNZSL and FNZSL need to be transformed.

########################
## Box-Cox Transformation - to determine optimal transformation for each of the three datasets
########################
## Graphical representation of lambda and Log-Likelihood
dev.off()
par(mfrow=c(3,1))

boxcox(lm(NovLeafArea~1, data = WSL)) #WSL_boxcox_LA
boxcox(lm(NovLeafArea~1, data = FSL)) #FSL_boxcox_LA
boxcox(lm(NovLeafArea~1, data = x))   #ALL_boxcox_LA
## Save each box-cox graph separately W = 517 H = 421

## Find lambda values for boxcox transformations
y <- boxcox(lm(NovLeafArea~1, data = WSL)) 
list(y) ## log-Likelihood: -59.72610, lambda: 0.18
y <- boxcox(lm(NovLeafArea~1, data = FSL)) 
list(y) ## -54.70156, -0.02
y <- boxcox(lm(NovLeafArea~1, data = x)) 
list(y) ## -538.8931,0.38

## Check normality after transformation
shapiro.test(FSL$log10_LA_1) ## This is 1+log10(LA)
# W = 0.98024, p-value = 0.4387 *
shapiro.test(log10(FSL$NovLeafArea)) ## This is log10(LA)
# W = 0.98129, p-value = 0.4856 *
## Both conform to normality now but which looks better:

dev.off()
par(mfrow=c(2,1))
hist(FSL$log10_LA_1)
hist(log10(FSL$NovLeafArea)) ## log10(LA) looked better than adding 1

## Check normality after transformation
shapiro.test(WSL$log10_LA_1)
# W = 0.97131, p-value = 0.1691 *
shapiro.test(log10(WSL$NovLeafArea))
# W = 0.97364, p-value = 0.2192 *
## Both conform to normality now but which looks better:

hist(WSL$log10_LA_1)
hist(log10(WSL$NovLeafArea))## log10 looked better than adding 1 for this pool too

## Check normality after transformation
shapiro.test(sqrt(x$NovLeafArea))
# W = 0.99534, p-value = 0.5156 *
hist(sqrt(x$NovLeafArea))

########################
## Correlation analysis for raw and transformed data. Find adjusted R squared values and correponding p-values for each dataset:
########################
summary(lm(SSS~NovLeafArea, data = WLL))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = WLL)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -86.603 -23.310  -7.068  26.606  93.010 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   74.575     22.225   3.355  0.00140 **
# NovLeafArea    4.645      1.504   3.088  0.00309 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 39.25 on 58 degrees of freedom
# Multiple R-squared:  0.1412,	Adjusted R-squared:  0.1264 
# F-statistic: 9.534 on 1 and 58 DF,  p-value: 0.003093

summary(lm(SSS~NovLeafArea, data = WSL))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = WSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -55.966 -20.298  -1.412  18.666  65.833 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  81.7817    10.7348   7.618 2.68e-10 ***
# NovLeafArea   3.9231     0.7223   5.431 1.15e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 27.83 on 58 degrees of freedom
# Multiple R-squared:  0.3371,	Adjusted R-squared:  0.3257 
# F-statistic:  29.5 on 1 and 58 DF,  p-value: 1.154e-06

summary(lm(SSS~NovLeafArea, data = WUS))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = WUS)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -68.043 -27.775  -5.318  27.761  65.679 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  68.5431    14.5482   4.711 1.71e-05 ***
# NovLeafArea   3.6533     0.8027   4.551 2.99e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 34.33 on 55 degrees of freedom
# Multiple R-squared:  0.2736,	Adjusted R-squared:  0.2604 
# F-statistic: 20.71 on 1 and 55 DF,  p-value: 2.991e-05

summary(lm(SSS~NovLeafArea, data = FLL))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = FLL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -76.078 -22.997   1.075  21.959  60.996 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  108.010     15.159   7.125 1.81e-09 ***
# NovLeafArea    2.155      0.811   2.657   0.0102 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 31.98 on 58 degrees of freedom
# Multiple R-squared:  0.1085,	Adjusted R-squared:  0.09315 
# F-statistic:  7.06 on 1 and 58 DF,  p-value: 0.01016

summary(lm(SSS~NovLeafArea, data = FSL))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = FSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -93.975 -38.997  -4.174  40.239 123.186 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 111.5707    21.2200   5.258 2.19e-06 ***
# NovLeafArea   0.6161     1.7501   0.352    0.726    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 51.98 on 58 degrees of freedom
# Multiple R-squared:  0.002132,	Adjusted R-squared:  -0.01507 
# F-statistic: 0.1239 on 1 and 58 DF,  p-value: 0.7261

summary(lm(SSS~NovLeafArea, data = x))
# Call:
#   lm(formula = SSS ~ NovLeafArea, data = x)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -107.126  -25.035   -2.038   25.962  117.594 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  89.8091     6.8208  13.167  < 2e-16 ***
# NovLeafArea   3.0156     0.4297   7.019 1.54e-11 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 38.58 on 295 degrees of freedom
# Multiple R-squared:  0.1431,	Adjusted R-squared:  0.1402 
# F-statistic: 49.26 on 1 and 295 DF,  p-value: 1.542e-11

## Transformed data:
summary(lm(SSS~log10(NovLeafArea), data = WSL))
# Call:
#   lm(formula = SSS ~ log10(NovLeafArea), data = WSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -54.357 -16.950  -2.418  16.822  76.167 
# 
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -10.58      26.27  -0.403    0.689    
# log10(NovLeafArea)   131.50      23.24   5.658 4.95e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 27.44 on 58 degrees of freedom
# Multiple R-squared:  0.3557,	Adjusted R-squared:  0.3446 
# F-statistic: 32.02 on 1 and 58 DF,  p-value: 4.951e-07

summary(lm(SSS~log10(NovLeafArea), data = FSL))
# Call:
#   lm(formula = SSS ~ log10(NovLeafArea), data = FSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -95.059 -39.750  -2.653  40.929 122.680 
# 
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           93.97      50.33   1.867   0.0669 .
# log10(NovLeafArea)    23.78      48.04   0.495   0.6226  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 51.92 on 58 degrees of freedom
# Multiple R-squared:  0.004205,	Adjusted R-squared:  -0.01296 
# F-statistic: 0.2449 on 1 and 58 DF,  p-value: 0.6226

summary(lm(SSS~sqrt(NovLeafArea), data = x))
# Call:
#   lm(formula = SSS ~ sqrt(NovLeafArea), data = x)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -108.478  -26.102   -1.122   24.713  117.299 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         43.127     12.889   3.346 0.000926 ***
# sqrt(NovLeafArea)   24.095      3.328   7.239 3.93e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 38.4 on 295 degrees of freedom
# Multiple R-squared:  0.1509,	Adjusted R-squared:  0.148 
# F-statistic: 52.41 on 1 and 295 DF,  p-value: 3.931e-12

########################
## Plot correlation analysis graphs for raw and transformed data. Use the above adjusted R squared values and correponding p-values for each dataset.
########################
WLL_gg <- ggscatter(WLL, x = "NovLeafArea", y = "SSS", title = expression(bold(A)~~"WNZLL"),
                    add = "reg.line", fullrange = TRUE, conf.int = FALSE, shape = 21, fill = "#e6194b", size = 3, ylim = c(0, 250), xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=30, y=38, label="atop(r ^ 2 == 0.13, italic(p) == 0.0031)", parse = TRUE)+
  theme_bw()

WSL_gg <- ggscatter(WSL, x = "NovLeafArea", y = "SSS",  title = expression(bold(B)~~"WNZSL"),
                    add = "reg.line", conf.int = FALSE, shape = 22, fill = "#ffe119", size = 3, ylim = c(0, 250),xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=30, y=38, label="atop(r ^ 2 == 0.33, italic(p) == 1.2e-06)", parse = TRUE)+
  theme_bw()

WUS_gg <- ggscatter(WUS, x = "NovLeafArea", y = "SSS",  title = expression(bold(C)~~"WUSLL"),
                    add = "reg.line", conf.int = FALSE, shape = 23, fill = "#3cb44b", size = 3, ylim = c(0, 250),xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)')) +
  annotate("text", x=30, y=38, label="atop(r ^ 2 == 0.26, italic(p) == 3.0e-05)", parse = TRUE)+
  theme_bw()

FLL_gg <- ggscatter(FLL, x = "NovLeafArea", y = "SSS",  title = expression(bold(D)~~"FNZLL"),
                    add = "reg.line", conf.int = FALSE, shape = 25, fill = "#4363d8", size = 3, ylim = c(0, 250),xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)')) +
  annotate("text", x=30, y=38, label="atop(r ^ 2 == 0.09, italic(p) == 0.01)", parse = TRUE)+
  theme_bw()

FSL_gg <- ggscatter(FSL, x = "NovLeafArea", y = "SSS", title = expression(bold(E)~~"FNZSL"),
                    add = "reg.line", conf.int = FALSE, shape = 24, fill = "#42D4F4", size = 3, ylim = c(0, 250),xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=30, y=38, label="atop(r ^ 2 == -0.015, italic(p) == 0.73)", parse = TRUE)+
  theme_bw()

ALL_gg <- ggscatter(x, x = "NovLeafArea", y = "SSS",  title = expression(bold(F)~~"All pools"),
                    add = "reg.line", conf.int = FALSE, shape = 8, size = 3, ylim = c(0, 250),xlim = c(4 , 36),
                    #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                    xlab = expression("Leaf area" ~~ (cm^2)), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=30, y=38, label="atop(r ^ 2 == 0.14, italic(p) == 1.5e-11)", parse = TRUE)+
  theme_bw()

## Transformed data:
WSL_T_gg <- ggscatter(WSL, x = "log10_LA", y = "SSS",  title = expression(bold(G)~~"WNZSL - Transformed"),
                      add = "reg.line", conf.int = FALSE, shape = 22, fill = "#ffe119", size = 3, ylim = c(0, 250),xlim = c(0.7, 1.6),
                      #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                      xlab = expression(paste(log[10] (leaf~area~cm^2))), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=1.45, y=38, label="atop(r ^ 2 == 0.34, italic(p) == 5.0e-07)", parse = TRUE)+
  theme_bw()

FSL_T_gg <- ggscatter(FSL, x = "log10_LA", y = "SSS",  title = expression(bold(H)~~"FNZSL - Transformed"),
                      add = "reg.line", conf.int = FALSE, shape = 24, fill = "#42D4F4", size = 3, ylim = c(0, 250),xlim = c(0.7, 1.6),
                      #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                      xlab = expression(paste(log[10] (leaf~area~cm^2))), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=1.45, y=38, label="atop(r ^ 2 == -0.013, italic(p) == 0.62)", parse = TRUE)+
  theme_bw()

ALL_T_gg <- ggscatter(x, x = "sqrt_LA", y = "SSS",  title = expression(bold(I)~~"All pools - Transformed"),
                      add = "reg.line", conf.int = FALSE, shape = 8, size = 3, ylim = c(0, 250),
                      #cor.coef = TRUE, cor.method = "pearson",cor.coeff.args = list(label.sep = "\n"),
                      xlab = expression(paste(sqrt((leaf~area~cm^2)))), ylab = expression("SSS  (g kg"^"-1"*' DM)'))+
  annotate("text", x=5.3, y=38, label="atop(r ^ 2 == 0.15, italic(p) == 3.9e-12)", parse = TRUE)+
  theme_bw()

## Plot all graphs onto a 3x3 grid
ggarrange(WLL_gg,FLL_gg, WSL_T_gg,WSL_gg,FSL_gg,FSL_T_gg,WUS_gg,ALL_gg,ALL_T_gg, ncol = 3, nrow = 3, font.label = list(family = "sans"))
## Saved image as "Correlation_combined_plots_SSS_g_kg_DM.TIFF" W = 883 H = 796 (Supplementary Figure 4)