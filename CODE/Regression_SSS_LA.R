library(readr)
library(ggplot2)
library(ggpubr)

########################
## Input data
########################
WLL <- read.table("WNZLL-CORR_pop.txt", header = TRUE)
WSL <- read.table("WNZSL-CORR_pop.txt", header = TRUE)
WUS <- read.table("WUSLL-CORR_pop.txt", header = TRUE)
FLL <- read.table("FNZLL-CORR_pop.txt", header = TRUE)
FSL <- read.table("FNZSL-CORR_pop.txt", header = TRUE)
## Data is split into each pool. 

########################
## Use Parent generation as the baseline for comparison for regression
########################
WLL$Generation = factor(WLL$Generation,
                        levels=unique(WLL$Generation))
WLL$Population = factor(WLL$Population,
                        levels=unique(WLL$Population))
WLL$Population_base_parent = factor(WLL$Population,levels(WLL$Population)[c(3,1,2,4,5)]) ## Makes parent Population baseline for comparison

WSL$Generation = factor(WSL$Generation,
                        levels=unique(WSL$Generation))
WSL$Population = factor(WSL$Population,
                        levels=unique(WSL$Population))
WSL$Population_base_parent = factor(WSL$Population,levels(WSL$Population)[c(3,1,2,4,5)])

WUS$Generation = factor(WUS$Generation,
                        levels=unique(WUS$Generation))
WUS$Population = factor(WUS$Population,
                        levels=unique(WUS$Population))
WUS$Population_base_parent = factor(WUS$Population,levels(WUS$Population)[c(3,1,2,4,5)])

FLL$Generation = factor(FLL$Generation,
                        levels=unique(FLL$Generation))
FLL$Population = factor(FLL$Population,
                        levels=unique(FLL$Population))
FLL$Population_base_parent = factor(FLL$Population,levels(FLL$Population)[c(3,1,2,4,5)])

FSL$Generation = factor(FSL$Generation,
                        levels=unique(FSL$Generation))
FSL$Population = factor(FSL$Population,
                        levels=unique(FSL$Population))
FSL$Population_base_parent = factor(FSL$Population,levels(FSL$Population)[c(3,1,2,4,5)])

########################
## Regression analysis
########################
summary(lm(SSS~0+Population_base_parent*NovLeafArea, data = WLL))
# Call:
#   lm(formula = SSS ~ 0 + Population * NovLeafArea, data = WLL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -4.9396 -1.5113 -0.5288  1.8166  5.5749
# 
# Coefficients:
#                                Estimate Std. Error t value Pr(>|t|)
# PopulationParent               16.11830    2.45889   6.555 2.98e-08 ***
# PopulationLow-End              15.37634    6.44806   2.385  0.02093 *
# PopulationLow-Mid               7.64951    3.96107   1.931  0.05914 .
# PopulationHigh-Mid             12.15628    4.26132   2.853  0.00629 **
# PopulationHigh-End             21.80625    4.23823   5.145 4.49e-06 ***
# NovLeafArea                    -0.19243    0.17091  -1.126  0.26560
# PopulationLow-End:NovLeafArea  -0.31102    0.55459  -0.561  0.57743
# PopulationLow-Mid:NovLeafArea   0.53055    0.35140   1.510  0.13739
# PopulationHigh-Mid:NovLeafArea  0.45328    0.29442   1.540  0.12997
# PopulationHigh-End:NovLeafArea  0.08466    0.31194   0.271  0.78719
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2.433 on 50 degrees of freedom
# Multiple R-squared:  0.9773,	Adjusted R-squared:  0.9728
# F-statistic: 215.3 on 10 and 50 DF,  p-value: < 2.2e-16

par(mfrow=c(2,2))
plot(lm(SSS~0+Population_base_parent*NovLeafArea, data = WLL)) ## Save as "WLL_residuals.tiff" W = 883 H = 796
dev.off()
plot(SSS~NovLeafArea,data=WLL,pch=16,col=WLL$Population_base_parent)+
  abline(16.11830,-0.19243,col=1)+
  abline(15.37634,(-0.19243+-0.31102),col=2)+
  abline(7.64951,(-0.19243+0.53055),col=3)+
  abline(12.15628,(-0.19243+0.45328),col=4)+
  abline(21.80625,(-0.19243+0.08466),col=5)

summary(lm(SSS~0+Population_base_parent*NovLeafArea, data = WSL))
# Call:
#   lm(formula = SSS ~ 0 + Population * NovLeafArea, data = WSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.9632 -1.1883 -0.0016  0.9157  6.0826 
# 
# Coefficients:
#                                  Estimate Std. Error t value Pr(>|t|)    
#   PopulationParent               10.22798    1.66786   6.132 1.36e-07 ***
#   PopulationLow-End               9.35455    2.80269   3.338 0.001601 ** 
#   PopulationLow-Mid               4.51093    2.82921   1.594 0.117146    
#   PopulationHigh-Mid             12.56957    2.21432   5.676 6.94e-07 ***
#   PopulationHigh-End             18.17037    4.85781   3.740 0.000474 ***
#   NovLeafArea                     0.15870    0.13420   1.183 0.242563    
#   PopulationLow-End:NovLeafArea   0.04609    0.25677   0.180 0.858269    
#   PopulationLow-Mid:NovLeafArea   0.44488    0.25685   1.732 0.089425 .  
#   PopulationHigh-Mid:NovLeafArea  0.03563    0.17467   0.204 0.839177    
#   PopulationHigh-End:NovLeafArea -0.18582    0.31380  -0.592 0.556407    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2.335 on 50 degrees of freedom
# Multiple R-squared:  0.9771,	Adjusted R-squared:  0.9725 
# F-statistic: 213.2 on 10 and 50 DF,  p-value: < 2.2e-16

par(mfrow=c(2,2))
plot(lm(SSS~0+Population_base_parent*NovLeafArea, data = WSL)) ## Save as "WSL_residuals.tiff" W = 883 H = 796
dev.off()
plot(SSS~NovLeafArea,data=WSL,pch=16,col=WSL$Population_base_parent)+
  abline(10.22798,0.15870,col=1)+
  abline(9.35455,(0.15870+0.04609),col=2)+
  abline(4.51093,(0.15870+0.44488),col=3)+
  abline(12.56957,(0.15870+0.03563),col=4)+
  abline(18.17037,(0.15870+-0.18582),col=5)

summary(lm(SSS~0+Population_base_parent*NovLeafArea, data = WUS))
# Call:
#   lm(formula = SSS ~ 0 + Population * NovLeafArea, data = WUS)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -8.741 -1.422 -0.054  2.157  6.214 
# 
# Coefficients:
#                                  Estimate Std. Error t value Pr(>|t|)    
#   PopulationParent                6.84153    2.25114   3.039 0.003868 ** 
#   PopulationLow-End               2.51739    3.29967   0.763 0.449323    
#   PopulationLow-Mid               6.01122    4.03068   1.491 0.142547    
#   PopulationHigh-Mid             14.35131    3.81415   3.763 0.000466 ***
#   PopulationHigh-End             10.04162    3.39238   2.960 0.004807 ** 
#   NovLeafArea                     0.36622    0.12099   3.027 0.004000 ** 
#   PopulationLow-End:NovLeafArea   0.19323    0.25563   0.756 0.453474    
#   PopulationLow-Mid:NovLeafArea  -0.07915    0.26480  -0.299 0.766320    
#   PopulationHigh-Mid:NovLeafArea -0.33649    0.22372  -1.504 0.139260    
#   PopulationHigh-End:NovLeafArea -0.05133    0.22165  -0.232 0.817868    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.199 on 47 degrees of freedom
# Multiple R-squared:  0.9552,	Adjusted R-squared:  0.9457 
# F-statistic: 100.2 on 10 and 47 DF,  p-value: < 2.2e-16

par(mfrow=c(2,2))
plot(lm(SSS~0+Population_base_parent*NovLeafArea, data = WUS)) ## Save as "WUS_residuals.tiff" W = 883 H = 796
dev.off()
plot(SSS~NovLeafArea,data=WUS,pch=16,col=WUS$Population_base_parent)+
  abline(6.84153,0.36622,col=1)+
  abline(2.51739,(0.36622+0.19323),col=2)+
  abline(6.01122,(0.36622+-0.07915),col=3)+
  abline(14.35131,(0.36622+-0.33649),col=4)+
  abline(10.04162,(0.36622+-0.05133),col=5)

summary(lm(SSS~0+Population_base_parent*NovLeafArea, data = FLL))
# Call:
#   lm(formula = SSS ~ 0 + Population * NovLeafArea, data = FLL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.8514 -1.6785 -0.1038  2.1082  5.7416 
# 
# Coefficients:
#                                  Estimate Std. Error t value Pr(>|t|)    
#   PopulationParent                18.6911     3.2488   5.753 5.28e-07 ***
#   PopulationLow-End                7.1071     3.6163   1.965  0.05494 .  
#   PopulationLow-Mid                9.9789     3.3738   2.958  0.00472 ** 
#   PopulationHigh-Mid              13.4225     4.3751   3.068  0.00348 ** 
#   PopulationHigh-End              13.5778     2.7537   4.931 9.44e-06 ***
#   NovLeafArea                     -0.2251     0.1627  -1.383  0.17278    
#   PopulationLow-End:NovLeafArea    0.6714     0.3132   2.143  0.03696 *  
#   PopulationLow-Mid:NovLeafArea    0.3564     0.2568   1.388  0.17132    
#   PopulationHigh-Mid:NovLeafArea   0.3895     0.2847   1.368  0.17752    
#   PopulationHigh-End:NovLeafArea   0.4363     0.2058   2.120  0.03902 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2.71 on 50 degrees of freedom
# Multiple R-squared:  0.973,	Adjusted R-squared:  0.9676 
# F-statistic:   180 on 10 and 50 DF,  p-value: < 2.2e-16

par(mfrow=c(2,2))
plot(lm(SSS~0+Population_base_parent*NovLeafArea, data = FLL)) ## Save as "FLL_residuals.tiff" W = 883 H = 796
dev.off()
plot(SSS~NovLeafArea,data=FLL,pch=16,col=FLL$Population_base_parent)+
  abline(18.6911,-0.2251,col=1)+
  abline(7.1071,(-0.2251+0.6714),col=2)+
  abline(9.9789,(-0.2251+0.3564),col=3)+
  abline(13.4225,(-0.2251+0.3895),col=4)+
  abline(13.5778,(-0.2251+0.4363),col=5)

summary(lm(SSS~0+Population_base_parent*NovLeafArea, data = FSL))
# Call:
#   lm(formula = SSS ~ 0 + Population * NovLeafArea, data = FSL)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6.0664 -1.8003  0.0879  1.7507  5.4921 
# 
# Coefficients:
#                                  Estimate Std. Error t value Pr(>|t|)    
#   PopulationParent                 7.1690     1.8647   3.845 0.000342 ***
#   PopulationLow-End                6.6783     3.2428   2.059 0.044679 *  
#   PopulationLow-Mid                8.9937     6.8359   1.316 0.194288    
#   PopulationHigh-Mid              11.5395     4.5589   2.531 0.014558 *  
#   PopulationHigh-End              21.1143     3.9724   5.315 2.48e-06 ***
#   NovLeafArea                      0.2340     0.1367   1.712 0.093069 .  
#   PopulationLow-End:NovLeafArea   -0.1651     0.3456  -0.478 0.634938    
#   PopulationLow-Mid:NovLeafArea   -0.2527     0.5593  -0.452 0.653438    
#   PopulationHigh-Mid:NovLeafArea   0.1401     0.4164   0.336 0.737975    
#   PopulationHigh-End:NovLeafArea  -0.4339     0.3991  -1.087 0.282242    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.153 on 50 degrees of freedom
# Multiple R-squared:  0.9504,	Adjusted R-squared:  0.9405 
# F-statistic: 95.77 on 10 and 50 DF,  p-value: < 2.2e-16

par(mfrow=c(2,2))
plot(lm(SSS~0+Population_base_parent*NovLeafArea, data = FSL)) ## Save as "FSL_residuals.tiff" W = 883 H = 796
dev.off()
plot(SSS~NovLeafArea,data=FSL,pch=16,col=FSL$Population_base_parent)+
  abline(7.1690,0.2340,col=1)+
  abline(6.6783,(0.2340+-0.1651),col=2)+
  abline(8.9937,(0.2340+-0.2527),col=3)+
  abline(11.5395,(0.2340+0.1401),col=4)+
  abline(21.1143,(0.2340+-0.4339),col=5)

## Output above compiled and used for Supplementary Table 3

########################
## Plot regression lines for each Population per pool. Use same colouring as PCAdapt score plots
########################
## order of Populations = Low-End, Low-Mid, Parent, High-Mid and High-End
#scale_color_manual(values = c("#80cdc1","#018571","#3e3e3d","#a6611a","#c2a5cf"),aesthetics = c("colour","fill")) ## WLL
#scale_color_manual(values = c("#a6dba0","#008837","#3e3e3d","#d01c8b","#f1b6da"),aesthetics = c("colour","fill")) ## WSL
#scale_color_manual(values = c("#b2abd2","#5e3c99","#3e3e3d","#e66101","#fdb863"),aesthetics = c("colour","fill")) ## WUS
#scale_color_manual(values = c("#92c5de","#0571b0","#3e3e3d","#ca0020","#f4a582"),aesthetics = c("colour","fill")) ## FLL
#scale_color_manual(values = c("#b8e186","#4dac26","#3e3e3d","#7b3294","#c2a5cf"),aesthetics = c("colour","fill")) ## FSL

WLL_gen <- ggplot(WLL, aes(NovLeafArea, SSS, shape = Population, colour = Population))+
  geom_smooth(method = "lm", fullrange = T, se = F)+
  geom_point(size = 3)+
  theme_bw()+
  xlab(expression(paste(Leaf~area~~(cm^2))))+
  ylab(expression("SSS  (g kg"^"-1"*' DM)'))+
  ggtitle(expression(bold(A)~~"WNZLL"))+
  expand_limits(y = c(0, 250))+
  expand_limits(x = c(0, 35))+
  scale_shape_manual(values = c(17,18,3,19,15))+
  scale_color_manual(values = c("#80cdc1","#018571","#3e3e3d","#a6611a","#dfc27d"),aesthetics = c("colour","fill"))+
  theme(legend.key.size = unit(.5, "cm"))#+
#geom_abline(intercept = 0, slope = 1,linetype= "dashed")

WSL_gen <- ggplot(WSL, aes(NovLeafArea, SSS, shape = Population, colour = Population))+
  geom_smooth(method = "lm", fullrange = T, se = F)+
  geom_point(size = 3)+
  theme_bw()+
  xlab(expression(paste(Leaf~area~~(cm^2))))+
  ylab(expression("SSS  (g kg"^"-1"*' DM)'))+
  ggtitle(expression(bold(B)~~"WNZSL"))+
  expand_limits(y = c(0, 250))+
  expand_limits(x = c(0, 35))+
  scale_y_continuous(name = waiver(), breaks = c(0,50,100,150,200,250))+
  scale_shape_manual(values = c(17,18,3,19,15))+
  scale_color_manual(values = c("#a6dba0","#008837","#3e3e3d","#d01c8b","#f1b6da"),aesthetics = c("colour","fill"))+
  theme(legend.key.size = unit(.5, "cm"))#+
#geom_abline(intercept = 0, slope = 1,linetype= "dashed")

WUS_gen <- ggplot(WUS, aes(NovLeafArea, SSS, shape = Population, colour = Population))+
  geom_smooth(method = "lm", fullrange = T, se = F)+
  geom_point(size = 3)+
  theme_bw()+
  xlab(expression(paste(Leaf~area~~(cm^2))))+
  ylab(expression("SSS  (g kg"^"-1"*' DM)'))+
  ggtitle(expression(bold(C)~~"WUSLL"))+
  expand_limits(y = c(0, 250))+
  expand_limits(x = c(0, 35))+
  scale_shape_manual(values = c(17,18,3,19,15))+
  scale_color_manual(values = c("#b2abd2","#5e3c99","#3e3e3d","#e66101","#fdb863"),aesthetics = c("colour","fill"))+
  theme(legend.key.size = unit(.5, "cm"))#+
#geom_abline(intercept = 0, slope = 1,linetype= "dashed")

FLL_gen <- ggplot(FLL, aes(NovLeafArea, SSS, shape = Population, colour = Population))+
  geom_smooth(method = "lm", fullrange = T, se = F)+
  geom_point(size = 3)+
  theme_bw()+
  xlab(expression(paste(Leaf~area~~(cm^2))))+
  ylab(expression("SSS  (g kg"^"-1"*' DM)'))+
  ggtitle(expression(bold(D)~~"FNZLL"))+
  expand_limits(y = c(0, 250))+
  expand_limits(x = c(0, 35))+
  scale_shape_manual(values = c(17,18,3,19,15))+
  scale_color_manual(values = c("#92c5de","#0571b0","#3e3e3d","#ca0020","#f4a582"),aesthetics = c("colour","fill"))+
  theme(legend.key.size = unit(.5, "cm"))#+
#geom_abline(intercept = 0, slope = 1,linetype= "dashed")

FSL_gen <- ggplot(FSL, aes(NovLeafArea, SSS, shape = Population, colour = Population))+
  geom_smooth(method = "lm", fullrange = T, se = F)+
  geom_point(size = 3)+
  theme_bw()+
  xlab(expression(paste(Leaf~area~~(cm^2))))+
  ylab(expression("SSS  (g kg"^"-1"*' DM)'))+
  ggtitle(expression(bold(E)~~"FNZSL"))+
  expand_limits(y = c(0, 250))+
  expand_limits(x = c(0, 35))+
  scale_shape_manual(values = c(17,18,3,19,15))+
  scale_color_manual(values = c("#b8e186","#4dac26","#3e3e3d","#7b3294","#c2a5cf"),aesthetics = c("colour","fill"))+
  theme(legend.key.size = unit(.5, "cm"))#+
#geom_abline(intercept = 0, slope = 1,linetype= "dashed")

ggarrange(WLL_gen, FLL_gen, WSL_gen,FSL_gen,WUS_gen, ncol = 2, nrow = 3,font.label = list(family = "sans"))
## Save as "pool_regression_45_line_pop_SSS_g_kg_DM.tiff" W = 800 H = 900 (Supplementary Figure 5)