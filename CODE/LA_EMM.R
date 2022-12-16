library(here)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(emmeans)
library(stringr)
library(tibble)
library(lme4)
library(lmerTest)
library(pbkrtest)

########################
## Import data
########################
LA_REML <- read.table(file = here::here("Leaf_area_1788_samples.txt"), header = TRUE) %>%
  separate(Generation, into = c("Trt", "Time"), "-") %>%
  mutate(Trt = factor(Trt, levels = c("Parent", "Low", "High")),
         Time = factor(ifelse(Trt == "Parent", "Start", Time), 
                       levels = c("Start", "Mid", "End")),
         EID = ifelse(Pool %in% c("WNZLL", "WNZSL", "WUSLL"), "W", "F"),
         Parent = factor(ifelse(Trt == "Parent", "Y", "N"),
                         levels = c("Y", "N"))) %>%
  transform(TrtC = interaction(Pool, Time, Trt))

########################
## Explore data
########################
LA_REML %>%
  group_by(Pool, Trt, Time) %>% 
  summarise(mean_LA = mean(Leaf_Area)) %>%
  ggplot(aes(Time, mean_LA, col = Trt)) +
  geom_point() + facet_wrap(~Pool)

## Three way interactions hard to interpret, so it can be helpful to do exploratory plots and models on subsets of the data.
## Plot means for each group. Low WSC Treatment leaf areas are on average smaller than High WSC for all Pools except FNZSL. 
## The pattern compared to Parent changes for each Pool. The time pattern changes for each Pool. 
## Evidence that a 3 way interaction is necessary. Remember these are raw means that have not been adjusted for row, column, block etc effects.

########################
## Remove the Row:Column effects and average the leaf area for each plant:
## Model on mean leaf area for each plant
########################
LA_REML_means <- LA_REML %>% 
  group_by(Pool, Trt, Time, TrtC, Block, Row, Column) %>% 
  summarise(Leaf_Area = mean(Leaf_Area)) %>% ungroup()
mod <- lmer(sqrt(Leaf_Area) ~ TrtC + 
              (1 | Block) + (1 | Block:Column) + (1 | Block:Row),
            # average leaf area within each plant, 
            # now 1 measurement per plot, so don't 
            # need GenotypeID or Row:Column random effects
            data = LA_REML_means)

summary(mod)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: sqrt(Leaf_Area) ~ TrtC + (1 | Block) + (1 | Block:Column) + (1 |      Block:Row)
#    Data: LA_REML_means
# 
# REML criterion at convergence: 703.4
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.3066 -0.5454 -0.0176  0.5685  3.2741 
# 
# Random effects:
#  Groups       Name        Variance Std.Dev.
#  Block:Row    (Intercept) 0.027533 0.16593 
#  Block:Column (Intercept) 0.013245 0.11509 
#  Block        (Intercept) 0.009913 0.09956 
#  Residual                 0.229759 0.47933 
# Number of obs: 447, groups:  Block:Row, 90; Block:Column, 15; Block, 3
# 
# Fixed effects:
#                         Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)              4.31029    0.11289  12.80769  38.180 1.41e-14 ***
# TrtCFNZSL.Start.Parent  -0.95316    0.12857 395.52312  -7.414 7.54e-13 ***
# TrtCWNZLL.Start.Parent  -0.68342    0.12823 392.22825  -5.329 1.67e-07 ***
# TrtCWNZSL.Start.Parent  -1.00238    0.12803 389.56185  -7.829 4.69e-14 ***
# TrtCWUSLL.Start.Parent  -0.17842    0.12838 393.40050  -1.390 0.165369    
# TrtCFNZLL.Mid.Low       -0.42566    0.15753 397.07738  -2.702 0.007184 ** 
# TrtCFNZSL.Mid.Low       -0.95339    0.15626 384.46401  -6.101 2.57e-09 ***
# TrtCWNZLL.Mid.Low       -0.60259    0.15753 397.17151  -3.825 0.000152 ***
# TrtCWNZSL.Mid.Low       -0.73929    0.15652 386.83130  -4.723 3.25e-06 ***
# TrtCWUSLL.Mid.Low       -0.38341    0.16081 392.12423  -2.384 0.017591 *  
# TrtCFNZLL.End.Low       -0.67551    0.15752 397.04904  -4.288 2.26e-05 ***
# TrtCFNZSL.End.Low       -1.18387    0.15624 383.93768  -7.577 2.67e-13 ***
# TrtCWNZLL.End.Low       -0.80054    0.15752 397.08036  -5.082 5.76e-07 ***
# TrtCWNZSL.End.Low       -0.91906    0.15751 396.91629  -5.835 1.12e-08 ***
# TrtCWUSLL.End.Low       -0.75177    0.16545 396.89495  -4.544 7.34e-06 ***
# TrtCFNZLL.Mid.High       0.08335    0.15596 381.01392   0.534 0.593351    
# TrtCFNZSL.Mid.High      -1.06664    0.15751 396.91430  -6.772 4.60e-11 ***
# TrtCWNZLL.Mid.High      -0.32374    0.15806 400.41611  -2.048 0.041183 *  
# TrtCWNZSL.Mid.High      -0.20816    0.15753 397.12337  -1.321 0.187114    
# TrtCWUSLL.Mid.High      -0.04268    0.15752 397.01668  -0.271 0.786557    
# TrtCFNZLL.End.High       0.13654    0.15807 400.61523   0.864 0.388221    
# TrtCFNZSL.End.High      -1.22773    0.15752 397.09239  -7.794 5.75e-14 ***
# TrtCWNZLL.End.High      -0.30202    0.15708 392.04266  -1.923 0.055239 .  
# TrtCWNZSL.End.High      -0.25839    0.15653 386.96433  -1.651 0.099598 .  
# TrtCWUSLL.End.High      -0.28528    0.15833 401.93142  -1.802 0.072321 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation matrix not shown by default, as p = 25 > 12.
# Use print(x, correlation=TRUE)  or
#     vcov(x)        if you need it
## Some evidence of differing results for each Pool, so final model should have 3 way interaction

########################
## View residuals
########################
LA_REML_means <- LA_REML_means %>%
  mutate(resids_pearson = resid(mod, type = "pearson"),
         fitteds = fitted(mod))
p <- list()

p[[length(p) + 1]] <- ggplot(LA_REML_means,aes(fitteds, resids_pearson)) + 
  geom_point(aes(color = Trt), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c("#ffe119", "#4363d8", "#3cb44b"),name = "Treatment", labels = c("None", "Low", "High"))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Treatment"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(fitteds, resids_pearson, col = Time)) + 
  geom_point(aes(color = Time), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c("#ffe119", "#4363d8", "#3cb44b"),name = "Generation", labels = c("Parent", "Mid", "End"))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Generation"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- ggplot(LA_REML_means, aes(fitteds, resids_pearson)) + 
  geom_point(aes(color = Pool), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c("#4363d8", "#42d4f4", "#e6194b","#ffe119", "#3cb44b"))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Pool"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(fitteds, resids_pearson, col = paste0(Trt, Time))) + 
  geom_point(aes(color = paste0(Trt, Time)), size = 2, alpha = 0.8) + 
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Treatment and Generation"))+
  scale_color_manual(values = c("#4363d8", "#42d4f4", "#e6194b","#ffe119", "#3cb44b"), name = "Trt & Gen", labels = c("High & End", "High & Mid", "Low & End", "Low & Mid", "None & Parent"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

grid.arrange(grobs = p, ncol = 1)
## Save plot as "Treatments_residual_Leaf_area_Gen_edited.tiff" W = 650 H = 850 (Supplementary Figure 1, right)

########################
## Extract random effects from model
########################
block_ranefs <-  ranef(mod)$Block %>% rownames_to_column("Block") %>% as_tibble() %>%
  rename("re_value_block" = "(Intercept)")
blockCol_ranefs <-  ranef(mod)$`Block:Column` %>% rownames_to_column("Block:Column") %>% as_tibble() %>%
  rename("re_value_column" = "(Intercept)") %>%
  separate('Block:Column', into = c("Block", "Column"))
blockRow_ranefs <-  ranef(mod)$`Block:Row` %>% rownames_to_column("Block:Row") %>% as_tibble() %>%
  rename("re_value_row" = "(Intercept)") %>%
  separate('Block:Row', into = c("Block", "Row"))

########################
## Add random effects values into the main data
########################
LA_REML_means <- LA_REML_means %>%
  mutate(resids_pearson = resid(mod, type = "pearson"),
         fitteds = fitted(mod)) %>%
  left_join(block_ranefs %>% mutate_all(~as.numeric(.))) %>%
  left_join(blockCol_ranefs %>% mutate_all(~as.numeric(.))) %>%
  left_join(blockRow_ranefs %>% mutate_all(~as.numeric(.))) %>%
  mutate(total_re = re_value_block + re_value_row + re_value_column )

########################
## Plot
########################
p <- list()

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_block)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Block"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_row)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Row"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_column)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Column"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = total_re)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Total random effects"))

grid.arrange(grobs = p, ncol = 2)
## Save plot as "Random_effects_residual_Leaf_area_edited.tiff" W = 800 H = 800

########################
## Plot spatial residuals
########################
p <- list()

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = resids_pearson)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Residuals"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = paste0(Trt, Time))) + 
  geom_tile() +
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Treatment and Time"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = Trt)) + 
  geom_tile()+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Treatment"))

p[[length(p) + 1]] <- LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = Time)) + 
  geom_tile() +
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Time"))

grid.arrange(grobs = p, ncol = 2)
## Save plot as "Treatments_residual_spatial_Leaf_area_edited.tiff" W = 800 H = 800

## Just need the first plot as there is no pattern for the residuals
LA_REML_means %>% 
  ggplot(aes(Column, Row, fill = resids_pearson)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Residuals for Leaf Area"))
## Save plot as "Residual_spatial_LA_edited.tiff" W = 700 H = 800

########################
## Building the model to produce estimated LA means per pop and comparisons to the Parent baseline
########################
em <- emmeans(mod, "TrtC", type = "response") ## This gives us the non-transformed values at the end
plot(em)
em_levels <- levels(em)$TrtC
adj_type <- "mvt" ## Multiple adjustment when combining all constrasts
adj_type_sub <- "none" ## Multiple adjustment for sub-constrasts
Time_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
  apply(1, function(i) paste0(i, collapse = ".")) ## Contrast comparing Time-Trt averaging over pools. Create vector with the correct Time-Trt labels by removing the pool labels from the em object.
## You can make a manual vector, but this is risky because it could give mismatched vectors if you change the model and the order of the TrtC labels changes. 
## It is much better to extract the information from the em object as this guarantees the labels will be correct.
em_TimeTrt_ <- add_grouping(em,
                            newname = "Time_Trt",
                            refname = "TrtC",
                            Time_Trt_values) ## Group the TrtC levels by Time_Trt_values. This effectively replaces the TrtC values with the Time_Trt_values.
em_TimeTrt <- emmeans(em_TimeTrt_, trt.vs.ctrl ~ Time_Trt, ref = "Start.Parent", adjust = adj_type_sub, type = "response") # Do the comparison- each Time-Trt combination is tested only against Parent
## You can change this to pairwise comparisons if they are of interest, however, this is far more comparisons and the multiple testing penalty will be much higher.
## So only do this if you are interested in the results.

summary(em_TimeTrt, type = "response", adj = adj_type_sub) ## plots results on transformed scale. p-values wrong.
# $emmeans
# Time_Trt     response    SE   df lower.CL upper.CL
# End.High         15.4 0.691 4.77     13.6     17.2
# End.Low          11.9 0.611 4.91     10.3     13.5
# Mid.High         16.0 0.705 4.79     14.2     17.9
# Mid.Low          13.6 0.652 4.84     12.0     15.4
# Start.Parent     14.0 0.587 2.99     12.2     16.0
# 
# Results are averaged over the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Intervals are back-transformed from the sqrt scale 
# 
# $contrasts
# contrast                estimate     SE  df t.ratio p.value
# End.High - Start.Parent   0.1761 0.0708 398  2.489  0.0132 
# End.Low - Start.Parent   -0.3027 0.0717 402 -4.219  <.0001 
# Mid.High - Start.Parent   0.2519 0.0710 403  3.545  0.0004 
# Mid.Low - Start.Parent   -0.0574 0.0711 397 -0.807  0.4202 
# 
# Results are averaged over the levels of: TrtC 
# Note: contrasts are still on the sqrt scale 
# Degrees-of-freedom method: kenward-roger 

plot(em_TimeTrt, type = "response")
Time_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
  apply(1, function(i) paste0(i, collapse = "."))
## Get list of all Pool IDs
pool_ids <- str_split(em_levels, "\\.", simplify = T)[, 1] %>% unique()
tsts <- lapply(pool_ids, function(pool_id) { # loop over results from each pool
  # browser()
  # print(pool_id)
  i <- which(grepl(pool_id, em_levels))
  # print(i)
  # print(em[i])
  res <- em[i] %>% 
    add_grouping(.,
                 newname = "Time_Trt",
                 refname = "TrtC",
                 paste0(pool_id, ".",
                        c("Start.Parent", "Mid.Low", "End.Low", 
                          "Mid.High", "End.High")))
  # print(res)
  res <- res %>%
    emmeans(., trt.vs.ctrl ~ Time_Trt, ref = paste0(pool_id, ".", "Start.Parent"), adjust = adj_type_sub, type = "response")
  # print(res)
  return(res)
})

## Extract out results by looping over the lists and pulling out the appropriate data
comparisons_TimeTrt_withinPool <- 
  lapply(tsts, function(i) i$contrasts) %>% ## Contrasts holds the comparisons
  do.call("rbind", .) %>% 
  update(adj = adj_type_sub)
em_TimeTrt_withinPool <-
  lapply(tsts, function(i) i$emmeans) %>% ## emmeans holds the estimated means values
  do.call("rbind", .) %>% 
  update(adjust = adj_type_sub) 

print(em_TimeTrt_withinPool) ## Use these emmeans and SE values (Supplementary Table 2)
# Time_Trt           response    SE   df lower.CL upper.CL
# FNZLL.End.High        19.77 1.288 33.0    17.24     22.5
# FNZLL.End.Low         13.21 1.052 33.0    11.16     15.4
# FNZLL.Mid.High        19.30 1.272 33.0    16.80     22.0
# FNZLL.Mid.Low         15.09 1.125 33.0    12.89     17.5
# FNZLL.Start.Parent    18.58 0.975 12.8    16.53     20.7
# FNZSL.End.High         9.50 0.893 33.0     7.77     11.4
# FNZSL.End.Low          9.77 0.905 33.0     8.02     11.7
# FNZSL.Mid.High        10.52 0.939 33.0     8.70     12.5
# FNZSL.Mid.Low         11.27 0.972 33.0     9.38     13.3
# FNZSL.Start.Parent    11.27 0.757 12.6     9.69     13.0
# WNZLL.End.High        16.07 1.161 33.0    13.79     18.5
# WNZLL.End.Low         12.32 1.016 33.0    10.34     14.5
# WNZLL.Mid.High        15.89 1.154 33.0    13.63     18.3
# WNZLL.Mid.Low         13.75 1.074 33.0    11.65     16.0
# WNZLL.Start.Parent    13.15 0.816 12.5    11.44     15.0
# WNZSL.End.High        16.42 1.173 33.0    14.12     18.9
# WNZSL.End.Low         11.50 0.982 33.0     9.59     13.6
# WNZSL.Mid.High        16.83 1.188 33.0    14.50     19.3
# WNZSL.Mid.Low         12.75 1.034 33.0    10.74     14.9
# WNZSL.Start.Parent    10.94 0.749 12.8     9.38     12.6
# WUSLL.End.High        16.20 1.165 33.0    13.92     18.7
# WUSLL.End.Low         12.66 1.093 41.0    10.55     15.0
# WUSLL.Mid.High        18.21 1.236 33.0    15.79     20.8
# WUSLL.Mid.Low         15.42 1.170 36.6    13.14     17.9
# WUSLL.Start.Parent    17.07 0.935 12.8    15.11     19.2
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Intervals are back-transformed from the sqrt scale 

print(comparisons_TimeTrt_withinPool) ## Don't use these p-values
# contrast                             estimate    SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.Start.Parent  0.136540 0.159 400  0.861  0.3900 
# FNZLL.End.Low - FNZLL.Start.Parent  -0.675505 0.158 397 -4.274  <.0001 
# FNZLL.Mid.High - FNZLL.Start.Parent  0.083352 0.156 380  0.533  0.5944 
# FNZLL.Mid.Low - FNZLL.Start.Parent  -0.425663 0.158 397 -2.693  0.0074 
# FNZSL.End.High - FNZSL.Start.Parent -0.274563 0.157 387 -1.748  0.0812 
# FNZSL.End.Low - FNZSL.Start.Parent  -0.230706 0.157 385 -1.472  0.1419 
# FNZSL.Mid.High - FNZSL.Start.Parent -0.113480 0.158 395 -0.718  0.4730 
# FNZSL.Mid.Low - FNZSL.Start.Parent  -0.000226 0.157 391 -0.001  0.9989 
# WNZLL.End.High - WNZLL.Start.Parent  0.381395 0.157 383  2.436  0.0153 
# WNZLL.End.Low - WNZLL.Start.Parent  -0.117119 0.157 386 -0.746  0.4559 
# WNZLL.Mid.High - WNZLL.Start.Parent  0.359674 0.157 383  2.297  0.0221 
# WNZLL.Mid.Low - WNZLL.Start.Parent   0.080830 0.158 392  0.513  0.6081 
# WNZSL.End.High - WNZSL.Start.Parent  0.743992 0.157 384  4.748  <.0001 
# WNZSL.End.Low - WNZSL.Start.Parent   0.083325 0.158 397  0.527  0.5984 
# WNZSL.Mid.High - WNZSL.Start.Parent  0.794219 0.158 397  5.024  <.0001 
# WNZSL.Mid.Low - WNZSL.Start.Parent   0.263097 0.158 392  1.669  0.0959 
# WUSLL.End.High - WUSLL.Start.Parent -0.106856 0.158 397 -0.676  0.4994 
# WUSLL.End.Low - WUSLL.Start.Parent  -0.573351 0.167 403 -3.431  0.0007 
# WUSLL.Mid.High - WUSLL.Start.Parent  0.135741 0.158 397  0.859  0.3910 
# WUSLL.Mid.Low - WUSLL.Start.Parent  -0.204988 0.160 383 -1.278  0.2020 
# 
# Note: contrasts are still on the sqrt scale 
# Degrees-of-freedom method: kenward-roger

plot(em_TimeTrt_withinPool)+
  labs(x = "Estimated marginal mean (LA)", y = "Population") ## Blue lines are confidence intervals for the EMMs

########################
## Get CIs
########################
res_overall <- rbind(
  comparisons_TimeTrt_withinPool, 
  em_TimeTrt$contrasts, 
  # em_Trt$contrasts,
  # em_Time$contrasts, 
  # em_Pool$contrasts,
  # em_ExperimentID$contrasts, 
  adjust = adj_type)

print(res_overall) ## Use these p-values (Table 1)
# contrast                             estimate     SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.Start.Parent  0.136540 0.1587 400  0.861  0.9999 
# FNZLL.End.Low - FNZLL.Start.Parent  -0.675505 0.1581 397 -4.274  0.0006 
# FNZLL.Mid.High - FNZLL.Start.Parent  0.083352 0.1564 380  0.533  1.0000 
# FNZLL.Mid.Low - FNZLL.Start.Parent  -0.425663 0.1581 397 -2.693  0.1473 
# FNZSL.End.High - FNZSL.Start.Parent -0.274563 0.1570 387 -1.748  0.8024 
# FNZSL.End.Low - FNZSL.Start.Parent  -0.230706 0.1567 385 -1.472  0.9434 
# FNZSL.Mid.High - FNZSL.Start.Parent -0.113479 0.1580 395 -0.718  1.0000 
# FNZSL.Mid.Low - FNZSL.Start.Parent  -0.000227 0.1574 391 -0.001  1.0000 
# WNZLL.End.High - WNZLL.Start.Parent  0.381395 0.1566 383  2.436  0.2744 
# WNZLL.End.Low - WNZLL.Start.Parent  -0.117118 0.1569 386 -0.746  1.0000 
# WNZLL.Mid.High - WNZLL.Start.Parent  0.359674 0.1566 383  2.297  0.3675 
# WNZLL.Mid.Low - WNZLL.Start.Parent   0.080830 0.1575 392  0.513  1.0000 
# WNZSL.End.High - WNZSL.Start.Parent  0.743992 0.1567 384  4.748  0.0001 
# WNZSL.End.Low - WNZSL.Start.Parent   0.083325 0.1581 397  0.527  1.0000 
# WNZSL.Mid.High - WNZSL.Start.Parent  0.794218 0.1581 397  5.024  <.0001 
# WNZSL.Mid.Low - WNZSL.Start.Parent   0.263098 0.1576 392  1.669  0.8530 
# WUSLL.End.High - WUSLL.Start.Parent -0.106857 0.1581 397 -0.676  1.0000 
# WUSLL.End.Low - WUSLL.Start.Parent  -0.573351 0.1671 403 -3.431  0.0153 
# WUSLL.Mid.High - WUSLL.Start.Parent  0.135741 0.1581 397  0.859  0.9999 
# WUSLL.Mid.Low - WUSLL.Start.Parent  -0.204988 0.1604 383 -1.278  0.9850 
# End.High - Start.Parent              0.176102 0.0708 398  2.489  0.2440 
# End.Low - Start.Parent              -0.302671 0.0717 402 -4.219  0.0007 
# Mid.High - Start.Parent              0.251901 0.0710 403  3.545  0.0102 
# Mid.Low - Start.Parent              -0.057390 0.0711 397 -0.807  1.0000 
# 
# Results are averaged over some or all of the levels of: TrtC 
# Note: contrasts are still on the sqrt scale 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: mvt method for 24 tests 

## For the differences between the means (Parent to other populations) had to use the backtransformed population means differences.
## Unable to backtransform the contrasts from the sqrt scale. The WSC data showed that the comparisons and differences between means are identical, 
## therefore can just use the differences between the Parent means and the others.
res_overall_ci <- rbind(
  em,
  em_TimeTrt$emmeans,
  # em_Trt$emmeans,
  # em_Time$emmeans,
  # em_Pool$emmeans,
  # em_ExperimentID$emmeans,
  adjust = adj_type)

print(res_overall_ci) ## This takes a while, Use these CIs. The mean and SE are same as first output so can just use these ones below (Supplementary Table 2)
# TrtC               Time_Trt     response    SE    df lower.CL upper.CL
# FNZLL.Start.Parent .               18.58 0.975 12.80    15.16     22.3
# FNZSL.Start.Parent .               11.27 0.757 12.61     8.62     14.3
# WNZLL.Start.Parent .               13.15 0.816 12.52    10.28     16.4
# WNZSL.Start.Parent .               10.94 0.749 12.81     8.36     13.9
# WUSLL.Start.Parent .               17.07 0.935 12.81    13.80     20.7
# FNZLL.Mid.Low      .               15.09 1.125 33.01    11.60     19.0
# FNZSL.Mid.Low      .               11.27 0.972 33.02     8.28     14.7
# WNZLL.Mid.Low      .               13.75 1.074 33.02    10.43     17.5
# WNZSL.Mid.Low      .               12.75 1.034 33.00     9.56     16.4
# WUSLL.Mid.Low      .               15.42 1.170 36.64    11.82     19.5
# FNZLL.End.Low      .               13.21 1.052 33.00     9.96     16.9
# FNZSL.End.Low      .                9.77 0.905 33.00     7.01     13.0
# WNZLL.End.Low      .               12.32 1.016 33.01     9.19     15.9
# WNZSL.End.Low      .               11.50 0.982 33.00     8.48     15.0
# WUSLL.End.Low      .               12.66 1.093 41.05     9.35     16.5
# FNZLL.Mid.High     .               19.30 1.272 32.99    15.32     23.7
# FNZSL.Mid.High     .               10.52 0.939 33.01     7.65     13.9
# WNZLL.Mid.High     .               15.89 1.154 32.99    12.31     19.9
# WNZSL.Mid.High     .               16.83 1.188 33.01    13.13     21.0
# WUSLL.Mid.High     .               18.21 1.236 33.01    14.36     22.5
# FNZLL.End.High     .               19.77 1.288 33.01    15.75     24.3
# FNZSL.End.High     .                9.50 0.893 33.00     6.78     12.7
# WNZLL.End.High     .               16.07 1.161 33.00    12.46     20.1
# WNZSL.End.High     .               16.42 1.173 33.01    12.76     20.5
# WUSLL.End.High     .               16.20 1.165 32.98    12.58     20.3
# .                  End.High        15.39 0.691  4.77    12.20     19.0
# .                  End.Low         11.86 0.611  4.91     9.07     15.0
# .                  Mid.High        15.99 0.705  4.79    12.73     19.6
# .                  Mid.Low         13.61 0.652  4.84    10.61     17.0
# .                  Start.Parent    14.04 0.587  2.99    10.40     18.2
# 
# Results are averaged over some or all of the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: mvt method for 30 estimates 
# Intervals are back-transformed from the sqrt scale 

plot(em_TimeTrt_withinPool, comparisons = TRUE)
## Blue bars are the confidence intervals for the EMMs and the red arrows are the comparisons among them. 
## If an arrow from one mean overlaps an arrow from another group, the difference is not significant at an alpha of 0.05
plot(em_TimeTrt_withinPool, comparisons = TRUE, type = "response")+
  theme_bw()+
  theme_minimal()

########################
## Create text file with means, SE and CI for each population (include Pool, pop and generation info) (LA_EMMs_pop.txt) and plot using the code below for figure       
########################
library(agricolae)
library(Rmisc)
library(doBy)
library(psych)
library(car)
library(readr)
library(ggpubr)
library(units)
library(ggplot2)
library(MASS)
library(lm.beta)
library(lmerTest)

LA_SE <- read.table(file = "LA_EMMs_pop.txt", header = TRUE)
LA_SE$Generation = factor(LA_SE$Generation,
                          levels=unique(LA_SE$Generation))
LA_SE$Pool = factor(LA_SE$Pool,
                    levels=unique(LA_SE$Pool))
LA_SE$Population = factor(LA_SE$Population,
                          levels=unique(LA_SE$Population))
pd <- position_dodge(.2) ## How much to jitter the points on the plot

ggplot(LA_SE, aes(x = Population, y = LA, colour = Pool, group = Pool))+
  geom_errorbar(aes(ymin=LA-SE, ymax=LA+SE), colour="black", width=.3, position=pd, size = 1) +
  geom_line(position=pd, size = 1) +
  geom_point(position = pd, aes(shape=Pool, fill=Pool), size=8) +
  scale_y_continuous(expression("Leaf Area" ~~ (cm^2)), limits = c(8, 22), breaks=seq(8, 22, by = 2)) +
  scale_shape_manual(values=c(21,22,23,25,24)) +
  scale_fill_manual(values=c("#e6194b","#ffe119","#3cb44b","#4363d8","#42d4f4")) + 
  ylab("LA  (cm2)") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.key.size = unit(1.2, "cm"), legend.key.width = unit(0.75,"cm")) #
# save as "LA_EMM_SE_by_Population_pop.TIFF" W = 1034 H = 842