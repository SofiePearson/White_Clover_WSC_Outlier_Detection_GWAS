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

########################
## Import data
########################
SSS_REML <- read.table(file = here::here("SSS_600_samples.txt"), header = TRUE) %>%
  separate(Generation, into = c("Trt", "Generation"), "-") %>%
  mutate(Trt = factor(Trt, levels = c("Parent", "Low", "High")),
         Generation = factor(ifelse(Trt == "Parent", "Start", Generation), 
                             levels = c("Start", "Mid", "End")),
         EID = ifelse(Pool %in% c("WNZLL", "WNZSL", "WUSLL"), "W", "F"),
         Parent = factor(ifelse(Trt == "Parent", "Y", "N"),
                         levels = c("Y", "N"))) %>%
  transform(TrtC = interaction(Pool, Generation, Trt))

########################
## Explore data
########################
SSS_REML %>%
  group_by(Pool, Trt, Generation) %>% 
  summarise(mean_SSS = mean(SSS)) %>%
  ggplot(aes(Generation, mean_SSS, col = Trt)) +
  geom_point() + facet_wrap(~Pool)

## Baseline is Start:Parent, so comparisons are to this group
pool_mods <- lapply(unique(SSS_REML$Pool), function(p) {
  lmer(SSS ~ TrtC + 
         (1 | Block) + (1 | Block:Column) + (1 | Block:Row),
       data = SSS_REML %>% filter(Pool == p) %>%
         transform(TrtC = interaction(Generation, Trt)))
})
lapply(pool_mods, summary)
## Some evidence of differing results for each Pool, so final model should have 3 way interaction

########################
## Starting model
########################
mod <- lmer(SSS ~ TrtC + 
              (1 | Block) + (1 | Block:Column) + (1 | Block:Row),
            data = SSS_REML)
summary(mod)

########################
## View residuals
########################
SSS_REML_means <- SSS_REML %>%
  mutate(resids_pearson = resid(mod, type = "pearson"),
         fitteds = fitted(mod))

p <- list()

p[[length(p) + 1]] <- ggplot(SSS_REML_means, aes(fitteds, resids_pearson))+
  geom_point(aes(color = Trt), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#ffe119", "#4363d8", "#3cb44b"),name = "Treatment", labels = c("None", "Low", "High"))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Treatment"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(fitteds, resids_pearson, col = Generation)) + 
  geom_point(aes(color = Generation), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c("#ffe119", "#4363d8", "#3cb44b"),name = "Generation", labels = c("Parent", "Mid", "End"))+
  geom_hline(yintercept = 0)+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Generation"))+
  # scale_color_discrete(name = "Generation", labels = c("Parent", "Mid", "End"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- ggplot(SSS_REML_means, aes(fitteds, resids_pearson)) + 
  geom_point(aes(color = Pool), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#4363d8", "#42d4f4", "#e6194b","#ffe119", "#3cb44b"))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Pool"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(fitteds, resids_pearson, col = paste0(Trt, Generation))) + 
  geom_point(aes(color = paste0(Trt, Generation)), size = 2, alpha = 0.8) + 
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Treatment and Generation"))+
  scale_color_manual(values = c("#4363d8", "#42d4f4", "#e6194b","#ffe119", "#3cb44b"), name = "Trt & Gen", labels = c("High & End", "High & Mid", "Low & End", "Low & Mid", "None & Parent"))+
  scale_x_continuous(name="Fitted values") +
  scale_y_continuous(name="Pearson residuals")

grid.arrange(grobs = p, ncol = 1)
## Save plot as "Treatments_residual_SSS_edited.tiff" W = 650 H = 850 (Supplementary Figure 1, left)

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
SSS_REML_means <- SSS_REML_means %>%
  mutate(resids_pearson = resid(mod, type = "pearson"),
         fitteds = fitted(mod)) %>%
  left_join(block_ranefs %>% mutate_all(~as.numeric(.))) %>%
  left_join(blockCol_ranefs %>% mutate_all(~as.numeric(.))) %>%
  left_join(blockRow_ranefs %>% mutate_all(~as.numeric(.))) %>%
  mutate(total_re = re_value_block + re_value_row + re_value_column )

p <- list()

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_block)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Block"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_row)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Row"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = re_value_column)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Column"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = total_re)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Total random effects"))

grid.arrange(grobs = p, ncol = 2)
## Save plot as "Random_effects_residual_SSS_edited.tiff" W = 800 H = 800

########################
## Plot spatial residuals
########################
p <- list()

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = resids_pearson)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Residuals"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = paste0(Trt, Generation))) + 
  geom_tile() +
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(C)~~"Treatment and Generation"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = Trt)) + 
  geom_tile()+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(B)~~"Treatment"))

p[[length(p) + 1]] <- SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = Generation)) + 
  geom_tile() +
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(D)~~"Generation"))
grid.arrange(grobs = p, ncol = 2)
## Save plot as "Treatments_residual_spatial_SSS_edited.tiff" W = 800 H = 800

## Just need the first plot as there is no pattern for the residuals
SSS_REML_means %>% 
  ggplot(aes(Column, Row, fill = resids_pearson)) + 
  geom_tile() +
  scale_fill_gradient2(name = "Residuals")+
  theme_bw()+
  theme_minimal()+
  ggtitle(expression(bold(A)~~"Residuals for SSS"))
## Save plot as "Residual_spatial_SSS_edited.tiff" W = 700 H = 800

########################
## Building the model to produce estimated SSS means per pop and comparisons to the Parent baseline
########################
em <- emmeans(mod, "TrtC", type = "response")
plot(em)
em_levels <- levels(em)$TrtC
adj_type <- "mvt" 
adj_type_sub <- "none" 
Generation_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
  apply(1, function(i) paste0(i, collapse = "."))
em_GenerationTrt_ <- add_grouping(em,
                                  newname = "Generation_Trt",
                                  refname = "TrtC",
                                  Generation_Trt_values)
em_GenerationTrt <- emmeans(em_GenerationTrt_, trt.vs.ctrl ~ Generation_Trt, ref = "Start.Parent", adjust = adj_type_sub)
summary(em_GenerationTrt, type = "response", adj = adj_type_sub) 

# $emmeans
# Generation_Trt emmean   SE   df lower.CL upper.CL
# End.High          180 3.26 7.40    172.5      188
# End.Low           102 3.26 7.39     94.2      109
# Mid.High          164 3.25 7.34    156.0      171
# Mid.Low           110 3.26 7.44    102.5      118
# Start.Parent      127 2.56 2.84    119.0      136
# 
# Results are averaged over the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate   SE  df t.ratio p.value
# End.High - Start.Parent     52.7 3.48 542 15.160  <.0001 
# End.Low - Start.Parent     -25.6 3.47 538 -7.366  <.0001 
# Mid.High - Start.Parent     36.3 3.46 530 10.467  <.0001 
# Mid.Low - Start.Parent     -17.2 3.49 547 -4.941  <.0001 
# 
# Results are averaged over the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 

plot(em_GenerationTrt, type = "response")
Generation_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
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
                 newname = "Generation_Trt",
                 refname = "TrtC",
                 paste0(pool_id, ".",
                        c("Start.Parent", "Mid.Low", "End.Low", 
                          "Mid.High", "End.High")))
  # print(res)
  res <- res %>%
    emmeans(., trt.vs.ctrl ~ Generation_Trt, ref = paste0(pool_id, ".", "Start.Parent"), adjust = adj_type_sub)
  # print(res)
  return(res)
})

comparisons_GenerationTrt_withinPool <- 
  lapply(tsts, function(i) i$contrasts) %>% ## contrasts holds the comparisons
  do.call("rbind", .) %>% 
  update(adj = adj_type_sub)
em_GenerationTrt_withinPool <-
  lapply(tsts, function(i) i$emmeans) %>% ## emmeans holds the estimated means values
  do.call("rbind", .) %>% 
  update(adjust = adj_type_sub) 

print(em_GenerationTrt_withinPool) ## Use these emmeans and SE values (Supplementary Table 1)
# Generation_Trt     emmean   SE    df lower.CL upper.CL
# FNZLL.End.High      184.6 6.55 102.0    171.6    197.5
# FNZLL.End.Low       129.8 6.55 102.0    116.8    142.8
# FNZLL.Mid.High      168.4 6.55 102.0    155.4    181.3
# FNZLL.Mid.Low       120.3 6.55 102.0    107.4    133.3
# FNZLL.Start.Parent  143.2 4.77  32.7    133.5    152.9
# FNZSL.End.High      173.5 6.55 102.0    160.5    186.5
# FNZSL.End.Low        74.8 6.55 102.0     61.8     87.7
# FNZSL.Mid.High      159.6 6.55 102.0    146.7    172.6
# FNZSL.Mid.Low        83.2 6.55 102.0     70.2     96.2
# FNZSL.Start.Parent  101.9 4.76  32.4     92.2    111.6
# WNZLL.End.High      191.6 6.55 102.0    178.6    204.6
# WNZLL.End.Low        92.0 6.55 102.0     79.0    105.0
# WNZLL.Mid.High      170.5 6.55 102.0    157.5    183.5
# WNZLL.Mid.Low       123.3 6.55 102.0    110.3    136.3
# WNZLL.Start.Parent  134.2 4.76  32.4    124.5    143.9
# WNZSL.End.High      188.9 6.55 102.0    175.9    201.9
# WNZSL.End.Low       119.1 6.55 102.0    106.2    132.1
# WNZSL.Mid.High      167.9 6.55 102.0    154.9    180.9
# WNZSL.Mid.Low       119.7 6.55 102.0    106.7    132.7
# WNZSL.Start.Parent  119.2 4.77  32.7    109.5    128.9
# WUSLL.End.High      162.1 6.55 102.0    149.1    175.1
# WUSLL.End.Low        93.4 6.55 102.0     80.4    106.4
# WUSLL.Mid.High      151.9 6.55 102.0    138.9    164.9
# WUSLL.Mid.Low       104.3 6.55 102.0     91.3    117.2
# WUSLL.Start.Parent  138.5 4.77  32.7    128.8    148.2
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 

print(comparisons_GenerationTrt_withinPool) ## Don't use these p-values ## Use these estimates and SE values (Table 1)
# contrast                            estimate   SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.Start.Parent  41.3677 7.79 545  5.311  <.0001 
# FNZLL.End.Low - FNZLL.Start.Parent  -13.4124 7.77 539 -1.726  0.0850 
# FNZLL.Mid.High - FNZLL.Start.Parent  25.1781 7.78 541  3.238  0.0013 
# FNZLL.Mid.Low - FNZLL.Start.Parent  -22.8424 7.79 545 -2.932  0.0035 
# FNZSL.End.High - FNZSL.Start.Parent  71.6089 7.77 540  9.211  <.0001 
# FNZSL.End.Low - FNZSL.Start.Parent  -27.1083 7.76 534 -3.495  0.0005 
# FNZSL.Mid.High - FNZSL.Start.Parent  57.7625 7.78 543  7.422  <.0001 
# FNZSL.Mid.Low - FNZSL.Start.Parent  -18.6531 7.77 537 -2.402  0.0166 
# WNZLL.End.High - WNZLL.Start.Parent  57.4637 7.77 537  7.400  <.0001 
# WNZLL.End.Low - WNZLL.Start.Parent  -42.1846 7.77 539 -5.430  <.0001 
# WNZLL.Mid.High - WNZLL.Start.Parent  36.3558 7.77 539  4.679  <.0001 
# WNZLL.Mid.Low - WNZLL.Start.Parent  -10.8535 7.78 542 -1.395  0.1635 
# WNZSL.End.High - WNZSL.Start.Parent  69.6310 7.78 541  8.954  <.0001 
# WNZSL.End.Low - WNZSL.Start.Parent   -0.0982 7.79 545 -0.013  0.9899 
# WNZSL.Mid.High - WNZSL.Start.Parent  48.6744 7.76 536  6.270  <.0001 
# WNZSL.Mid.Low - WNZSL.Start.Parent    0.4428 7.78 542  0.057  0.9546 
# WUSLL.End.High - WUSLL.Start.Parent  23.6055 7.76 536  3.041  0.0025 
# WUSLL.End.Low - WUSLL.Start.Parent  -45.1185 7.80 547 -5.786  <.0001 
# WUSLL.Mid.High - WUSLL.Start.Parent  13.3676 7.76 536  1.722  0.0857 
# WUSLL.Mid.Low - WUSLL.Start.Parent  -34.2444 7.77 539 -4.406  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 

plot(em_GenerationTrt_withinPool)+
  labs(x = "Estimated marginal mean (SSS)", y = "Population") ## Blue lines are confidence intervals for the EMMs

########################
## Get CIs
########################         
res_overall <- rbind(
  comparisons_GenerationTrt_withinPool, 
  em_GenerationTrt$contrasts, 
  # em_Trt$contrasts,
  # em_Generation$contrasts, 
  # em_Pool$contrasts,
  # em_ExperimentID$contrasts, 
  adjust = adj_type)

print(res_overall) ## Use these p-values (Table 1)
# contrast                            estimate   SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.Start.Parent  41.3677 7.79 545  5.311  <.0001 
# FNZLL.End.Low - FNZLL.Start.Parent  -13.4124 7.77 539 -1.726  0.8187 
# FNZLL.Mid.High - FNZLL.Start.Parent  25.1781 7.78 541  3.238  0.0285 
# FNZLL.Mid.Low - FNZLL.Start.Parent  -22.8424 7.79 545 -2.932  0.0745 
# FNZSL.End.High - FNZSL.Start.Parent  71.6089 7.77 540  9.211  <.0001 
# FNZSL.End.Low - FNZSL.Start.Parent  -27.1083 7.76 534 -3.495  0.0120 
# FNZSL.Mid.High - FNZSL.Start.Parent  57.7625 7.78 543  7.422  <.0001 
# FNZSL.Mid.Low - FNZSL.Start.Parent  -18.6531 7.77 537 -2.402  0.2943 
# WNZLL.End.High - WNZLL.Start.Parent  57.4637 7.77 537  7.400  <.0001 
# WNZLL.End.Low - WNZLL.Start.Parent  -42.1846 7.77 539 -5.430  <.0001 
# WNZLL.Mid.High - WNZLL.Start.Parent  36.3558 7.77 539  4.679  0.0001 
# WNZLL.Mid.Low - WNZLL.Start.Parent  -10.8535 7.78 542 -1.395  0.9651 
# WNZSL.End.High - WNZSL.Start.Parent  69.6310 7.78 541  8.954  <.0001 
# WNZSL.End.Low - WNZSL.Start.Parent   -0.0982 7.79 545 -0.013  1.0000 
# WNZSL.Mid.High - WNZSL.Start.Parent  48.6744 7.76 536  6.270  <.0001 
# WNZSL.Mid.Low - WNZSL.Start.Parent    0.4428 7.78 542  0.057  1.0000 
# WUSLL.End.High - WUSLL.Start.Parent  23.6055 7.76 536  3.041  0.0542 
# WUSLL.End.Low - WUSLL.Start.Parent  -45.1185 7.80 547 -5.786  <.0001 
# WUSLL.Mid.High - WUSLL.Start.Parent  13.3676 7.76 536  1.722  0.8213 
# WUSLL.Mid.Low - WUSLL.Start.Parent  -34.2444 7.77 539 -4.406  0.0003 
# End.High - Start.Parent              52.7353 3.48 542 15.160  <.0001 
# End.Low - Start.Parent              -25.5844 3.47 538 -7.366  <.0001 
# Mid.High - Start.Parent              36.2677 3.46 530 10.467  <.0001 
# Mid.Low - Start.Parent              -17.2301 3.49 547 -4.941  <.0001 
# 
# Results are averaged over some or all of the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: mvt method for 24 tests 

res_overall_ci <- rbind(
  em,
  em_GenerationTrt$emmeans,
  # em_Trt$emmeans,
  # em_Generation$emmeans,
  # em_Pool$emmeans,
  # em_ExperimentID$emmeans,
  adjust = adj_type)

print(res_overall_ci) ## This takes a while, Use these CIs. The mean and SE are same as first output so can just use these ones below (Supplementary Table 1)
# TrtC               Generation_Trt emmean   SE     df lower.CL upper.CL
# FNZLL.Start.Parent .               143.2 4.77  32.67    127.1    159.3
# FNZSL.Start.Parent .               101.9 4.76  32.42     85.8    117.9
# WNZLL.Start.Parent .               134.2 4.76  32.42    118.1    150.2
# WNZSL.Start.Parent .               119.2 4.77  32.67    103.2    135.3
# WUSLL.Start.Parent .               138.5 4.77  32.67    122.4    154.6
# FNZLL.Mid.Low      .               120.3 6.55 102.00     99.4    141.3
# FNZSL.Mid.Low      .                83.2 6.55 101.99     62.3    104.1
# WNZLL.Mid.Low      .               123.3 6.55 102.00    102.4    144.2
# WNZSL.Mid.Low      .               119.7 6.55 102.00     98.8    140.6
# WUSLL.Mid.Low      .               104.3 6.55 101.99     83.4    125.2
# FNZLL.End.Low      .               129.8 6.55 101.99    108.8    150.7
# FNZSL.End.Low      .                74.8 6.55 101.99     53.8     95.7
# WNZLL.End.Low      .                92.0 6.55 101.99     71.0    112.9
# WNZSL.End.Low      .               119.1 6.55 102.00     98.2    140.0
# WUSLL.End.Low      .                93.4 6.55 101.99     72.4    114.3
# FNZLL.Mid.High     .               168.4 6.55 101.99    147.4    189.3
# FNZSL.Mid.High     .               159.6 6.55 101.99    138.7    180.6
# WNZLL.Mid.High     .               170.5 6.55 101.99    149.6    191.4
# WNZSL.Mid.High     .               167.9 6.55 101.99    147.0    188.8
# WUSLL.Mid.High     .               151.9 6.55 101.99    131.0    172.8
# FNZLL.End.High     .               184.6 6.55 102.00    163.6    205.5
# FNZSL.End.High     .               173.5 6.55 101.99    152.6    194.4
# WNZLL.End.High     .               191.6 6.55 101.99    170.7    212.5
# WNZSL.End.High     .               188.9 6.55 101.99    167.9    209.8
# WUSLL.End.High     .               162.1 6.55 101.99    141.2    183.0
# .                  End.High        180.1 3.26   7.40    165.8    194.5
# .                  End.Low         101.8 3.26   7.39     87.5    116.2
# .                  Mid.High        163.7 3.25   7.34    149.3    178.0
# .                  Mid.Low         110.2 3.26   7.44     95.8    124.5
# .                  Start.Parent    127.4 2.56   2.84    109.8    145.0
# 
# Results are averaged over some or all of the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: mvt method for 30 estimates 

plot(em_GenerationTrt_withinPool, comparisons = TRUE)
## Blue bars are the confidence intervals for the EMMs and the red arrows are the comparisons among them. 
## If an arrow from one mean overlaps an arrow from another group, the difference is not significant at an alpha of 0.05
plot(em_GenerationTrt_withinPool, comparisons = TRUE, type = "response")+
  theme_bw()+
  theme_minimal()

########################
## Create text file with means, SE and CI for each population (include Pool, pop and generation info) (SSS_EMMs_pop.txt) and plot using the code below for Figure 2       
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

SSS_SE <- read.table(file = "SSS_EMMs_pop.txt", header = TRUE)
SSS_SE$Generation = factor(SSS_SE$Generation,
                           levels=unique(SSS_SE$Generation))
SSS_SE$Pool = factor(SSS_SE$Pool,
                     levels=unique(SSS_SE$Pool))
SSS_SE$Population = factor(SSS_SE$Population,
                           levels=unique(SSS_SE$Population))
pd <- position_dodge(.2) ## How much to jitter the points on the plot

ggplot(SSS_SE, aes(x = Population, y = SSS, colour = Pool, group = Pool))+
  geom_errorbar(aes(ymin=SSS-SE, ymax=SSS+SE), colour="black", width=.3, position=pd, size = 1) +
  geom_line(position=pd, size = 1) +
  geom_point(position = pd, aes(shape=Pool, fill=Pool), size=8) +
  scale_y_continuous(expression("SSS  (g kg"^"-1"*' DM)'), limits = c(65, 200), breaks=seq(25, 200, by = 25)) +
  scale_shape_manual(values=c(21,22,23,25,24)) +
  scale_fill_manual(values=c("#e6194b","#ffe119","#3cb44b","#4363d8","#42d4f4")) + 
  ylab("SSS  (%DM)") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.key.size = unit(1.2, "cm"), legend.key.width = unit(0.75,"cm")) #
## Save as "SSS_EMM_SE_by_Population_pop.TIFF" W = 1034 H = 842 (Figure 2)

########################
## Test differences between just High-End and Low-End:
## Make the base line Low-End
########################
em <- emmeans(mod, "TrtC", type = "response")
plot(em)
em_levels <- levels(em)$TrtC
adj_type <- "mvt" 
adj_type_sub <- "none" 
Generation_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
  apply(1, function(i) paste0(i, collapse = "."))
em_GenerationTrt_ <- add_grouping(em,
                                  newname = "Generation_Trt",
                                  refname = "TrtC",
                                  Generation_Trt_values)
em_GenerationTrt <- emmeans(em_GenerationTrt_, trt.vs.ctrl ~ Generation_Trt, ref = "End.Low", adjust = adj_type_sub)

summary(em_GenerationTrt, type = "response", adj = adj_type_sub) ## Check the new baseline worked - yup it did
# $emmeans
# Generation_Trt emmean   SE   df lower.CL upper.CL
# End.High          180 3.26 7.40    172.5      188
# End.Low           102 3.26 7.39     94.2      109
# Mid.High          164 3.25 7.34    156.0      171
# Mid.Low           110 3.26 7.44    102.5      118
# Start.Parent      127 2.56 2.84    119.0      136
# 
# Results are averaged over the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate   SE  df t.ratio p.value
# End.High - End.Low        78.32 4.02 543 19.491  <.0001 
# Mid.High - End.Low        61.85 4.01 539 15.416  <.0001 
# Mid.Low - End.Low          8.35 4.03 550  2.071  0.0388 
# Start.Parent - End.Low    25.58 3.47 538  7.366  <.0001
# 
# Results are averaged over the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 

plot(em_GenerationTrt, type = "response")
Generation_Trt_values <- str_split(em_levels, "\\.", simplify = T)[, -1] %>%
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
                 newname = "Generation_Trt",
                 refname = "TrtC",
                 paste0(pool_id, ".",
                        c("Start.Parent", "Mid.Low", "End.Low", 
                          "Mid.High", "End.High")))
  # print(res)
  res <- res %>%
    emmeans(., trt.vs.ctrl ~ Generation_Trt, ref = paste0(pool_id, ".", "End.Low"), adjust = adj_type_sub)
  # print(res)
  return(res)
})

comparisons_GenerationTrt_withinPool <- 
  lapply(tsts, function(i) i$contrasts) %>% ## Contrasts holds the comparisons
  do.call("rbind", .) %>% 
  update(adj = adj_type_sub)
em_GenerationTrt_withinPool <-
  lapply(tsts, function(i) i$emmeans) %>% ## Emmeans holds the estimated means values
  do.call("rbind", .) %>% 
  update(adjust = adj_type_sub) 

print(em_GenerationTrt_withinPool)
# Generation_Trt     emmean   SE    df lower.CL upper.CL
# FNZLL.End.High      184.6 6.55 102.0    171.6    197.5
# FNZLL.End.Low       129.8 6.55 102.0    116.8    142.8
# FNZLL.Mid.High      168.4 6.55 102.0    155.4    181.3
# FNZLL.Mid.Low       120.3 6.55 102.0    107.4    133.3
# FNZLL.Start.Parent  143.2 4.77  32.7    133.5    152.9
# FNZSL.End.High      173.5 6.55 102.0    160.5    186.5
# FNZSL.End.Low        74.8 6.55 102.0     61.8     87.7
# FNZSL.Mid.High      159.6 6.55 102.0    146.7    172.6
# FNZSL.Mid.Low        83.2 6.55 102.0     70.2     96.2
# FNZSL.Start.Parent  101.9 4.76  32.4     92.2    111.6
# WNZLL.End.High      191.6 6.55 102.0    178.6    204.6
# WNZLL.End.Low        92.0 6.55 102.0     79.0    105.0
# WNZLL.Mid.High      170.5 6.55 102.0    157.5    183.5
# WNZLL.Mid.Low       123.3 6.55 102.0    110.3    136.3
# WNZLL.Start.Parent  134.2 4.76  32.4    124.5    143.9
# WNZSL.End.High      188.9 6.55 102.0    175.9    201.9
# WNZSL.End.Low       119.1 6.55 102.0    106.2    132.1
# WNZSL.Mid.High      167.9 6.55 102.0    154.9    180.9
# WNZSL.Mid.Low       119.7 6.55 102.0    106.7    132.7
# WNZSL.Start.Parent  119.2 4.77  32.7    109.5    128.9
# WUSLL.End.High      162.1 6.55 102.0    149.1    175.1
# WUSLL.End.Low        93.4 6.55 102.0     80.4    106.4
# WUSLL.Mid.High      151.9 6.55 102.0    138.9    164.9
# WUSLL.Mid.Low       104.3 6.55 102.0     91.3    117.2
# WUSLL.Start.Parent  138.5 4.77  32.7    128.8    148.2
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 

print(comparisons_GenerationTrt_withinPool) ## Don't use these p-values
# contrast                           estimate   SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.End.Low      54.7801 8.98 541  6.101  <.0001 
# FNZLL.Mid.High - FNZLL.End.Low      38.5905 8.98 541  4.298  <.0001 
# FNZLL.Mid.Low - FNZLL.End.Low       -9.4300 8.98 541 -1.050  0.2941 
# FNZLL.Start.Parent - FNZLL.End.Low  13.4124 7.77 539  1.726  0.0850 
# FNZSL.End.High - FNZSL.End.Low      98.7172 8.98 541 10.994  <.0001 
# FNZSL.Mid.High - FNZSL.End.Low      84.8707 8.97 539  9.460  <.0001 
# FNZSL.Mid.Low - FNZSL.End.Low        8.4551 8.99 543  0.941  0.3472 
# FNZSL.Start.Parent - FNZSL.End.Low  27.1083 7.76 534  3.495  0.0005 
# WNZLL.End.High - WNZLL.End.Low      99.6483 8.99 545 11.079  <.0001 
# WNZLL.Mid.High - WNZLL.End.Low      78.5405 8.99 543  8.740  <.0001 
# WNZLL.Mid.Low - WNZLL.End.Low       31.3311 8.99 543  3.486  0.0005 
# WNZLL.Start.Parent - WNZLL.End.Low  42.1846 7.77 539  5.430  <.0001 
# WNZSL.End.High - WNZSL.End.Low      69.7292 8.98 541  7.766  <.0001 
# WNZSL.Mid.High - WNZSL.End.Low      48.7726 8.99 543  5.427  <.0001 
# WNZSL.Mid.Low - WNZSL.End.Low        0.5410 8.99 543  0.060  0.9520 
# WNZSL.Start.Parent - WNZSL.End.Low   0.0982 7.79 545  0.013  0.9899 
# WUSLL.End.High - WUSLL.End.Low      68.7240 8.98 541  7.654  <.0001 
# WUSLL.Mid.High - WUSLL.End.Low      58.4860 8.98 541  6.514  <.0001 
# WUSLL.Mid.Low - WUSLL.End.Low       10.8741 8.98 541  1.211  0.2264 
# WUSLL.Start.Parent - WUSLL.End.Low  45.1185 7.80 547  5.786  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 

res_overall <- rbind(
  comparisons_GenerationTrt_withinPool, 
  em_GenerationTrt$contrasts, 
  # em_Trt$contrasts,
  # em_Generation$contrasts, 
  # em_Pool$contrasts,
  # em_ExperimentID$contrasts, 
  adjust = adj_type)
print(res_overall) ## Use these p-values and the differences between the Low-End and High-End values for each Pool

# contrast                           estimate   SE  df t.ratio p.value
# FNZLL.End.High - FNZLL.End.Low      54.7801 8.98 541  6.101  <.0001 
# FNZLL.Mid.High - FNZLL.End.Low      38.5905 8.98 541  4.298  0.0005 
# FNZLL.Mid.Low - FNZLL.End.Low       -9.4300 8.98 541 -1.050  0.9962 
# FNZLL.Start.Parent - FNZLL.End.Low  13.4124 7.77 539  1.726  0.7689 
# FNZSL.End.High - FNZSL.End.Low      98.7172 8.98 541 10.994  <.0001 
# FNZSL.Mid.High - FNZSL.End.Low      84.8707 8.97 539  9.460  <.0001 
# FNZSL.Mid.Low - FNZSL.End.Low        8.4551 8.99 543  0.941  0.9990 
# FNZSL.Start.Parent - FNZSL.End.Low  27.1083 7.76 534  3.495  0.0117 
# WNZLL.End.High - WNZLL.End.Low      99.6483 8.99 545 11.079  <.0001 
# WNZLL.Mid.High - WNZLL.End.Low      78.5405 8.99 543  8.740  <.0001 
# WNZLL.Mid.Low - WNZLL.End.Low       31.3311 8.99 543  3.486  0.0116 
# WNZLL.Start.Parent - WNZLL.End.Low  42.1846 7.77 539  5.430  <.0001 
# WNZSL.End.High - WNZSL.End.Low      69.7292 8.98 541  7.766  <.0001 
# WNZSL.Mid.High - WNZSL.End.Low      48.7726 8.99 543  5.427  <.0001 
# WNZSL.Mid.Low - WNZSL.End.Low        0.5410 8.99 543  0.060  1.0000 
# WNZSL.Start.Parent - WNZSL.End.Low   0.0982 7.79 545  0.013  1.0000 
# WUSLL.End.High - WUSLL.End.Low      68.7240 8.98 541  7.654  <.0001 
# WUSLL.Mid.High - WUSLL.End.Low      58.4860 8.98 541  6.514  <.0001 
# WUSLL.Mid.Low - WUSLL.End.Low       10.8741 8.98 541  1.211  0.9832 
# WUSLL.Start.Parent - WUSLL.End.Low  45.1185 7.80 547  5.786  <.0001 
# End.High - End.Low                  78.3197 4.02 543 19.491  <.0001 
# Mid.High - End.Low                  61.8521 4.01 539 15.416  <.0001 
# Mid.Low - End.Low                    8.3543 4.03 550  2.071  0.4980 
# Start.Parent - End.Low              25.5844 3.47 538  7.366  <.0001 
# 
# Results are averaged over some or all of the levels of: TrtC 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: mvt method for 24 tests  

## Use the last output for differences between the Low-End and High-End populations for each pool (Table 1 - last five rows)
