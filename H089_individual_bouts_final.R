################################################################################
#H089 individual bout analysis 

#make sure to have H089n data available (Run analysis script)
################################################################################

library(tidyverse)
library(reshape2)
library(irr)
library(broom)
library(modifiedmk)
library (ggrepel)




#Import data
df.t <- read.csv( "./data/CSF_H089_final.csv", header = TRUE, check.names = FALSE)

#store in working object:
H089<- df.t



# Apply sqrt transformation centering and scaling
met.sc <- scale(sqrt(H089[,6:ncol(H089)]), center = TRUE, scale = TRUE)

rownames(met.sc) <- H089[,3] #Making row names the samples_nr

#Centered and scaled working object:
H089n <- met.sc


#colnames(H089n)
head(H089n[,1:8])

#bind to get tb_mean
H089n <- cbind(tb_mean = H089$tb_mean, H089n)
H089n <- cbind(sample_nr = H089$sample_nr, H089n)

H089n <- as_tibble(H089n)

#Remove qc vars
H089n <- H089n%>%filter(sample_nr >0)

head(H089n)



# mk plot to get overview
Met <- "Adenosine"

t <- cbind(df.t[3], met.sc) %>%
  melt(id = "sample_nr") %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(variable %in% Met) %>%
  filter(sample_nr > 0)

ggplot(t, aes(x = as.numeric(sample_nr) * 2, y = value)) +
  geom_ribbon(
    data = H089n, inherit.aes = FALSE,
    aes(
      x = sample_nr * 2,
      ymin = -2,
      ymax = (tb_mean / 9) - 2.5 + 0.1
    ),
    fill = "steelblue4",
    alpha = 0.4
  ) +
  geom_point(color = "gray15", size = 1) +
  geom_line(na.rm = TRUE, color = "gray15", linewidth = 0.5) +
  scale_x_continuous(
    name = "Time (hours)",
    breaks = seq(0, 220, 20)
  ) +
  theme_bw()


################################################################################
# need outlier removal before any between-bout variance analysis can do much.
# metabolites covary with Tb cyclically, so fit abundance ~ tb_mean and flag
# points more than 3 SD off. Replace with mean of nearest non-outlier neighbors.
################################################################################

# Convert to long format for outlier detection
H089n_long <- H089n %>%
  pivot_longer(
    cols = -c(sample_nr, tb_mean),
    names_to = "Metabolite",
    values_to = "Abundance"
  )

length(unique(H089n_long$Metabolite))


Met <- "Adenosine"

t <- cbind(df.t[5], met.sc) %>%
  melt(id = "tb_mean") %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(variable %in% Met)

fit <- lm(value ~ tb_mean, data = t)
t$fitted <- fitted(fit)
rsd <- sigma(fit)

ggplot(t, aes(x = tb_mean, y = value)) +
  geom_ribbon(aes(ymin = fitted - 3*rsd, ymax = fitted +
                    3*rsd),
              fill = "steelblue", alpha = 0.15) +
  geom_point(color = "gray15", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color =
                "steelblue") +
  theme_bw()

# That looks good, let build to loop across all. 


detect_outliers_tb <- function(abundance, tb, k = 3) {
  valid_idx <- !is.na(abundance) & !is.na(tb)
  if(sum(valid_idx) < 4) {
    return(rep(FALSE, length(abundance)))
  }
  
  fit <- lm(abundance[valid_idx] ~ tb[valid_idx])
  residuals <- rep(NA, length(abundance))
  residuals[valid_idx] <- resid(fit)
  
  res_sd <- sigma(fit)
  
  is_outlier <- abs(residuals) > k * res_sd
  is_outlier[is.na(is_outlier)] <- FALSE
  
  return(is_outlier)
}

# Detect outliers per metabolite (based on deviation from Tb relationship)
H089n_long <- H089n_long %>%
  group_by(Metabolite) %>%
  mutate(
    is_outlier = detect_outliers_tb(Abundance, tb_mean, k = 3)
  ) %>%
  ungroup()

# let see a summary of outliers
outlier_summary <- H089n_long %>%
  group_by(Metabolite) %>%
  summarise(
    n_outliers = sum(is_outlier, na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  )

head(outlier_summary)
summary(outlier_summary)


# it seems many have one or 2 outliers which is fine. No one with more than 6 outliers so that is good.

# To preserve time series lets replace removed outliers with mean of nearest neighbors:
replace_outliers_with_neighbors <- function(abundance, is_outlier, sample_idx) {
  # Order by sample index
  ord <- order(sample_idx)
  abundance_ord <- abundance[ord]
  outlier_ord <- is_outlier[ord]
  
  corrected <- abundance_ord
  n <- length(abundance_ord)
  
  for(i in which(outlier_ord)) {
    # Find nearest non-outlier neighbor before
    before_val <- NA
    for(j in (i-1):1) {
      if(j < 1) break
      if(!outlier_ord[j] && !is.na(abundance_ord[j])) {
        before_val <- abundance_ord[j]
        break
      }
    }
    
    # Find nearest non-outlier neighbor after
    after_val <- NA
    for(j in (i+1):n) {
      if(j > n) break
      if(!outlier_ord[j] && !is.na(abundance_ord[j])) {
        after_val <- abundance_ord[j]
        break
      }
    }
    
    # Take mean of available neighbors
    neighbors <- c(before_val, after_val)
    neighbors <- neighbors[!is.na(neighbors)]
    
    if(length(neighbors) > 0) {
      corrected[i] <- mean(neighbors)
    }
  }
  
  # Return in original order
  corrected[order(ord)]
}

H089n_long <- H089n_long %>%
  group_by(Metabolite) %>%
  mutate(
    Abundance_corrected = replace_outliers_with_neighbors(Abundance, is_outlier, sample_nr),
    Abundance_original = Abundance,
    Abundance = Abundance_corrected
  ) %>%
  ungroup()

# Convert back to wide format
H089n <- H089n_long %>%
  select(sample_nr, tb_mean, Metabolite, Abundance) %>%
  pivot_wider(
    names_from = Metabolite,
    values_from = Abundance
  )


# Plot check: visualize outlier correction for Adenosine with Tb overlay
H089n_long %>%
  filter(Metabolite == "Adenosine") %>% #"298.1757/10.202" is the one with most
  ggplot(aes(x = sample_nr * 2)) +
  geom_ribbon(
    aes(ymin = -2, ymax = (tb_mean / 9) - 2.5 + 0.1),
    fill = "steelblue4",
    alpha = 0.3
  ) +
  # Show original values as open circles for outliers
  geom_point(aes(y = Abundance_original, shape = is_outlier, color = is_outlier),
             size = 2, na.rm = TRUE) +
  # Show corrected line
  geom_line(aes(y = Abundance), color = "gray30", na.rm = TRUE) +
  geom_point(aes(y = Abundance), color = "gray30", size = 1, na.rm = TRUE) +
  scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "firebrick"),
                     labels = c("Normal", "Outlier (original)")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1),
                     labels = c("Normal", "Outlier (original)")) +
  labs(
    title = "Outlier Correction Check - Adenosine",
    subtitle = "Red open circles = original outliers; Black = corrected (neighbor mean); Blue shading = Tb",
    x = "Time (hours)",
    y = "Normalized Abundance",
    color = "",
    shape = ""
  ) +
  scale_x_continuous(breaks = seq(0, 220, 20)) +
  theme_bw() +
  theme(legend.position = "bottom")





################################################################################
# OOK, lets convert the bout hours to percentages from 0-100% bout completion and overlap.
################################################################################

#Split by time and assign Tb ~30-30 as 0 to 100% bout progression, Ignore IBE for now. 
bout1 <- H089n%>%
  dplyr::filter(sample_nr > 1 & sample_nr < 33)%>%
  mutate(Percent = seq(-4, 117, length.out = n()))%>%
  relocate(Percent)

bout2 <- H089n%>%
  dplyr::filter(sample_nr > 42 & sample_nr < 86)%>%
  mutate(Percent = seq(-2, 110, length.out = n()))%>%
  relocate(Percent)

#check they overlap well
ggplot()+
  geom_point(data= bout1, aes(x = Percent, y=tb_mean), color= "firebrick")+
  geom_line(data= bout1, aes(x = Percent, y=tb_mean))+
  geom_point(data= bout2, aes(x = Percent, y=tb_mean), color= "steelblue")+
  geom_line(data= bout2, aes(x = Percent, y=tb_mean), color = "gray20")+
  geom_vline(xintercept = 100)


#Looks good.


################################################################################
# Aligning bouts to standard % grid for direct comparison
################################################################################

# Convert to long format
bout1_long <- bout1 %>%
  select(-sample_nr) %>%
  pivot_longer(
    cols = -c(Percent, tb_mean),
    names_to = "Metabolite",
    values_to = "Abundance"
  ) %>%
  mutate(Bout = "Bout1")

bout2_long <- bout2 %>%
  select(-sample_nr) %>%
  pivot_longer(
    cols = -c(Percent, tb_mean),
    names_to = "Metabolite",
    values_to = "Abundance"
  ) %>%
  mutate(Bout = "Bout2")

# Combine bouts
bouts_combined <- bind_rows(bout1_long, bout2_long)

# Interpolate to standard % grid (every 2.5%) So that measurments snap to same overlapping percent
standard_percent <- seq(-5, 120, by = 2.5)

bouts_aligned <- bouts_combined %>%
  group_by(Metabolite, Bout) %>%
  arrange(Percent) %>%
  summarise(
    Aligned_abundance = list(approx(Percent, Abundance, xout = standard_percent)$y),
    Aligned_tb = list(approx(Percent, tb_mean, xout = standard_percent)$y),
    .groups = "drop"
  ) %>%
  unnest(c(Aligned_abundance, Aligned_tb)) %>%
  mutate(Percent = rep(standard_percent, n() / length(standard_percent)))

#check plot:

# Filter data for the metabolite of interest
plot_data <- bouts_aligned %>%
  filter(Metabolite == "278.17305/7.026") # Top accumulator: 185.2144/6.52

# Calculate dynamic ranges for flexible scaling
tb_range <- range(plot_data$Aligned_tb, na.rm = TRUE)
abund_range <- range(plot_data$Aligned_abundance, na.rm = TRUE)

# Calculate scaling factor and offset to map abundance to tb scale
# Formula: tb_scale = abund * scale_factor + offset
scale_factor <- diff(tb_range) / diff(abund_range)
offset <- tb_range[1] - abund_range[1] * scale_factor

# Add some padding to limits (5% on each side)
tb_padding <- diff(tb_range) * 0.05
tb_limits <- c(tb_range[1] - tb_padding, tb_range[2] + tb_padding)

plot_data %>%
  ggplot(aes(x = Percent, y = Aligned_tb, color = Bout)) +
  geom_vline(xintercept = 100, alpha = 0.6, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = tb_limits[1], ymax = Aligned_tb), alpha = 0.2) +
  geom_line(aes(y = Aligned_abundance * scale_factor + offset)) +
  geom_point(aes(y = Aligned_abundance * scale_factor + offset)) +
  scale_x_continuous(limits = c(-4, 120), breaks = seq(-5, 120, 5)) +
  scale_y_continuous(
    limits = tb_limits,
    name = "Aligned_tb",
    sec.axis = sec_axis(
      ~ (. - offset) / scale_factor,
      name = "Aligned_abundance"
    )
  )


# Looks like good staring point to assess accumulating or depleting factors
# The overlap is not perfect but that seems to be due to different dynamics over time
# at the thermal inflection points.


################################################################################
# Check between-bout agreement using the CCF-significant metabolites.
# Then take the median of the two bouts per metabolite, Mann-Kendall
# (https://www.statisticshowto.com/mann-kendall-trend-test/), Benjamini-Hochberg
# for multiple testing at FDR 5%.
################################################################################

#the ccf results:
ccf_results <- read.csv("./data/cross_correlation_permutation_results.csv")

#How many significant?
ccf_results%>%filter(significant)%>%nrow()

#Subset the data
ccf_results_filtered <- ccf_results%>% filter(significant)


# Filter to only significant CCF metabolites
sig_metabolites <- ccf_results_filtered$Metabolite

bouts_aligned_sigs <- bouts_aligned %>%
  filter(Metabolite %in% sig_metabolites) #comment out to include all mets

#no spelling errors causing droppouts?
length(unique(bouts_aligned_sigs$Metabolite)) #No, good


################################################################################
# ICC for between-bout agreement.
# Koo & Li (2016) PMC4913118. ICC(3,1): two-way mixed, single measurement, consistency.
# Lest have a sligly forgiving interperatin thershold values:
# <0.4 poor, 0.4-0.65 moderate, 0.65-0.9 good, >0.9 excellent.
################################################################################

# Prepare data in wide format for ICC calculation
# Each row = one timepoint (Percent), columns = Bout1 and Bout2 abundance
icc_data_wide <- bouts_aligned_sigs  %>%
  select(Metabolite, Percent, Bout, Aligned_abundance) %>%
  pivot_wider(
    names_from = Bout,
    values_from = Aligned_abundance
  )

# Calculate ICC for each metabolite
calculate_icc <- function(data) {
  # Need at least 3 paired observations
  data_clean <- data %>%
    filter(!is.na(Bout1) & !is.na(Bout2))
  
  if(nrow(data_clean) < 3) {
    return(tibble(
      ICC = NA_real_,
      ICC_lower = NA_real_,
      ICC_upper = NA_real_,
      ICC_category = NA_character_
    ))
  }
  
  # Prepare matrix for ICC: rows = subjects (timepoints), cols = raters (bouts)
  icc_matrix <- data_clean %>%
    select(Bout1, Bout2) %>%
    as.matrix()
  
  # Calculate ICC(3,1) - two-way mixed, single measurement, consistency
  tryCatch({
    icc_result <- icc(icc_matrix, model = "twoway", type = "consistency", unit = "single")
    
    # Categorize ICC following Koo & Li (2016) guidelines
    icc_val <- icc_result$value
    category <- case_when(
      icc_val < 0.4 ~ "Poor",
      icc_val < 0.65 ~ "Moderate",
      icc_val < 0.9 ~ "Good",
      TRUE ~ "Excellent"
    )
    
    tibble(
      ICC = icc_result$value,
      ICC_lower = icc_result$lbound,
      ICC_upper = icc_result$ubound,
      ICC_category = category
    )
  }, error = function(e) {
    tibble(
      ICC = NA_real_,
      ICC_lower = NA_real_,
      ICC_upper = NA_real_,
      ICC_category = NA_character_
    )
  })
}

# Apply ICC calculation to each metabolite
icc_results <- icc_data_wide %>%
  group_by(Metabolite) %>%
  group_modify(~ calculate_icc(.x)) %>%
  ungroup()


# Visualization: ICC distribution
p_icc <- ggplot(icc_results %>% filter(!is.na(ICC)), aes(x = ICC)) +
  geom_histogram(binwidth = 0.025, fill = "steelblue4", color = "gray20", alpha = 0.7) +
  geom_vline(xintercept = c(0.4, 0.65, 0.9), linetype = "dashed", color = "firebrick", alpha = 1) +
  annotate("text", x = 0.25, y = Inf, label = "Poor", vjust = 2, size = 3) +
  annotate("text", x = 0.47, y = Inf, label = "Moderate", vjust = 2, size = 3) +
  annotate("text", x = 0.7, y = Inf, label = "Good", vjust = 2, size = 3) +
  annotate("text", x = 0.97, y = Inf, label = "Excellent", vjust = 2, size = 3) +
  labs(
    title = "Between-Bout Agreement: ICC Distribution",
    subtitle = "ICC(3,1) - Two-way mixed, consistency",
    x = "Intraclass Correlation Coefficient (ICC)",
    y = "Count"
  ) +
  theme_classic() 

p_icc

ggsave("./figures/ICC_distribution_final.svg", plot = p_icc, width = 7, height = 6)


#interesting... a gradual falloff of ICC.. We might restrict top picks to ICC >0.4
# we can check to see how good/bad between bout performers do by using plot strting at line 330 

tibble(threshold = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) %>%
  mutate(count = map_int(threshold, ~sum(icc_results$ICC > ., na.rm = TRUE)))

# Sample 16 metabolites evenly across the full ICC range to show what different values look like
icc_range_examples <- icc_results %>%
  filter(!is.na(ICC)) %>%
  arrange(ICC) %>%
  mutate(bin = cut(ICC, breaks = c(-Inf, 0.0, 0.1, 0.2,0.3, 0.4, seq(0.45, 0.90, by = 0.05), Inf))) %>%
  group_by(bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  arrange(ICC) %>%
  mutate(label = factor(paste0(Metabolite, "\nICC = ", round(ICC, 2)), levels = paste0(Metabolite, "\nICC = ", round(ICC, 2))))

examples_icc <- bouts_aligned_sigs %>%
  filter(Metabolite %in% icc_range_examples$Metabolite) %>%
  left_join(icc_range_examples %>% select(Metabolite, ICC, label), by = "Metabolite") %>%
  ggplot(aes(x = Percent, y = Aligned_abundance, color = Bout)) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 0.7, na.rm = TRUE) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Bout1" = "firebrick", "Bout2" = "steelblue4")) +
  scale_x_continuous(breaks = seq(0, 100, 50)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance", color = NULL) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))

#save for sub fig
ggsave("./figures/examples_icc.svg", plot = examples_icc, width = 7, height = 6)


# Sampling through a few times it looks like it might be worth keeping metabolites above ICC 0.2-0.3, 
# maybe above >0.3
# But the best is consistetly above 0.5. Keep this in mind if filter further along.. 



################################################################################
# For the trend analysis we apperntly need the median value
# Calculate median abundance across both bouts for each metabolite at each % 
# (https://www.statisticshowto.com/mann-kendall-trend-test/)
################################################################################

bouts_aligned_sigs_median <- bouts_aligned_sigs  %>%
  group_by(Metabolite, Percent) %>%
  summarise(
    Median_abundance = median(Aligned_abundance, na.rm = TRUE),
    Bout1_abundance = Aligned_abundance[Bout == "Bout1"][1],
    Bout2_abundance = Aligned_abundance[Bout == "Bout2"][1],
    Median_tb = median(Aligned_tb, na.rm = TRUE),
    Bout1_tb = Aligned_tb[Bout == "Bout1"][1],
    Bout2_tb = Aligned_tb[Bout == "Bout2"][1],
    .groups = "drop"
  )

# Check: top 6 by CCF score, individual bouts + median trace, Tb as ribbon
top_ccf_met <- ccf_results_filtered %>% arrange(desc(Max_Correlation)) %>% slice_head(n = 6) %>% pull(Metabolite)

bouts_aligned_sigs_median %>%
  filter(Metabolite %in% top_ccf_met) %>%
  ggplot(aes(x = Percent)) +
  geom_ribbon(aes(ymin = -2.5, ymax = (Median_tb / 10) - 3.5), fill = "steelblue4", alpha = 0.2) +
  geom_line(aes(y = Bout1_abundance, color = "Bout 1"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance, color = "Bout 2"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Median_abundance, color = "Median"), linewidth = 1, na.rm = TRUE) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Bout 1" = "firebrick", "Bout 2" = "steelblue4", "Median" = "gray15")) +
  scale_x_continuous(breaks = seq(0, 100, 50)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# Middle CCF
mid_ccf_met <- ccf_results_filtered %>%
  arrange(desc(Max_Correlation)) %>%
  slice(round(n()/2) + (-2:3)) %>%
  pull(Metabolite)

bouts_aligned_sigs_median %>%
  filter(Metabolite %in% mid_ccf_met) %>%
  ggplot(aes(x = Percent)) +
  geom_ribbon(aes(ymin = -2.5, ymax = (Median_tb / 10) - 3.5), fill = "steelblue4", alpha = 0.2) +
  geom_line(aes(y = Bout1_abundance, color = "Bout 1"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance, color = "Bout 2"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Median_abundance, color = "Median"), linewidth = 1, na.rm = TRUE) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Bout 1" = "firebrick", "Bout 2" = "steelblue4", "Median" = "gray15")) +
  scale_x_continuous(breaks = seq(0, 100, 50)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance", color = NULL, title = "Middle CCF") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# Lowest CCF
low_ccf_met <- ccf_results_filtered %>% arrange(Max_Correlation) %>% slice_head(n = 6) %>% pull(Metabolite)

bouts_aligned_sigs_median %>%
  filter(Metabolite %in% low_ccf_met) %>%
  ggplot(aes(x = Percent)) +
  geom_ribbon(aes(ymin = -2.5, ymax = (Median_tb / 10) - 3.5), fill = "steelblue4", alpha = 0.2) +
  geom_line(aes(y = Bout1_abundance, color = "Bout 1"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance, color = "Bout 2"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Median_abundance, color = "Median"), linewidth = 1, na.rm = TRUE) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Bout 1" = "firebrick", "Bout 2" = "steelblue4", "Median" = "gray15")) +
  scale_x_continuous(breaks = seq(0, 100, 50)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance", color = NULL, title = "Lowest CCF") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")



################################################################################
# Autocorrelation check, we ended up using Hamed-Rao modified MK.
# So picked up this because though perhaps too many were significant... Autocorrelation 
# inflates test-statistic (but not coefficients), so we overestimate n sigs. 
# so we need to chack for autocorrelation, that is if an observation is
# correlating with it self in a lagged position (like cross-correlation).
# This is very typical for timeseries data. The observations are likely to not be 
# independent, but "history dependent". The points adjacent are likely to predict the obsevariton point. 

# Standard MK assumes independent obs. good for us packange for handling this library(modifiedmk) 

# The two plots below check that on the
# median-aggregated 0-100 % series. If lag-1 residual ACF sits outside the
# white-noise band the naive MK p-values are likely inflated/overestimated, and the
# Hamed-Rao (1998) variance inflation inside mmkh() handles it.
################################################################################

compute_acf_residuals <- function(data, lags = 1:10) {
  d <- data %>% filter(!is.na(Median_abundance)) %>% arrange(Percent)
  if(nrow(d) < 6) return(tibble(lag = lags, acf_val = NA_real_))
  resids <- resid(lm(Median_abundance ~ Percent + I(Percent^2), data = d))
  ac <- as.numeric(acf(resids, lag.max = max(lags), plot = FALSE)$acf)
  tibble(lag = lags, acf_val = ac[lags + 1])
}

acf_by_met <- bouts_aligned_sigs_median %>%
  filter(Percent >= 0, Percent <= 100) %>%
  group_by(Metabolite) %>%
  group_modify(~ compute_acf_residuals(.x)) %>%
  ungroup()

lag1_acf <- acf_by_met %>% filter(lag == 1)
n_grid   <- bouts_aligned_sigs_median %>%
  filter(Percent >= 0, Percent <= 100) %>% pull(Percent) %>% unique() %>% length()
ac_ci    <- 2 / sqrt(n_grid)   # ~95 % white-noise band

# Pre-compute N/N* per metabolite for the effective-n plot
nn_star_only <- function(data) {
  d <- data %>% filter(!is.na(Median_abundance)) %>% arrange(Percent)
  if(nrow(d) < 4) return(tibble(n_over_n_star = NA_real_))
  res <- tryCatch(mmkh(d$Median_abundance), error = function(e) NULL)
  if(is.null(res)) return(tibble(n_over_n_star = NA_real_))
  tibble(n_over_n_star = unname(res["N/N*"]))
}

ac_diagnostic <- bouts_aligned_sigs_median %>%
  filter(Percent >= 0, Percent <= 100) %>%
  group_by(Metabolite) %>%
  group_modify(~ nn_star_only(.x)) %>%
  ungroup()

# Plot 1: lag-1 residual ACF distribution across CCF-significant metabolites
ggplot(lag1_acf, aes(x = acf_val)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = c(-ac_ci, ac_ci), linetype = "dotted", colour = "firebrick") +
  geom_vline(xintercept = median(lag1_acf$acf_val, na.rm = TRUE), linetype = "dashed") +
  labs(
    title = "Lag-1 residual autocorrelation across CCF-significant metabolites",
    subtitle = sprintf("Residuals after OLS linear+quadratic detrend; red dotted = white-noise 95%% band (±%.2f)", ac_ci),
    x = "Lag-1 ACF", y = "Metabolite count"
  ) +
  theme_bw()

# Plot 2: Hamed-Rao effective-n ratio N/N*
ggplot(ac_diagnostic, aes(x = n_over_n_star)) +
  geom_histogram(binwidth = 0.2, fill = "darkorange", colour = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 6)) +
  labs(
    title = "Hamed-Rao effective-n ratio (N / N*)",
    subtitle = ">1 means nominal n overstates information; MK variance is inflated",
    x = "N / N*", y = "Metabolite count"
  ) +
  theme_bw()

################################################################################
# Mann-Kendall trend test (non-parametric).
# Hamed-Rao (1998) AR-corrected MK via modifiedmk::mmkh. Senss slope as effect size.
# 
# One-directional restiction: Hamed-Rao is meant to PENALISE positive autocorrelation
# (inflate the variance, raise the p-value). On weak/noisy median series a couple
# of negative-lag rank autocorrelations can breach the white-noise band by chance,
# which instead DEFLATES the variance (when n/n < 1) and manufactures significance
# (Zc balloons, in the extreme the variance goes negative and p = NaN,
# e.g Np-020400 and 273.15769/3.877). We therefore only let the correction raise
# the p-value, never lower it: when the Hamed-Rao p is below the standard MK p
# (equivalently n/n < 1) or is NaN, we fall back to the standard MK p-value.
# Raw Hamed-Rao p is retained as p_value_mk_hr and the fallback flagged in
# hr_reverted for transparency.
################################################################################

mann_kendall_test <- function(data) {
  data_clean <- data %>%
    filter(!is.na(Median_abundance)) %>%
    arrange(Percent)
  
  abundance <- data_clean$Median_abundance
  n <- length(abundance)
  
  empty <- tibble(
    tau             = NA_real_,
    sen_slope       = NA_real_,
    p_value_mk      = NA_real_,
    p_value_mk_hr   = NA_real_,
    p_value_mk_orig = NA_real_,
    n_over_n_star   = NA_real_,
    hr_reverted     = NA
  )
  
  if(n < 4) return(empty)
  
  res <- tryCatch(mmkh(abundance), error = function(e) NULL)
  if(is.null(res)) return(empty)
  
  p_hr   <- unname(res["new P-value"])   # Hamed-Rao AR-corrected
  p_orig <- unname(res["old P.value"])   # standard MK
  
  # Clamp: never allow the correction to lower the p-value (see header note).
  reverted <- is.na(p_hr) || (!is.na(p_orig) && p_hr < p_orig)
  p_use    <- if (reverted) p_orig else p_hr
  
  tibble(
    tau             = unname(res["Tau"]),
    sen_slope       = unname(res["Sen's slope"]),
    p_value_mk      = p_use,
    p_value_mk_hr   = p_hr,
    p_value_mk_orig = p_orig,
    n_over_n_star   = unname(res["N/N*"]),
    hr_reverted     = reverted
  )
}

# Apply Mann-Kendall to each metabolite (0-100% only)
mk_results <- bouts_aligned_sigs_median %>%
  filter(Percent >= 0, Percent <= 100) %>%
  group_by(Metabolite) %>%
  group_modify(~ mann_kendall_test(.x)) %>%
  ungroup()

# Check plot: top 6 metabolites by absolute tau
mk_top_slope_met <- mk_results %>% arrange(desc(abs(tau))) %>% slice_head(n = 6) %>% pull(Metabolite)

bouts_aligned_sigs_median %>%
  filter(Metabolite %in% mk_top_slope_met) %>%
  ggplot(aes(x = Percent, y = Median_abundance)) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  geom_smooth(data = . %>% filter(Percent >= 0, Percent <= 100),
              method = "loess", se = TRUE, color = "firebrick", linewidth = 0.8) +
  # raw per-bout aligned values behind the median
  geom_point(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(colour = "Median"), size = 1, alpha = 0.75) +
  
  facet_wrap(~ Metabolite, scales = "free_y") +
  scale_colour_manual(
    name   = NULL,
    values = c("Bout 1" = "forestgreen", "Bout 2" = "firebrick2", "Median" = "gray20")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")







################################################################################
# Combine results and apply Benjamini-Hochberg FDR correction
################################################################################
trend_results <- mk_results %>%
  mutate(
    p_adj_mk   = p.adjust(p_value_mk, method = "BH"),
    sig_mmk_adjusted     = !is.na(p_adj_mk) & p_adj_mk < 0.05,
    trend_type = case_when(
      sig_mmk_adjusted & sen_slope > 0 ~ "Monotonic-up",
      sig_mmk_adjusted & sen_slope < 0 ~ "Monotonic-down",
      TRUE                   ~ "Flat"
    )
  )

str(trend_results)

# Add CCF and ICC info
trend_results <- trend_results %>%
  left_join(
    ccf_results_filtered %>% 
      select(Metabolite, spear_rho = Max_Correlation, Lag),
    by = "Metabolite"
  ) %>%
  left_join(
    icc_results %>% 
      select(Metabolite, ICC, ICC_category),
    by = "Metabolite"
  ) %>%
  relocate(Metabolite, tau, sen_slope, p_adj_mk, ICC, spear_rho) %>%
  arrange(desc(abs(tau)))


################################################################################
# Overview: distributions of the per-metabolite trend metrics
################################################################################
trend_results %>%
  mutate(
    abs_tau = abs(tau),
    abs_sen_slope = abs(sen_slope)
  ) %>%
  select(Metabolite, abs_tau, abs_sen_slope, p_adj_mk, ICC, spear_rho) %>%
  pivot_longer(-Metabolite, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("abs_tau", "abs_sen_slope", "p_adj_mk", "ICC", "spear_rho"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 40, fill = "steelblue4", colour = "white", alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  labs(
    x = NULL,
    y = "Metabolite count",
    title = "Distributions of trend metrics across metabolites"
  ) +
  theme_bw(base_size = 10)

# Defining filter thershold operations. 
# p-value BH adjusted: 0.05

# Tau: No good sources, various say around >.3 is good to consider or make similar to Persons rank correlation levels >.3 moderate
# tau dist suggest .3 as good threshold filtering ca 50% of data:
median(abs(trend_results$tau)) # jupp, lets use this level 

# ICC we already considered 0.3 as a forgiving filtering point
median(abs(trend_results$ICC)) #.22 as median but we already eval .3 as conservative so we use that. 

# Effect size, sen's slope, we might want to filter on this as well.lets look at dist after we filer on above: 
trend_results %>%
  mutate(
    abs_tau = abs(tau),
    abs_sen_slope = abs(sen_slope)
  ) %>%
  filter(p_adj_mk < 0.05, abs_tau > 0.3, ICC > 0.4)%>%
  ggplot(aes(x = abs_sen_slope)) +
  geom_histogram(bins = 40, fill = "steelblue4", colour = "white", alpha = 0.85, na.rm = TRUE)

trend_results %>%
  mutate(
    abs_tau = abs(tau),
    abs_sen_slope = abs(sen_slope)
  ) %>%
  filter(p_adj_mk < 0.05, abs_tau > 0.3, ICC > 0.4) %>%
  pull(abs_sen_slope) %>%
  median()

#and without filter:
trend_results %>%
  mutate(
    abs_tau = abs(tau),
    abs_sen_slope = abs(sen_slope)
  )%>%
  pull(abs_sen_slope) %>%
  median()

# Lets filter on the approx median at 0.01. it means very little since this is change in zNorm values per 2.5 percent 
# on a range -2 to 2 changes at best. so a 0.01 means we have a 0.4 change over 0-100% which is a generous filter level
# for change over the torpor bout to pass. i.e if metabolite dont pass abs(sen_slope)>0.01 its not going to be a 
# meaningful candidate to consider. 

# So whats the percetange left if we implement all filters: 
trend_results_top_perf <- trend_results %>%
  filter(p_adj_mk < 0.05, abs(tau) > 0.3, ICC > 0.4, abs(sen_slope)>0.01)


nrow(trend_results_top_perf)/ nrow(trend_results) # the top 15% remains!! thats a good filter level.  



################################################################################
# Monotonic-shape descriptor: "elbow-early" vs "elbow-late" vs "linear".
# Describes where the monotonic change happens. All front-loaded:
#   change_center = bout % at which the cumulative |change| reaches its halfway
#                   point (where the bulk of the change happens),
#   concentration = fraction of total |change| packed into the steepest 25 % of
#                   the axis (how abrupt; 0.25 = perfectly even). Reported only.
# Classify (thresholds on change_center; late_edge sets the front & rear bands):
#   elbow-early   : < 20 %                       -> sharp, very front-loaded corner.
#   elbow-late    : 20 % to < late_edge (~30 %)  -> front-loaded elbow, a bit later.
#   near-linear   : late_edge to 100-late_edge   -> centred change (line or sigmoid).
#   elbow-reverse : > 100-late_edge (~70 %)      -> rapid LATE change (reverse elbow).
# Tune with `early_edge` / `late_edge` (the band boundaries in %).
# NB: intended for monotonic trends: a curve with a reversal (peak/dip)
# still gets a label but it is not meaningful there.
################################################################################

shape_profile <- function(data, span = 0.75, early_edge = 100/5, late_edge = 100/3.3) {
  d <- data %>% filter(!is.na(Median_abundance)) %>% arrange(Percent)
  empty <- tibble(change_center = NA_real_, concentration = NA_real_, shape_class = NA_character_)
  if (nrow(d) < 6) return(empty)
  
  fit <- tryCatch(loess(Median_abundance ~ Percent, data = d, span = span),
                  error = function(e) NULL)
  if (is.null(fit)) return(empty)
  
  xg <- seq(min(d$Percent), max(d$Percent), length.out = 200)
  yg <- as.numeric(predict(fit, xg))
  dy <- diff(yg); ad <- abs(dy); tot <- sum(ad)
  if (!is.finite(tot) || tot == 0)
    return(tibble(change_center = 50, concentration = 0.25, shape_class = "linear"))
  
  cumf   <- cumsum(ad) / tot
  xm     <- (xg[-1] + xg[-length(xg)]) / 2
  center <- approx(cumf, xm, xout = 0.5, ties = "ordered")$y
  center <- (center - min(xg)) / (max(xg) - min(xg)) * 100         # -> 0-100 %
  k      <- ceiling(0.25 * length(ad))
  conc   <- sum(sort(ad, decreasing = TRUE)[1:k]) / tot            # concentration/abruptness (reported)
  
  # front elbow split early/late; last third = reverse elbow; middle = near-linear
  cls <- if (center < early_edge)           "elbow-early"
  else if (center < late_edge)       "elbow-late"
  else if (center > 100 - late_edge) "elbow-reverse"
  else                               "near-linear"
  
  tibble(change_center = center, concentration = conc, shape_class = cls)
}

shape_results <- bouts_aligned_sigs_median %>%
  filter(Percent >= 0, Percent <= 100) %>%
  group_by(Metabolite) %>%
  group_modify(~ shape_profile(.x)) %>%
  ungroup()

table(shape_results$shape_class)   # elbow-early ~299, elbow-late ~162, near-linear ~284, elbow-reverse ~16

# attach shape descriptor + MSI tag to trend_results
metabolite_tags <- read.csv("./data/metabolite_tags.csv", check.names = FALSE, stringsAsFactors = FALSE) %>%
  distinct(species, .keep_all = TRUE)

trend_results_msi <- trend_results %>%
  left_join(shape_results %>% select(Metabolite, change_center, shape_class),
            by = "Metabolite") %>%
  left_join(metabolite_tags, by = c("Metabolite" = "species"))%>%
  relocate(tags, Metabolite, tau, sen_slope, p_adj_mk, ICC, spear_rho, shape_class, change_center) %>%
  arrange(desc(abs(tau)))%>%
  mutate(tags = ifelse(is.na(tags), "F", tags))

# How many are MSI A are there?
trend_results_msi %>%
  filter(tags == "A")%>%
  nrow()


################################################################################
# Plotting the single top performer of each shape and direction
# (elbow-earl,linear accumulating or depleating) 
################################################################################
shape_pick <- trend_results_msi %>%
  filter(!is.na(shape_class), 
         p_adj_mk < 0.05, abs(tau) > 0.3, ICC > 0.4, abs(sen_slope)>0.01) %>%
  mutate(direction = if_else(tau > 0, "Accumulating", "Depleting")) %>%
  select(Metabolite, tau, ICC, shape_class, direction) %>%
  mutate(shape_class = factor(shape_class, levels = c("elbow-early","elbow-late", "near-linear")),
         direction   = factor(direction,   levels = c("Accumulating", "Depleting"))) %>%
  group_by(shape_class, direction) %>%
  mutate(perf_rank = rank(-abs(tau)) + rank(-ICC)) %>%   # top performer = strongest & most reproducible
  slice_min(perf_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(direction, shape_class) %>%
  mutate( panel = factor(
    sprintf("%s | %s\n%s\ntau = %.2f | ICC = %.2f", shape_class, direction, Metabolite, tau, ICC),
    levels = sprintf("%s | %s\n%s\ntau = %.2f | ICC = %.2f", shape_class, direction, Metabolite, tau, ICC)
  ))

Plot_shapes <- bouts_aligned_sigs_median %>%
  filter(Metabolite %in% shape_pick$Metabolite) %>%
  left_join(shape_pick %>% select(Metabolite, panel), by = "Metabolite") %>%
  ggplot(aes(x = Percent, y = Median_abundance)) +
  geom_line(aes(y = (Median_tb / 10) - 2.5), colour = "steelblue4", linewidth = 0.4, na.rm = TRUE) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  geom_smooth(data = . %>% filter(Percent >= 0, Percent <= 100),
              method = "loess", se = TRUE, color = "firebrick4", linewidth = 0.8) +
  geom_point(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout1_abundance,  colour = "Bout 1"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance,  colour = "Bout 2"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(colour = "Median"), size = 1, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_colour_manual(
    name   = NULL,
    values = c("Bout 1" = "#AA4499", "Bout 2" = "#117733", "Median" = "gray20")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

Plot_shapes

ggsave("./figures/Elbows_linear_top_final.svg", plot = Plot_shapes, width = 10, height = 8)


################################################################################
# Top performers of each shape class, as three separate plots. To show a
# representative direction mix (shape x direction is skewed), take the top 4
# ACCUMULATING + top 4 DEPLETING per shape (row 1 = accumulating, row 2 = depleting).
# "Top" = best combined rank of trend strength (|tau|) + reproducibility (ICC),
# within the top-performer gate (p<0.05 & |tau|>0.3 & ICC>0.4 & |sen|>0.01).
################################################################################

# top-4-per-direction picker for one shape class (8 total: 4 accum + 4 deplete)
top8_for <- function(sc) {
  trend_results_msi %>%
    filter(shape_class == sc,
           p_adj_mk < 0.05, abs(tau) > 0.3, ICC > 0.4, abs(sen_slope) > 0.01) %>%
    mutate(direction = if_else(tau > 0, "Accumulating", "Depleting")) %>%
    group_by(direction) %>%
    mutate(perf_rank = rank(-abs(tau)) + rank(-ICC)) %>%
    slice_min(perf_rank, n = 4, with_ties = FALSE) %>%      # best 4 in each direction
    ungroup() %>%
    mutate(direction = factor(direction, levels = c("Accumulating", "Depleting"))) %>%
    arrange(direction, perf_rank) %>%                       # accum row then deplete row
    mutate(panel = factor(
      sprintf("%s\n%s | tau=%.2f | ICC=%.2f", Metabolite, direction, tau, ICC),
      levels = sprintf("%s\n%s | tau=%.2f | ICC=%.2f", Metabolite, direction, tau, ICC)))
}

# plot builder for one shape class
plot_top8 <- function(sc, ncol = 4) {
  pick <- top8_for(sc)
  bouts_aligned_sigs_median %>%
    filter(Metabolite %in% pick$Metabolite) %>%
    left_join(pick %>% select(Metabolite, panel), by = "Metabolite") %>%
    ggplot(aes(x = Percent, y = Median_abundance)) +
    geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
    geom_smooth(data = . %>% filter(Percent >= 0, Percent <= 100),
                method = "loess", se = TRUE, color = "firebrick", linewidth = 0.8) +
    geom_point(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
    geom_line(aes(y = Bout1_abundance,  colour = "Bout 1"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
    geom_point(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
    geom_line(aes(y = Bout2_abundance,  colour = "Bout 2"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
    geom_point(aes(colour = "Median"), size = 1, alpha = 0.75, na.rm = TRUE) +
    facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
    scale_colour_manual(
      name   = NULL,
      values = c("Bout 1" = "forestgreen", "Bout 2" = "firebrick2", "Median" = "gray20")
    ) +
    scale_x_continuous(breaks = seq(0, 100, 25)) +
    labs(x = "Bout progression (%)", y = "Normalised abundance",
         title = paste0("Top 4 accumulating + 4 depleting: ", sc)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
}

Plot_early  <- plot_top8("elbow-early")
Plot_late   <- plot_top8("elbow-late")
Plot_linear <- plot_top8("near-linear")

Plot_early
Plot_late
Plot_linear


################################################################################
# Reverse-elbow check done in the function so we can check if there are any. 
################################################################################
reverse_pick <- trend_results_msi %>%
  filter(shape_class == "elbow-reverse", p_adj_mk < 0.05, abs(tau) > 0.3, ICC > 0.4, abs(sen_slope)>0.01)

head(reverse_pick)
#no


################################################################################
# Non-monotonic check: are there any clear U-shapes / humps among the reproducible,
# strongly-trending metabolites? Localised confirmation (no columns added to
# trend_results_msi / CSV / app). reversal_prominence() loesses the median over
# 0-100 %
################################################################################

reversal_prominence <- function(data, span = 0.75) {
  d <- data %>% filter(!is.na(Median_abundance), Percent >= 0, Percent <= 100) %>% arrange(Percent)
  if (nrow(d) < 6) return(NA_real_)
  fit <- tryCatch(loess(Median_abundance ~ Percent, data = d, span = span), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  xg  <- seq(min(d$Percent), max(d$Percent), length.out = 200)
  yg  <- as.numeric(predict(fit, xg))
  rng <- diff(range(yg)); if (!is.finite(rng) || rng == 0) return(0)
  s <- sign(diff(yg)); s[s == 0] <- NA
  for (i in seq_along(s)) if (is.na(s[i]) && i > 1) s[i] <- s[i - 1]   # carry sign forward over flats
  turn <- which(diff(s) != 0) + 1                                     # interior turning points
  if (length(turn) == 0) return(0)
  ext <- c(1, turn, length(yg)); ey <- yg[ext]; ji <- 2:(length(ext) - 1)
  max(pmin(abs(ey[ji] - ey[ji - 1]), abs(ey[ji] - ey[ji + 1])) / rng)  # dominant reversal / range
}

nonmono_amp_tbl <- bouts_aligned_sigs_median %>%
  group_by(Metabolite) %>%
  summarise(nonmono_amp = reversal_prominence(pick(everything())), .groups = "drop")

nonmono_check <- trend_results_msi %>%
  left_join(nonmono_amp_tbl, by = "Metabolite") %>%
  filter(p_adj_mk < 0.05, abs(tau) > 0.3, ICC > 0.4, nonmono_amp >= 0.25) %>%
  arrange(desc(nonmono_amp)) %>%
  select(Metabolite, nonmono_amp, tau, ICC, p_adj_mk, shape_class)

head(nonmono_check)


# Plot the non-monotonic metabolites that pass the filter (nonmono_check)
nonmono_panel <- nonmono_check %>%
  mutate(panel = factor(
    sprintf("%s\ntau = %.2f | ICC = %.2f | amp = %.2f", Metabolite, tau, ICC, nonmono_amp),
    levels = sprintf("%s\ntau = %.2f | ICC = %.2f | amp = %.2f", Metabolite, tau, ICC, nonmono_amp)))

bouts_aligned_sigs_median %>%
  filter(Metabolite %in% nonmono_panel$Metabolite) %>%
  left_join(nonmono_panel %>% select(Metabolite, panel), by = "Metabolite") %>%
  ggplot(aes(x = Percent, y = Median_abundance)) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  geom_smooth(data = . %>% filter(Percent >= 0, Percent <= 100),
              method = "loess", se = TRUE, color = "firebrick", linewidth = 0.8) +
  geom_point(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout1_abundance,  colour = "Bout 1"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance,  colour = "Bout 2"), linewidth = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_point(aes(colour = "Median"), size = 1, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_colour_manual(
    name   = NULL,
    values = c("Bout 1" = "forestgreen", "Bout 2" = "firebrick2", "Median" = "gray20")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance",
       title = "Non-monotonic mets") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

#all of these do not look worth considering further. 



################################################################################
################################################################################

#Write to disk
write.csv(bouts_aligned_sigs_median, "./data/bouts_aligned_sigs_median.csv", row.names = FALSE)
write.csv(trend_results_msi, "./data/trend_results_msi.csv", row.names = FALSE)

################################################################################
################################################################################
colnames(trend_results_msi)



#What are some numbers: 
trend_results_msi %>%
  summarise(
    total = n(),
    sig_any = sum(sig_mmk_adjusted),
    sig_accumulating = sum(sig_mmk_adjusted & sen_slope > 0),
    sig_depleting = sum(sig_mmk_adjusted & sen_slope < 0),
    pct_sig = round(sig_any / total * 100, 1),
    n_early_elbow = sum(sig_mmk_adjusted & shape_class == "elbow-early"),
    n_late_elbow = sum(sig_mmk_adjusted & shape_class == "elbow-late"),
    n_near_linear = sum(sig_mmk_adjusted & shape_class == "near-linear")
  )


# Cross-tabulate the candidate filters: ICC, |tau|, |sen's slope| and sig_mmk_adjusted.
# Uses the top-performer gate thresholds (ICC>0.4, |tau|>0.3, |sen|>0.01, p_adj<0.05),
# so `all_four` == the top-performer set.
icc_sig_overlap <- trend_results_msi %>%
  mutate(
    sig      = !is.na(sig_mmk_adjusted) & sig_mmk_adjusted,     # p_adj_mk < 0.05
    icc_pass = !is.na(ICC)       & ICC > 0.4,
    tau_pass = !is.na(tau)       & abs(tau) > 0.3,
    sen_pass = !is.na(sen_slope) & abs(sen_slope) > 0.01
  ) %>%
  summarise(
    total_tested  = n(),
    sig           = sum(sig),
    icc_over_0.4  = sum(icc_pass),
    tau_over_0.3  = sum(tau_pass),
    sen_over_0.01 = sum(sen_pass),
    sig_and_icc   = sum(sig & icc_pass),                        # reproducible AND significant
    all_four      = sum(sig & icc_pass & tau_pass & sen_pass),  # top performers
    pct_all_four  = round(100 * sum(sig & icc_pass & tau_pass & sen_pass) / n(), 1)
  )

icc_sig_overlap

# The high-confidence shortlist itself (reproducible AND significantly trending)
trend_results_sig_shortlist <- trend_results_msi %>%
  #filter(!is.na(ICC), ICC > 0.4, sig_mmk_adjusted) %>%
  arrange(desc(abs(tau))) %>%
  select(tags, Metabolite, tau, sen_slope, p_adj_mk, spear_rho, ICC, ICC_category, trend_type, shape_class)


trend_results_sig_shortlist%>%
  head()


trend_results%>%
  head()


write.csv(trend_results_sig_shortlist, "./data/trend_results_sig_msi_shortlist.csv", row.names = FALSE)
writexl::write_xlsx(trend_results_sig_shortlist, "./data/trend_results_sig_msi_shortlist.xlsx")


################################################################################
# Visualising
################################################################################

# Volcano plot: Sen's slope vs -log10(adjusted p-value)
p1 <- trend_results_msi %>%
  mutate(log_p = -log10(p_adj_mk)) %>%
  ggplot(aes(x = sen_slope, y = log_p, text = Metabolite)) +
  geom_point(aes(color = trend_type), alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0,            linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c(
    "Flat"                = "grey60",
    "Monotonic-up"        = "firebrick",
    "Monotonic-down"      = "steelblue3"
  )) +
  labs(
    title = "Within-torpor trend: Sen's slope vs significance",
    x     = "Sen's slope (z per 2.5% step)",
    y     = "-log10(BH-adjusted p-value)",
    color = "Trend type"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

plotly::ggplotly(p1, tooltip = "text")

#Lets lable the top 20-30:

top_mets <- trend_results_msi %>%
  filter(sig_mmk_adjusted) %>%
  mutate(rank_score = rank(p_adj_mk) + rank(-abs(sen_slope))) %>%
  arrange(rank_score) %>%
  slice_head(n = 30)

p1 <- trend_results_msi %>%
  mutate(log_p = -log10(p_adj_mk)) %>%
  ggplot(aes(x = sen_slope, y = log_p)) +
  geom_point(aes(color = trend_type), alpha = 0.6, size = 1.8) +
  geom_text_repel(
    data = top_mets %>% mutate(log_p = -log10(p_adj_mk)),
    aes(label = Metabolite),
    size          = 2.8,
    max.overlaps  = 20,
    box.padding   = 0.4,
    segment.color = "grey40",
    segment.size  = 0.3
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0,            linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c(
    "Flat"               = "grey60",
    "Monotonic-up"       = "firebrick",
    "Monotonic-down"     = "steelblue3"
  )) +
  labs(
    title = "Within-torpor trend: Sen's slope vs significance \n
                                                               Depleating - Accumulating",
    x     = "Sen's slope (z per 2.5% step)",
    y     = "-log10(BH-adjusted p-value)",
    color = "Trend type"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

p1


head(trend_results_msi)



################################################################################
# Volcano coloured by between-bout reproducibility (ICC), shape = elbow/shape class.
################################################################################

volcano_icc_df <- trend_results_msi %>%
  mutate(
    log_p     = -log10(p_adj_mk),
    highlight = !is.na(ICC) & ICC > 0.4 & !is.na(tau) & abs(tau) > 0.3,   # reproducible AND trending
    ICC_bin   = if_else(
      highlight,
      as.character(cut(ICC, breaks = c(0.40, 0.75, 0.85, Inf),
                       labels = c("0.40-0.65", "0.65-0.9", ">= 0.9"), right = FALSE)),
      "gray"
    ),
    ICC_bin   = factor(ICC_bin, levels = c("gray", "0.40-0.65", "0.75-0.9", ">= 0.9")),
    shape_grp = factor(coalesce(shape_class, "unclassified"),
                       levels = c("elbow-early", "elbow-late", "near-linear",
                                  "elbow-reverse", "unclassified"))
  ) %>%
  arrange(highlight)   # draw highlighted points last (on top)

icc_cols <- c(
  "gray"      = "grey75",
  "0.40-0.74" = "darkolivegreen4",
  "0.75-0.84" = "darkorange2",
  ">= 0.85"   = "purple3"
)

shape_vals <- c("elbow-early" = 17, "elbow-late" = 15, "near-linear" = 16,
                "elbow-reverse" = 18, "unclassified" = 4)



volcano_icc_labels <- volcano_icc_df %>%
  filter(highlight) %>%
  mutate(rank_score = rank(p_adj_mk) + rank(-abs(sen_slope))) %>%
  arrange(rank_score) %>%
  slice_head(n = 30)

p_volcano_icc <- ggplot(volcano_icc_df,
                        aes(x = sen_slope, y = log_p,
                            colour = ICC_bin, shape = shape_grp, alpha = highlight)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0,            linetype = "dashed", alpha = 0.5) +
  geom_point(size = 2) +
  geom_text_repel(
    data          = volcano_icc_labels,
    aes(label = Metabolite),
    inherit.aes   = TRUE,
    alpha         = 1,
    size          = 2.8,
    max.overlaps  = 20,
    box.padding   = 0.4,
    segment.color = "grey40",
    segment.size  = 0.3,
    show.legend   = FALSE
  ) +
  scale_colour_manual(values = icc_cols, name = "ICC (between-bout)", drop = FALSE) +
  scale_shape_manual(values = shape_vals, name = "Shape (monotonic shape)") +
  scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3), guide = "none") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 16)),
         shape  = guide_legend(override.aes = list(alpha = 1, colour = "gray30"))) +
  labs(
    title    = "Within-torpor trend: reproducible (ICC>0.4) & trending (|tau|>0.3) highlighted",
    subtitle = "Grey = ICC <= 0.4 or |tau| < 0.3; colour = ICC bin; shape = monotonic shape class",
    x = "Sen's slope (z per 2.5% step)",
    y = "-log10(BH-adjusted MK p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

p_volcano_icc


ggsave("./figures/volcano_final.svg", plot = p_volcano_icc, width = 10, height = 9)



################################################################################


str(bouts_aligned_sigs_median)

Met <- c("Adenosine", "Tyrosine", "Phenylalanine", "Cytosine")

# strip label with rho (CCF) and ICC pulled from trend_results_msi
fig_lab <- trend_results_msi %>%
  filter(Metabolite %in% Met) %>%
  mutate(Metabolite = factor(Metabolite, levels = Met)) %>%
  arrange(Metabolite) %>%
  transmute(Metabolite = as.character(Metabolite),
            panel = sprintf("%s\nrho = %.2f | ICC = %.2f", Metabolite, spear_rho, ICC)) %>%
  mutate(panel = factor(panel, levels = panel))

fig <- bouts_aligned_sigs_median %>%
  filter(Metabolite %in% Met) %>%
  left_join(fig_lab, by = "Metabolite") %>%
  ggplot(aes(x = Percent, y = Median_abundance)) +
  geom_line(aes(y = (Median_tb / 10) - 2.5), colour = "steelblue4", linewidth = 0.4, na.rm = TRUE) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", alpha = 0.4) +
  geom_smooth(data = . %>% filter(Percent >= 0, Percent <= 100),
              method = "loess", se = TRUE, color = "firebrick4", linewidth = 0.8) +
  #geom_line(linewidth = 1.0, alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = Bout1_abundance, colour = "Bout 1"), size = 0.5, alpha = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Bout1_abundance,  colour = "Bout 1"), linewidth = 0.5, alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = Bout2_abundance, colour = "Bout 2"), size = 0.5, alpha = 0.5, na.rm = TRUE) +
  geom_line(aes(y = Bout2_abundance,  colour = "Bout 2"), linewidth = 0.5, alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(colour = "Median"), size = 1, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(
    name   = NULL,
    values = c("Bout 1" = "#AA4499", "Bout 2" = "#117733", "Median" = "gray20")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(x = "Bout progression (%)", y = "Normalised abundance") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

fig

ggsave("./figures/traces_as_tyr_pln_cyt.svg", plot = fig3d, width = 15, height = 5)



################################################################################
# If you made it this far, Well Done :) 
################################################################################
################################################################################



