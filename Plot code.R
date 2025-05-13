

#example correlation test for each metabolite with tb: 
cor.test(H089n$tb_mean, H089n$`Dl-Glutamine`, method = "spearman")



# Calculate correlation of each metabolite with tb_mean
cor_results <- cor(H089n[,-1], use = "complete.obs", method = "spearman")

head(H089n[,1:6])

library(purrr)
cor_tb <- H089n %>%
  select(-sample_nr) %>%
  pivot_longer(cols = -tb_mean, names_to = "Metabolite", values_to = "value") %>%
  group_by(Metabolite) %>%
  summarise(
    test = list(cor.test(tb_mean, value, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    correlation = map_dbl(test, ~ .x$estimate %||% NA_real_),
    p_value = map_dbl(test, ~ .x$p.value %||% NA_real_)
  ) %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  arrange(desc(abs(correlation))) %>%
  select(Metabolite, correlation, p_value)


head(cor_tb)

# Sort by correlation strength
cor_tb <- cor_tb %>%
  arrange(desc(abs(correlation)))


# View top 10 strongest correlations
head(cor_tb, 10)  

#Fancy ggpot of mets over .5 in r correlation with tb
cor_tb%>%
  filter(correlation< -0.5 | correlation >0.5)%>%
  ggplot( aes(x = reorder(Metabolite, correlation), y = correlation, fill = correlation)) +
  geom_col(color="black", size=.05) +
  coord_flip() +
  scale_fill_gradient2(low = "forestgreen", mid = "white", high = "purple", midpoint = 0) +
  scale_color_gradient2(low = "forestgreen", mid = "white", high = "purple", midpoint = 0) +
  theme_classic() +
  labs(
    x = "Metabolite",
    y = "Correlation with Tb",
    title = "Correlation of Metabolites with Body Temperature"
  )


#Filter H089n columns to contain only mets with p<0.05
mets <- cor_tb %>%
  pull(Metabolite)

H089n<- H089n %>%
  select(sample_nr, tb_mean, all_of(mets))

cor_tb %>%
  filter(abs(correlation) > 0.5)%>%
  pull(Metabolite)%>%
  length()


###############################################################################




# Plot faceted correlation for all metabolite with correlation over 0.5
H089n_long <- H089n %>%
  tidyr::pivot_longer(cols = -c(sample_nr, tb_mean), 
                      names_to = "Metabolite", 
                      values_to = "Abundance")

H089n_long%>%
  filter(Metabolite %in% c("Dl-Glutamine"))%>%
  ggplot( aes(x = tb_mean, y = Abundance)) +
  geom_point(size = 2, color = "#2A9D8F") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  #facet_wrap(~ Metabolite, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = "Body Temperature (tb_mean)",
    y = "Metabolite Abundance",
    title = "Correlation between Body Temp and Metabolites"
  )



#Lets group metabolites base on lag
#Lag < 0  metabolite precedes Tb change (potential signal/metabolic trigger)
#Lag = 0  synchronous
#Lag > 0  metabolite follows Tb change (potential response/metabolic consequence)

ccf_results_filtered <- ccf_results_filtered  %>%
  mutate(Timing = case_when(
    Lag < 0 ~ "Precedes Tb",
    Lag == 0 ~ "Synchronous",
    Lag > 0 ~ "Follows Tb"
  ))

head(ccf_results_filtered)

(ggplot(ccf_results_filtered, aes(x = Lag * 2, y = Max_Correlation, color = Timing)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    #scale_x_continuous(limits= c(-6,6),breaks = seq(-6,6,2))+
    # Add metabolite names only for points where abs(Lag) > 1
    theme_bw() +
    labs(
      title = "Correlation vs Lag for Metabolites",
      x = "Lag in Hours",
      y = "Spearman r"
    ) +
    scale_color_manual(values = c(
      "Precedes Tb" = "purple",
      "Synchronous" = "gray20",
      "Follows Tb" = "forestgreen"
    )))

pretb <- ccf_results_filtered%>%
  filter(Lag<0)%>%
  pull(Metabolite)

posttb <-ccf_results_filtered%>%
  filter(Lag>0)%>%
  pull(Metabolite)

synchtb <- ccf_results_filtered%>%
  filter(Lag==0)%>%
  pull(Metabolite)



ccf_all_df <- map_dfr(
  mets,
  ~{
    met_name <- .x
    ccf_out <- ccf(H089n[[met_name]], H089n$tb_mean, lag.max = 30, plot = FALSE)
    tibble(
      Metabolite = met_name,
      lag = ccf_out$lag * 2,  # in hours
      correlation = as.vector(ccf_out$acf)  # flatten array
    )
  }
)

ccf_all_df

# Add CI threshold (same across all)
ci_threshold <- 1.96 / sqrt(nrow(H089n))

ccf_all_df %>%
  filter(Metabolite %in% pretb)%>%
  ggplot(aes(x = lag, y = correlation)) +
  geom_col(fill = "steelblue4", width = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_hline(yintercept = c(ci_threshold, -ci_threshold), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_bw(base_size = 10)+
  scale_x_continuous(breaks = seq(-16, 16, 2), limits = c(-16, 16)) +
  labs(
    title = "Cross-Correlation: Preceeding & Tb",
    x = "Lag in Hours",
    y = "Correlation (r)"
  ) +
  annotate(
    "text",
    x = 6, y = 0.5,
    label = "Lag = 0: Tb & metabolite in synch\nLag < 0: Metabolite precedes Tb\nLag > 0: Metabolite follows Tb",
    hjust = 0, vjust = 1,
    size = 1.5, color = "gray20"
  ) +
  annotate(
    "text",
    x = max(ccf_df$lag), y = 0.23,
    label = "95% CI",
    hjust = 1, vjust = 1,
    size = 1, color = "gray20"
  ) +
  theme_bw(base_size = 6)



ccf_all_df %>%
  filter(Metabolite %in% posttb)%>%
  ggplot(aes(x = lag, y = correlation)) +
  geom_col(fill = "forestgreen", width = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_hline(yintercept = c(ci_threshold, -ci_threshold), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_bw(base_size = 10)+
  scale_x_continuous(breaks = seq(-16, 16, 2), limits = c(-16, 16)) +
  labs(
    title = "Cross-Correlation: Following Tb",
    x = "Lag in Hours",
    y = "Correlation (r)"
  ) +
  annotate(
    "text",
    x = 6, y = 0.5,
    label = "Lag = 0: Tb & metabolite in synch\nLag < 0: Metabolite precedes Tb\nLag > 0: Metabolite follows Tb",
    hjust = 0, vjust = 1,
    size = 1.5, color = "gray20"
  ) +
  annotate(
    "text",
    x = max(ccf_df$lag), y = 0.23,
    label = "95% CI",
    hjust = 1, vjust = 1,
    size = 1, color = "gray20"
  ) +
  theme_bw(base_size = 6)

metfollowstb <- H089n %>%
  melt(id.vars = c("sample_nr", "tb_mean")) %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(as.numeric(sample_nr) >= 0) %>%
  filter(variable %in% posttb)

#write.csv(metfollowstb, file= "./data/shiny/metabolitesfollowstb.csv", row.names = FALSE)


metpretb <- H089n %>%
  melt(id.vars = c("sample_nr", "tb_mean")) %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(as.numeric(sample_nr) >= 0) %>%
  filter(variable %in% pretb)

#write.csv(metpretb, file= "./data/shiny/metabolitesPREtb.csv", row.names = FALSE)


metsynchtb <- H089n %>%
  melt(id.vars = c("sample_nr", "tb_mean")) %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(as.numeric(sample_nr) >= 0) %>%
  filter(variable %in% synchtb)


#write.csv(metsynchtb, file= "./data/shiny/metabolitesSYNCHtb.csv", row.names = FALSE)






ccf_all_df %>%
  filter(Metabolite %in% synchtb)%>%
  ggplot(aes(x = lag, y = correlation)) +
  geom_col(fill = "gray20", width = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_hline(yintercept = c(ci_threshold, -ci_threshold), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_bw(base_size = 10)+
  scale_x_continuous(breaks = seq(-16, 16, 2), limits = c(-16, 16)) +
  labs(
    title = "Cross-Correlation: Following Tb",
    x = "Lag in Hours",
    y = "Correlation (r)"
  ) +
  annotate(
    "text",
    x = 6, y = 0.5,
    label = "Lag = 0: Tb & metabolite in synch\nLag < 0: Metabolite precedes Tb\nLag > 0: Metabolite follows Tb",
    hjust = 0, vjust = 1,
    size = 1.5, color = "gray20"
  ) +
  annotate(
    "text",
    x = max(ccf_df$lag), y = 0.23,
    label = "95% CI",
    hjust = 1, vjust = 1,
    size = 1, color = "gray20"
  ) +
  theme_bw(base_size = 6)



#See individual metabolite dynamicss for more. 
#https://60afyc-fredrik0markussen.shinyapps.io/MetaboliteDynamicsvsBodyTemperature/


#we need to define the start and end of each cycle. 

plotly::ggplotly(ggplot(H089n, aes( x = sample_nr, y = tb_mean))+
                   geom_point())

#cycle starts and ends 3 samples before the first drop in tb

# Define sample ranges
cycle1_range <- 1:41
cycle2_range <- 42:94

#so then we can split the data by cycle and ad relative phase
# Cycle 1: sample 1 to 41 (41 total) (torpor end nr 27)
# Cycle 2: sample 42 to 94 (53 total) (torpor end nr 81)

cycle1_df <- H089n %>%
  filter(sample_nr %in% cycle1_range) %>%
  mutate(
    rel_phase_deg =NA_real_,
    Cycle = 1
  )%>%
  relocate(rel_phase_deg, Cycle, .after = sample_nr)

cycle2_df <- H089n %>%
  filter(sample_nr %in% cycle2_range) %>%
  mutate(
    rel_phase_deg =NA_real_,
    Cycle = 2
  )%>%
  relocate(rel_phase_deg, Cycle, .after = sample_nr)


#cycle 1 phase degree step per sample
deg_step_c1 <- 360/(41-1) 
#cycle 2 phase degree step per sample
deg_step_c2 <- 360/(94-42) 

#Fill in rel_phase_deg from deg_step_c1 and deg_step_c2 in respectve dfs

cycle1_df <- cycle1_df %>%
  mutate(rel_phase_deg = seq(0, by = deg_step_c1, length.out = n()))


cycle2_df <- cycle2_df %>%
  mutate(rel_phase_deg = seq(0, by = deg_step_c2, length.out = n()))


#plot for sanity check
ggplot() +
  geom_point(data = cycle1_df, aes(x = rel_phase_deg, y = tb_mean), color = "black") +
  geom_point(data = cycle2_df, aes(x = rel_phase_deg, y = tb_mean), color = "firebrick") +
  labs(title = "Tb_mean vs. rel_phase_deg",
       x = "Relative Phase (°)", y = "Tb_mean") +
  theme_minimal()


################################################################################
#cycle 1

#Set metabolite
Met<- c("Dodecanedioic Acid")


ccf_result <- ccf_permutation_ci(cycle1_df[[Met]], cycle1_df$tb_mean, lag.max = 30, level = 0.99)


ggplot(ccf_result, aes(x = lag *2, y = correlation)) +
  geom_vline(xintercept = 0, size= 0.5,  linetype = "dashed")+
  geom_hline(yintercept = 0, color = "black") +
  geom_line(aes(y = ci_upper), size = 0.5, color= 'firebrick',  linetype = "dashed") + 
  geom_line(aes(y = ci_lower), size = 0.5,color= 'firebrick',  linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "gray80", alpha = 0.4) +
  geom_col(fill = "black", width = 0.5) +
  labs(
    title = paste("Cross-correlation with", Met),
    x = "Lag in hours)",
    y = "Spearman correlation"
  ) +
  annotate(
    "text",
    x = 10, y = 0.9,
    label = "Lag = 0: Tb & metabolite change in synch\nLag < 0: Metabolite precedes Tb\nLag > 0: Metabolite follows Tb",
    hjust = 0, vjust = 1,
    size = 4, color = "gray20"
  )+
  theme_bw()



extract_max_ccf <- function(met_name, df, lag.max = 15, level = 0.99, n_perm = 1000) {
  ci_df <- ccf_permutation_ci(df[[met_name]], df$tb_mean, lag.max = lag.max, level = level, n_perm = n_perm)
  
  # Compute top 5 correlations by absolute value
  top5_idx <- order(abs(ci_df$correlation), decreasing = TRUE)[1:5]
  top5_corrs <- ci_df$correlation[top5_idx]
  
  # Keep only those with same sign as the max *positive* correlation
  pos_corrs <- ci_df$correlation
  pos_idx <- which.max(pos_corrs)  # max *positive* correlation
  
  # Handle case where no positive correlations exist
  if (pos_corrs[pos_idx] <= 0) {
    max_corr <- NA_real_
    lag_at_max <- NA_integer_
    mean_top5 <- NA_real_
    min_top5 <- NA_real_
    significant <- FALSE
  } else {
    max_corr <- pos_corrs[pos_idx]
    lag_at_max <- ci_df$lag[pos_idx]
    
    # Recalculate top 5 correlations using only positive values
    top5_pos <- sort(pos_corrs[pos_corrs > 0], decreasing = TRUE)[1:min(5, sum(pos_corrs > 0))]
    
    mean_top5 <- mean(top5_pos)
    min_top5 <- top5_pos[which.min(abs(top5_pos))]
    
    # CI bounds (same across lags)
    upper_bound <- quantile(ci_df$ci_upper, probs = 0.5)
    lower_bound <- quantile(ci_df$ci_lower, probs = 0.5)
    
    significant <- (min_top5 > upper_bound)
  }
  
  tibble(
    Metabolite = met_name,
    Max_Correlation = max_corr,
    mean_top5 = mean_top5,
    min_top5 = min_top5,
    Lag = lag_at_max,
    ci_lower = quantile(ci_df$ci_lower, probs = 0.5),
    ci_upper = quantile(ci_df$ci_upper, probs = 0.5),
    significant = significant
  )
}


# Loop over all metabolites using map_dfr 

################## !!! WARNING!!! HEAVY ON COMPUTATION!!! ##################### Cycle 1
ccf_results <- map_dfr(
  colnames(cycle1_df)[5:ncol(cycle1_df)],
  extract_max_ccf,
  df = cycle1_df,
  lag.max = 21,
  level = 0.99,
  n_perm = 1000  # Reduce for speed if needed, but needs to mach CI level!
) %>%
  arrange(desc(abs(Max_Correlation)))

summary(ccf_results)

ccf_results%>%filter(significant)%>%nrow()


ccf_results_filtered <- ccf_results%>% filter(significant)
ccf_results_filtered <- ccf_results_filtered %>%
  mutate(Lag_abs_hours = abs(Lag * 2))


ccf_results_filtered <- ccf_results_filtered  %>%
  mutate(Timing = case_when(
    Lag < 0 ~ "Precedes Tb",
    Lag == 0 ~ "Synchronous",
    Lag > 0 ~ "Follows Tb"
  ))




(ggplot(ccf_results_filtered, aes(x = Lag * 2, y = Max_Correlation, color = Timing)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    #scale_x_continuous(limits= c(-6,6),breaks = seq(-6,6,2))+
    # Add metabolite names only for points where abs(Lag) > 1
    theme_bw() +
    labs(
      title = "Correlation vs Lag for Metabolites",
      x = "Lag in Hours",
      y = "Spearman r"
    ) +
    scale_color_manual(values = c(
      "Precedes Tb" = "purple",
      "Synchronous" = "gray20",
      "Follows Tb" = "forestgreen"
    )))


cycle_period <- 41  # samples, adjust   o your TA-cycle estimate

ccf_results_filtered <- ccf_results_filtered %>%
  mutate(
    Phase_frac = Lag / cycle_period,  # keep the fractional value
    Phase_deg = (Phase_frac * 360) %% 360,
    Phase_rad = (Phase_frac * 2 * pi) %% (2 * pi),
    
    Timing = case_when(
      Phase_deg >= 345 | Phase_deg <= 15 ~ "Synchronous",
      Phase_deg < 180 ~ "Precedes Tb",
      TRUE ~ "Follows Tb"
    )
  )




ggplot(ccf_results_filtered %>% filter(significant), 
       aes(x = Phase_rad, y = abs(Max_Correlation), color = Timing)) +
  geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 2*pi, by = pi/2),
                     labels = c("0°", "90°", "180°", "270°", "360°")) +
  coord_polar(theta = "x", start = 0) +
  scale_color_manual(values = c("Precedes Tb" = "skyblue", 
                                "Synchronous" = "forestgreen", 
                                "Follows Tb" = "firebrick")) +
  theme_minimal(base_size = 12) +
  labs(title = "Phase Distribution of Significant Metabolites",
       x = NULL, y = "Abs Max Correlation")

head(ccf_results_filtered)








################################################################################



