#load libraries
library(dplyr)
library(gridExtra)
library(reshape2)
library(ca)
library(vegan)
library(circlize)
library(RColorBrewer)
library(ggalt)
library(ComplexHeatmap) 
library(stringr)
library(PCAtools)
library(tidyr)


#Import data
df.t <- read.csv( "./data/CSF_H040_final.csv", header = TRUE, check.names = FALSE)

H040<- df.t


names(df.t)  
ncol(df.t)

#total metabolites
metabolites <- as.data.frame(colnames(H040)[6:ncol(H040)] , check.names = F)

colnames(metabolites) <- c("Feature")

#how many just mz/rt peaks?
metabolites_mzpeaks <- metabolites %>%
  filter(grepl("^\\d+\\.\\d+/\\d+\\.\\d+$", Feature))

#how many with name?
metabolites <- metabolites %>%
  filter(!grepl("^\\d+\\.\\d+/\\d+\\.\\d+$", Feature))

#save as csv
#write.csv(metabolites, file = "./data/metabolites_H040.csv", row.names = FALSE)





#Normalization process:
#normalization against the QC samples using the raw values. 
#aims to correct for batch effects, instrumental drift, systematic variations -> systemic variations are minimized

#transformation to QC-normalized data.
#transformation stabilizes the variance across the range of data, making the data more suitable for subsequent analysis.
#It also helps to make the data more normally distributed, which is an assumption of many statistical techniques.

#Here I choose sqrt() transform because it is the same transformation done when creating vulcano plots.
#(Vulcano LogFC calculation needs the data to be centered around 1 becasue negative values are impossible in log calcs.
#Thus, here, I am keeping the sqrt transfom the same as the vulcano plot to keep the data consistent, but I chose to center the data
#around 0, as it is a common practice in metabolomics to center the data around 0.

#to start normalize against QC samples first makes sense beacuseit directly addresses the technical variability introduced during the sample processing and acquisition.
#thus I have ensured that any subsequent transformations or scaling are applied to data that have already been adjusted for technical artifacts,
#which increases reliabillity and interpretabillity of the results. 

#Scale and center the data:
#Helps issues where some variables dominate solely due to their scale.
#Centering (subtracting the mean) can be particularly useful for techniques like PCA, where the variation around the mean is of interest. 


#Extract the QC samples
qc <- H040 %>%
  filter(str_detect(sample, "QC")) %>%  
  select(6:ncol(.)) 

met.sc <- H040
# Apply sqrt transformation centering and scaling
met.sc <- scale(sqrt(met.sc[,6:ncol(met.sc)]), center = T, scale = T)





#view(met.sc)

#Making 4 plots showing effects of normalization
m.sc.plot <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, variable))+
  geom_boxplot(fill='gray20', color= 'gray50') +
  labs(y='Compound',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )


m.sc.dens <- cbind(df.t[2], met.sc)%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density(linewidth=1, color ='grey20')+
  labs(title = 'After Normalization', 
       y='Density',
       x='Normalized Peak Intensity') + 
  theme_classic()+
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )


df.t.unnorm.box <-
  cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value, variable))+
  geom_boxplot(fill='gray20', color= 'gray50') +
  labs(y='Compound',
       x='Peak Intensity') + 
  theme_classic()+
  theme(
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )

df.t.unnorm.dens <-cbind(df.t[2], df.t[4:ncol(df.t)])%>%
  melt(., id= "sample")%>%
  mutate(value =as.numeric(unlist(value)))%>%
  ggplot(aes(value))+
  geom_density(linewidth=1, color ='grey20')+
  labs(title = 'Before Normalization', 
       y='Density',
       x=' Peak Intensity') + 
  theme_classic()+
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    #plot.margin = margin(1, 1, 1, 1, "cm")
  )



#Drawing the 4 plots together
grid.arrange(df.t.unnorm.dens,m.sc.dens, df.t.unnorm.box, m.sc.plot, ncol=2)
###############################################################################
#Initial exploration using heatmap
###############################################################################
rownames(met.sc) <- H040[,3] #Making row names the samples_nr



h.df<-as.matrix(met.sc)




##Heatmap preamble

c <- rev(brewer.pal(5, "PRGn"))  # Reverse order of Brewer's RdBu palette
col_fun = colorRamp2(c(-1,-1,0,1,1.5),c) 


#any NAs in data h.df?
sum(is.na(h.df))

#what are they?
which(is.na(h.df), arr.ind = TRUE)

colnames(h.df)

length(rownames(h.df))

custom_order <- setdiff(seq(-12, 40, by = 1), 0)
custom_order <- as.character(custom_order) 

rownames(h.df) <- as.character(rownames(h.df))  



# Define color mapping for temperature (blue = cold, red = hot)
temp_col_fun <- colorRamp2(c(5, 20, 35), c("steelblue4", "white", "firebrick"))

# Create top annotation
top_anno <- HeatmapAnnotation(
  Temperature = anno_simple(
    H040$tb_mean,  # Use tb_mean as the annotation value
    col = temp_col_fun,  # Apply the color function
    border = TRUE  # Optional: Add border
  ),
  annotation_name_gp = gpar(fontsize = 10)  # Customize text size
)

# Create heatmap with the top annotation
heat.ordered <- Heatmap(
  t(h.df),
  show_row_names = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_method_rows = "complete",
  column_dend_height = unit(2, "cm"), 
  row_dend_width = unit(2, "cm"),
  row_km = 5,
  col = col_fun,
  column_order = custom_order,
  top_annotation = top_anno,  
  heatmap_legend_param = list(
    title = "Relative abundance",
    legend_height = unit(3, "cm"),
    title_position = "leftcenter-rot"
  )
)

# Print heatmap
ht <-draw(heat.ordered)




#Not very clear patterne here. Probeplavepent suspisious... 

###############################################################################
#correlation analysis to make row order for heat map 
###############################################################################

#Run this if keeping all metabolites
H040n <- met.sc


#extract the meabolites of interest (only named ones) Skip if keeping all 1312 metabolites
#H040n <- met.sc[, !grepl("^\\d+\\.\\d+/\\d+\\.\\d+$", colnames(met.sc))]




#colnames(H040n)

head(H040n[,1:8])
#bind to get tb_mean
H040n <- cbind(tb_mean = H040$tb_mean, H040n)
H040n <- cbind(sample_nr = H040$sample_nr, H040n)

H040n <- as_tibble(H040n)

#Remove qc vars
H040n <- H040n%>%filter(sample_nr >0)



#do not detrend data for now


################################################################################
################################################################################
#Run correlation lag analyis
#in this order: Metabolite, Tb  positive lag indicates that change in Metabolite after change in body temperature
#lag = 0	Metabolite and Tb change together (synchronous)
#lag < 0	metabolite changes before Tb
#lag > 0	metabolite changes after Tb

#we can also just look at positive correlation at a larger lag span to indicate phase relationship if we use appropriate cylce. 

#to correlate this is function we use, example:
ccf_result <- ccf(H040n$`Dodecanedioic Acid`, 
                  H040n$tb_mean, 
                  lag.max = 23, 
                  plot = TRUE)

#the estimate (fishers?) significance threshold is:
1.96 / sqrt(length(H040$tb_mean))

#need to create a more formal calculation of the 95% CI since 1.96 / sqrt(length(H040$tb_mean) is only an approximation 
#need to be done for every lag calculation. We can permute Tb to thest correlation angainst thid null distributiion at each lagged position. 
#probably computer heavy
library(stats)
library(purrr)
library(dplyr)
library(tibble)

#leans heavly on "Covariate-Adjusted Spearman’s Rank Correlation with Probability-Scale Residuals" doi: 10.1111/biom.12812
#and "A robust Spearman correlation coefficient permutation test." doi:10.1080/03610926.2022.2121144.


ccf_permutation_ci <- function(x, y, lag.max = 23, level = 0.99, n_perm = 1000, seed = 42) {
  set.seed(seed)
  n <- length(x)
  lags <- -lag.max:lag.max
  
  # Function to compute Spearman correlation at a given lag
  spearman_ccf <- function(x, y, lag) {
    if (lag < 0) {
      cor(x[1:(n + lag)], y[(1 - lag):n], method = "spearman", use = "complete.obs")
    } else if (lag > 0) {
      cor(x[(1 + lag):n], y[1:(n - lag)], method = "spearman", use = "complete.obs")
    } else {
      cor(x, y, method = "spearman", use = "complete.obs")
    }
  }
  
  # Observed Spearman correlations
  obs_ccf <- sapply(lags, function(lag) spearman_ccf(x, y, lag))
  
  # Null distribution via permutation
  perm_ccf <- replicate(n_perm, {
    y_perm <- sample(y)
    sapply(lags, function(lag) spearman_ccf(x, y_perm, lag))
  })
  
  # Confidence interval bounds from permutations
  alpha <- 1 - level
  ci_lower <- apply(perm_ccf, 1, quantile, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- apply(perm_ccf, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  
  # Output
  tibble(
    lag = lags,
    correlation = obs_ccf,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant = (obs_ccf < ci_lower) | (obs_ccf > ci_upper)
  )
}



#Set metabolite
Met<- c("Creatine")

#Run for metabolite of interest:
ccf_result <- ccf_permutation_ci(H040n[[Met]], H040n$tb_mean, lag.max = 30, level = 0.99)


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
  theme_bw(base_size = 16)


#Looks good, So now we create a function to go over all 1312 metabolites
head(H040n)


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


# Loop over all metabolites using

################## !!! WARNING!!! HEAVY ON COMPUTATION!!! #####################
#check number of cores available:
parallel::detectCores()

#We can use 12 cores for this but lets limit to 8
library(furrr)
library(future)

#set up cores used
plan(multisession, workers = 10)

ccf_results <- future_map_dfr(
  colnames(H040n)[3:ncol(H040n)],
  extract_max_ccf,
  df = H040n,
  lag.max = 28,
  level = 0.99,
  n_perm = 1000, #reduce for speed if needed but higher shuffeling instanses is needed for proper estimation of null distribution of spearmanns rho. 
  .progress = TRUE
) %>%
  arrange(desc(abs(Max_Correlation)))

summary(ccf_results)



#How many significant?
ccf_results%>%filter(significant)%>%nrow()

#Subset the data
ccf_results_filtered <- ccf_results%>% filter(significant)



#Phase estimation, roguh estimate over the two cycles
#becomes a sort of "psudophase" over the two cycles since the period of the two phases ar asymmetric (unequal) 
#But is acctualy a good estimate for the second cycle. 


cycle_period <- 104  #  period in hours, we have to take Lag *2 because 1 lag is 2 hours  

ccf_results_filtered <- ccf_results_filtered %>%
  dplyr::mutate(
    Phase_deg = (Lag * 2) %% cycle_period / cycle_period * 360,    
    Phase_rad = (Lag * 2) %% cycle_period / cycle_period * 2 * pi   
  )

colnames(ccf_results_filtered)


#Assign before and after
ccf_results_filtered <- ccf_results_filtered %>%
  mutate(
    Timing = case_when(
      abs(Phase_deg - 0) <= 15 | abs(Phase_deg - 360) <= 15 ~ "Synchronous",
      Phase_deg < 180 ~ "After Tb",
      TRUE ~ "Precedes Tb"
    )
  )


#linear plot
(ggplot(ccf_results_filtered, 
        aes(x = Lag*2, y = abs(Max_Correlation))) +
    geom_point(size = 2) +
    #scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 2*pi, by = pi/2),
    #labels = c("0°", "90°", "180°", "270°", "360°")) +
    scale_x_continuous( breaks = seq(-60, 60, by = 5))+
    #labels = c("0°", "90°", "180°", "270°", "360°"))
    # coord_polar(theta = "x", start = 0) +
    theme_bw(base_size = 12) +
    labs(title = "Phase Distribution of Significant Metabolites",
         x = NULL, y = "Max Correlation"))

head(ccf_results_filtered)
#scale_color_manual(values = c("Precedes Tb" = "skyblue", 
#                             "Synchronous" = "forestgreen", 
#                            "After Tb" = "firebrick")) +


#one lag is 2 hours, which in turn corresponds to 6.9 degrees 
#so we define antiphase as 180 +- 6.9

ggplot(ccf_results_filtered, 
       aes(x = Phase_deg, y = abs(Max_Correlation))) +
  geom_point(size = 2) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = seq(0, 360, by = 90),
    labels = c("0°", "90°", "180°", "270°", "360°")
  ) +
  coord_polar( start = 0) +
  geom_vline(xintercept = 173)+
  geom_vline(xintercept = 187)+
  geom_vline(xintercept = 353)+
  geom_vline(xintercept = 7)+
  theme_bw(base_size = 12) +
  labs(
    title = "Phase Distribution of Significant Metabolites",
    x = NULL,
    y = "Max Correlation"
  )  




head(ccf_results_filtered)
#scale_color_manual(values = c("Precedes Tb" = "skyblue", 
#                             "Synchronous" = "forestgreen", 
#                            "After Tb" = "firebrick")) +


#####################################################################################
#heatmap of significant Metabolites
#####################################################################################
#Cluster analysis with heat map
rownames(met.sc) <- H040[,3]

#filter colnames to cor_tb significant metabolites
h.df <- met.sc %>%
  as.data.frame() %>%
  select(any_of(c(ccf_results_filtered$Metabolite))) %>%
  .[as.numeric(rownames(met.sc)) > 0, ] %>%  # Filter rows by rownames > 0
  as.matrix()

nrow(h.df)
nrow(H040n)

rownames(h.df) <-H040n %>% pull(1)#Making row names the samples_nr

h.df<-as.matrix(h.df)

rownames(h.df)


##Heatmap preamble

c <- rev(brewer.pal(5, "PRGn"))  # Reverse order of Brewer's RdBu palette
col_fun = colorRamp2(c(-1,-1,0,1,1.5),c) 


#any NAs in data h.df?
sum(is.na(h.df))

#what are they?
which(is.na(h.df), arr.ind = TRUE)

colnames(h.df)
rownames(h.df)


# set the column order:
custom_order <- setdiff(seq(1, 40, by = 1), 0)
custom_order <- as.character(custom_order) 
rownames(h.df) <- as.character(rownames(h.df))  


# Set the row order by lag
ccf_results_filtered$Lag


met_order_by_phase <- ccf_results_filtered %>%
  filter(Metabolite %in% colnames(h.df)) %>%
  arrange(Phase_deg, desc(Max_Correlation)) %>%
  pull(Metabolite)

# Reorder h.df columns (i.e., metabolite order in heatmap rows after transpose)
h.df <- h.df[, met_order_by_phase]



# Define color mapping for temperature (blue = cold, red = hot)
temp_col_fun <- colorRamp2(c(5, 20, 35), c("steelblue4", "white", "firebrick"))

# Create top annotation
top_anno <- HeatmapAnnotation(
  Temperature = anno_simple(
    H040n$tb_mean,  # Use tb_mean as the annotation value
    col = temp_col_fun,  # Apply the color function
    border = TRUE  # Optional: Add border
  ),
  annotation_name_gp = gpar(fontsize = 10)  # Customize text size
)


#create row annotation
row_phase <- ccf_results_filtered %>%
  filter(Metabolite %in% colnames(h.df)) %>%
  arrange(Phase_deg, desc(Max_Correlation)) %>%
  pull(Phase_deg)

row_anno <- rowAnnotation(
  Phase = anno_text(sprintf("%.1f°", row_phase), gp = gpar(fontsize = 7)),
  annotation_name_gp = gpar(fontsize = 5)
)


# Create Lag vector matching ordered metabolites
row_lags <- ccf_results_filtered %>%
  filter(Metabolite %in% colnames(h.df)) %>%
  arrange(Lag, desc(abs(Max_Correlation))) %>%
  pull(Lag)

# Factor to retain lag order
lag_factor <- factor(row_lags, levels = sort(unique(row_lags)))



# Create heatmap with the top annotation orderd by lag and time
heat.ordered <- Heatmap(
  t(h.df),
  row_names_side = "left",
  show_column_names =T,
  show_row_names = F,
  cluster_rows = FALSE,
  col = col_fun,
  column_order = custom_order,
  right_annotation = row_anno,
  top_annotation = top_anno,
  row_names_gp = gpar(fontsize = 6),       
  column_names_gp = gpar(fontsize = 6),  
  
  heatmap_legend_param = list(
    title = "Relative abundance",
    legend_height = unit(3, "cm"),
    title_position = "leftcenter-rot"
  )
)

ht <-draw(heat.ordered)

#save as svg
#svg("./figures/H040_TorporCSF_heatmap_Tb_LAGG_by_phase_correlated_all_metabolites_Spearmanns.svg", width = 8, height = 12)
#draw(heat.ordered)
#dev.off()

################################################################################


################################################################################
#Individual metabolite plots - Adjust as needed
################################################################################
library(reshape2)
library(dplyr)
library(ggplot2)

ccf_results_filtered$Metabolite

Met <-c("180.07631/7.976"  ,                                                                               
         "N,N-Dimethylaniline"  ,                                                                           
         "236.11376/3.671" ,                                                                                
         "Dethiobiotin" ,                                                                                   
        "239.13695/2.807" ,                                                                                
         "Leu-Gln" ,                                                                                        
         "266.13666/3.421" ,                                                                                
         "Methdilazine")


# Melt and filter
t <- cbind(df.t[3], met.sc) %>%
  melt(id = "sample_nr") %>%
  mutate(value = as.numeric(unlist(value))) %>%
  filter(variable %in% Met) %>%
  filter(sample_nr > 0)

# Plot all in one with facet
plot <- ggplot(t, aes(x = as.numeric(sample_nr) * 2, y = value)) +
  geom_ribbon(
    data = H040n, inherit.aes = FALSE,
    aes(
      x = sample_nr * 2,
      ymin = -2,
      ymax = (tb_mean / 9) - 2.5 + 0.1
    ),
    fill = "gray50",
    alpha = 0.4
  ) +
  geom_point(color = "gray15", size = 1) +
  geom_line(na.rm = TRUE, color = "gray15", linewidth = 0.5) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y") +
  scale_y_continuous( limits = c(-2,3),
                      name = "Cluster Mean Abundance",
                      sec.axis = sec_axis(~ . * 9 + 22, name = "Body Temperature (°C)")
  ) +
  scale_x_continuous(
    name = "Time (hours)",
    breaks = seq(0, 220, 20)
  ) +
  labs(title = NULL) +
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45))

plot

#filename <- paste0("./figures/", gsub(" ", "_", as.character(Met)), "_Tb_trace.svg")

#ggsave(filename, plot, width = 10, height =8)



################################################################################
#export data.
################################################################################
head(H040n[,1:6])
ccf_results_filtered

amplitude_df <- H040n %>%
  select(-sample_nr, -tb_mean) %>%
  summarise(across(everything(), ~ (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)) / 2)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Amplitude")

ccf_results_filtered <- ccf_results_filtered %>%
  left_join(amplitude_df, by = "Metabolite")


#from wrangling find metabolite tags and add tags to ccf_results_filtered
head(metabolite_tags)
head(ccf_results_filtered)


metabolite_tags_1 <- metabolite_tags%>%
  dplyr::filter(species %in% ccf_results_filtered$Metabolite)%>%distinct()%>%
  mutate(tag_rank = match(tags, LETTERS)) %>%  # A = 1, F = 6,
  group_by(species) %>%
  slice_min(order_by = tag_rank, n = 1) %>%    # Keep the row with highest priority tag
  ungroup() %>%
  select(-tag_rank)%>%
  rename(Metabolite=species, tag=tags)



ccf_results_filtered_tag<-full_join(ccf_results_filtered, metabolite_tags_1)%>%relocate(tag)

ccf_results_tag<-full_join(ccf_results, metabolite_tags_1)%>%relocate(tag)

ccf_results_tag%>% filter(tag %in% c("F"))%>%
  nrow()




#export:
#write.csv(ccf_results_filtered_tag, "./data/ccf_results_filtered_all_tag.csv")


#Correlation with tb not entierly clear on this one. Probe placement perhaps not ideal. Exclude?  





