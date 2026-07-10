#Clean script FM
library(dplyr)
library(arrow, warn.conflicts = T)
library(stringr)
library(readxl)
library(tidyverse)



######## Read dataframe: 
library(readr)



df<- read.csv("data/250129_CSF-hamsters_Compounds_Fresdrik-Markussen.csv",fileEncoding = "ISO-8859-2" ,header = T, check.names = F)

str(df) 

d <- df$Name #check for unique metabolites

#fill in all spaces in column headings with underscores
colnames(df) <- gsub(" ", "_", colnames(df))

#remove all . and : and "[""]" from column headings
colnames(df) <- gsub("\\.", "", colnames(df))
colnames(df) <- gsub(":", "", colnames(df))
colnames(df) <- gsub("\\[", "", colnames(df))
colnames(df) <- gsub("]", "", colnames(df))

str(df) 

df <-df %>%
mutate(Tags = ifelse(Annot_Source_mzVault_Search == "Full match", "A", Tags))

df <- df %>%
  filter(Annot_DeltaMass_ppm  > -5 & Annot_DeltaMass_ppm  < 5) # This makes sure that we are not working with bad data from the MS/MS instrument. 



colnamesList <- as.data.frame(colnames(df))



## Select name, CalcMW, RT, norm_area and peak rating, Select out group area and group CV:
df_sel <- df %>% 
  dplyr::select(species=Name, tags=Tags, CalcMW=Calc_MW , RTmin=RT_min, contains(c("Norm_Area_", "Peak")), -contains(c("RSD")))


#colnames(df_sel) <- ifelse(grepl("Area", colnames(df_sel)), 
 #                          paste0("Norm_", colnames(df_sel)), 
  #                         colnames(df_sel))

colnamesList_sel <- as.data.frame(colnames(df_sel))
colnamesList_sel


clean.names <- function(name) {
  name <- str_trim(name)  # Trim whitespace
  name <- str_replace_all(name, "˛|Ľ|Ł|´|Ą|î|ž|~|ł|Â|ą|â|\u0088|\u0092", "")  # Remove listed symbols
  name <- str_replace_all(name, "[\\{\\}\\[\\]?�]", "")  # Remove existing problematic symbols
  name <- str_replace_all(name, "Â", "")  # Remove existing problematic symbols
  name <- str_replace_all(name, "4�\u17e", "")  # Remove existing problematic symbols
  name <- str_replace_all(name, "�\u008e", "")  # Remove existing problematic symbols
  name <- str_replace_all(name, "-�", "")
  name <- str_to_title(name)  # Capitalize each word
  
  return(name)
}



#Apply the function to the species column
df_sel$species<- sapply(df_sel$species, clean.names)



#add column that have CalcMW/RT, is unique 
df_sel$CalcMWRT<-paste(df_sel$CalcMW ,df_sel$RTmin ,sep="/")


## reorder so CalcMWRT is first:
df_sel<- df_sel %>% 
  relocate(tags, CalcMWRT)


##Fill in Species based on CalcMWRT
df_sel$species[is.na(df_sel$species)] <- "" #convert NA to empty string
sum(df_sel$species == "") #how many unidentified metabolites are there

df_sel$species <- ifelse(df_sel$species == "",  df_sel$CalcMWRT, df_sel$species ) #if else statement filling in MW/RT if metabolite ID is missing

df_sel$tags[is.na(df_sel$tags)] <- "" #convert NA to empty string
sum(df_sel$tags == "") #how many tag classes

#assigning level 4 class (un-IDable) as "F"
df_sel$tags <- ifelse(df_sel$tags == "",  'F', df_sel$tags ) #if else statement filling in MW/RT if metabolite ID is missing

#check to see distribution of confidence categories are. To keep in mind for later filtering. 
barplot(table(df_sel$tags), xlab='conf level category', ylab='count')


str(df_sel)


#write.csv(colnamesList_sel,"./data/colnameslist.csv")


#renaming sample location in 96-well plate to match Animal Individual ID

# Load rename mapping
rename_normA <- read.csv("./data/250129_CSF-sample-renaming.csv", header = TRUE)
head(rename_normA)
rename_normA <- rename_normA[, 1:2]  # Select only relevant columns


rename_vector <- setNames(rename_normA$sample, rename_normA$original)

rename_vector <- rename_vector[names(rename_vector) %in% names(df_sel)]


df_sel <- df_sel %>%
  rename_with(~ ifelse(.x %in% names(rename_vector), rename_vector[.x], .x), .cols = names(df_sel))

#clean \t form column names
colnames(df_sel) <- gsub("\t", "", colnames(df_sel))

# Print updated column names
print(colnames(df_sel))







#find index range to perform peak filtering on
str_which(names(df_sel), "PR")
ncol(df_sel)


# filter out only metabolites that have at least 4 samples with > 5 in peakrating
df_PR <- df_sel%>%
  filter(rowSums(.[,162:ncol(.)] >=5, na.rm = T) >=4) #rowsums switches to counts if we put a condition behind the selection (>=5). Then we filter on the basis if that count is >=4


#inspect structure:
str(df_PR)
ncol(df_PR)

# Write file to have it. 
# write.csv(df_PR, "./data/df_PF_CSF_RPLC_peakrated_FM.csv")

#Consider filtering on tags at this point. 
#Filering to keep MSI level 1-3 (A-E). 

sum(df_PR$tags %in% c("E"))

df_PR <- df_PR %>%
  dplyr::filter(tags %in% c("A", "B", "C", "D", "E", "F"))

df_PR_A <- df_PR %>%
  dplyr::filter(tags %in% c("A"))
df_PR_B <- df_PR %>%
  dplyr::filter(tags %in% c("B"))
df_PR_C <- df_PR %>%
  dplyr::filter(tags %in% c("C"))
df_PR_D <- df_PR %>%
  dplyr::filter(tags %in% c("D"))
df_PR_E <- df_PR %>%
  dplyr::filter(tags %in% c("E"))

#write.csv(unique(df_PR_A$species), "./data/df_HibExp3_MSI-A-RPLC_list.csv", row.names = FALSE)
#write.csv(unique(df_PR_B$species), "./data/df_HibExp3_MSI-B-RPLC_list.csv", row.names = FALSE)
#write.csv(unique(df_PR_C$species), "./data/df_HibExp3_MSI-C-RPLC_list.csv", row.names = FALSE)
#write.csv(unique(df_PR_D$species), "./data/df_HibExp3_MSI-D-RPLC_list.csv", row.names = FALSE)
#write.csv(unique(df_PR_E$species), "./data/df_HibExp3_MSI-E-RPLC_list.csv", row.names = FALSE)


metabolite_tags <- df_PR%>%
  select(species, tags)

write.csv(metabolite_tags, "./data/metabolite_tags.csv", row.names = F)




#H089 ---#######################################################################
#H089
#Isolate data frames to ID because we are doing time series analysis

df_89 <- df_PR %>%
  dplyr::select( species, contains(c("89", "QC", "blank" )), -contains(c("PR", "RSD")))

colnames(df_89) 

tb89 <- read.csv("./data/H089_temps_sampling_start.csv", header = TRUE)
tb89$datetime <- as.POSIXct(tb89$datetime, format = "%d/%m/%Y %H:%M", tz= "GMT")

tb89 %>%
  ggplot(aes(x = datetime, y = Tb)) +
  geom_line() 

head(tb89)

#need to fill in sample numbers for each sample, every 4 hours so sample row 1-8 are sample 
#1 and 9-16 are sample 2 etc.
tb89$Sample_nr <- rep(seq_len(ceiling(nrow(tb89) / 8)), each = 8, length.out = nrow(tb89))


tb89 <- tb89 %>% 
  select(datetime, Tb, ID, Sample_nr)




# Convert df_89 to long format 
df_89_long <- df_89%>%
  dplyr::select(species, contains(c("89","QC","blank")), -contains(c("PR", "RSD"))) %>%
  tidyr::pivot_longer(cols = -c(species), names_to = "sample", values_to = "norm_area")

# Add sample number to df_89_long
df_89_long <- df_89_long %>%
  mutate(Sample_nr = as.numeric(str_extract(sample, "\\d{2,3}$"))) %>%
  arrange(Sample_nr)

summary(df_89_long)
head(df_89_long)

# Merge df_89_long with tb89 by Sample_nr
df_89_long <- df_89_long %>%
  full_join(tb89, by = "Sample_nr", relationship = "many-to-many")%>%
  dplyr::select(ID,sample,Sample_nr, datetime, Tb, species,  norm_area)

summary(df_89_long)
head(df_89_long)


# Rename columns 
names(df_89_long) <- c("ID","sample","sample_nr","datetime","tb","species", "norm_area")

# Need to round time stamps to nearest hour to be able to match to corresponding changes in metabolites
df_89_long <- df_89_long %>%
  group_by(sample_nr) %>% 
  mutate(tb_mean = round(mean(tb, na.rm = TRUE), 1)) %>% 
  ungroup()%>%
  select(-tb)

head(df_89_long)


df_89_long <- df_89_long %>%
  group_by(sample_nr) %>%
  mutate(datetime = rep(datetime[seq(1, n(), by = 8)], each = 8, length.out = n())) %>%
  ungroup()

# Plot it to sanity check
df_89_long %>%
  select(ID, datetime, tb_mean)%>%
  distinct() %>%
  filter(ID == "H089")%>%
  group_by(datetime) %>%
  ggplot(aes(x = datetime, y = tb_mean))+
  geom_point()


Metabolites <- as.data.frame(unique(df_89_long$species))

#Add layer of metabolite to see correlations
df_89_long %>%
  dplyr::filter(species== "Adenosine")%>%
  dplyr::filter(ID == "H089")%>%
  group_by(datetime) %>%
  ggplot(aes(x = sample_nr, y = norm_area))+
  geom_line()+
  geom_point(aes(x=sample_nr, y=tb_mean*100000), color = "red")
  


head(df_89_long)

# Cast df to wide format using pivot_wider
H089_final <- df_89_long %>%
  group_by(ID, sample, sample_nr, datetime, species) %>%
  pivot_wider(
    names_from = species,
    values_from = norm_area,
    values_fill = list(norm_area = 0),
    values_fn = list(norm_area = mean) 
  ) 



#reorder rows to have them in order of datetime
H089_final <- H089_final %>%
  arrange(datetime)

#fill in sample_nr with NAs. QC1-9: -1to -9, blank_start1:-10, blank_start2:-11, blank_end:-12
H089_final <- H089_final %>%
  mutate(sample_nr = case_when(
    grepl("^QC[1-9]$", sample) ~ as.numeric(sub("QC", "-", sample)),
    sample == "blank_start1" ~ -10,
    sample == "blank_start2" ~ -11,
    sample == "blank_end" ~ -12,
    is.na(sample_nr) ~ NA_real_, 
    TRUE ~ sample_nr  
  ))

H089_final <- H089_final %>%
  mutate(ID = case_when(
    grepl("^QC[1-9]$", sample) ~ "QC",        
    sample %in% c("blank_start1", "blank_start2", "blank_end") ~ "blank", 
    !is.na(sample) ~ as.character(sample), 
    TRUE ~ NA_character_  
  ))



#need to remove rows with missing observations in "sample" column
H089_final <-H089_final[!is.na(H089_final$sample), ]


#paste datetime NAs: fixed datetime 2023-05-28 10:00:00
H089_final <- H089_final %>%
  mutate(datetime = case_when(
    is.na(datetime) ~ as.POSIXct("2023-05-28 10:00:00", format = "%Y-%m-%d %H:%M:%S"),
    TRUE ~ datetime
  ))

#paste tb_mean NaNs: fixed tb_mean tp 15
H089_final <- H089_final %>%
  mutate(tb_mean = case_when(
    is.na(tb_mean) ~ 15,
    TRUE ~ tb_mean
  ))

#remove columns that contains NAs
H089_final <- H089_final %>%
  select(where(~ !any(is.na(.)))) 

#remove column names NA
H089_final <- H089_final %>%
  select(-"NA")



#Write to file to have it.
#write.csv(H089_final, "./data/CSF_H089_final.csv", row.names = FALSE)

#clean env
#rm(list = ls())
gc()

#H040##############################################################################################################
#H040 NO GO - No to little correlation with Tb - Bad probe placement?

#Isolate dataframes to ID because we are doing timeseries analysis

df_40 <- df_PR %>%
  dplyr::select( species, contains(c("40", "QC", "blank" )), -contains(c("PR","89")))

colnames(df_40) 

tb40 <- read.csv("./data/H040_temps_sampling_start.csv", header = TRUE)
head(tb40)


tb40$datetime <- as.POSIXct(tb40$datetime, format = "%d/%m/%Y %H:%M", tz= "GMT")

tb40 %>%
  ggplot(aes(x = datetime, y = Tb)) +
  geom_line() 



#need to fill in sample numbers for each sample, every 4 hours so sample row 1-8 are sample 
#1 and 9-16 are sample 2 etc.
tb40$Sample_nr <- rep(seq_len(ceiling(nrow(tb40) / 8)), each = 8, length.out = nrow(tb40))


tb40 <- tb40 %>% 
  select(datetime, Tb, ID, Sample_nr)

summary(tb40)
#remove all rows with sample_nr 41-76
tb40 <- tb40 %>%
  filter(Sample_nr <= 40)



# Convert df_40 to long format 
df_40_long <- df_40%>%
  dplyr::select(species, contains(c("40","QC","blank")), -contains(c("PR", "RSD"))) %>%
  tidyr::pivot_longer(cols = -c(species), names_to = "sample", values_to = "norm_area")

unique(df_40_long$sample)

df_40_long <- df_40_long %>%
  mutate(Sample_nr = case_when(
    str_detect(sample, "QC|blank") ~ NA_real_,  
    TRUE ~ as.numeric(str_extract(sample, "\\d{1,3}$"))  
  )) %>%
  arrange(Sample_nr) 

#how many of each unique sample_nr
table(df_40_long$Sample_nr)
summary(df_40_long)
head(df_40_long)



# Merge df_40_long with tb40 by Sample_nr
df_40_long <- df_40_long %>%
  full_join(tb40, by = "Sample_nr", relationship = "many-to-many")%>%
  dplyr::select(ID,sample, Sample_nr, datetime, Tb, species,  norm_area)

table(df_40_long$Sample_nr)

summary(df_40_long)
head(df_40_long)


# Rename columns 
names(df_40_long) <- c("ID","sample","sample_nr","datetime","tb","species", "norm_area")

df_40_long <- df_40_long %>%
  group_by(sample_nr) %>% 
  mutate(tb_mean = round(mean(tb, na.rm = TRUE), 1)) %>% 
  ungroup()%>%
  select(-tb)

head(df_40_long)

df_40_long <- df_40_long %>%
  group_by(sample_nr) %>%
  mutate(datetime = rep(datetime[seq(1, n(), by = 8)], each = 8, length.out = n())) %>%
  ungroup()



df_40_long %>%
  select(ID, datetime, tb_mean)%>%
  distinct() %>%
  filter(ID == "H040")%>%
  group_by(datetime) %>%
  ggplot(aes(x = datetime, y = tb_mean))+
  geom_point()

Metabolites <- as.data.frame(unique(df_40_long$species))

df_40_long %>%
  dplyr::filter(species== "Tryptophan")%>%
  dplyr::filter(ID == "H040")%>%
  group_by(datetime) %>%
  ggplot(aes(x = sample_nr, y = norm_area))+
  geom_line()+
  geom_point(aes(x=sample_nr, y=tb_mean*100), color = "red")



#df_40_long$species <-iconv(df_40_long$species, "UTF-8", "ASCII", sub = "")

head(df_40_long)

# Cast df to wide format using pivot_wider
H040_final <- df_40_long %>%
  group_by(ID, sample, sample_nr, datetime, species) %>%
  pivot_wider(
    names_from = species,
    values_from = norm_area,
    values_fill = list(norm_area = 0),
    values_fn = list(norm_area = mean) 
  ) 



#reorder rows to have them in order of datetime
H040_final <- H040_final %>%
  arrange(datetime)

#fill in sample_nr with NAs. QC1-9: -1to -9, blank_start1:-10, blank_start2:-11, blank_end:-12
H040_final <- H040_final %>%
  mutate(sample_nr = case_when(
    grepl("^QC[1-9]$", sample) ~ as.numeric(sub("QC", "-", sample)),
    sample == "blank_start1" ~ -10,
    sample == "blank_start2" ~ -11,
    sample == "blank_end" ~ -12,
    is.na(sample_nr) ~ NA_real_, 
    TRUE ~ sample_nr  
  ))

H040_final <- H040_final %>%
  mutate(ID = case_when(
    grepl("^QC[1-9]$", sample) ~ "QC",        
    sample %in% c("blank_start1", "blank_start2", "blank_end") ~ "blank", 
    !is.na(sample) ~ as.character(sample), 
    TRUE ~ NA_character_  
  ))




#need to remove rows with missing observations in "sample" column
H040_final <-H040_final[!is.na(H040_final$sample), ]


#paste datetime NAs: fixed datetime 2023-05-28 10:00:00
H040_final <- H040_final %>%
  mutate(datetime = case_when(
    is.na(datetime) ~ as.POSIXct("2023-06-08 18:15:00", format = "%Y-%m-%d %H:%M:%S"),
    TRUE ~ datetime
  ))

#paste tb_mean NaNs: fixed tb_mean tp 15
H040_final <- H040_final %>%
  mutate(tb_mean = case_when(
    is.na(tb_mean) ~ 15,
    TRUE ~ tb_mean
  ))
#remove columns that contains NAs
H040_final <- H040_final %>%
  select(where(~ !any(is.na(.)))) 

#remove column names NA
H040_final <- H040_final %>%
  select(-"NA")



#Write to file to have it.
write.csv(H040_final, "./data/CSF_H040_final.csv", row.names = FALSE)

#clean env
rm(list = ls())
gc()



