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

unique(df$Name) #check for unique metabolites

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
  if (is.na(name) || name == "" || str_detect(name, "Similar to")) {
    return(NA)  # Return NA for missing values or those containing "Similar to..."
  }
  
  name <- str_trim(name)  # Trim whitespace
  name <- str_replace_all(name, "[˛ĽŁ´Ąîž~łÂąâ\u0088\u0092]", "")  # Remove listed symbols
  name <- str_replace_all(name, "[\\{\\}\\[\\]?�]", "")  # Remove existing problematic symbols
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


barplot(table(df_PR$tags), xlab='conf level category', ylab='count')

metabolite_tags <- df_PR%>%
  select(species, tags)




#H089 ---#######################################################################
#H089
#Isolate dataframes to ID because we are doing timeseries analysis

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



df_89_long %>%
  select(ID, datetime, tb_mean)%>%
  distinct() %>%
  filter(ID == "H089")%>%
  group_by(datetime) %>%
  ggplot(aes(x = datetime, y = tb_mean))+
  geom_point()

Metabolites <- as.data.frame(unique(df_89_long$species))

df_89_long %>%
  dplyr::filter(species== "Pyridoxal")%>%
  dplyr::filter(ID == "H089")%>%
  group_by(datetime) %>%
  ggplot(aes(x = sample_nr, y = norm_area))+
  geom_line()+
  geom_point(aes(x=sample_nr, y=tb_mean*100000), color = "red")
  


#df_89_long$species <-iconv(df_89_long$species, "UTF-8", "ASCII", sub = "")

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
#H040

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
#rm(list = ls())
#gc()





#################################################################################################################
#Desciding to keep D compound that we can identify for the time being. and run separate analysis for compund we cannot get a hit on a DB.
#Thus we focus on analysis on infromation we most likley can follow up on,
#avoiding getting stuck with compunds we cannot say much about their biological place or relevance.

#using Metaboanalyst compound identifier tool. 
#Matching to DBs using metaboanalys tool (https://www.metaboanalyst.ca/faces/upload/MSI_Match.jsp)
#this yielded:
#Imporing DB.match files
A.DB.match <- readr::read_delim("./data/MSI-A-RPLC_feature_DB_match.txt", delim = ",")
B.DB.match <- readr::read_delim("./data/MSI-B-RPLC_feature_DB_match.txt", delim = ",")
C.DB.match <- readr::read_delim("./data/MSI-C-RPLC_feature_DB_match.txt", delim = ",")
D.DB.match <- readr::read_delim("./data/MSI-D-RPLC_feature_DB_match.txt", delim = ",")

#merge all DB.match datasets and assign tag value to each dataset according to the original
DB.match <- rbind(A.DB.match, B.DB.match, C.DB.match, D.DB.match)
DB.match$tag <- c(rep("A", nrow(A.DB.match)), rep("B", nrow(B.DB.match)), rep("C", nrow(C.DB.match)), rep("D", nrow(D.DB.match)))


#for each tag class, calculate the percentage of HMDB matches compared to the total number of unique Querys and plot in bargraph
DB.match%>%
  group_by(tag)%>%
  summarise(HMDB_matches = length(na.omit(unique(HMDB))), total_Querys = length(unique(Query )))%>%
  mutate(percentage = HMDB_matches/total_Querys*100)%>%
  ggplot(aes(x = tag, y = percentage, fill = tag))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(title = "Percentage of HMDB matches per tag class", x = "Tag class", y = "Percentage of HMDB matches")


head(DB.match)

head(RPLC.final[,1:6])


#filter RPLC.final colums to keep only the ones that have a match in HMDB or PubChem or KEGG columns of DB.match
#Filter out rows in DB.match that have NA in HMDB, PubChem, and KEGG to store in a new tibble called DB.no.match
DB.no.match <- DB.match %>%
  dplyr::filter(is.na(HMDB) & is.na(PubChem) & is.na(KEGG) & is.na(ChEBI) & is.na(METLIN))

#remove all rows in DB.match with tag D that have NA in HMDB, PubChem, and KEGG (keep unmatched Tag A B C )
DB.match.1 <- DB.match %>%
  dplyr::filter(!(tag == "D" & is.na(HMDB) & is.na(PubChem) & is.na(KEGG)))


ids <- DB.match.1$Query


rplc.colnames <- names(RPLC.final[,3:ncol(RPLC.final)])
df.RPLC.colnames <- as.data.frame(rplc.colnames)

matching.colnames <- rplc.colnames[rplc.colnames %in% ids]


#some mismatch in lengths here. Need check: 
n.mia <- length(ids) - length(matching.colnames)
n.mia

#which of the ids are missing in the matching.colnames
mia <- setdiff(ids, rplc.colnames)
mia

#mysterious... duplicatess?

#Check for duplicates in 'ids'
print(length(ids))
print(length(unique(ids)))

length(names(RPLC.final))
length(unique(names(RPLC.final)))

#Check for duplicates in 'matching.colnames'
print(length(matching.colnames))
print(length(unique(matching.colnames)))

#Check for duplicates in hilic.colnames
print(length(rplc.colnames))
print(length(unique(rplc.colnames)))

#which of the ids are duplicated
duplicated.ids <- ids[duplicated(ids)]
duplicated.ids

#filter DB.match to show only the duplicated ids
DB.match.duplicated <- DB.match %>%
  dplyr::filter(Query %in% duplicated.ids)%>%
  arrange(Query)

DB.match.duplicated

#inspection of hmdb seems that its ok. we can keep the first one and remove the rest thus keep .final as is. 
####################################
####################################
#What about the unmatched metabolites
no.ids <- DB.no.match$Query 
matching.colnames.nomatch <- rplc.colnames[rplc.colnames %in% no.ids]

#some mismatch in lengths here. Need check: 
n.mia <- length(no.ids) - length(matching.colnames.nomatch)
n.mia

#which of the ids are missing in the matching.colnames
mia <- setdiff(no.ids,matching.colnames.nomatch)
mia


#seems there is some spelling issues.the missing of the mia of the no mach compund dont have a "-" in the end of the name. fix this
#add "-" to the end of the no.ids that are in the mia list
DB.no.match$Query <- ifelse(DB.no.match$Query %in% mia, paste0(DB.no.match$Query, "-"), DB.no.match$Query)

#check it worked
no.ids <- DB.no.match$Query 
matching.colnames.nomatch <- rplc.colnames[rplc.colnames %in% no.ids]
mia <- setdiff(no.ids,matching.colnames.nomatch)
mia
#found 3 wrong spellings:
#"6-Thio-5â€²-Xanthylic Acid-"                                                  
#"Glycine, N-(3Î±,5Î²,7Î±,12Î±)-3,12-Dihydroxy-24-Oxo-7-(Sulfooxy)Cholan-24-Yl-"
#"6,6-Bis(Ethylsulfanyl)-5-Methoxy-1,2,3,4-Hexanetetro-" 

#correct Spelling:
#6-Thio-5-Xanthylic Acid
#Glycine, N-(3,5,7,12)-3,12-Dihydroxy-24-Oxo-7-(Sulfooxy)Cholan-24-Yl-	
#6,6-Bis(Ethylsulfanyl)-5-Methoxy-1,2,3,4-Hexanetetrol

DB.no.match$Query <- ifelse(DB.no.match$Query == "6-Thio-5â€²-Xanthylic Acid-", "6-Thio-5-Xanthylic Acid", DB.no.match$Query)
DB.no.match$Query <- ifelse(DB.no.match$Query == "Glycine, N-(3Î±,5Î²,7Î±,12Î±)-3,12-Dihydroxy-24-Oxo-7-(Sulfooxy)Cholan-24-Yl-", "Glycine, N-(3,5,7,12)-3,12-Dihydroxy-24-Oxo-7-(Sulfooxy)Cholan-24-Yl-", DB.no.match$Query)
DB.no.match$Query <- ifelse(DB.no.match$Query == "6,6-Bis(Ethylsulfanyl)-5-Methoxy-1,2,3,4-Hexanetetro-", "6,6-Bis(Ethylsulfanyl)-5-Methoxy-1,2,3,4-Hexanetetrol", DB.no.match$Query)

#check it worked
no.ids <- DB.no.match$Query 
matching.colnames.nomatch <- rplc.colnames[rplc.colnames %in% no.ids]
mia <- setdiff(no.ids,matching.colnames.nomatch)
mia

#good

###################################
#filter RPLC.final to keep only the columns that DO NOT have a match in any DB
RPLC.nomatch.final <- RPLC.final %>%
  dplyr::select(sample, group, all_of(no.ids))



#filter RPLC.final to keep only the columns that DO have a match in DB.match and all unmatched MSI A, B, C
RPLC.final <- RPLC.final %>%
  dplyr::select(sample, group, all_of(ids))



# View the first few rows of the filtered tibble
head(RPLC.final)
head(RPLC.nomatch.final)



#write to disk
#write.csv(RPLC.final, "./data/RPLC_pos_final_FM.csv")
#write.csv(RPLC.nomatch.final, "./data/RPLC_pos_NO-DB-MATCH_FM.csv")


#clean env to only keep RPLC.final and DB.no.match
rm(list=setdiff(ls(), c("RPLC.final", "DB.match.1")))
