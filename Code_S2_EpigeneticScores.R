################################################### 
#####  Code S2 DNAm epigenetic scores       #######     
################################################### 


#Beta data = beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples
#Breakdown of the long label:  
#beta = betafile
#dropCpG = dropped CpG that did not pass QC
#IncCrossReac = including crossreactive probes
#noob = Noob normalized
#nodupl = collapsing EpicV2 probes that are on the same CpG site (to have lables in line with EpicV1 which is necessary for clock computation)
#dropSamples = dropping samples that did not pass QC

load("~/4_Processed Data/beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.rda")

#pheno data
load("~/4_Processed Data/pheno.dropSamples.df.rda") 
load("~/4_Processed Data/pheno.dropsamples.sw.df.rda") 

#### DNAm measures ####

# Here we create our DNAm measures of interest and add it to the phenotype data

#### Epigenetic-G with cross-reactive probes ####

# Scale methylation data across rows
beta_dataZ <- t(scale(t(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples)))

# Identify CpG probes present in both beta data and EpigeneticGProbes
EpigeneticGProbes<- fread("~/4_Processed Data/g_processed_Mean_Beta_PIP.txt") #see earlier code for reference

postqc_EpiGprobes <- Reduce(intersect, list(rownames(beta_dataZ), EpigeneticGProbes$CpG))
cat("Number of matching CpGs:", length(postqc_EpiGprobes), "\n")  # 678765 CpGs (so 566 more than the 678199  CpGs of the non-cross reactive probes score below)

# Match and subset EpigeneticGProbes data to the common CpGs
EpigeneticGProbes <- EpigeneticGProbes[match(postqc_EpiGprobes, EpigeneticGProbes$CpG), ]
EpigeneticGProbes$Coefficient <- as.numeric(as.character(EpigeneticGProbes$Mean_Beta))  # Convert coefficients to numeric

# Subset and handle missing values in beta data
beta_EpiGprobes <- beta_dataZ[match(postqc_EpiGprobes, rownames(beta_dataZ)), ]
beta_EpiGprobes[is.na(beta_EpiGprobes)] <- 0  # Impute missing values with 0

# Compute the epigenetic score for each sample
epi_g_sa_crrea <- unlist(lapply(1:ncol(beta_dataZ), function(x) {
  sum(EpigeneticGProbes$Coefficient * beta_EpiGprobes[, x])  
}))


# Before merging the epigenetic score to the phenotype data
methylation_samples <- colnames(beta_dataZ)
phenotypic_samples <- pheno.dropSamples.df$BaseID

# They have the same order: Are sample IDs identical and in the same order?  TRUE 
all_identical <- identical(methylation_samples, phenotypic_samples)
cat("Are sample IDs identical and in the same order? ", all_identical, "\n")


#### Add the epigenetic score to the phenotypic data ####
# merge it to the old phenotype data without swop, so the order of CpGs is the same
pheno.dropSamples.df$epi_g_sa_crrea <- epi_g_sa_crrea

epi_g_sa_crrea_data <- pheno.dropSamples.df[,c("BaseID",
                                               "epi_g_sa_crrea")]


save(epi_g_sa_crrea_data, file = file.path("~/4_Processed Data/epi_g_sa_crrea_data.rda"))
load("~/4_Processed Data/epi_g_sa_crrea_data.rda")

# add it to the phenotype data with the swopped samples

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  left_join(epi_g_sa_crrea_data, by = "BaseID")


#### Epigenetic-G with cross-reactive probes, separately for mothers and children ####

#For epigenetic-g we scale methylation data across rows. We discovered that doing this across mothers and children simulaniously affects the overall mean it should be 0 for children and 0 for mothers, but creating it across creates 0 across mothers and children)
# We therefore create epigentic-g here separatly for mothers and children, and add them later again as one variable to the phenotype data

#Load file with Betas for CpGs if not loaded in already
#load("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.rda")


#select only betas for children
# Subset BaseIDs in the phenotypic data where the column value is "C"
selected_base_ids <- pheno.dropsamples.sw.df$BaseID[pheno.dropsamples.sw.df$Code.C..child..M..mother. == "C"]

# Subset the large matrix based on the selected BaseIDs
beta_children <- beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples[, 
                                                                   colnames(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples) %in% selected_base_ids]


# Scale methylation data across rows
beta_dataZ <- t(scale(t(beta_children)))

# Identify CpG probes present in both beta data and EpigeneticGProbes
EpigeneticGProbes<- fread("~/4_Processed Data/g_processed_Mean_Beta_PIP.txt") #see earlier code for reference

postqc_EpiGprobes <- Reduce(intersect, list(rownames(beta_dataZ), EpigeneticGProbes$CpG))
cat("Number of matching CpGs:", length(postqc_EpiGprobes), "\n")  # 678765 CpGs 

# Match and subset EpigeneticGProbes data to the common CpGs
EpigeneticGProbes <- EpigeneticGProbes[match(postqc_EpiGprobes, EpigeneticGProbes$CpG), ]

EpigeneticGProbes$Coefficient <- as.numeric(as.character(EpigeneticGProbes$Mean_Beta))  # Convert coefficients to numeric

# Subset and handle missing values in beta data
beta_EpiGprobes <- beta_dataZ[match(postqc_EpiGprobes, rownames(beta_dataZ)), ]
beta_EpiGprobes[is.na(beta_EpiGprobes)] <- 0  # Impute missing values with 0

# Compute the epigenetic score for each sample
epi_g_sa_crrea_ChMumSep <- unlist(lapply(1:ncol(beta_dataZ), function(x) {
  sum(EpigeneticGProbes$Coefficient * beta_EpiGprobes[, x])  # No intercept added
}))


# Before merging the epigenetic score to the phenotype data

#create phenotype data with only kids
pheno.dropsamples.sw.df.Kidsonly <- pheno.dropsamples.sw.df[pheno.dropsamples.sw.df$Code.C..child..M..mother. == "C", ]

#check if methylation and phenotypic data is in same order
methylation_samples <- colnames(beta_dataZ)
phenotypic_samples <- pheno.dropsamples.sw.df.Kidsonly$BaseID

# They have the same order: Are sample IDs identical and in the same order?  TRUE 
all_identical <- identical(methylation_samples, phenotypic_samples)
cat("Are sample IDs identical and in the same order? ", all_identical, "\n")

#they are the same order so no need to re-arrange
sample_ids_beta <- colnames(beta_dataZ)
sample_ids_pheno <- pheno.dropsamples.sw.df.Kidsonly$BaseID
pheno.dropsamples.sw.df.Kidsonly  <- pheno.dropsamples.sw.df.Kidsonly[match(sample_ids_beta, sample_ids_pheno), ]

#### Add the epigenetic score to the phenotypic kids data

pheno.dropsamples.sw.df.Kidsonly$epi_g_sa_crrea_ChMumSep <- epi_g_sa_crrea_ChMumSep

describe(pheno.dropsamples.sw.df.Kidsonly$epi_g_sa_crrea_ChMumSep)


epi_g_sa_crrea_ChMumSep_Ch <- pheno.dropsamples.sw.df.Kidsonly[,c("BaseID",
                                                                  "epi_g_sa_crrea_ChMumSep")]


### do the same as above but then for mothers 

#select only betas for mothers
# Subset BaseIDs in the phenotypic data where the column value is "M"
selected_base_ids <- pheno.dropsamples.sw.df$BaseID[pheno.dropsamples.sw.df$Code.C..child..M..mother. == "M"]

# Subset the large matrix based on the selected BaseIDs
beta_mothers <- beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples[, 
                                                                  colnames(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples) %in% selected_base_ids]


# Scale methylation data across rows
beta_dataZ <- t(scale(t(beta_mothers)))

# Identify CpG probes present in both beta data and EpigeneticGProbes
EpigeneticGProbes<- fread("~/4_Processed Data/g_processed_Mean_Beta_PIP.txt") #see earlier code for reference

postqc_EpiGprobes <- Reduce(intersect, list(rownames(beta_dataZ), EpigeneticGProbes$CpG))
cat("Number of matching CpGs:", length(postqc_EpiGprobes), "\n")  # 678765 CpGs 

# Match and subset EpigeneticGProbes data to the common CpGs
EpigeneticGProbes <- EpigeneticGProbes[match(postqc_EpiGprobes, EpigeneticGProbes$CpG), ]

EpigeneticGProbes$Coefficient <- as.numeric(as.character(EpigeneticGProbes$Mean_Beta))  # Convert coefficients to numeric

# Subset and handle missing values in beta data
beta_EpiGprobes <- beta_dataZ[match(postqc_EpiGprobes, rownames(beta_dataZ)), ]
beta_EpiGprobes[is.na(beta_EpiGprobes)] <- 0  # Impute missing values with 0

# Compute the epigenetic score for each sample
epi_g_sa_crrea_ChMumSep <- unlist(lapply(1:ncol(beta_dataZ), function(x) {
  sum(EpigeneticGProbes$Coefficient * beta_EpiGprobes[, x])  # No intercept added
}))


# Before merging the epigenetic score to the phenotype data

#create phenotype data with only mothers
pheno.dropsamples.sw.df.Mumssonly <- pheno.dropsamples.sw.df[pheno.dropsamples.sw.df$Code.C..child..M..mother. == "M", ]

#check if methylation and phenotypic data is in same order
methylation_samples <- colnames(beta_dataZ)
phenotypic_samples <- pheno.dropsamples.sw.df.Mumssonly$BaseID

# They have the same order: Are sample IDs identical and in the same order?  TRUE 
all_identical <- identical(methylation_samples, phenotypic_samples)
cat("Are sample IDs identical and in the same order? ", all_identical, "\n")

#they are the same order so no need to re-arrange
sample_ids_beta <- colnames(beta_dataZ)
sample_ids_pheno <- pheno.dropsamples.sw.df.Mumssonly$BaseID
pheno.dropsamples.sw.df.Mumssonly  <- pheno.dropsamples.sw.df.Mumssonly[match(sample_ids_beta, sample_ids_pheno), ]

#### Add the epigenetic score to the phenotypic kids data

pheno.dropsamples.sw.df.Mumssonly$epi_g_sa_crrea_ChMumSep <- epi_g_sa_crrea_ChMumSep

describe(pheno.dropsamples.sw.df.Mumssonly$epi_g_sa_crrea_ChMumSep)


epi_g_sa_crrea_ChMumSep_Mum <- pheno.dropsamples.sw.df.Mumssonly[,c("BaseID",
                                                                    "epi_g_sa_crrea_ChMumSep")]


#combine the epigenetic-g data for mothers and children

epi_g_sa_crrea_ChMumSep_combined <- rbind(epi_g_sa_crrea_ChMumSep_Mum, epi_g_sa_crrea_ChMumSep_Ch)

describe(epi_g_sa_crrea_ChMumSep_combined$epi_g_sa_crrea_ChMumSep)

#save the data
save(epi_g_sa_crrea_ChMumSep_combined, file = file.path("~/4_Processed Data/epi_g_sa_crrea_ChMumSep_combined.rda"))
load("~/4_Processed Data/epi_g_sa_crrea_ChMumSep_combined.rda")

#add it to the full phenotype data

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  left_join(epi_g_sa_crrea_ChMumSep_combined, by = "BaseID")

#check means of mothers and children, mean should be 0 as it is standardized within the group

# For Code.C..child..M..mother. == "C"
describe_C <- describe(pheno.dropsamples.sw.df$epi_g_sa_crrea_ChMumSep[
  pheno.dropsamples.sw.df$Code.C..child..M..mother. == "C"
])

# For Code.C..child..M..mother. == "M"
describe_M <- describe(pheno.dropsamples.sw.df$epi_g_sa_crrea_ChMumSep[
  pheno.dropsamples.sw.df$Code.C..child..M..mother. == "M"
])

#### DunedinPACE ####

#DunedinPACE in chunks cause of problems with memory. 
# We tested whether it makes a difference to compute DunedinPACE for mothers and children together or seperately, and it shows it doesnt matter. This does matter for epigenetic-G as you standardize across rows before computing the score
# I recreated DunedinPACE based on 5 chunks per time, instead of 20 as a sanity check. Results are the same, so chunk size does not influence Dunedinpace score

devtools::install_github("danbelsky/DunedinPACE", build_vignettes = FALSE)
library(DunedinPACE)

#load in CpGs 
load("--/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.dropSamples.rda")

# Initialize a list to store the outputs
pace_outputs <- list()

# Get the total number of columns in the dataset
total_samples <- ncol(beta.dropCpG.IncCrossReac.noob.dropSamples)

# Split the data into chunks of 20 columns
chunk_size <- 20
num_chunks <- ceiling(total_samples / chunk_size)

for (i in 1:num_chunks) {
  # Determine column indices for the current chunk
  start_col <- (i - 1) * chunk_size + 1
  end_col <- min(i * chunk_size, total_samples)
  
  # Subset the data for the current chunk
  beta_subset <- beta.dropCpG.IncCrossReac.noob.dropSamples[, start_col:end_col]
  
  # Run the PACEProjector function
  pace_output <- PACEProjector(beta_subset, 0.7)
  
  # Convert to data frame and store in the list
  pace_outputs[[paste0("chunk_", i)]] <- as.data.frame(pace_output)
  
  # Optionally, write the output to a file (e.g., CSV)
  output_filename <- paste0("pace_output_chunk_", i, ".csv")
  write.csv(as.data.frame(pace_output), file = output_filename, row.names = FALSE)
}

# Combine all outputs into one data frame if needed
combined_output <- do.call(rbind, pace_outputs)

# Save the combined output to a file
write.csv(combined_output, file = "combined_pace_output.csv", row.names = FALSE)

# Print the combined output or the list of outputs
print(pace_outputs)

# Move rownames to a new column called BaseID
combined_output$BaseID <- rownames(combined_output)

# Reset the rownames to default integers
rownames(combined_output) <- NULL

# Remove "chunk_X." from the BaseID column
combined_output$BaseID <- gsub("chunk_\\d+\\.", "", combined_output$BaseID)

# Trim any leading or trailing spaces that might remain
combined_output$BaseID <- trimws(combined_output$BaseID)

# View the modified data
print(head(combined_output))

DunedinPACE_chunks_incCrossReact <- combined_output

save(DunedinPACE_chunks_incCrossReact, file = "/DunedinPACE_chunks_incCrossReact.rda")


load("--/data/4_Processed Data/DunedinPACE_chunks_incCrossReact.rda")
colnames(DunedinPACE_chunks_incCrossReact)

colnames(DunedinPACE_chunks_incCrossReact)[colnames(DunedinPACE_chunks_incCrossReact) == "DunedinPACE"] <- "DunedinPACE_crrea"


pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  left_join(DunedinPACE_chunks_incCrossReact, by = "BaseID")

colnames(pheno.dropsamples.sw.df)


##### Creating DNAm age based clocks with MethylClock ####

#to compute skin-horvath to test age in the sample 
#this we computed to see if "kids" samples are really kids, and "mother" samples were really mums. 

BiocManager::install("methylclockData")
library(methylclockData)

#check missing CpGs for clocks we are interested in
cpgs.missing <- checkClocks(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples)

#this shows that when computing it including cross-reactive probes, there are no differences
#so we can stick with our earlier data!

#compute clocks
#so here we use package methylclock to compute "Horvath", "Hannum", "Levine", "skinHorvath", "PedBE", "Wu",
#we do use here the CpG file without duplicates, otherwise it does not recognize the CpG labels
bioage <- DNAmAge(beta.dropCpG.noob.nodupl.dropSamples) 
bioage

head(bioage)
colnames(bioage)

colnames(bioage)[colnames(bioage) == "id"] <- "BaseID"

save(bioage, file="~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/bioage.rda") 
load("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/bioage.rda")

#### Add DNAm aging clocks to pheno 

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  mutate(BaseID = as.character(BaseID))

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  left_join(bioage, by = "BaseID")

colnames(pheno.dropsamples.sw.df)


#### Data with all DNAm before residualizing them

pheno.dropsamples.sw.df.DNAmmeasures <- pheno.dropsamples.sw.df

save(pheno.dropsamples.sw.df.DNAmmeasures, file="--/data/4_Processed Data/pheno.dropsamples.sw.df.DNAmmeasures.rda") 
load("--/data/4_Processed Data/pheno.dropsamples.sw.df.DNAmmeasures.rda")

colnames(pheno.dropsamples.sw.df.DNAmmeasures)

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df.DNAmmeasures

#### Create residualized DNAm measures ####

#We create residualized DNAm measures, seperately for mothers and children

# Filter data so we create dataset for mothers and children 
mothers_data <- pheno.dropsamples.sw.df[pheno.dropsamples.sw.df$Mum0Child1 == 0 , ]
children_data <- pheno.dropsamples.sw.df[pheno.dropsamples.sw.df$Mum0Child1 == 1 , ]

# skinHorvath Mothers
reg_skinHorvath_mum <- lm(skinHorvath ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                            as.factor(plate) + as.factor(Slide),
                          data = mothers_data)
mothers_data$skinHorvath_res_cell_tech <- residuals(reg_skinHorvath_mum)
mothers_data$skinHorvath_res_cell_tech_Z <- as.numeric(scale(residuals(reg_skinHorvath_mum)))

# skinHorvath Children
reg_skinHorvath_child <- lm(skinHorvath ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                              as.factor(plate) + as.factor(Slide),
                            data = children_data)
children_data$skinHorvath_res_cell_tech <- residuals(reg_skinHorvath_child)
children_data$skinHorvath_res_cell_tech_Z <- as.numeric(scale(residuals(reg_skinHorvath_child)))


# epi_g_sa_crrea Mothers
reg_epi_g_sa_crrea_mum <- lm(epi_g_sa_crrea ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                               as.factor(plate) + as.factor(Slide),
                             data = mothers_data)
mothers_data$epi_g_sa_crrea_res_cell_tech <- residuals(reg_epi_g_sa_crrea_mum)
mothers_data$epi_g_sa_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_epi_g_sa_crrea_mum)))

# epi_g_sa_crrea Children
reg_epi_g_sa_crrea_child <- lm(epi_g_sa_crrea ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                                 as.factor(plate) + as.factor(Slide),
                               data = children_data)
children_data$epi_g_sa_crrea_res_cell_tech <- residuals(reg_epi_g_sa_crrea_child)
children_data$epi_g_sa_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_epi_g_sa_crrea_child)))


# epi_g_sa_crrea_ChMumSep for Mothers
reg_epi_g_sa_crrea_ChMumSep_mum <- lm(epi_g_sa_crrea_ChMumSep ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                                        as.factor(plate) + as.factor(Slide),
                                      data = mothers_data)
mothers_data$epi_g_sa_crrea_ChMumSep_res_cell_tech <- residuals(reg_epi_g_sa_crrea_ChMumSep_mum)
mothers_data$epi_g_sa_crrea_ChMumSep_res_cell_tech_Z <- as.numeric(scale(residuals(reg_epi_g_sa_crrea_ChMumSep_mum)))

# epi_g_sa_crrea_ChMumSep for Children
reg_epi_g_sa_crrea_ChMumSep_child <- lm(epi_g_sa_crrea_ChMumSep ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                                          as.factor(plate) + as.factor(Slide),
                                        data = children_data)
children_data$epi_g_sa_crrea_ChMumSep_res_cell_tech <- residuals(reg_epi_g_sa_crrea_ChMumSep_child)
children_data$epi_g_sa_crrea_ChMumSep_res_cell_tech_Z <- as.numeric(scale(residuals(reg_epi_g_sa_crrea_ChMumSep_child)))

# DunedinPACE_crrea Mothers
reg_DunedinPACE_crrea_mum <- lm(DunedinPACE_crrea ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                                  as.factor(plate) + as.factor(Slide),
                                data = mothers_data)
mothers_data$DunedinPACE_crrea_res_cell_tech <- residuals(reg_DunedinPACE_crrea_mum)
mothers_data$DunedinPACE_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_DunedinPACE_crrea_mum)))

# DunedinPACE_crrea Children
reg_DunedinPACE_crrea_child <- lm(DunedinPACE_crrea ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                                    as.factor(plate) + as.factor(Slide),
                                  data = children_data)
children_data$DunedinPACE_crrea_res_cell_tech <- residuals(reg_DunedinPACE_crrea_child)
children_data$DunedinPACE_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_DunedinPACE_crrea_child)))


summary_stats <- function(x) {
  c(
    Mean = round(mean(x, na.rm = TRUE), 3),
    SD = round(sd(x, na.rm = TRUE), 3),
    Min = round(min(x, na.rm = TRUE), 3),
    Max = round(max(x, na.rm = TRUE), 3)
  )
}


colnames(mothers_data)

#double check if Mean=0, SD=1
summary_stats(mothers_data$epi_g_sa_crrea_res_cell_tech_Z)
summary_stats(children_data$epi_g_sa_crrea_res_cell_tech_Z)

summary_stats(mothers_data$epi_g_sa_crrea_ChMumSep_res_cell_tech_Z)
summary_stats(children_data$epi_g_sa_crrea_ChMumSep_res_cell_tech_Z)

summary_stats(mothers_data$DunedinPACE_crrea_res_cell_tech_Z)
summary_stats(children_data$DunedinPACE_crrea_res_cell_tech_Z)

#check specifically also for epigenetic-g standardized seperatly for mothers and children

summary_stats(mothers_data$epi_g_sa_crrea_ChMumSep)  #M=0,00 and SD =0.20 , so that is correct!
summary_stats(children_data$epi_g_sa_crrea_ChMumSep) #M=0,00 and SD =0.22


## re-merge children and mothers data

# Compare column names
mothers_cols <- colnames(mothers_data)
children_cols <- colnames(children_data)

# Check if column names match but may be in a different order
setequal(mothers_cols, children_cols) # TRUE if all names are the same, regardless of order

# Identify mismatched or out-of-order columns
which(mothers_cols != children_cols) # Identify positions where column names differ

# Combine the two datasets row-wise
combined_data <- rbind(mothers_data, children_data)



combined_data <- combined_data[c("BaseID",
                                 "skinHorvath_res_cell_tech",
                                 "skinHorvath_res_cell_tech_Z",
                                 "epi_g_sa_crrea_res_cell_tech",
                                 "epi_g_sa_crrea_res_cell_tech_Z",
                                 "epi_g_sa_crrea_ChMumSep_res_cell_tech",
                                 "epi_g_sa_crrea_ChMumSep_res_cell_tech_Z",
                                 "DunedinPACE_crrea_res_cell_tech",
                                 "DunedinPACE_crrea_res_cell_tech_Z")]

pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df %>%
  left_join(combined_data , by = "BaseID")

colnames(pheno.dropsamples.sw.df)

pheno.dropsamples.sw.df_full <- pheno.dropsamples.sw.df
pheno.dropsamples.sw.df <- pheno.dropsamples.sw.df_full 

save(pheno.dropsamples.sw.df_full, file="--/data/4_Processed Data/pheno.dropsamples.sw.df_full.rda") 
load("--/data/4_Processed Data/pheno.dropsamples.sw.df_full.rda")

pheno.dropsamples.sw.df_full_14Jan <- pheno.dropsamples.sw.df_full

save(pheno.dropsamples.sw.df_full_14Jan, file="--/data/4_Processed Data/pheno.dropsamples.sw.df_full_14Jan.rda") 
load("--/data/4_Processed Data/pheno.dropsamples.sw.df_full_14Jan.rda")

### look at mean differences Children and Mothers on DNAm aging measures ####

# generate the plot
pdf("--/data analysis/00_QC/Plots_full/MeanDiff Mothers Children.pdf") #open PDF device to save the plot


# Filter data where Techreplicate == 0
filtered_data <- pheno.dropsamples.sw.df[pheno.dropsamples.sw.df$Techreplicate == 0, ]

# Create a grouping variable based on Code.C..child..M..mother.
filtered_data$Group <- ifelse(filtered_data$Code.C..child..M..mother. == "C", "Children", "Mothers")

# List of variables to create boxplots for
variables <- c("skinHorvath", "epi_g_sa_crrea","epi_g_sa_crrea_ChMumSep", "DunedinPACE_crrea")

# Loop through each variable and create a boxplot
for (var in variables) {
  boxplot(
    filtered_data[[var]] ~ filtered_data$Group,
    main = paste("Comparison of", var, "Means"),
    ylab = var,
    xlab = "Group",
    col = c("lightblue", "lightgreen"), # Optional colors
    notch = TRUE # Notches indicate confidence intervals for medians
  )
  
  # Calculate means for the groups
  means <- aggregate(filtered_data[[var]] ~ filtered_data$Group, FUN = mean)
  
  # Add mean points to the boxplot
  points(1:2, means[, 2], col = "red", pch = 19)
  
  # Pause to display each plot in an interactive session
  Sys.sleep(1) # Adjust delay if needed
}


dev.off() #closes the PDF device, ensuring that the file is properly saved.

###### Pheno With Epi-G, Aging Clocks, DunedinPACE, Pheno

pheno.dropsamples.sw.df.DNAm_14Jan <- pheno.dropsamples.sw.df

save(pheno.dropsamples.sw.df.DNAm_14Jan, file="--/4_Processed Data/pheno.dropsamples.sw.df.DNAm_14Jan.rda") 
load("--/4_Processed Data/pheno.dropsamples.sw.df.DNAm_14Jan.rda")



###### Create overall dataset ####

load("--/data/4_Processed Data/pheno.dropsamples.sw.df.DNAm_14Jan.rda")

pheno.dropsamples.sw.df.DNAm <- pheno.dropsamples.sw.df.DNAm_14Jan

#Do a couple of cross-checks if this is the correct data
#Check colnames, are all clocks in there?
colnames(pheno.dropsamples.sw.df.DNAm)

#Check if the first samples are starting with 208226070014 (cause those are on plate 1)
head(pheno.dropsamples.sw.df.DNAm)

#Check if in this dataset, the BaseIDs were indeed swopped
load("--/data/4_Processed Data/swoppedIDs.rda")

#Check if in the dataset, skinHorvath for kids does not exceed 11 and for mothers is not lower than 11

# Subset for mothers (Code = "M")
mothers_stats <- pheno.dropsamples.sw.df.DNAm %>%
  filter(Code.C..child..M..mother. == "M") %>%
  summarise(
    Mean = mean(skinHorvath, na.rm = TRUE),
    Min = min(skinHorvath, na.rm = TRUE),
    Max = max(skinHorvath, na.rm = TRUE)
  )

# Subset for children (Code = "C")
children_stats <- pheno.dropsamples.sw.df.DNAm %>%
  filter(Code.C..child..M..mother. == "C") %>%
  summarise(
    Mean = mean(skinHorvath, na.rm = TRUE),
    Min = min(skinHorvath, na.rm = TRUE),
    Max = max(skinHorvath, na.rm = TRUE)
  )


print(mothers_stats) 
print(children_stats)



# Create dataset for First pre-registration paper

data_DNAmScores_BFY <- pheno.dropsamples.sw.df.DNAm


# Select variables of interest
colnames(data_DNAmScores_BFY)

data_DNAmScores_BFY_InclTechRep_14Jan <- data_DNAmScores_BFY

save(data_DNAmScores_BFY_InclTechRep_14Jan, file="--/data_DNAmScores_BFY_InclTechRep_14Jan.rda") 
load("--/data_DNAmScores_BFY_InclTechRep_14Jan.rda")

#also save as .csv file
write.csv(data_DNAmScores_BFY_InclTechRep_14Jan, 
          file = "--/data_DNAmScores_BFY_InclTechRep_14Jan.csv", 
          row.names = FALSE)


colnames(data_DNAmScores_BFY_InclTechRep_14Jan)


# Omit technical replicates
data_DNAmScores_BFY_ExclTechRep_14Jan <- data_DNAmScores_BFY[data_DNAmScores_BFY$Techreplicate == 0, ]

save(data_DNAmScores_BFY_ExclTechRep_14Jan, file="--/data_DNAmScores_BFY_ExclTechRep_14Jan.rda") 
load("--/data_DNAmScores_BFY_ExclTechRep_14Jan.rda")

#also save as csv file
write.csv(data_DNAmScores_BFY_ExclTechRep_14Jan, 
          file = "--/data_DNAmScores_BFY_ExclTechRep_14Jan.csv", 
          row.names = FALSE)


# Check descriptives
table(data_DNAmScores_BFY_ExclTechRep_14Jan$Mum0Child1) #735 children, #777 mothers
table(data_DNAmScores_BFY_ExclTechRep_14Jan$Mum0Child1, data_DNAmScores_BFY_ExclTechRep_14Jan$sex) #in kids, 369 boys, 366 girls

#count number of mother-child pairs
library(dplyr)

# Count the occurrences of each block
block_counts <- data_DNAmScores_BFY_ExclTechRep_14Jan %>%
  group_by(block) %>%
  summarise(
    count = n(),
    mothers = sum(Mum0Child1 == 0, na.rm = TRUE),
    children = sum(Mum0Child1 == 1, na.rm = TRUE)
  )

# Count family pairs (block IDs that appear exactly twice)
num_family_pairs <- block_counts %>%
  filter(count == 2) %>%
  nrow()

# Count unpaired mothers (block IDs where only mothers are present)
unpaired_mothers <- block_counts %>%
  filter(count != 2, mothers > 0) %>%
  summarise(total_unpaired_mothers = sum(mothers)) %>%
  pull(total_unpaired_mothers)

# Count unpaired children (block IDs where only children are present)
unpaired_children <- block_counts %>%
  filter(count != 2, children > 0) %>%
  summarise(total_unpaired_children = sum(children)) %>%
  pull(total_unpaired_children)

# Output results
cat("Number of family pairs:", num_family_pairs, "\n")
cat("Number of unpaired mothers:", unpaired_mothers, "\n")
cat("Number of unpaired children:", unpaired_children, "\n")


### On SILO servers ####

# the phenotype data with age cannot be transfered, so PhenoAge and GrimAge have been calculated on the secured server.


#### Create Pheno Age Acceleration and Grim Age Acceleration ####

#load packages
library("dplyr")
library("tibble")
library("tidyr")
library("glmnet")
library("devtools")


#Import phenotype data
data_DNAmScores_BFY_age <- read_csv("--/clean_data/data_DNAmScores_BFY_age_16Jan.csv")
colnames(data_DNAmScores_BFY_age)

describe(data_DNAmScores_BFY_age$age4_mothers_inyears)

#Import phenotype data
data_DNAmScores_BFY_age_2 <- read_csv("--/clean_data/03_data_DNAmScores_BFY_InclTechRep_14Jan.csv")
colnames(data_DNAmScores_BFY_age_2)

#Import Betas
load("--/raw_data/beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.rda")

betas <- as.data.frame(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples)
pheno <- data_DNAmScores_BFY_age

#the original beta file still contains technical replicates, which we remove here
betas_bl <- betas
pheno_bl <- pheno

#transpose betas
betas_blT = t(betas_bl)

#check if order of samples is the same
head(pheno_bl$BaseID)   
head(rownames(betas_blT))
sample_ids_beta <-rownames(betas_blT)
sample_ids_pheno <- pheno_bl$BaseID

pheno_bl <- pheno_bl[match(sample_ids_beta, sample_ids_pheno), ]

head(betas_bl)

#load in files from PC clocks
load("/data4/project/irp4/BabyFirstYear/R/4.0.3/RPackages and Scripts DNAm clocks/clock data updated/CalcAllPCClocks.RData")

betas_blT <- as.data.frame(betas_blT)
length(colnames(betas_blT))

#here you set the missing CpGs of pheno
missingCpGs <- CpGs[!(CpGs %in% colnames(betas_blT))] 

betas_blT[,missingCpGs] <- NA

for(i in 1:length(missingCpGs)){
  betas_blT[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
}

#Prepare methylation data for calculation of PC Clocks 
betas_blT <- betas_blT[,CpGs]
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
betas_blT <- apply(betas_blT,2,meanimpute)

# If you did not calculate the original CpG-based clocks, initialize a data frame for PC clocks
pheno_bl= as.data.frame(pheno_bl) #it was already a phenotype
DNAmAge_bl <- pheno_bl   ##replace DNAmAge with DNAmAge_bl

#Double check that order of samples is still correct
sum(DNAmAge_bl$BaseID == pheno_bl$BaseID)
head(DNAmAge_bl$BaseID);  head (rownames(betas_blT))

#Loading PCPhenoAge model
load("/data4/project/irp4/BabyFirstYear/R/4.0.3/RPackages and Scripts DNAm clocks/clock data updated/CalcPCPhenoAge.RData")

#This line loads the model for PCPhenoAge that was previously saved as an RData file.
datPCPhenoAge = predict(PCPhenoAge, betas_blT)

#Predicting
library(glmnet)             # Load the package

#the predict function uses the loaded PCPhenoAge model to estimate the PCPhenoAge values using your methylation data (betas_epic_buT)
DNAmAge_bl$PCPhenoAge = as.numeric(predict.glmnet(fit, datPCPhenoAge, s = cv$lambda.min))

library(psych)

describe(DNAmAge_bl$PCPhenoAge)


#### Create GrimAge - way that Sepideh & Laurel do
pheno_bl$age4_mothers_children_inyears

#converting age and sex to numeric types
pheno_bl$sex = as.numeric(pheno_bl$sex)
pheno_bl$age4_mothers_children_inyears = as.numeric(pheno_bl$age4_mothers_children_inyears)

#load in grimage
load("/data4/project/irp4/BabyFirstYear/R/4.0.3/RPackages and Scripts DNAm clocks/clock data updated/CalcPCGrimAge.RData")

#prediction
datPCGrimAge = predict(PCGrimAge, betas_blT) [,-3935] #remove the duplicated one


#add the 'Female' and 'Age' variables to the datPCGrimAge dataframe.
datPCGrimAge = cbind(datPCGrimAge, Sex = pheno_bl$sex, Age = pheno_bl$age4_mothers_children_inyears)
#This loop applies the predict.glmnet function to each of the GrimAge components (except the first one), using the respective fits and lambda values. 
#The resulting values are added as new columns to the DNAmAge_bl dataframe.

for(i in 2:9){
  DNAmAge_bl[,GrimAgeComponents[i]] = as.numeric(predict.glmnet(fit[[i]], datPCGrimAge, s = cv[[i]]$lambda.min))
}

#create a new matrix that contains only the GrimAge components, 'Age', and 'Female' columns from the DNAmAge_bl dataframe.
datGrimAge = data.matrix(DNAmAge_bl)[,c(GrimAgeComponents[2:9],"age4_mothers_children_inyears","sex")]

#This line applies the predict.glmnet function to the 10th GrimAge component and adds the resulting values as a new column
DNAmAge_bl[,GrimAgeComponents[10]] = as.numeric(predict.glmnet(fit[[10]], newx = as(datGrimAge, "dgCMatrix"), s = cv[[10]]$lambda.min))

pheno_bl$PCPhenoAge_bl <- DNAmAge_bl$PCPhenoAge
pheno_bl$PCGrimAge_bl <- DNAmAge_bl$PCGrimAge

describe(pheno_bl$PCPhenoAge_bl) 
describe(pheno_bl$PCGrimAge_bl)  


##### Below we generate pheno and grim age for mothers and children 

#Omit technical replicates and only select those with age available

Mum_data0 <- Mum_data
Mum_data <- Mum_data0

Mum_data <- Mum_data[Mum_data$Techreplicate==0,]
Mum_data <- Mum_data[!is.na(Mum_data$age4_mothers_inyear),]

Child_data <- Child_data[Child_data$Techreplicate==0,]
Mum_data <- Mum_data[!is.na(Mum_data$age4_mothers_inyear),]

#We use the one based on the original script, so PCPhenoAge_bl and PCGrimAge_bl

#Create residualized values for mothers
reg_Pheno_mum <- lm(PCPhenoAge_bl ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                      as.factor(plate) + as.factor(Slide) + age4_mothers_inyears,
                    data=Mum_data)
Mum_data$PC_PhenoAge_Accel_crrea_res_cell_tech <- residuals(reg_Pheno_mum)
Mum_data$PC_PhenoAge_Accel_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_Pheno_mum)))

reg_Grim_mum <- lm(PCGrimAge_bl ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                     as.factor(plate) + as.factor(Slide) + age4_mothers_inyears,
                   data=Mum_data)
Mum_data$PC_GrimAge_Accel_crrea_res_cell_tech <- residuals(reg_Grim_mum)
Mum_data$PC_GrimAge_Accel_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_Grim_mum)))


#Create residualized values forchildren
reg_Pheno_child <- lm(PCPhenoAge_bl ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                        as.factor(plate) + as.factor(Slide) + age4_children_inyears,
                      data=Child_data)
Child_data$PC_PhenoAge_Accel_crrea_res_cell_tech <- residuals(reg_Pheno_child)
Child_data$PC_PhenoAge_Accel_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_Pheno_child)))

reg_Grim_child <- lm(PCGrimAge_bl ~ celltype_1 + celltype_2 + celltype_3 + celltype_4 + celltype_5 +
                       as.factor(plate) + as.factor(Slide) + age4_children_inyears,
                     data=Child_data)
Child_data$PC_GrimAge_Accel_crrea_res_cell_tech <- residuals(reg_Grim_child)
Child_data$PC_GrimAge_Accel_crrea_res_cell_tech_Z <- as.numeric(scale(residuals(reg_Grim_child)))

#check if with all M=0, sd=1
describe(Mum_data$PC_GrimAge_Accel_crrea_res_cell_tech_Z)
describe(Mum_data$PC_PhenoAge_Accel_crrea_res_cell_tech_Z)
describe(Child_data$PC_GrimAge_Accel_crrea_res_cell_tech_Z)
describe(Child_data$PC_PhenoAge_Accel_crrea_res_cell_tech_Z)

#check overall descriptives
describe(Mum_data$PC_GrimAge_Accel_crrea_res_cell_tech)
describe(Child_data$PC_GrimAge_Accel_crrea_res_cell_tech)

describe(Mum_data$PC_PhenoAge_Accel_crrea_res_cell_tech)
describe(Child_data$PC_PhenoAge_Accel_crrea_res_cell_tech)

describe(Mum_data$PCGrimAge_bl)
describe(Child_data$PCGrimAge_bl)

describe(Mum_data$PCPhenoAge_bl)
describe(Child_data$PCPhenoAge_bl)

cor.test(Child_data$DunedinPACE_crrea_res_cell_tech_Z, Child_data$PC_PhenoAge_Accel_crrea_res_cell_tech_Z)

#Merge mother and child data
mothers_cols <- colnames(Mum_data)
children_cols <- colnames(Child_data)
setequal(mothers_cols, children_cols) #colnames are the same

combined_data <- rbind(Mum_data, Child_data)

colnames(combined_data)
colnames(combined_data)[colnames(combined_data) == "PCGrimAge_bl"] <- "PC_GrimAge_crrea"
colnames(combined_data)[colnames(combined_data) == "PCPhenoAge_bl"] <- "PC_PhenoAge_crrea"



PhenoGrim_BFY_17jan2025 <- combined_data[,c("BaseID",
                                            "Sample_Name",
                                            "sampleId",
                                            "block",
                                            "Proben_Name",
                                            "PC_PhenoAge_crrea",
                                            "PC_PhenoAge_Accel_crrea_res_cell_tech",
                                            "PC_PhenoAge_Accel_crrea_res_cell_tech_Z",
                                            "PC_GrimAge_crrea",
                                            "PC_GrimAge_Accel_crrea_res_cell_tech",
                                            "PC_GrimAge_Accel_crrea_res_cell_tech_Z")]


write_csv(PhenoGrim_BFY_17jan2025, "--/clean_data/PhenoGrim_BFY_17jan2025.csv")



################# The End ##############################################                                              







