#!/usr/bin/Rscript

# PART ONE ----

# Set working directory
setwd("/home/nyamarim/Marion/bgeecall_exercice")

# install downloader from CRAN
install.packages("downloader")

library(downloader)
#download("https://github.com/wch/downloader/zipball/master",
#         "downloader.zip", mode = "wb")

# Specify URL where file is stored
#url <- "http://ftp.ensemblgenomes.org/pub/release-51/metazoa/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.51.chr.gtf.gz"
#url1 <- "http://ftp.ensemblgenomes.org/pub/release-51/metazoa/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"

# Specify destination where file should be saved
#destfile <- "C:/Users/Marion/Desktop/BgeeCall-assignment/Drosophila_melanogaster.BDGP6.32.51.chr.gtf.gz"
#destfile1 <- "C:/Users/Marion/Desktop/BgeeCall-assignment/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"

# Download files
#download.file(url, destfile)
#download.file(url1, destfile1)

# Install BgeeCall

#BiocManager::install("BgeeCall")

library(BgeeCall)


# List of all intergenic releases 
{
  list_intergenic_release()
}

# Verify which species are available for the current Bgee intergenic release ()
{
  bgee <- new("BgeeMetadata")
  list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)
  list_bgee_ref_intergenic_species(release = '1.0')
}

# Species that belong to the community. 
{
  list_community_ref_intergenic_species()
}

# Creating an object of the KallistoMetadata class

kallisto <- new("KallistoMetadata", download_kallisto = TRUE, cutoff=0.05)

# Creating a userMetadata object 

Metadata <- new("UserMetadata", species_id = '7227', reads_size = 732.3)

# set path to RNASeq library
{
  path <- "/home/nyamarim/Marion/bgeecall_exercice/SRX109278"
  Metadata <- setRNASeqLibPath(Metadata, path)
}

# set path to transcriptome
{
  transcriptome_path <- "/home/nyamarim/Marion/Drosophila_melanogaster.BDGP6.32.cdna.all.fa"
  Metadata <- setTranscriptomeFromFile(Metadata, transcriptome_path,
                                       "Drosophila_melanogaster")
}

# set path to gtf annotation file
{
  annotation_file <- "/home/nyamarim/Marion/Drosophila_melanogaster.BDGP6.32.51.chr.gtf"
  Metadata <- setAnnotationFromFile(Metadata, annotation_file,
                                    "Drosophila_melanogaster.BDGP6.32.51.chr.gtf.gz")
}

# set path to output directory
{
dir.create("/home/nyamarim/Marion/results")
results <- "/home/nyamarim/Marion/results"
Metadata <- setOutputDir(Metadata, results)
## Abundance file is missing : C:/Users/Marion/Desktop/BgeeCall-assignment/bgeecall_exercice/abundance.tsv.
}

# set path to working directory
Metadata <- setWorkingPath(Metadata, getwd())

## Generate  present and absent calls for the library SRX109278 
# generate transcriptome index

{
  myBgeeMetadata <- new("BgeeMetadata")
  create_kallisto_index(kallisto, myBgeeMetadata, Metadata, transcriptome_path)
}

# run generate_calls_workflow
calls_output <- generate_calls_workflow(abundanceMetadata = kallisto, bgeeMetadata = myBgeeMetadata, userMetadata = Metadata)


## Error in density.default(log2(counts$abundance[selected_intergenic] +  : 
## need at least 2 points to select a bandwidth automatically

# Plot the frequency of p-values for the correspondent library

#PART TWO ----

# Set working directory
setwd("/home/nyamarim/Marion/bgeecall_exercice")

library(BgeeCall)

#import TSV file into data frame
UserMetadataTemplate <- read.table(file = "/home/nyamarim/Marion/bgeecall_exercice/userMetadataTemplate.tsv.1", sep = '\t', header = TRUE)

# Remove columns I won't use
UserMetadataTemplate <- subset(UserMetadataTemplate, select = c('run_ids', 'output_directory', 'custom_intergenic_path'))


# Creating data frame 
t_path <- "/home/nyamarim/Marion/Drosophila_melanogaster.BDGP6.32.cdna.all.fa"
an_path <- "/home/nyamarim/Marion/Drosophila_melanogaster.BDGP6.32.51.chr.gtf"
wk_path <- "/home/nyamarim/Marion/bgeecall_excercice"



SRX109279 <- data.frame("7227", "581.5", "/home/nyamarim/Marion/bgeecall_exercice/SRX109279", t_path, an_path,
                        wk_path)      
SRX109278 <- data.frame("7227", "732.3", "/home/nyamarim/Marion/bgeecall_exercice/SRX109278", t_path, an_path,
                        wk_path)
SRX1720957 <- data.frame("7227", "786.3", "/home/nyamarim/Marion/bgeecall_exercice/SRX1720957", t_path, an_path,
                         wk_path)
SRX493950 <- data.frame("7227", "604.2", "/home/nyamarim/Marion/bgeecall_exercice/SRX493950", t_path, an_path,
                        wk_path)
SRX1720958 <- data.frame("7227", "868.2", "/home/nyamarim/Marion/bgeecall_exercice/SRX1720958", t_path, an_path,
                         wk_path)
SRX493999 <- data.frame("7227", "254.1", "/home/nyamarim/Marion/bgeecall_exercice/SRX493999", t_path, an_path,
                        wk_path)

#Naming the Data Frame - Step 2  
names(SRX109279) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  
names(SRX109278) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  
names(SRX1720957) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  
names(SRX493950) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  
names(SRX1720958) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  
names(SRX493999) <- c("species_id", "reads_size", "rnaseq_lib_path", "transcriptome_path", "annotation_path", "working_path")  

#Using rbind() function to insert above observation 
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX109279)
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX109278)
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX1720957)
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX493950)
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX1720958)
UserMetadataTemplate <- rbind(UserMetadataTemplate, SRX493999)


# writing a tsv file
write.table(UserMetadataTemplate, file = "/home/nyamarim/Marion/userMetadata.tsv", row.names=FALSE, sep="\t")

# Define path to userfile 
uf <- "/home/nyamarim/Marion/userMetadata.tsv"

# run generation of present and absent calls from the user file with default values

calls_output <- generate_calls_workflow(userFile = uf)

#gene_cutoff_info <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109279/gene_cutoff_info_file.tsv")
#head(gene_cutoff_info)

#gene_abundance <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109278/abundance.tsv")
#head(gene_abundance)

#calls_tsv <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109278/gene_level_abundance+calls.tsv")
#head(calls_tsv)

#PART THREE ----

#import TSV file into data frame
data <- read.table(file = "/home/nyamarim/Marion/userMetadata.tsv", sep = '\t', header = TRUE)
head(data)

# assigning the column name to a new name
#colnames(data)[1] <- "specied_id"

# writing a new tsv file
write.table(data, file = "/home/nyamarim/Marion/userMetadata_final.tsv", row.names=FALSE, sep="\t")

# Combine multiple libraries per species using the merging_libraries() function
mergingLibraries_pValue <- merging_libraries(userFile ="/home/nyamarim/Marion/userMetadata_final.tsv", approach = "BH", condition = "species_id", outDir = "/home/nyamarim/Marion/")

#Using 0 libraries for condition: species_id devStage  with values: 7227 FBdv:00007079 
#Error in approachesMerging(allFiles = allFiles, approach = approach, cutoff = cutoff) : 
#  Select the appropriated method for your quantitative metric: BH for p-values OR fdr_inverse for q-values

# Append condition 'devStage' to all libraries

{
  data$devStage <- c("FBdv:00007079", "FBdv:00007079", "UBERON:0000066", "FBdv:00007079", "UBERON:0000066", "FBdv:00007079")
  
  write.table(data, file = "/home/nyamarim/Marion/userMetadata_final1.tsv", row.names=FALSE, sep="\t")
}

# Combine libraries per species (species_id, devStage) 
mergingLibraries <- merging_libraries(userFile = "/home/nyamarim/Marion/userMetadata_final1.tsv", approach = "BH", condition = c("species_id", "devStage"), cutoff = 0.01, outDir = "/home/nyamarim/Marion/")

#Using 0 libraries for condition: species_id devStage  with values: 7227 FBdv:00007079 
#Error in approachesMerging(allFiles = allFiles, approach = approach, cutoff = cutoff) : 
#  Select the appropriated method for your quantitative metric: BH for p-values OR fdr_inverse for q-values

# Get summary stats of all libraries by using get_summary_stats() function.

stats <- read.table(calls_output$S4_slots_summary, header = TRUE, sep = "\t")
#Error in read.table(calls_output$S4_slots_summary, header = TRUE, sep = "\t") : 
#  'file' must be a character string or connection

# read summary stats per library
{
  summ <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109279/S4_slots_summary.tsv", sep = "\t")
  #head(summ)
  
  summ1 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX109278/S4_slots_summary.tsv", sep = "\t")
  
  summ2 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX493950/S4_slots_summary.tsv", sep = "\t")
  
  summ3 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX493999/S4_slots_summary.tsv", sep = "\t")
  
  summ4 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX1720957/S4_slots_summary.tsv", sep = "\t")
  
  summ5 <- read.table("/tmp/Rtmp5p73sC/R.INSTALLdb773187b3a/BgeeCall/intergenic_1.0/all_results/SRX1720958/S4_slots_summary.tsv", sep = "\t")
}

# Combine the summary statistics for all libraries
summ_stats <- cbind(summ, summ1, summ2, summ3, summ4, summ5)
head(summ_stats)

# Plot the proportion of protein coding genes of all libraries for each p-value = 0.05
read.table(calls_output$cutoff_info_file_path, stringsAsFactors = FALSE, header = TRUE, blank.lines.skip = TRUE)
# Error in read.table(calls_output$cutoff_info_file_path, stringsAsFactors = FALSE,  : 
# 'file' must be a character string or connection

# Read each the cutoff_info for each library and subset 'proportionCodingPresent'

{
  cut_info <- read.table(calls_output[[1]][["cutoff_info_file_path"]]) #SRX109279
  #head(cut_info)
  sub_1 <- cut_info[6, ]
  #proportionCodingPresent   79.5889135572585
  
  cut_info1 <- read.table(calls_output[[2]][["cutoff_info_file_path"]]) #SRX109278
  #head(cut_info1)
  sub_2 <- cut_info1[6, ]
  #proportionCodingPresent 71.8327007090167
  
  cut_info2 <- read.table(calls_output[[3]][["cutoff_info_file_path"]]) #SRX1720957
  sub_3 <- cut_info2[6, ]
  
  cut_info3 <- read.table(calls_output[[4]][["cutoff_info_file_path"]]) #SRX493950
  sub_4 <- cut_info3[6, ]
  
  cut_info4 <- read.table(calls_output[[5]][["cutoff_info_file_path"]]) #SRX1720958
  sub_5 <- cut_info4[6, ]
  
  cut_info5 <- read.table(calls_output[[6]][["cutoff_info_file_path"]]) #SRX493999
  sub_6 <- cut_info5[6, ]
}

#create data frame with 0 rows and 2 columns
prot_cod <- data.frame(matrix(ncol = 2, nrow = 0))

# append 'proportionCodingPresent' values to dataframe
prot_cod <- rbind(sub_1, sub_2, sub_3, sub_4, sub_5, sub_6)

# Plot the proportion of protein coding genes
library(ggplot2)

dev.off()

ggplot(prot_cod,aes(x = V1 ,fill = V2)) + 
  geom_bar(position = "fill") + scale_fill_brewer(palette="Dark1")

