#! R

library("DEP")
library("tidyverse")
library("limma")
library("SummarizedExperiment")

##source functions
#source("https://raw.githubusercontent.com/brucemoran/R/master/functions/expression/DEP/DEP_LFQ.func.R")
source("/Users/brucemoran/Google Drive/scripts/R/github/functions/expression/DEP/DEP_LFQ.func.R")

##read and parse XLSX input 
prot_xlsx <- dir(path = "proteomics", pattern = "xlsx", full.names = TRUE)
prot_tb <- readxl::read_xlsx(path = prot_xlsx)

prot_de_tb <- dplyr::mutate(.data = prot_tb, ProteinGenesid = paste0(Protein, "_", `Gene names`, "_", id) )%>%
              dplyr::select("ProteinGenesid", tidyselect::starts_with("Intensity "))
colnames(prot_de_tb) <- gsub("-", "_", gsub("Intensity ", "", colnames(prot_de_tb)))

##create sample map from colnames
samps <- grep("_", colnames(prot_de_tb), value = TRUE)
samp_map <- tibble::tibble(samples = samps,
                           cell_line = gsub("ph", "", unlist(lapply(samps, function(f){
                             stringi::stri_split(regex = "_", f)[[1]][1]
                             }))),
                           metformin = ifelse(unlist(lapply(samps, function(f){
                             stringi::stri_split(regex = "_", f)[[1]][2]
                             })) == "TR", "YES", "NO"),
                           tech_rep = unlist(lapply(samps, function(f){
                             paste(stringi::stri_split(regex = "_", f)[[1]][1:3], collapse = "_")
                             }))
                           )

sample_map_tb <- samp_map
sample_ID <- "samples"
group <- "cell_line"
prot_data_tb <- prot_de_tb
row_data <- "ProteinGenesid"
tech_reps <- "tech_rep"
exclude_group = NULL; convert_group = NULL; 

##
dir.create("DEP_output")
GLC_met <- run_expt_group(project_name = "GLC_metformin", 
                          sample_map_tb = samp_map,
                          sample_ID = "samples", 
                          group = c("cell_line", "metformin"), 
                          prot_data_tb = prot_de_tb, 
                          row_data = "ProteinGenesid",
                          tech_reps = "tech_rep",
                          exclude_group = "NLC", 
                          convert_group = NULL, 
                          p_adj_val = 0.1,
                          out_dir = "DEP_output")

GLC_NLC <- run_expt_group(project_name = "GLC_NLC", 
                          sample_map_tb = samp_map,
                          sample_ID = "samples", 
                          group = c("cell_line", "metformin"), 
                          prot_data_tb = prot_de_tb, 
                          row_data = "ProteinGenesid",
                          tech_reps = "tech_rep",
                          exclude_group = "YES", 
                          convert_group = NULL, 
                          p_adj_val = 0.1,
                          out_dir = "DEP_output")

##save that
readr::write_csv(prot_de_tb, file = "proteomics/parse_clean.csv")
