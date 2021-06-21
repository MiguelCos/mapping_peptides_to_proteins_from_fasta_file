#set working directory to current folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Preceding and following residues extraction from fasta file ----
## Miguel Cosenza v0.3

## This script would: ----

# Add columns to the input tabular file from DIA-NN containing the following infos:
# (similar to MQ peptides.txt) to be extracted from a fasta file
# - following 10 residues
# - preceding 10 residues
# - protein name; uniprot ID

# Load packages ----
library(tidyverse)
library(seqinr)

# Load data ----

expr_tab <- read.delim("data/sample/sample_tabular.tsv", stringsAsFactors = FALSE)

fasta <- read.fasta("data/sample/sample_sequences.fasta",
  seqtype = "AA", as.string = TRUE)

## Wrangle data ----
## Correct header for expression matrix ----

head_names <- str_remove_all(names(expr_tab), pattern = ".*\\\\") %>% 
                    str_remove_all(string = ., pattern = "\\.raw.dia") 

head_names[1] <- "Peptide"

print(head_names)

colnames(expr_tab) <- head_names


## Define variales that will go into the function  

fasta <- fasta

## Load the annotate_peptides functions into the R environment ----

source(file = "annotate_peptides.R")

#### Execution of the function -----

# Execute the function ----
annotated_peptides <- annotate_peptides(expr_mat = expr_tab, 
                                        fasta = fasta,
                                        decoy_tag = "^rev_")

# Join the original data table with the newly generated data table including the peptide annotation
peptides_w_annotation <- left_join(expr_tab,
  annotated_peptides, by = "Peptide")

write_delim(x = peptides_w_annotation,
            file = "output/sample/tabular_file_w_peptide_annotation.tsv",
            delim = "\t")
