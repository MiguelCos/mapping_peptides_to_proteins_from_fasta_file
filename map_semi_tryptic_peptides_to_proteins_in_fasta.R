## Test script for the usability of the annotate_peptides function 
## Miguel Cosenza v1.4

# Load packages ----
library(tidyverse)
library(seqinr)
library(here)

# Load data ----

expr_tab <- read.delim(here("data/sample/sample_tabular.tsv"), 
                       stringsAsFactors = FALSE)

fasta <- read.fasta(here("data/sample/sample_sequences.fasta"),
  seqtype = "AA", as.string = TRUE)

## Wrangle data ----
## Correct header for expression matrix ----

head_names <- str_remove_all(names(expr_tab), pattern = ".*\\\\") %>% 
                    str_remove_all(string = ., pattern = "\\.raw.dia") 

head_names[1] <- "Peptide"

#print(head_names)

colnames(expr_tab) <- head_names


## Define variables that will go into the function  

fasta <- fasta

## Load the annotate_peptides functions into the R environment ----

source(file = here("annotate_peptides.R"))

#### Execution of the function -----

# Execute the function ----
annotated_peptides <- annotate_peptides(expr_mat = expr_tab, 
                                        fasta = fasta,
                                        decoy_tag = "^rev_")

# Join the original data table with the newly generated data table including the peptide annotation
peptides_w_annotation <- left_join(expr_tab,
  annotated_peptides, by = "Peptide")

write_delim(x = peptides_w_annotation,
            file = here("output/sample/tabular_file_w_peptide_annotation.tsv"),
            delim = "\t")
