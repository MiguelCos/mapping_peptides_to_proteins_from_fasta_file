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

## Define variables that will go into the function  

fasta <- fasta

## Load the annotate_peptides functions into the R environment ----

source(file = here("annotate_peptides.R"))

#### Execution of the function -----

# Execute the function ----
annotated_peptides <- annotate_peptides(expr_mat = expr_tab, 
                                        fasta = fasta,
                                        decoy_tag = "^rev_",
                                        specificity = "R|K")


# Save the annotation table ----

write_delim(x = annotated_peptides,
            file = here("output/sample/annotated_peptides.tsv"),
            delim = "\t")

# Join the original data table with the newly generated data table including the peptide annotation
peptides_w_annotation <- left_join(expr_tab,
  annotated_peptides, by = "Peptide")

# Save the abundance matrix joined with N-ter annotation table
write_delim(x = peptides_w_annotation,
            file = here("output/sample/tabular_file_w_peptide_annotation.tsv"),
            delim = "\t")
