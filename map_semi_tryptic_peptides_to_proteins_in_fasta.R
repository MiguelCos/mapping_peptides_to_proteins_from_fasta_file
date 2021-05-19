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

peptides <- read.delim("data/sample/sample_tabular.tsv")

fasta <- read.fasta("data/sample/sample_sequences.fasta",
  seqtype = "AA", as.string = TRUE)

## Wrangle data ----

## Create a column 'peptide_seq' including the peptide sequence ----
#this is just of the sake of 'left_join' mergin later on (line 103)

peptides <- mutate(peptides,
                   peptide_seq = Stripped.Sequence)

## Load the annotate_peptides functions into the R environment ----

source(file = "annotate_peptides.R")

#### Execution of the function -----

# Extract a vector of peptide sequences
sequences <- peptides$peptide_seq

# Execute the function ----
annotated_peptides <- annotate_peptides(peptide_sequence = sequences, 
                                        fasta = fasta)

# Join the original data table with the newly generated data table including the peptide annotation
peptides_w_annotation <- left_join(peptides,
  annotated_peptides, by = "peptide_seq")

write_delim(x = peptides_w_annotation,
            file = "output/sample/tabular_file_w_peptide_annotation.tsv",
            delim = "\t")
