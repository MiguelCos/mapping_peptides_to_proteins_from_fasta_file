# Annotate peptides from proteins sequences in fasta file  

## The function `annotate_peptides.R`  

The script `annotate_peptides` contains the function `annotate_peptides` that takes two argumens:

- `peptide_sequence`: character vector containing peptide sequences.
- `fasta`: list containing the fasta file that was used by the search engine to identify the peptides in `peptide_sequence`. The list of fasta sequences should be loaded into R with the `read.fasta` function from the `seqinr` package.  

The output would be a `data.frame` containing the columns:

- protein_id
- protein_description
- peptide_seq: sequence of the input peptide.
- start_position: numeric start position of the peptide in the protein sequence.
- stop_position: numeric stop position of the peptide in the protein sequence.
- last_aa: last amino acid.
- aa_after: amino acid after input peptide
- aa_before: amino acid before input peptide
- following_10_resid: following 10 amino acids after input peptide
- previous_10_resid: previous 10 amino acids after input peptide

This can be used to join with your tabular data containing quantitative information.

__Note__: each mapping between a peptide sequence and the correspondent protein sequence in the fasta file can take ~5-30 seconds depending on the size of the fasta file so the execution can take a while depending on how many peptide sequences you have as an input and the size of the fasta file. The function has a counter to keep track of the execution.  

## The script `map_semi_tryptic_peptides_to_proteins_in_fasta.R`  

This script contains an example usage of the `annotate_peptides` function using a DIA-NN peptide tabular file as an input.


