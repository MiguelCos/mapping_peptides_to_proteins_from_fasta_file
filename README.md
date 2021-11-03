# Annotate peptides from proteins sequences in fasta file  

## Updates on version 1.0:

The function `annotate_peptides` now requires a data frame that contain the mapping between peptide sequences to their corresponding Uniprot protein ID. This was necessary because searching the fasta file by matching the peptide sequence of interest against all the protein sequences in the fasta file was too time intensive. An older version of this `annotate_peptides` function allows to match the peptide sequence against all protein sequences in the fasta file. 

## Updates on version 1.1:

The `annotate_peptides` function now forces the input columns for the annotation be characters. In older versions of R, `read.delim` loads strings as factors, making this input incompatible with `str_*` functions within `annotate_peptides`.

## Updates on version 1.2:

The `annotate_peptides` function now creates two additional columns: `semi_type` and `specificity`. Look function description for details.

## Update on version 1.3:

The `annotate_peptides` function now creates an additional column: `previous_all_resid`.  

## Update on version 1.4: 

The `annotate_peptides` function now creates an additional column: `protein_length`.  

## Update on version 1.5: 

The `annotate_peptides` function now includes the argument `specificity`. With that, the user can select which protease you are working with. Now the labelling changed from `semi_tryptic` to `semi_specific` to make it general in case the user is using a protease different than trypsin.

## The function `annotate_peptides`  

The script `annotate_peptides.R` contains the function `annotate_peptides` that takes three argumens:

- `expr_mat`: a data frame containing __at least__ two columns with the next _exact_ names: `Peptide` (containing the peptide sequence that you want to annotatate) `Genes` (containing the Uniprot ID that can be used to look for the interesting protein sequence in the fasta file).
- `fasta`: list containing the fasta file that was used by the search engine to identify the peptides in `peptide_sequence`. The list of fasta sequences should be loaded into R with the `read.fasta` function from the `seqinr` package.  
- `decoy_tag`: a string defining how to differentiate the decoy sequences in the fasta sequences file for their exclusion. 
- `specificity` = a string defining the enzyme specificity of your experiment it is set to `"R|K"` by default for trypsin. You can set it to `"R"` for Argc specificity, for example.

The output would be a `data.frame` containing the columns:

- `protein_id`
- `protein_description`
- `protein_length`
- `peptide_seq`: sequence of the input peptide.
- `start_position`: numeric start position of the peptide in the protein sequence.
- `stop_position`: numeric stop position of the peptide in the protein sequence.
- `last_aa`: last amino acid.
- `aa_after`: amino acid after input peptide
- `aa_before`: amino acid before input peptide
- `following_10_resid`: following 10 amino acids after input peptide
- `previous_10_resid`: previous 10 amino acids after input peptide
- `previous_all_resid`: all the previous residues up to the protein N-termini as annotated in the fasta file. This can be regarded as an 'artificial N-termini'. When the non-tryptic cleavage of the peptide is at the C-termini, this artificial one can be used as a proxy to get potential N-terminal peptides arousing from protease activity. 
- `specificity`: tryptic or semi-tryptic peptide.
- `semi_type`: if semi-tryptic; C-term or N-term semi-tryptic.

This can be used to join with your tabular data containing quantitative information.

__Note on the older version of this function__: each mapping between a peptide sequence and the correspondent protein sequence in the fasta file can take ~5-30 seconds depending on the size of the fasta file so the execution can take a while depending on how many peptide sequences you have as an input and the size of the fasta file. The function has a counter to keep track of the execution.  

Our test with the new version show that it takes ~10 minutes of execution to annotate ~27800 peptides against a fasta file with ~20K sequences.

## The script `map_semi_tryptic_peptides_to_proteins_in_fasta.R`  

This script contains an example usage of the `annotate_peptides` function using an abundance matrix with two columns `Peptide` and `Genes`. This can

Small sample files, both for the quantitative matrix and the fasta file are provided in the `data/sample` folder to test the usage of this function and script.

The input file is an quantitative matrix of peptides in which the first two columns correspond to peptide sequence and protein ID.

In line `46-49`, you can see how you can join this newly generated information with your general quantitative matrix.


