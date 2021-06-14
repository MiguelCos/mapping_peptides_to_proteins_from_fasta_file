# Annotate peptides function ----
# Miguel Cosenza v 1.0 

annotate_peptides <- function(expr_mat, fasta,
                              decoy_tag = "^rev"){ # input should be a vector of peptide sequences and the fasta file
   
   pept2prot <- expr_mat %>% 
      select(Peptide, Genes)
   
   peptide_sequence <- pept2prot$Peptide
   
   
   if (!is.null(decoy_tag)){
      fasta <- fasta %>% 
         purrr::discard(str_detect(names(.), decoy_tag))
   }
   
   
   annotation <- data.frame()
   
   for (i in 1:length(peptide_sequence)){
      
      # extract the list element that matches the peptide sequence
      list_elem <- fasta %>%
         purrr::keep(str_detect(names(.), pept2prot$Genes[i]))
      
      
      if(is_empty(list_elem)){
         protein_id <- NA
         protein_description <- NA
         peptide_seq <- peptide_sequence[i]
         peptide_position <- NA
         start_position <- NA
         end_position <- NA
         following_10_resid <- NA
         previous_10_resid <- NA
         aa_after <- NA
         aa_before <- NA
         last_aa <- NA
         
         temprow <- cbind(protein_id,
                          protein_description,
                          peptide_seq, 
                          start_position,
                          end_position, 
                          last_aa, aa_after,
                          aa_before, 
                          following_10_resid, 
                          previous_10_resid)
         
         annotation <- rbind(annotation,
                             temprow)
         
         row.names(annotation) <- NULL
         
      } else {
         
         # extract the annotation info stored in the fasta header, check for an Uniprot fasta pattern and extract protein ID and description
         if(str_detect(attr(list_elem[[1]], "Annot"), 
                       pattern = "^\\>sp\\|") | str_detect(attr(list_elem[[1]], "Annot"), 
                                                           pattern = "^\\>tr\\|")){
            split1 <- str_split_fixed(attr(list_elem[[1]], "Annot"), 
                                      pattern = "\\|", 
                                      n = 3)
            protein_id <- split1[,2] # get protein id for 'uniprot' IDs
            split2 <- str_split_fixed(split1[,3], 
                                      pattern = "\\s", 
                                      n = 2)
            protein_description <- str_extract(split2[,2], 
                                               pattern = ".+?(?=OS=)") %>% 
               str_trim() # get protein description for non-uniprot IDs
         } else {
            split1 <- str_split_fixed(attr(list_elem[[1]], 
                                           "Annot"), 
                                      pattern = "\\s", 
                                      n = 2)
            
            protein_id  <- split1[,1] # get protein ID
            
            protein_description <- split1[,2] # get protein description for non Uniprot IDs
         }
         
         # extract protein sequence from the extracted list element
         
         protein_seq <- list_elem[[1]][1] # get the matched protein sequence
         
         peptide_seq <- peptide_sequence[i] # get the sequence of the peptide
         
         peptide_position <- str_locate(protein_seq, 
                                        pattern = peptide_seq) # get the position of the peptide within the matched protein sequence
         
         start_position <- peptide_position[,1] # start position of the peptide within the protein
         
         end_position <- peptide_position[,2] # end positionj of the peptide within the protein
         
         following_10_resid <- str_sub(protein_seq, 
                                       start = peptide_position[,2]+1, 
                                       end = peptide_position[,2]+10) # 10 residues after the end of the peptide sequence
         
         previous_10_resid <- str_sub(protein_seq, 
                                      start = peptide_position[,1]-10, 
                                      end = peptide_position[,1]-1) # 10 residues before the start of the peptide sequence
         
         aa_after <- str_sub(protein_seq, 
                             start = peptide_position[,2]+1, 
                             end = peptide_position[,2]+1) # amino acid after
         
         aa_before <- str_sub(protein_seq, 
                              start = peptide_position[,1]-1, 
                              end = peptide_position[,1]-1) # amino acid before
         
         last_aa <- str_sub(peptide_seq, 
                            start = str_count(peptide_seq), 
                            end = str_count(peptide_seq)) # last amino acid of the identified peptide
         
         temprow <- cbind(protein_id,
                          protein_description,
                          peptide_seq, 
                          start_position,
                          end_position, 
                          last_aa, 
                          aa_after,
                          aa_before, 
                          following_10_resid, 
                          previous_10_resid)
         
         annotation <- rbind(annotation,
                             temprow)
         
         row.names(annotation) <- NULL
         
         print(paste(i, peptide_sequence[i], "out of", length(peptide_sequence)))
      }
   }
   
   return(annotation)
}