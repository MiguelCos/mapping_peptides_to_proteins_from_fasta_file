# Annotate peptides from argc specificity function ----
# Miguel Cosenza v 0.1

annotate_peptides_argc <- function(expr_mat, fasta,
                              decoy_tag = "^rev"){ # input should be a vector of peptide sequences and the fasta file
   
   pept2prot <- expr_mat %>% 
      dplyr::select(Peptide, Genes) %>% 
      mutate(Peptide = as.character(Peptide),
             Genes = as.character(Genes))
   
   peptide_sequence <- pept2prot$Peptide
   
   
   if (!is.null(decoy_tag)){
      fasta <- fasta %>% 
         purrr::discard(str_detect(names(.), decoy_tag))
   }
   
   
   annotation <- data.frame()
   
   for (i in 1:length(peptide_sequence)){
      
      # extract the list element that matches the peptide sequence
      list_elem <- fasta %>%
         purrr::keep(str_detect(names(.), fixed(pept2prot$Genes[i])))
      
      
      if(is_empty(list_elem)){
         protein_id <- NA
         protein_description <- NA
         protein_length <- NA
         Peptide <- peptide_sequence[i]
         peptide_position <- NA
         start_position <- NA
         end_position <- NA
         following_10_resid <- NA
         previous_10_resid <- NA
         previous_all_resid <- NA
         aa_after <- NA
         aa_before <- NA
         last_aa <- NA
         
         temprow <- cbind(protein_id,
                          protein_description,
                          protein_length,
                          Peptide, 
                          start_position,
                          end_position, 
                          last_aa, aa_after,
                          aa_before, 
                          following_10_resid, 
                          previous_10_resid,
                          previous_all_resid)
         
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
         
         protein_length <- str_length(protein_seq) # get the matched protein length  
         
         Peptide <- peptide_sequence[i] # get the sequence of the peptide
         
         peptide_position <- str_locate(protein_seq, 
                                        pattern = Peptide) # get the position of the peptide within the matched protein sequence
         
         start_position <- peptide_position[,1] # start position of the peptide within the protein
         
         end_position <- peptide_position[,2] # end positionj of the peptide within the protein
         
         following_10_resid <- str_sub(protein_seq, 
                                       start = peptide_position[,2]+1, 
                                       end = peptide_position[,2]+10) # 10 residues after the end of the peptide sequence
         
         previous_10_resid <- str_sub(protein_seq,
                                      start = ifelse(peptide_position[,1]-10 > 0, peptide_position[,1]-10, 1),
                                      end = ifelse(peptide_position[,1]-1 > 0, peptide_position[,1]-1, 0)) # 10 residues before the start of the peptide sequence
         
         
         previous_all_resid <- str_sub(protein_seq, 
                                       start = 2, 
                                       end = peptide_position[,1]-1) # all residues before up to annotated protein start
         
         aa_after <- str_sub(protein_seq, 
                             start = peptide_position[,2]+1, 
                             end = peptide_position[,2]+1) # amino acid after
         
         aa_before <- str_sub(protein_seq, 
                              start = peptide_position[,1]-1, 
                              end = peptide_position[,1]-1) # amino acid before
         
         last_aa <- str_sub(Peptide, 
                            start = str_count(Peptide), 
                            end = str_count(Peptide)) # last amino acid of the identified peptide
         
         temprow <- cbind(protein_id,
                          protein_description,
                          protein_length,
                          Peptide, 
                          start_position,
                          end_position, 
                          last_aa, 
                          aa_after,
                          aa_before, 
                          following_10_resid, 
                          previous_10_resid,
                          previous_all_resid)
         
         annotation <- rbind(annotation,
                             temprow)
         
         row.names(annotation) <- NULL
         
         print(paste(i, peptide_sequence[i], "out of", length(peptide_sequence)))
      }
   }
   
   annotation <-  annotation %>%
      mutate(semi_type = case_when(str_detect(last_aa, "R") & str_detect(aa_before, "R", negate = TRUE) ~ "semi_Nterm",
                                   str_detect(last_aa, "R", negate = TRUE) & str_detect(aa_before, "R") ~ "semi_Cterm",
                                   str_detect(last_aa, "R", negate = TRUE) & str_detect(aa_before, "R", negate = TRUE) ~ "unspecific",
                                   TRUE ~ "specific")) %>%
      mutate(specificity = ifelse(semi_type == "specific",
                                  yes = "specific",
                                  no = "semi_specific"))
   
   return(annotation)
}
