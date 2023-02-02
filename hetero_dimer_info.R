### Please keep your data in ths format
#This file contains information about your heterodimer
ProteinA_ProteinB.info <- list(pdb    = "none", 
                    resolution  = 0.0,
                    pdb_local   = "none",
                    pdb_chain   = list(x = "ProteinA", y = "ProteinB"),
                    sequence_A  = "FASTA-data/ProteinA.fasta",
                    sequence_B  = "FASTA-data/ProteinB.fasta",
                    descrip_A   = "ProteinA",
                    descrip_B   = "ProteinB",
                    gene_names  = c("ProteinA","ProteinB"),
                    chain_A     = "ProteinA", 
                    chain_B     = "ProteinB", 
                    uniprot_A   = "none",
                    uniprot_B   = "none",
                    dca_signal  = "results/ProteinA_ProteinB/ProteinA_ProteinB.dca_scores",
                    cMSA        = "results/ProteinA_ProteinB/ProteinA_aligned_ProteinB_aligned_msa.fasta.nodups",
                    rand_dca    = "results/ProteinA_ProteinB/rand_dca/"
)

# ProteinA_ProteinC.info <- list(pdb    = "none", 
#                                resolution  = 0.0,
#                                pdb_local   = "none",
#                                pdb_chain   = list(x = "A", y = "B"),
#                                sequence_A  = "FASTA-data/ProteinA.fasta",
#                                sequence_B  = "FASTA-data/ProteinC.fasta",
#                                descrip_A   = "ProteinA",
#                                descrip_B   = "ProteinC",
#                                gene_names  = c("ProteinA","ProteinC"),
#                                chain_A     = "ProteinA", 
#                                chain_B     = "ProteinC", 
#                                uniprot_A   = "none",
#                                uniprot_B   = "none",
#                                dca_signal  = "results/ProteinA_ProteinC/ProteinA_ProteinC.dca_scores",
#                                cMSA        = "results/ProteinA_ProteinC/ProteinA_aligned_ProteinC_aligned_msa.fasta.nodups",
#                                rand_dca    = "results/ProteinA_ProteinC/rand_dca/"
# )