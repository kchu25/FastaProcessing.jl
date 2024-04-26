module FastaProcessing

using SeqShuffle, Random

export dna_data, 
       dna_data_labels,
       combined_dna_dataset, 
       multi_2_one_dna_dataset,       
       multi_2_multi_dna_dataset,
       functional_data,
       fullstack_dna_dataset,
       make_cross_linked_DNA_matrix

# Write your package code here.
include("structs_motifs.jl")
include("structs_functional.jl")
include("helper.jl")
include("crossinglinking.jl")

end
