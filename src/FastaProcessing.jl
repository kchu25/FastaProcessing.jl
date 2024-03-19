module FastaProcessing

using SeqShuffle, Random

export dna_data, 
       combined_dna_dataset, 
       fullstack_dna_dataset,
       make_cross_linked_DNA_matrix

# Write your package code here.
include("structs.jl")
include("helper.jl")
include("crossinglinking.jl")

end
