module FastaProcessing

using SeqShuffle, Random

export dna_data, 
       combined_dna_dataset, 
       fullstack_dna_dataset

# Write your package code here.
include("structs.jl")
include("helper.jl")

end
