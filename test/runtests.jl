using FastaProcessing
using Test

@testset "FastaProcessing.jl" begin
    # Write your tests here.
    test_1 = "jaspar/MA1337.1.sites"
    test_2 = "jaspar/MA1349.1.sites"
    test_3 = "jaspar/MA1351.1.sites"

    trim_len=30
    dna_datasets = dna_data.([test_1, test_2, test_3]; trim=trim_len)
    @test all([i.len for i in dna_datasets] .== 115-2*trim_len)
    f_dna_data = fullstack_dna_dataset{Float16}(dna_datasets)
end
