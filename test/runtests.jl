using FastaProcessing
using Test

@testset "FastaProcessing.jl" begin
    # Write your tests here.
    test_1 = "jaspar/MA1337.1.sites"
    test_2 = "jaspar/MA1349.1.sites"
    test_3 = "jaspar/MA1351.1.sites"

    trim_len=30

    dna_datasets = dna_data.([test_1]; trim=trim_len)
    f_dna_data = fullstack_dna_dataset{Float16}(dna_datasets)

    dna_datasets = dna_data.([test_1, test_2, test_3]; trim=trim_len)
    @test all([i.len for i in dna_datasets] .== 115-2*trim_len)
    f_dna_data = fullstack_dna_dataset{Float16}(dna_datasets)

    test_f1 = "jaspar/test_1_dna.fa"
    test_f2 = "jaspar/test_2_dna.fa"
    test_f1l = "jaspar/test_1_labels.txt"
    test_f2l = "jaspar/test_2_labels.txt"
    dna_datasets = dna_data.([test_f1, test_f2], trim=30)
    labels = dna_data_labels{Float16}(test_f1l)
    m2one = multi_2_one_dna_dataset{Float16}(dna_datasets, labels)
    m2one_f = functional_data{Float16}(m2one)
    # directly from paths (strings)
    m2one_f = functional_data{Float16}([test_f1, test_f2], test_f1l)

    labels = dna_data_labels{Float16}.([test_f1l, test_f2l])
    m2m = multi_2_multi_dna_dataset{Float16}(dna_datasets, labels)
    m2m_f = functional_data{Float16}(m2m)
    m2m_f = functional_data{Float16}([test_f1, test_f2],[test_f1l, test_f2l])


end
