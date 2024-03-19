function obtain_dna_heads_and_dna_reads(fastapath::String)
    f = open(fastapath)
    reads = read(f, String)
    close(f)
    dna_heads = Vector{String}();
    dna_reads = Vector{String}();
    for i in split(reads, ">")
        if !isempty(i)
            splits = split(i, "\n");
            this_read_head= splits[1];
            this_read = join(splits[2:end]);
            push!(dna_heads, this_read_head);
            push!(dna_reads, this_read);
        end
    end
    return dna_heads, dna_reads
end

"""
trim the dna reads on both sides by n bases
"""
function trim_dna_reads!(dna_reads::Vector{String}, dna_heads::Vector{String}, n::Int)
    keep = falses(length(dna_reads))
    @assert unique(length.(dna_reads)) |> length == 1 "All reads should have the same length"
    len = length(dna_reads[1]) - 2*n
    @assert n < length(dna_reads[1]) ÷ 2 "n should be less than the half the length of the reads"
    for i in 1:length(dna_reads)
        shortened_read = dna_reads[i][n+1:end-n]
        # so only get rid of that shortened-read if it has an n or N in it
        if !occursin("N", shortened_read) && !occursin("n", shortened_read)
            keep[i] = true
            dna_reads[i] =shortened_read
        end
    end
    dna_heads = @view dna_heads[keep]
    dna_reads = @view dna_reads[keep]
    return dna_heads, dna_reads, len
end

struct dna_data
    dna_heads::AbstractVector{String}
    dna_reads::AbstractVector{String}
    len::Int
    function dna_data(fastapath::String; trim::Int=10)
        dna_heads, dna_reads = obtain_dna_heads_and_dna_reads(fastapath)
        dna_heads, dna_reads, len = trim_dna_reads!(dna_reads, dna_heads, trim)
        new(dna_heads, dna_reads, len)
    end 
end

"""
combined the dna reads from all datasets
following the order of the input dna_datasets
useful for 
    1. referencing the different crosslinking signitures (read_range)
    2. generate a datamatrix that includes all of the sequences
"""
struct combined_dna_dataset{F}
    dna_datasets::Vector{dna_data}
    # data_matrix::Array{F, 3}
    read_range::Dict{UnitRange{Int}, Int}
    L::Int
    function combined_dna_dataset{F}(dna_datasets::Vector{dna_data}) where {F <: Real}
        lens = [dna_datasets[i].len for i in 1:length(dna_datasets)]
        @assert length(unique(lens)) == 1 "All datasets trimmed reads should have the same length"
        read_range = get_range_dataset(dna_datasets)
        new(dna_datasets, read_range, dna_datasets[1].len)
    end
end

#=
data_matrix: full dataset of dna reads
data_matrix_bg: full dataset of shuffled dna reads
data_matrix_shuffled: training set of dna reads (using shuffled indices of data_matrix)
data_matrix_shuffled_bg: training set of shuffled dna reads (using shuffled indices of data_matrix_bg)
data_matrix_shuffled_test: test set of dna reads (using shuffled indices of data_matrix)
data_matrix_shuffled_test_bg: test set of shuffled dna reads (using shuffled indices of data_matrix_bg)
=#
function make_data_matrices(cdna_datasets; train_test_ratio=0.9, k_freq=1)
    cat_reads = get_all_reads(cdna_datasets)
    train_set_inds, test_set_inds = 
        get_shuffled_train_test_inds(cat_reads; gamma=train_test_ratio)
    data_matrix = data_2_dummy(cat_reads);
    data_matrix_bg = data_2_dummy(seq_shuffle.(cat_reads, k=k_freq))
    data_matrix_shuffled_train      = @view data_matrix[:, :, train_set_inds]
    data_matrix_shuffled_train_bg   = @view data_matrix_bg[:, :, train_set_inds]
    data_matrix_shuffled_test       = @view data_matrix[:, :, test_set_inds]
    data_matrix_shuffled_test_bg    = @view data_matrix_bg[:, :, test_set_inds]
    return data_matrix, data_matrix_bg, 
           train_set_inds, test_set_inds,
           data_matrix_shuffled_train, data_matrix_shuffled_train_bg, 
           data_matrix_shuffled_test, data_matrix_shuffled_test_bg
end

struct fullstack_dna_dataset{F}
    cdna_datasets::combined_dna_dataset{F}
    data_matrix_full::Array{F, 3}
    data_matrix_bg::Array{F, 3}
    train_set_inds::Vector{Int}
    test_set_inds::Vector{Int}
    data_matrix_shuffled_train::AbstractArray{F, 3}
    data_matrix_shuffled_train_bg::AbstractArray{F, 3}
    data_matrix_shuffled_test::AbstractArray{F, 3}
    data_matrix_shuffled_test_bg::AbstractArray{F, 3}
    function fullstack_dna_dataset{F}(
        dna_datasets::Vector{dna_data},
        train_test_ratio=0.9, kmer_bg_freq=1
        ) where {F <: Real}
        cdna_datasets = combined_dna_dataset{F}(dna_datasets)
        data_matrix, data_matrix_bg, 
        train_set_inds, test_set_inds,
        data_matrix_shuffled_train, data_matrix_shuffled_train_bg, 
        data_matrix_shuffled_test, data_matrix_shuffled_test_bg = 
            make_data_matrices(cdna_datasets; 
                train_test_ratio=train_test_ratio, k_freq=kmer_bg_freq)
        new(cdna_datasets, 
            data_matrix, data_matrix_bg,
            train_set_inds, test_set_inds,
            data_matrix_shuffled_train, data_matrix_shuffled_train_bg, 
            data_matrix_shuffled_test, data_matrix_shuffled_test_bg)
    end
end