function dna2dummy(dna_string::String, dummy::Dict; F=Float32)
    v = Array{F,2}(undef, (4, length(dna_string)));
    @inbounds for (index, alphabet) in enumerate(dna_string)
        # start = (index-1)*4+1;
        # v[start:start+3] = dummy[uppercase(alphabet)];
        v[:,index] = dummy[uppercase(alphabet)];
    end
    return v
end

function data_2_dummy(dna_strings; F = Float32) 
    dummy = Dict('A'=>Array{F}([1, 0, 0, 0]), 
                 'C'=>Array{F}([0, 1, 0, 0]),
                 'G'=>Array{F}([0, 0, 1, 0]), 
                 'T'=>Array{F}([0, 0, 0, 1]));
    how_many_strings = length(dna_strings);
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    _len_ = length(dna_strings[1]); # length of each dna string in data    
    _S_ = Array{F, 4}(undef, (4, _len_, 1, how_many_strings));
    @inbounds for i = 1:how_many_strings
        @assert length(dna_strings[i]) == _len_ "DNA strings must be same length in the dataset"
        _S_[:, :, 1, i] = dna2dummy(dna_strings[i], dummy; F=F)
    end
    return _S_
end

"""
example:
dataset A has 500 reads
dataset B has 1000 reads
return Dict(1=>1:500, 2=>501:1500)
"""
function get_range_dataset(dna_datasets)
    cum_len = [length(dna_datasets[i].dna_reads) 
                    for i in eachindex(dna_datasets)] |> cumsum
    read_range = Dict(1:cum_len[1]=>1)
    for i in 2:length(dna_datasets)
        read_range[cum_len[i-1]+1:cum_len[i]] = i
    end
    return read_range
end



"""
    get_all_reads
    
"""
function get_all_reads(cdna_datasets::combined_dna_dataset) 
    cat_reads = reduce(vcat, 
        [i.dna_reads for i in cdna_datasets.dna_datasets])
    return cat_reads
end

function get_shuffled_train_test_inds(cat_reads; gamma=0.9)
    how_many_reads = length(cat_reads)
    how_may_entries_in_test = 
        Int(floor((1-gamma) * how_many_reads));
    shuffled_inds = randperm(how_many_reads)
    test_set_inds  = @view shuffled_inds[1:how_may_entries_in_test]
    train_set_inds = @view shuffled_inds[how_may_entries_in_test+1:end]
    return train_set_inds, test_set_inds
end
