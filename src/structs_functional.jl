
struct dna_data_labels{F}
    vals::Union{F, AbstractVector{F}}
    function dna_data_labels{F}(label_path::String) where {F <: Real}
        f = open(label_path)
        reads = read(f, String)
        close(f)
        vals = [parse(F,i) for i in split(reads, "\n")]
        new(vals)
    end
end


"""
multi_2_one_dna_dataset
    multiple datasets of dna seuqneces, e.g. 
        dataset_1: [s₁, s₂,...,sₙ]
        dataset_2: [t₁, t₂,...,tₙ]
    and one label set
        labels: [l₁, l₂,...,lₙ]

    note that all datasets and the label set must have the same cardinality
    
    This can happen in splicing datasets, where we only want to 
    model multiple (non-endogeneous) intronic regions and how 
    they affect gene-expression
"""
struct multi_2_one_dna_dataset{F}
    dna_datasets::Vector{dna_data}
    labels::dna_data_labels{F}
    function multi_2_one_dna_dataset{F}(
        dna_datasets::Vector{dna_data}, labels::dna_data_labels{F}) where F <: Real
        num_seq_each = [length(i.dna_reads) for i in dna_datasets]
        @assert length(unique(num_seq_each)) == 1 "all datasets must have the same number of sequences"
        @assert num_seq_each[1] == length(labels.vals) "number of labels and seqs must be the same"
        new(dna_datasets, labels)
    end
end

"""
multi_2_multi_dna_dataset
    multiple datasets of dna seuqneces, e.g. 
        dataset_1: [s₁, s₂,...,sₙ]
        dataset_2: [t₁, t₂,...,tₘ]
    and multiple label sets
        labels_1: [l₁, l₂,...,lₙ]
        labels_2: [l₁, l₂,...,lₘ]

    note that each dataset and its corresponding label set 
    must have the same cardinality
    
"""
struct multi_2_multi_dna_dataset{F}
    dna_datasets::Vector{dna_data}
    labels::Vector{dna_data_labels{F}}
    function multi_2_multi_dna_dataset{F}(
        dna_datasets::Vector{dna_data}, labels::Vector{dna_data_labels{F}}
    ) where F <: Real
        num_seq_each = [length(i.dna_reads) for i in dna_datasets]
        num_label_each = [length(i.vals) for i in labels]
        cond = all(length.(num_seq_each) .== length.(num_label_each))
        @assert cond "number of seqs and labels for each dataset must match"
        new(dna_datasets, labels)
    end
end



# TODO define train_test_ratio as constant here

struct functional_data{F}
    data::Union{multi_2_one_dna_dataset, multi_2_multi_dna_dataset}
    data_matrices_full::Vector{Array{F, 4}}
    train_set_inds::Vector{Int}
    test_set_inds::Vector{Int}
    data_matrices_shuffled_train::Vector{AbstractArray{F, 4}}
    labels_shuffled_train::Union{AbstractArray{F,1}, Vector{Vector{F}}}
    data_matrices_shuffled_test::Vector{AbstractArray{F, 4}}
    labels_shuffled_test::Union{AbstractArray{F,1}, Vector{Vector{F}}}
    function functional_data{F}(
        m2one::multi_2_one_dna_dataset
        ; train_test_ratio=0.9) where F <: Real
        data_matrices_full = [data_2_dummy(i.dna_reads) for i in m2one.dna_datasets];
        train_set_inds, test_set_inds = 
            get_shuffled_train_test_inds(m2one.dna_datasets[1].dna_reads; gamma=train_test_ratio)
        data_matrices_shuffled_train = [@view i[:, :, :, train_set_inds] for i in data_matrices_full]
        data_matrices_shuffled_test  = [@view i[:, :, :, test_set_inds] for i in data_matrices_full]
        labels_shuffled_train = @view m2one.labels.vals[train_set_inds]
        labels_shuffled_test = @view m2one.labels.vals[test_set_inds]
        new(m2one, data_matrices_full, train_set_inds, test_set_inds,
            data_matrices_shuffled_train, labels_shuffled_train,
            data_matrices_shuffled_test, labels_shuffled_test)
    end
    function functional_data{F}(
        m2m::multi_2_multi_dna_dataset
        ; train_test_ratio=0.9) where F <: Real
        data_matrices_full = [data_2_dummy(i.dna_reads) for i in m2m.dna_datasets];
        train_set_inds, test_set_inds = 
            get_shuffled_train_test_inds(m2m.dna_datasets[1].dna_reads; gamma=train_test_ratio)
        data_matrices_shuffled_train = [@view i[:, :, :, train_set_inds] for i in data_matrices_full]
        data_matrices_shuffled_test  = [@view i[:, :, :, test_set_inds] for i in data_matrices_full]
        labels_shuffled_train = [i.vals[train_set_inds] for i in m2m.labels]
        labels_shuffled_test = [i.vals[test_set_inds] for i in m2m.labels]
        new(m2m, data_matrices_full, train_set_inds, test_set_inds,
            data_matrices_shuffled_train, labels_shuffled_train,
            data_matrices_shuffled_test, labels_shuffled_test)
    end
    # directly from multiple data paths
    function functional_data{F}(
            dna_data_paths::Vector{String}, data_labels::String; 
            trim=0, train_test_ratio=0.9) where F <: Real
        dna_datasets = dna_data.(dna_data_paths, trim=trim)
        labels = dna_data_labels{F}(data_labels)
        m2one = multi_2_one_dna_dataset{F}(dna_datasets, labels)
        data_matrices_full = [data_2_dummy(i.dna_reads) for i in m2one.dna_datasets];
        train_set_inds, test_set_inds = 
            get_shuffled_train_test_inds(m2one.dna_datasets[1].dna_reads; gamma=train_test_ratio)
        data_matrices_shuffled_train = [@view i[:, :, :, train_set_inds] for i in data_matrices_full]
        data_matrices_shuffled_test  = [@view i[:, :, :, test_set_inds] for i in data_matrices_full]
        labels_shuffled_train = @view m2one.labels.vals[train_set_inds]
        labels_shuffled_test = @view m2one.labels.vals[test_set_inds]
        new(m2one, data_matrices_full, train_set_inds, test_set_inds,
            data_matrices_shuffled_train, labels_shuffled_train,
            data_matrices_shuffled_test, labels_shuffled_test)
    end
    # directly from multiple data paths and multiple labels
    function functional_data{F}(
                dna_data_paths::Vector{String}, data_labels::Vector{String}; 
                trim=0, train_test_ratio=0.9) where F <: Real
            cond = length(dna_data_paths) == length(data_labels)
            @assert cond "length of data_paths and lables must be the same"
            dna_datasets = dna_data.(dna_data_paths, trim=trim)
            labels = dna_data_labels{Float16}.(data_labels)
            m2m = multi_2_multi_dna_dataset{Float16}(dna_datasets, labels)
            data_matrices_full = [data_2_dummy(i.dna_reads) for i in m2m.dna_datasets];
            train_set_inds, test_set_inds = 
                get_shuffled_train_test_inds(m2m.dna_datasets[1].dna_reads; gamma=train_test_ratio)
            data_matrices_shuffled_train = [@view i[:, :, :, train_set_inds] for i in data_matrices_full]
            data_matrices_shuffled_test  = [@view i[:, :, :, test_set_inds] for i in data_matrices_full]
            labels_shuffled_train = [i.vals[train_set_inds] for i in m2m.labels]
            labels_shuffled_test = [i.vals[test_set_inds] for i in m2m.labels]
            new(m2m, data_matrices_full, train_set_inds, test_set_inds,
                data_matrices_shuffled_train, labels_shuffled_train,
                data_matrices_shuffled_test, labels_shuffled_test)
    end
end

# TODO save the data
# only need to save the filepath, and the shuffle inds

