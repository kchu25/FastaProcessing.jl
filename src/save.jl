"""
the particular data that we save for future reference
    
"""

struct save_data
    paths::Vector{String}
    train_set_inds::Vector{Int}
    test_set_inds::Vector{Int}
end