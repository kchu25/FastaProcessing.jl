
function make_cross_linked_DNA_matrix(data_matrix; crosslink_on=false)
    L, _, N = size(data_matrix)
    data_matrix_reshaped = reshape(data_matrix, 4, L÷4, N);
    crosslinked_position = size(data_matrix_reshaped,2)÷2+1
    # create a one hot vector and concatenate it as the fifth row 
    # of each data_matrix wrt the first dimension
    one_hot = zeros(eltype(data_matrix), (5, L÷4, N));
    crosslink_on && (@inbounds one_hot[5, crosslinked_position, :] .= 1;)
    @inbounds one_hot[1:4,:,:] .= data_matrix_reshaped;
    return reshape(one_hot, (5*(L÷4), 1, N))
end