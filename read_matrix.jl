using DelimitedFiles

# ============ Setting =================================
DIRECTORY = joinpath(pwd(), "tmp_example");
NAME = "matrix2";

struct MTRX
    trace::Matrix{Float64}
    ind_i::Matrix{Int64}
    ind_j::Matrix{Int64}
    el_ij::Matrix{Float64}   
end

function get_mtrx(directory::String, name::String)
    file_trace = joinpath(directory, name*"_trace.txt")
    file_index_i = joinpath(directory, name*"_i.txt")
    file_index_j = joinpath(directory, name*"_j.txt")
    file_element_ij = joinpath(directory, name*"_a.txt")

    trace = DelimitedFiles.readdlm(file_trace, ' ', Float64)
    indexs_i = DelimitedFiles.readdlm(file_index_i, ' ', Int64)
    indexs_j = DelimitedFiles.readdlm(file_index_j, ' ', Int64)
    elements_ij = DelimitedFiles.readdlm(file_element_ij, ' ', Float64)
    
    mtrx = MTRX(trace, indexs_i, indexs_j, elements_ij)
    return mtrx
end

function mul(mtrx::MTRX, vec::Array{Float64})
    n = length(vec)
    iter = length(mtrx.el_ij)
    rez = Array{Float64}(undef, n, 1)
    for k in range(1, n)
        rez[k] = vec[k] * mtrx.trace[k]
    end
    for k in range(1, iter)  
        rez[mtrx.ind_i[k]] += mtrx.el_ij[k] * vec[mtrx.ind_j[k]]
        rez[mtrx.ind_j[k]] += mtrx.el_ij[k] * vec[mtrx.ind_i[k]]
    end
    return rez
end

mtrx = get_mtrx(DIRECTORY, NAME)
x = [1., 1., 1., 1., 1., 1.1, 1.1]
println(mul(mtrx, x))


