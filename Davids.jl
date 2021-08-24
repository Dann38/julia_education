using DelimitedFiles
DITECTORY = joinpath(pwd(), "tmp_example")
NAME = "matrix2"
file_trace = joinpath(DITECTORY, NAME*"_trace.txt")
file_index_i = joinpath(DITECTORY, NAME*"_i.txt")
file_index_j = joinpath(DITECTORY, NAME*"_j.txt")
file_element_ij = joinpath(DITECTORY, NAME*"_a.txt")
trace = DelimitedFiles.readdlm(file_trace, ' ', Float64)
indexs_i = DelimitedFiles.readdlm(file_index_i, ' ', Int64)
indexs_j = DelimitedFiles.readdlm(file_index_j, ' ', Int64)
elements_ij = DelimitedFiles.readdlm(file_element_ij, ' ', Float64)

N = size(trace)[2]
num_nodiag_el = size(elements_ij)[2]

H = zeros(Float64, N, N)
for i in range(1,N)
    H[i, i] = trace[i]
end

for k in range(1, num_nodiag_el)
    H[indexs_i[k], indexs_j[k]] = elements_ij[k]
    H[indexs_j[k], indexs_i[k]] = elements_ij[k]
end

H

rayleigh_quotient(x) = x' * H * x / (x' * x)

function step_davids(i::Int64, k::Int64, phi)
    rayl = rayleigh_quotient(phi)
    a_old = 0.0
    for j in range(1, k)
        phi_now = phi
        phi_now[i] = phi_now[i]+a_old
        rayl = rayleigh_quotient(phi_now)
        a_now = (H[i,:]'*phi-rayl*phi[i])/(rayl-H[1,1])
        a_old = a_now
    end
    return a_old
end

phi = zeros(N);

for t in range(1, N)
    phi[t] = eps()
end
phi[1] = 1.0+eps();
a =zeros(N);
for j in range(1, N)
   a[j] = step_davids(j, 100, phi)
end
phi_new = phi + a
lambda = rayleigh_quotient(phi_new)
println("max  Eigenvalues = ",lambda)







