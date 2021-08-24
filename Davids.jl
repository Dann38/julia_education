using DelimitedFiles

# ============ Setting =================================
DIRECTORY = joinpath(pwd(), "tmp_example")
NAME = "matrix2"


# ============Function =================================
function get_matrix(directory::String, name::String)
    file_trace = joinpath(directory, name*"_trace.txt")
    file_index_i = joinpath(directory, name*"_i.txt")
    file_index_j = joinpath(directory, name*"_j.txt")
    file_element_ij = joinpath(directory, name*"_a.txt")

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

    return H
end

rayleigh_quotient(x) = x' * H * x / (x' * x)

function step_davids(i::Int64, k::Int64, phi)
    """
    ...
    # Arguments
    - `i::Int64`: Координата по которой идёт приближение.
    - `k::Int64`: Кол-во итераций
    - `phi` : Вектор первого приближения
    ...
    # Return
    - Вернет шаг по i-й координате для приближения к собственному вектору
    """
    rayl = rayleigh_quotient(phi)
    a_old = 0.0
    for j in range(1, k)
        phi_new = phi
        phi_new[i] = phi_new[i]+a_old
        rayl = rayleigh_quotient(phi_new)
        #       (H*phi)_i - rho(phi+a*e_i)
        # a  = ----------------------------
        #        rho(phi+a*e_i) - H_ii

        a_new = (H[i,:]'*phi-rayl*phi[i])/(rayl-H[1,1])
        a_old = a_new
    end
    return a_old
end


# ============Program =================================
H = get_matrix(DIRECTORY, NAME)
N = size(H)[1]

# Приближение
# phi = [1.0..01,   0.0..01,   ...,   0.0..01]
phi = zeros(N);
for t in range(1, N) phi[t] = eps() end
phi[1] = 1.0+eps();

# Считаем вектор a для шага
a =zeros(N);
for j in range(1, N)
   a[j] = step_davids(j, 100, phi)
end

# Делаем шаг
phi_new = phi + a

# Считаем собственное значения из отношения Рэлея
lambda = rayleigh_quotient(phi_new)

# Результат
println("max  Eigenvalues = ",lambda)
