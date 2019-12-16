
function make_y(net :: Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                    start_pins::Array{Symbol}, end_pins::Array{Symbol}, s :: Complex)

    n = length(dict[:node_list])
    Y_matrix = zeros(Complex, n, n)
    for element in dict[:element_list]
        element = net.elements[element]
        Y = get_y(element, s)
        for (key₁, val₁) in element.pins, (key₂, val₂) in element.pins
            i = parse(Int, string(key₁)[1]) - 1 + parse(Int, string(key₁)[3])
            j = parse(Int, string(key₂)[1]) - 1 + parse(Int, string(key₂)[3])
            iₚ = findfirst(p -> p == val₁, dict[:node_list])
            (iₚ == nothing) ? continue : nothing
            jₚ = findfirst(p -> p == val₂, dict[:node_list])
            (jₚ == nothing) ? continue : nothing
            Y_matrix[iₚ, jₚ] += Y[i,j]
        end
    end

    # HRMA
    Λ = Diagonal(eigvals(Y_matrix))
    T = eigvecs(Y_matrix)

    Z = T * inv(Λ) * inv(T)

    I_matrix = zeros(n,n)
    indexes = []
    for node in [start_pins end_pins]
        i = findfirst(p -> p == node, dict[:node_list])
        push!(indexes, i)
        I_matrix[i,i] = 1
    end

    sol = (Z * I_matrix)[indexes, :]
    sol = sol[:, indexes]

    return convert(Array{Complex}, sol)
end
