# New version, only being used to generate the admittance matrix of the entire network.
# Function description
# Potentially makes sense to make this function general again, as initially intended.
# Call make_y with make_y_edge!
# Call make_y with make_y_node!
function make_y(net :: Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}}, s :: Complex)

    n = length(dict[:node_list])
    Y_matrix = zeros(Complex, n, n)
    for element in dict[:element_list]
        element = net.elements[element]
        if is_passive(element)
            Y = get_y(element, s)
            # Get number of phases
            phase = Int(length(element.pins) / 2)
            for (key₁, val₁) in element.pins, (key₂, val₂) in element.pins # key is the pin name, val is the node name
                i = (parse(Int, string(key₁)[1]) - 1) * phase + parse(Int, string(key₁)[3])
                j = (parse(Int, string(key₂)[1]) - 1) * phase + parse(Int, string(key₂)[3])
                iₚ = findfirst(p -> p == val₁, dict[:node_list]) # Find row index of node in the admittance matrix
                (iₚ === nothing) ? continue : nothing # If the node is not in the node list, skip it
                jₚ = findfirst(p -> p == val₂, dict[:node_list]) # Find column index of node in the admittance matrix
                (jₚ === nothing) ? continue : nothing
                Y_matrix[iₚ, jₚ] += Y[i,j]
            end
        end
    end

    return Y_matrix
end
