export make_abcd

"""
function make_abcd(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                    start_pins::Array{Symbol}, end_pins::Array{Symbol})

    Creates ABCD represntation of the network between start pins and end pins
    using data written in dictionary dict.
"""
function make_abcd(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                    start_pins::Array{Symbol}, end_pins::Array{Symbol}, s :: Complex)
    nᵥ = length(unique(dict[:node_list]))       # number of unknown node voltages
    nₒ = length(unique(dict[:output_list]))     # number of output nodes = grounds and end pin nodes
    nₙ = nᵥ + nₒ                                # number of nodes
    nᵢ = length(unique(start_pins))             # number of input pins
    mₚ = sum(nip_abcd(net.elements[element]) for element in dict[:element_list])
    mₛ = sum(np_abcd(net.elements[element]) for element in dict[:element_list])

    matrix = zeros(Complex, nₙ+mₚ, nᵥ+nᵢ+mₛ)
    output = zeros(Complex, nₙ+mₚ, 2nₒ)
    elim_rows = Int[]
    elim_cols = Int[]

    # fix input and output currents
    for j in 1:length(start_pins)
        k = findfirst(p -> p == start_pins[j], dict[:node_list])
        matrix[k, nᵥ+k] = -1
    end

    # indices for matrix elements
    nₑₚ =  nₙ   # element input position
    nₑₛ = nᵥ + nᵢ   # element output positions
    for element in dict[:element_list]
        if !isa(net.elements[element].element_value, MMC)
            # get ABCD parameters and split to matrices
            ABCD = get_abcd(net.elements[element], s)
            pₚ = Int(size(ABCD,1)/2)                                                # element input pins
            pₛ = Int(size(ABCD,2)/2)                                                # element output pins
            (a, b, c, d) = (ABCD[1:pₚ,1:pₛ], ABCD[1:pₚ,pₛ+1:end], ABCD[pₚ+1:end,1:pₛ], ABCD[pₚ+1:end, pₛ+1:end])

            I = convert(Array{Complex}, Diagonal([1 for dummy in 1:pₚ]))
            matrix[nₑₚ+1:nₑₚ+pₚ, nₑₛ+pₚ+1:nₑₛ+pₚ+pₛ] += b                               # BIₛ in element eq
            matrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, nₑₛ+1:nₑₛ+pₚ] = -I                                 # -Iₚ in element eq
            matrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, nₑₛ+pₚ+1:nₑₛ+pₚ+pₛ] += d                            # DIₛ in element eq

            for (pin, node_name) in pairs(net.elements[element].pins)
                i = findfirst(p -> p == node_name, dict[:node_list])                # position of the node in dict
                j = parse(Int,string(pin)[3:end])
                if (occursin("1.", string(pin)))   # input pin
                    if (i == nothing)
                        i = findfirst(p -> p == node_name, dict[:output_list])
                        if (i == nothing)
                            push!(elim_rows, nₑₚ+j)                                 # eliminate Vₚ[j] row
                            push!(elim_rows, nₑₚ+pₚ+j)                               # eliminate Iₚ[j] row
                            push!(elim_cols, nₑₛ+j)                                 # eliminate Iₚ[j] column
                            push!(elim_cols, nₑₛ+pₚ+j)                               # eliminate Iₛ[j] column
                            continue
                        end
                        matrix[nᵥ+i, nₑₛ+j] = 1                                     # +Iₚ[j] in node
                        output[nᵥ+i, 2i] = -1                                       # current value in node
                        output[nₑₚ+j, 2(i-1)+1] = 1                                 # +Vₚ[j] in element eq
                    else
                        matrix[i, nₑₛ+j] = 1                                        # +Iₚ[j] in node
                        matrix[nₑₚ+j, i] = -1                                       # -Vₚ[j] in element eq
                    end
                else    # output pin
                    if (i == nothing)
                        i = findfirst(p -> p == node_name, dict[:output_list])
                        (i == nothing) && continue
                        output[nᵥ+i, 2i] = 1                                        # current value in node
                        matrix[nᵥ+i,  nₑₛ+pₚ+j] = -1                                 # -Iₛ[j] in node
                        output[nₑₚ+1:nₑₚ+pₚ, 2(i-1)+1] += -a[1:end, j]               # -AVₛ[j] in element eq
                        output[nₑₚ+pₚ+1:nₑₚ+2pₚ, 2(i-1)+1] += -c[1:end, j]            # -CVₛ[j] in element eq
                    else
                        matrix[i, nₑₛ+pₚ+j] = -1                                     # -Iₛ[j] in node
                        matrix[nₑₚ+1:nₑₚ+pₚ, i] += a[1:end, j]                       # +AVₛ[j] in element eq
                        matrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, i] += c[1:end, j]                    # +CVₛ[j] in element eq
                    end
                end
            end

        else
            Y = get_abcd(net.elements[element], s)
            pₚ = nip_abcd(net.elements[element])  # element input pins
            I = convert(Array{Complex}, Diagonal([1 for dummy in 1:pₚ]))

            matrix[nₑₚ+1:nₑₚ+pₚ, nₑₛ+1:nₑₛ+pₚ] = -I                                 # -Iₚ in element eq
            for (pin, node_name) in pairs(net.elements[element].pins)
                n = findfirst(p -> p == node_name, dict[:node_list])  # position of the node in dict
                i = parse(Int,string(pin)[1])
                yi = (i-1) + parse(Int,string(pin)[3:end])
                if (n == nothing)
                    n = findfirst(p -> p == node_name, dict[:output_list])
                    (n == nothing) && continue
                    matrix[nᵥ+n,  nₑₛ+yi] = -(-1)^i
                    output[nᵥ+n, 2n] = -(-1)^i
                    output[nₑₚ+1:nₑₚ+pₚ, 2(n-1)+1] += -Y[1:end, yi]               # -Y[:, yi]
                else
                    matrix[n, nₑₛ+yi] = -(-1)^i
                    matrix[nₑₚ+1:nₑₚ+pₚ, n] += Y[1:end, yi]                       # +Y[:, yi]
                end
            end

        end
        nₑₚ += nip_abcd(net.elements[element]) # update input element position
        nₑₛ += np_abcd(net.elements[element])   # update output element position
    end

    sol = pinv(matrix) * output

    pᵢ = length(unique(start_pins))
    pₒ = length(unique(end_pins))
    abcd = zeros(Complex, 2pᵢ, 2pₒ)
    for i in 1:pᵢ, j in 1:pₒ
        k = findfirst(p -> p == start_pins[i], dict[:node_list])
        l = findfirst(p -> p == end_pins[j], dict[:output_list])
        abcd[i,j] = sol[k, end-2nₒ+2(l-1)+1]
        abcd[i,j+pₒ] = sol[k, end-2nₒ+2l]
        abcd[i+pᵢ,j] = sol[k + nᵥ, end-2nₒ+2(l-1)+1]
        abcd[i+pᵢ,j+pₒ] = sol[k + nᵥ, end-2nₒ+2l]
    end
    return abcd
end
