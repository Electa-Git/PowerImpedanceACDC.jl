export make_z

# include("../Components/tools/tools.jl")

"""
function make_z(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                    start_pins::Array{Symbol}, end_pins::Array{Symbol})

    Based on the individual component equations, this function obtains the impedance between the start and end pins.
"""
function make_z(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                    start_pins::Array{Symbol}, end_pins::Array{Symbol}, s :: Complex)
    nᵥ = length(unique(dict[:node_list]))       # number of unknown node voltages
    nₒ = length(unique(dict[:output_list]))     # number of output nodes = grounds and end pin nodes
    nₙ = nᵥ + nₒ                                # number of nodes
    nᵢ = length(unique(start_pins))             # number of input pins
    mₚ = sum(nip_abcd(net.elements[element]) for element in dict[:element_list])
    mₛ = sum(np_abcd(net.elements[element]) for element in dict[:element_list])

    sysmatrix = zeros(ComplexF64, nₙ+mₚ, nᵥ+nₒ+mₛ)      # Kirchoff plus component equations, the unknown variables are the node voltages, component currents and the output currents 
    inmatrix = zeros(ComplexF64, nₙ+mₚ, nᵢ)             # Kirchoff plus component equations, the known variables are the input currents

    # One element of the KCL equation for component inputs
    for j in 1:length(start_pins)
        k = findfirst(p -> p == start_pins[j], dict[:node_list])
        # sysmatrix[k, nᵥ+k] = -1
        inmatrix[k,j] = 1
    end

    # indices for matrix elements
    nₑₚ =  nₙ   # element input position
    nₑₛ = nᵥ + nₒ   # element output positions
    for element in dict[:element_list]

        if isa(net.elements[element].element_value, MMC)
            Y = get_abcd(net.elements[element], s)
            pₚ = nip_abcd(net.elements[element])  # element input pins
            I = convert(Array{ComplexF64}, Diagonal([1 for dummy in 1:pₚ]))
            
            sysmatrix[nₑₚ+1:nₑₚ+pₚ, nₑₛ+1:nₑₛ+pₚ] = -I                                 # -Iₚ in element eq
            for (pin, node_name) in pairs(net.elements[element].pins)
                n = findfirst(p -> p == node_name, dict[:node_list])  # position of the node in dict
                i = parse(Int,string(pin)[1])
                yi = (i-1) + parse(Int,string(pin)[3:end])
                if (n == nothing)
                    n = findfirst(p -> p == node_name, dict[:output_list])
                    (n == nothing) && continue
                    sysmatrix[nᵥ+n,  nₑₛ+yi] = -(-1)^i                  # Current direction (1.1, 2.1, 2.2)
                    sysmatrix[nᵥ+n, nᵥ+n] = (-1)^i
                    # inmatrix[nₑₚ+1:nₑₚ+pₚ, 2(n-1)+1] += -Y[1:end, yi] # All outputs are assumed to be grounded. This can be implemented in a future work.               # -Y[:, yi]
                else
                    # Original
                    sysmatrix[n, nₑₛ+yi] = -(-1)^i                          # The other part of the KCL
                    sysmatrix[nₑₚ+1:nₑₚ+pₚ, n] += Y[1:end, yi]               # +Y[:, yi] 
                end
            end
        else
            if isa(net.elements[element].element_value, SynchronousMachine)
                # TODO: Generalize this. Check instead if the component has an ABCD or a Y representation.
                pₚ = 2          # element input pins
                pₛ = 2          # element output pins
                
                Y = get_abcd(net.elements[element], s)
                Z = inv(Y)

                I = convert(Array{ComplexF64}, Diagonal([1 for dummy in 1:pₚ]))

                a = I
                b = Z
                c = convert(Array{ComplexF64}, zeros(2,2))
                d = I
            else
                ABCD = get_abcd(net.elements[element], s) 

                pₚ = Int(size(ABCD,1)/2)                                                # element input pins
                pₛ = Int(size(ABCD,2)/2)                                                # element output pins
                (a, b, c, d) = (ABCD[1:pₚ,1:pₛ], ABCD[1:pₚ,pₛ+1:end], ABCD[pₚ+1:end,1:pₛ], ABCD[pₚ+1:end, pₛ+1:end]) # TODO: Efficiency improvement possibility. Using the matrices without allocation.
            end

            I = convert(Array{ComplexF64}, Diagonal([1 for dummy in 1:pₚ]))
            sysmatrix[nₑₚ+1:nₑₚ+pₚ, nₑₛ+pₚ+1:nₑₛ+pₚ+pₛ] += b                             # BIₛ in element eq
            sysmatrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, nₑₛ+1:nₑₛ+pₚ] = -I                               # -Iₚ in element eq
            sysmatrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, nₑₛ+pₚ+1:nₑₛ+pₚ+pₛ] += d                          # DIₛ in element eq

            for (pin, node_name) in pairs(net.elements[element].pins)
                i = findfirst(p -> p == node_name, dict[:node_list])                # position of the node in dict
                j = parse(Int,string(pin)[3:end])
                if (occursin("1.", string(pin)))   # input pin
                    if (i == nothing)
                        i = findfirst(p -> p == node_name, dict[:output_list])
                        if (i == nothing)
                            continue
                        end
                        sysmatrix[nᵥ+i, nₑₛ+j] = 1                                     # +Iₚ[j] in node
                        sysmatrix[nᵥ+i, nᵥ+i] = 1                                       # current value in node
                        # sysmatrix[nₑₚ+j, 2(i-1)+1] = 1                                  # +Vₚ[j] in element eq # Output nodes are grounded and eliminated. TODO: generalization later
                    else
                        sysmatrix[i, nₑₛ+j] = 1                                        # +Iₚ[j] in node
                        sysmatrix[nₑₚ+j, i] = -1                                       # -Vₚ[j] in element eq
                    end
                else    # component output pin
                    if (i == nothing) # output node or input node
                        i = findfirst(p -> p == node_name, dict[:output_list])
                        (i == nothing) && continue
                        sysmatrix[nᵥ+i, nᵥ+i] = -1                                        # current value in node
                        sysmatrix[nᵥ+i,  nₑₛ+pₚ+j] = -1                                 # -Iₛ[j] in node
                        # The equations below are removed because the output nodes are grounded.
                        # sysmatrix[nₑₚ+1:nₑₚ+pₚ, 2(i-1)+1] += a[1:end, j]               # -AVₛ[j] in element eq
                        # sysmatrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, 2(i-1)+1] += c[1:end, j]            # -CVₛ[j] in element eq
                    else # regular node
                        sysmatrix[i, nₑₛ+pₚ+j] = -1                                     # -Iₛ[j] in node
                        sysmatrix[nₑₚ+1:nₑₚ+pₚ, i] += a[1:end, j]                       # +AVₛ[j] in element eq
                        sysmatrix[nₑₚ+pₚ+1:nₑₚ+2pₚ, i] += c[1:end, j]                    # +CVₛ[j] in element eq
                    end
                end
            end
        end
        nₑₚ += nip_abcd(net.elements[element]) # update input element position
        nₑₛ += np_abcd(net.elements[element])   # update output element position
    end

    sol = sysmatrix \ inmatrix
    Z = zeros(ComplexF64,size(inmatrix,2),size(inmatrix,2))
    for j in 1:size(inmatrix,2)
        k = findfirst(p -> p == start_pins[j], dict[:node_list])
        Z[k,:] = sol[k,:]
    end
    return Z
end
