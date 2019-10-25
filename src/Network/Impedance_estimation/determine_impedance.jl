export determine_impedance

"""
    function determine_impedance(network::Network; input_pins :: Array{Any},
        output_pins :: Array{Any}, elim_elements :: Array{Symbol},
        omega_range = (-3, 5, 100))

Estimation of the impedance visible from the port, between input and
output pins. Port pins can possibly be connected to some elements, which
should not be considered for the impedance estimation. Those elements
are listed as `elim_elements`.

Specification for the impedance estimation is given on the example of the following
network that consists of the DC voltage source and a cable.
```
net = @network begin
    vs = dc_source(voltage = 500e3)
    c = cable(length = 100e3, positions = [(0,1)], earth_parameters = (1,1,1),
    C1 = Conductor(rₒ = 24.25e-3, ρ = 1.72e-8), C2 = Conductor(rᵢ = 41.75e-3, rₒ = 46.25e-3, ρ = 22e-8),
    C3 = Conductor(rᵢ = 49.75e-3, rₒ = 60.55e-3, ρ = 18e-8, μᵣ = 10),
    I1 = Insulator(rᵢ = 24.25e-3, rₒ = 41.75e-3, ϵᵣ = 2.3),
    I2 = Insulator(rᵢ = 46.25e-3, rₒ = 49.75e-3, ϵᵣ = 2.3),
    I3 = Insulator(rᵢ = 60.55e-3, rₒ = 65.75e-3, ϵᵣ = 2.3))
    vs[1.1] ⟷ c[1.1] ⟷ Node1
    vs[2.1] ⟷ c[2.1] ⟷  gnd
end
```
To determine impedance visible from the voltage source `vs`, the following command
should be called:
```
imp, omega = determine_impedance(net, elim_elements = [:vs], input_pins = Any[:Node1], output_pins = Any[:gnd],
        omega_range = (-1,6,10000))
```
Impedance is determined inside network `net`, from the element `vs` and the port defined
with `input_pins` as array consisting of `Node1` and the `output_pins` containing
array with `gnd`. Impedance is estimated for the frequency in [rad/s] with the range
0.1 to 1000000 in 10000 points.

The function returns complex impedance map and two arrays: first is the impedance
array and the second one is frequency array.
"""
function determine_impedance(network::Network; input_pins :: Array{Any},
    output_pins :: Array{Any}, elim_elements :: Array{Symbol},
    omega_range = (-3, 5, 100))

    """
    function make_lists(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
        elim_elements::Array{Symbol}, start_pins::Array{Symbol}, end_pins::Array{Symbol})

        Makes dictionary entries containing list of elements, nodes and outputs
        between the start and end pins. It calls itself recursively.
        Dictionary entries are:
        :node_list = consists of nonground and nonoutput processed nodes
        :element_list = consists of processed elements
        :output_list = consists of processed ground and output (end) nodes

    """
    function make_lists(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
        elim_elements::Array{Symbol}, start_pins::Array{Symbol}, end_pins::Array{Symbol})

        for node_name in start_pins
            node = netfor!(net, node_name)

            # if it is the end of the port, add it to the output and return
            if in(node_name, end_pins)
                !in(node_name, dict[:output_list]) && push!(dict[:output_list], node_name)
                continue
            end
            if occursin("gnd", string(node_name))
                if in(node_name, start_pins)
                    push!(dict[:node_list], node_name)
                elseif !in(node_name, dict[:output_list])
                    push!(dict[:output_list], node_name)
                end
            else
                # add nodes to the node list
                !in(node_name, dict[:node_list]) && push!(dict[:node_list], node_name)
            end

            # find all elements inside the port connected to the node
            elements_pins = filter(p ->  !in(p[1], elim_elements) && !in(p[1], dict[:element_list]), node)

            for (element, pin) in elements_pins
                push!(dict[:element_list], element) # add element's symbol to the list

                other_nodes = get_nodes(net.elements[element], pin) # get the pins from the other side of element
                make_lists(net, dict, elim_elements, other_nodes, end_pins)
            end
        end
    end

    """
    function make_abcd(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                        start_pins::Array{Symbol}, end_pins::Array{Symbol})

        Creates ABCD represntation of the network between start pins and end pins
        using data written in dictionary dict.
    """
    function make_abcd(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
                        start_pins::Array{Symbol}, end_pins::Array{Symbol},
                        omega_range)
        (min_ω, max_ω, n_ω) = omega_range
        n = (max_ω - min_ω) / n_ω
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
        s = 1im * omegas
        nᵥ = length(dict[:node_list])       # number of unknown node voltages
        nₒ = length(dict[:output_list])     # number of output nodes = grounds and end pin nodes
        nₙ = nᵥ + nₒ                        # number of nodes
        nᵢ = length(start_pins)             # number of input pins
        mₚ = sum(nip_abcd(net.elements[element]) for element in dict[:element_list])
        mₛ = sum(nop_abcd(net.elements[element]) for element in dict[:element_list])

        p = length(output_pins)
        Zₜ = zeros(Complex, p, p)

        impedance = []
        omega = []

        for index in 1:length(omegas)
            matrix = zeros(Complex, nₙ+2mₚ, nᵥ+nᵢ+mₚ+mₛ)
            output = zeros(Complex, nₙ+2mₚ, 2nₒ)
            elim_rows = Int[]
            elim_cols = Int[]

            # fix input current
            for node_name in start_pins
                i = findfirst(p -> p == node_name, dict[:node_list])
                j = findfirst(p -> p == node_name, start_pins)
                if (i != nothing)
                    matrix[i, nᵥ+j] = -1
                else
                    i = findfirst(p -> p == node_name, dict[:output_list])
                    matrix[nᵥ+i, nᵥ+j] = -1
                end
            end

            for element in dict[:element_list]
                nₑ = findfirst(p -> p == element, dict[:element_list]) - 1              # element position in dict
                if nₑ > 0
                    nₑₚ = 2sum(nip_abcd(net.elements[element]) for element in
                        (dict[:element_list])[1:nₑ]) + nₙ                               # element input position
                    nₑₛ = sum(np_abcd(net.elements[element]) for element in
                        (dict[:element_list])[1:nₑ]) + nᵥ + nᵢ                          # element output position
                else
                    nₑₚ =  nₙ
                    nₑₛ = nᵥ + nᵢ
                end

                # get ABCD parameters and split to matrices
                ABCD = get_abcd(net.elements[element], s[index])
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
            end

            if det(matrix) == 0
                continue
                sol = pinv(matrix) * output
            else
                sol = matrix \ output
            end

            pᵢ = length(start_pins)
            pₒ = length(end_pins)
            abcd = zeros(Complex, 2pᵢ, 2pₒ)
            for iₚ in 1:pᵢ
                iₗ = findfirst(p -> p == start_pins[iₚ], dict[:node_list])
                for jₚ in 1:pₒ
                    jₗ = findfirst(p -> p == end_pins[jₚ], dict[:output_list])
                    abcd[iₚ,jₚ] = sol[iₗ, end-2nₒ+2(jₗ-1)+1]
                    abcd[iₚ,jₚ+pₒ] = sol[iₗ, end-2nₒ+2jₗ]
                    abcd[iₚ+pᵢ,jₚ] = sol[iₗ + nᵥ, end-2nₒ+2(jₗ-1)+1]
                    abcd[iₚ+pᵢ,jₚ+pₒ] = sol[iₗ + nᵥ, end-2nₒ+2jₗ]
                end
            end

            z = closing_impedance(abcd, Zₜ)
            push!(impedance, z)
            push!(omega, abs(s[index]))
        end
        return impedance, omega
    end


    isempty(input_pins) && throw(ArgumentError("Impedance cannot be determined from nonexistent input port."))
    isempty(output_pins) && throw(ArgumentError("Impedance cannot be determined from nonexistent output port."))
    for i in 1:length(input_pins)
        input_pins[i] = netname(network, input_pins[i])
    end
    input_pins = convert(Array{Symbol}, input_pins)
    for i in 1:length(output_pins)
        output_pins[i] = netname(network, output_pins[i])
    end
    output_pins = convert(Array{Symbol}, output_pins)

    dict = Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[],
        :output_list => Symbol[])
    make_lists(network, dict, elim_elements, input_pins, output_pins)

    # if end pin is not in dictionary
    if any(!in(end_pin, dict[:output_list]) for end_pin in output_pins)
        dict_o = Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[],
            :output_list => Symbol[])
        make_lists(network, dict_o, elim_elements, output_pins, input_pins)
        append!(dict[:output_list], setdiff(setdiff(dict_o[:node_list], dict[:output_list]), dict[:node_list]))
        if !isempty(setdiff(dict_o[:element_list], dict[:element_list]))
            append!(dict[:element_list], setdiff(dict_o[:element_list], dict[:element_list]))
        end
    end

    # println(dict)
    impedance, omega = make_abcd(network, dict, input_pins, output_pins, omega_range)
    return impedance, omega
end
