export make_y_matrix

function make_y_matrix(network :: Network; elim_elements :: Array{Symbol} = Array{Symbol}(undef,0,0), input_pins :: Array{Any}, omega_range = (-3, 5, 100))

    function make_lists(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}},
        elim_elements::Array{Symbol}, start_pins::Array{Symbol})

        for node_name in start_pins
            node = netfor!(net, node_name)

            # if it is the end of the port, add it to the output and return
            if occursin("gnd", string(node_name))
                if !in(node_name, dict[:output_list])
                    push!(dict[:output_list], node_name)
                    continue
                end
            else
                # add nodes to the node list
                !in(node_name, dict[:node_list]) && push!(dict[:node_list], node_name)
            end

            # find all elements inside the port connected to the node
            elements_pins = filter(p ->  !in(p[1], elim_elements) && !in(p[1], dict[:element_list]), node)

            for (element, pin) in elements_pins
                !in(element, dict[:element_list]) && push!(dict[:element_list], element) # add element's symbol to the list, only if the element has not been added before
                other_nodes = get_nodes(net.elements[element], pin) # get the pins from the other side of element
                make_lists(net, dict, elim_elements, other_nodes)
            end
        end
    end

    isempty(input_pins) && throw(ArgumentError("Input ports have to be specified to determine the admittance matrix."))
    for i in 1:length(input_pins)
        input_pins[i] = netname(network, input_pins[i])
    end
    input_pins = convert(Array{Symbol}, input_pins)

    dict = Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[],
        :output_list => Symbol[])
    make_lists(network, dict, elim_elements, unique(input_pins))

    # make frequency range
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    omegas= 2*pi* 10 .^range(min_ω, max_ω, length= n_ω) 

    Ybus=[] 

    dict[:node_list] = [input_pins; dict[:output_list]]
    for omega in omegas
        Y = make_y(network, dict, omega*1im)
        push!(Ybus, Y[1:end-length(dict[:output_list]),1:end-length(dict[:output_list])]) # Remove the output nodes that are at the end of the matrix.
    end  

    return Ybus
end
