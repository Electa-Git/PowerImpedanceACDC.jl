export make_y_matrix

function make_y_matrix(network :: Network; elim_elements :: Array{Symbol}, input_pins :: Array{Any}, omega_range = (-3, 5, 100), parameters_type = :ABCD)

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
    # function make_lists(net::Network, dict::Dict{Symbol, Array{Union{Symbol,Int}}}, start_pins::Array{Symbol})

    #     for node_name in start_pins
    #         node = netfor!(net, node_name)

    #         # if it is the end of the port, add it to the output and return
    #         if occursin("gnd", string(node_name))
    #             if !in(node_name, dict[:output_list])
    #                 push!(dict[:output_list], node_name)
    #                 continue
    #             end
    #         else
    #             # add nodes to the node list
    #             !in(node_name, dict[:node_list]) && push!(dict[:node_list], node_name)
    #         end

    #         # find all elements inside the port connected to the node
    #         elements_pins = filter(p ->  !in(p[1], dict[:element_list]), node)

    #         for (element, pin) in elements_pins
    #             !in(element, dict[:element_list]) && push!(dict[:element_list], element) # add element's symbol to the list, only if the element has not been added before
    #             other_nodes = get_nodes(net.elements[element], pin) # get the pins from the other side of element
    #             make_lists(net, dict, other_nodes)
    #         end
    #     end
    # end
    # TODO: Right now, if the components in elim_elements are not in the list, no error is given. This could be generalized.
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

    isempty(input_pins) && throw(ArgumentError("Impedance cannot be determined from nonexistent input port."))
    # isempty(output_pins) && throw(ArgumentError("Impedance cannot be determined from nonexistent output port."))
    for i in 1:length(input_pins)
        input_pins[i] = netname(network, input_pins[i])
    end
    input_pins = convert(Array{Symbol}, input_pins)
    # for i in 1:length(output_pins)
    #     output_pins[i] = netname(network, output_pins[i])
    # end
    # output_pins = convert(Array{Symbol}, output_pins)

    dict = Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[],
        :output_list => Symbol[])
    make_lists(network, dict, elim_elements, unique(input_pins))

    # reorder dictionary
    input_order = []
    output_order = []
    for input in unique(input_pins)
        i = findfirst(p -> p == input, dict[:node_list])
        push!(input_order, i)
    end
    dict[:node_list] = dict[:node_list][[input_order; setdiff(1:length(dict[:node_list]), input_order)]]

    # make frequency range
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    omegas= 2*pi* 10 .^range(min_ω, max_ω, length= n_ω) 

    Ybus=[] 

    # dict[:node_list] = [dict[:node_list]; dict[:output_list]]
    dict[:node_list] = [input_pins; dict[:output_list]]
        for omega in omegas
            z = make_y(network, dict, input_pins, omega*1im)
            z = z[1:end-length(dict[:output_list]),1:end-length(dict[:output_list])]
            # z = kron(z, Int[i for i in 1:length(input_pins)])
            push!(Ybus, z)
        end

    return Ybus
end
