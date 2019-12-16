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
function determine_impedance(network :: Network; input_pins :: Array{Any},
    output_pins :: Array{Any}, elim_elements :: Array{Symbol},
    omega_range = (-3, 5, 100), parameters_type = :ABCD)

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
                push!(dict[:element_list], element) # add element's symbol to the list
                other_nodes = get_nodes(net.elements[element], pin) # get the pins from the other side of element
                make_lists(net, dict, elim_elements, other_nodes, end_pins)
            end
        end
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

    # reorder dictionary
    input_order = []
    output_order = []
    for input in input_pins
        i = findfirst(p -> p == input, dict[:node_list])
        push!(input_order, i)
    end
    dict[:node_list] = dict[:node_list][[input_order; setdiff(1:length(dict[:node_list]), input_order)]]
    for output in output_pins
        i = findfirst(p -> p == output, dict[:output_list])
        push!(output_order, i)
    end
    dict[:output_list] = dict[:output_list][[output_order; setdiff(1:length(dict[:output_list]), output_order)]]

    # make frequency range
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]

    # make closing impedance
    p = length(output_pins)
    Zₜ = zeros(Complex, p, p)

    impedance = []

    if parameters_type == :ABCD
        for omega in omegas
            abcd = make_abcd(network, dict, input_pins, output_pins, omega*1im)
            z = closing_impedance(abcd, Zₜ)
            push!(impedance, z)
        end
    else
        dict[:node_list] = [dict[:node_list]; dict[:output_list]]
        delete!(dict, :output_list)
        for omega in omegas
            z = make_y(network, dict, input_pins, output_pins, omega*1im)
            z = kron(z, Int[i for i in 1:length(input_pins)])
            push!(impedance, z)
        end
    end

    return impedance, omegas
end
