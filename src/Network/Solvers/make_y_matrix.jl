export make_y_matrix
"""
    make_y_matrix(network::Network; elim_elements::Array{Symbol} = Symbol[], 
                  input_pins::Array{Any}, omega_range=(-3, 5, 100))

Computes the frequency-dependent admittance matrix (Y-matrix) for a given electrical network.  
It determines the nodal admittance representation based on the specified input pins and frequency range.

# Arguments
- `network::Network`: The electrical network model for which the admittance matrix is computed.
- `elim_elements::Array{Symbol}` (default `[]`): List of elements to be excluded from the computation.
- `input_pins::Array{Any}`: A list of input nodes (ports) that define the matrix rows/columns.
- `omega_range::Tuple{Real, Real, Int}` (default `(-3, 5, 100)`):  
  Defines the range of angular frequencies in logarithmic scale:
  - `min_ω`: Minimum exponent (base 10) for frequency.
  - `max_ω`: Maximum exponent (base 10) for frequency.
  - `n_ω`: Number of frequency points.

# Behavior
- Recursively explores the network from the input pins to determine the node and element connectivity.
- Constructs the Y-matrix by eliminating specified elements and non-essential nodes.
- Computes admittance matrices over a frequency range given by `omega_range`.
- Frequencies are logarithmically spaced and converted to radians per second (`ω = 2πf`).

# Output
- Returns `Ybus`, an array of admittance matrices (one per frequency point).
- The output matrices exclude ground-connected output nodes.

# Exceptions
- Throws `ArgumentError` if `input_pins` is empty.

# Example
```julia
net = create_network(...)  # Assume a valid network object
Y_matrices = make_y_matrix(net, input_pins=["N1", "N2"], omega_range=(-2, 4, 50))
```
"""
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
