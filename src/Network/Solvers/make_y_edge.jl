export make_y_edge
"""


"""

function make_y_edge(network::Network; freq_range = (1,1e3, 1000))

node_list= Symbol[]
element_list= Symbol[]

for node in keys(network.nets)
    
    if occursin("gnd", string(node))
        continue # Skip ground nodes
    end

    # Check whether node is already part of the node list of Yedge
    if in(node, node_list)
        continue # Skip nodes that are already in the list
    end

    # Check whether the current node is part of Yedge

    nets=netfor!(network, node) # Get all elements (designator, pin) connected to the node
    

    # Check whether the node is connected to a source
    # If it is, skip the node
    skip_node = false # Flag to skip the node if it is connected to a source
    for net in nets # Iterate over all elements connected to the node

        element=network.elements[net[1]] # Get the element, via the designator net[1]
        
        if is_source(element) 
            
            skip_node = true # Skip the node if it is connected to a source

        end

    end

    if skip_node
        continue # Skip the node if it is connected to a source
    end

    isDCnode=false # Flag to check whether the node is a DC node
    isACnode= false # Flag to check whether the node is an AC node


    # Identify DC or AC node by checking the pins .1 -->DC, .2 --> AC
    # If DC node elements connected to have only one pin connected
    # If AC node elements connected to have two pins connected
    for net in nets # Iterate over all elements connected to the node

       pin = string(net[2]) # Get the pin name

       if occursin(".2",pin)

        # AC port!

       else # DC or AC port further check base on 

            element=network.elements[net[1]]
            
            for (element_pin, node) in element.pins # Iterate over all pins of the element
 
                if replace(pin, ".1" => ".2") == string(element_pin) # Check if the pin is a AC pin
                   
                    # AC pin
                    push!(node_list, node) # Add node to the list
                    push!(element_list, net[1]) # Add element to the list
                    break # Break the loop as we found the DC pin
                else
                    #DC pin



                end

            

            end

       end



    end

    

# If yes order correctly [DC], ACd, ACq

    # If AC node, then get the other node --> pair
        # Push to node_list

    # If DC node 
        #Check whether MMC node or not
            # If MMC node, then get the triplet
                # Push to node_list
            # If not MMC node
                # Push to node_list










    push!(node_list, node.name) # Add node to the list
end





# Elements weed out all active elements, i.e. MMC, sources, SG, TLC


# Initialize dict to hold the nodes and elements for the admittance matrix
dict= Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[])

dict[:node_list]=node_list
dict[:element_list]=element_list

Yedge=[] # Preallocate the admittance matrix for each frequency

    for frequencies in freq_range

        # Get the admittance matrix for the current frequency
        Y = make_y(network, dict, frequencies*1im*2*pi)
        push!(Yedge,Y)

    end

    return Yedge

end