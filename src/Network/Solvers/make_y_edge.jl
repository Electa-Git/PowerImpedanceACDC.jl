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

    nets=PowerImpedanceACDC.netfor!(network, node) # Get all elements (designator, pin) connected to the node
    

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

    isDCnode= false # Flag to check whether the node is a DC node
    isACnode= false # Flag to check whether the node is an AC node

    # Identify DC or AC node by checking the pins connected to node:
    # If DC node: Only one pin for each port, e.g. 1.1, 2.1 etc.
    # If AC node: Two pins for each port, e.g. 1.1, 1.2, 2.1, 2.2 etc.
    
    # Get the pin of one of the element connected to the node
    pin = string(nets[1][2]) # Get the pin name

    if occursin(".2", pin) # Check whether the pin is a AC pin
        
        isACnode=true # Set the flag to true
    
    end
    if occursin(".1", pin) 
       
        pin=replace(pin, ".1" => ".2") # Replace .1 with .2 to search for the other AC pin if exists

        designator = nets[1][1] # Get the designator of the element
        pin=Symbol(pin)

        if PowerImpedanceACDC.netfor!(network, (designator,pin)) === Symbol("") # Check whether this net (designator, pin) exists
        
            isDCnode= true

        else
            
            isACnode=true
        end

    end


    if isACnode

        # Get the other AC node
        if occursin(".2", pin)
            pin=replace(pin, ".2" => ".1") # Replace .2 with .1 to search for the other AC pin if exists
        else
            pin=replace(pin, ".1" => ".2") # Replace .1 with .2 to search for the other AC pin if exists

        end

        designator = nets[1][1] # Get the designator of the element
        pin=Symbol(pin)
        node2=PowerImpedanceACDC.netfor!(network, (designator,pin))

        if occursin(".1",pin) # Node 2 is a d node
 
 
            push!(node_list, [node2, node]) # Add the node to the list

        else  # Node 2 is a q node

            push!(node_list, [node, node2]) # Add the node to the list

        end

    

    end


    isConverternode = false # Flag to check whether the node is a converter node
    if isDCnode
        

        # Check whether the DC node is a MMC node

        for net in nets # Iterate over all elements connected to the node

        element=network.elements[net[1]] # Get the element, via the designator net[1]
        
        if is_converter(element) 
            
            isConverternode=true
            break # Break the loop if a converter node is found

        end

        end

        if isConverternode # Proper node ordering required to match 3x3 admittance matrix of converter

            triplet=Array{Union{Symbol,Int}}(undef,3)
            triplet[1]=node # Add the node to the triplet
        
            for net in nets # Iterate over all elements connected to the node

            element=network.elements[net[1]] # Get the element, via the designator net[1]
        
            if is_converter(element) 
            
                for (element_pin,element_node) in element.pins

                    if element_node == node # Check whether the element pin is connected to the current node
                        continue # Skip the current node
                    else # AC node

                        if occursin(".1", string(element_pin)) # Check whether the pin is a d pin
                            triplet[2]=element_node # Add the other node to the triplet
                        end
                        if occursin(".2", string(element_pin)) # Check whether the pin is a q pin
                            triplet[3]=element_node # Add the other node to the triplet
                        end

                    end

                end

            end

        end
            
            
            triplet=[node,]


            # AC node pairs could be already part of node_list!

        else

            # If not a converter node, then it is a DC node
            # No restriction regarding ordering
            # Add the node to the list
            push!(node_list, node)


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