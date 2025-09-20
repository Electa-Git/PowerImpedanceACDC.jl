export make_y_node




function make_y_node(network::Network; nodelist = [], freq_range = (1,1e3, 1000))

node_list= Symbol[] #Node list to generate Ynode
element_list= Symbol[] #Element list to generate Ynode

if nodelist == []

        for elements in network.elements

            # Check whether element is an active element: converter, SG, source
            if is_passive(elements) 
                continue # Skip passives 
            end
            if is_source(elements) || is_converter(elements) || is_generator(elements)
                
                if is_converter(elements)

                    for (pin,node) in elements.pins # Iterate over all pins of the element


                        # Do AC DC pin ordering for MMC! See below copy paste

                         if occursin("2.", pin) 

                            if !in(node, node_list)
                                 push!(node_list,node) # Add the node to the list [ACd, ACq]
                            end
                            

                         end

                        end
                end

            end

        end

        # Create node list for Ynode
        for node in keys(network.nets)

        if occursin("gnd", string(node))
            continue # Skip ground nodes
        end

        # Check whether node is already part of the node list of Ynode
        if in(node, node_list)
            continue # Skip nodes that are already in the list
        end

        # Check whether the current node is part of Ynode

        nets=PowerImpedanceACDC.netfor!(network, node) # Get nets (designator, pin) connected to the node


        # Check whether the node is connected to a source
        # If it is, skip the node.
        # This assumes that AC & DC sources are only connected to an impedance 
        # TODO: Make this more general
        isSourceNode= false # Flag to skip the node if it is connected to a source
        for net in nets # Iterate over all elements connected to the node

            element=network.elements[net[1]] # Get the element, via the designator net[1]
            
            if is_source(element) 
                
                isSourceNode = true # Skip the node if it is connected to a source
                break
            end

        end

        if isSourceNode
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
            
            ACpin=replace(pin, ".1" => ".2") # Replace .1 with .2 to search for the other AC pin if existent

            designator = nets[1][1] # Get the designator of the element
            ACpin=Symbol(ACpin)

            if PowerImpedanceACDC.netname(network, (designator,ACpin)) === Symbol("") # Check whether this net (designator, pin) exists
            
                isDCnode= true

            else
                
                isACnode=true
            end

        end


        if isACnode 

            # Get the other AC node
            if occursin(".2", pin)
                ACpin=replace(pin, ".2" => ".1") # Replace .2 with .1 to search for the other AC pin 
            else
                ACpin=replace(pin, ".1" => ".2") # Replace .1 with .2 to search for the other AC pin 

            end

            designator = nets[1][1] # Get the designator of the element
            ACpin=Symbol(ACpin)
            node2=PowerImpedanceACDC.netname(network, (designator,ACpin))

            # Proper node ordering required to match 2x2 admittance matrix of elements [ACd, ACq]
            if occursin(".1", string(ACpin)) # Node 2 is a d node


                push!(node_list,node2) # Add the node to the list [ACd, ACq]
                push!(node_list,node)

            else  # Node 2 is a q node

                push!(node_list,node) # Add the node to the list [ACd, ACq]
                push!(node_list,node2)

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

            if isConverternode # Proper node ordering required to match 3x3 admittance matrix of converter [DC, ACd, ACq]

                triplet=Array{Union{Symbol}}(undef,3)
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
                        break

                    end
                end
                # Check whether the AC nodes are already in the node list
                # 
                index = findfirst(p -> p == triplet[2], node_list) 
                if index === nothing # AC nodes are not in the list add the entire triplet :) 
                    push!(node_list, triplet[1])
                    push!(node_list, triplet[2])
                    push!(node_list, triplet[3])
                else # AC nodes are in the list, place the DC node at the correction position
                    
                    insert!(node_list, index, triplet[1]) # Insert the DC node at the correct position
                end

            # If not a converter node, then it is a DC node
            # No restriction regarding ordering
            
            else
            
                # Add the node to the list
                push!(node_list, node)


            end


        end

        end

else

    # If the node list is given, use it
    node_list = nodelist

end



# 2. Element list creation
# Iterate over all active elements in the network
# Skip sources but include source impedance!
# Warning if not source connected to single impedance!

for (designator, element) in network.elements

    if is_passive(element) 
        continue # Skip passives 
    end

    # If not passive then SG, converter or source
    # If source, add (everything that is connected to the source to the element list)
    # If the connected element is a converter, then can skip it, as we only need the AC-side contribution,
    # when an ideal DC source is connected to the DC side
    if is_source(element) 
        for (pin,element_node) in element.pins # Search for the impedance 

            if occursin("gnd", string(element_node)) # Skip elements connected to ground
                continue
            end

            for net in netfor!(network,element_node)

                designator=net[1]
                if designator == element.symbol # Source itself
                    continue
                else
                    if !in(designator, element_list)

                        if !is_converter(network.elements[designator] ) # Only add to the element list when it is not a converter 
                            push!(element_list, designator) # Add the source to the list
                        end

                    end

                end
            end

        end

        
    else # Converters, SGs


        push!(element_list, designator) # Add the element to the list

    end


   
   

end

# Initialize dict to hold the nodes and elements for the admittance matrix
dict= Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[])
dict[:node_list]=node_list
dict[:element_list]=element_list


#TODO: Make sense of it
# make frequency range
(min_f, max_f, n_f) = freq_range
if !isa(n_f,Int)
    n_f = parse(Int, n_f) #Make Int to work with range (error when 1e4)
end
omegas= 2*pi* 10 .^range(log10(min_f), log10(max_f), length= n_f) 




Ynode=[] # Preallocate the admittance matrix for each frequency

for omega in omegas

    Y = make_y(network, dict, omega*1im)
    push!(Ynode,Y)


end

return Ynode,node_list


end