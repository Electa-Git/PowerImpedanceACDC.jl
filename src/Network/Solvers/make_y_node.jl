export make_y_node




function make_y_node(network::Network; freq_range = (1,1e3, 1000))


# 1. Node list creation
# Create node list for Ynode, similar to make_y_edge
# Either take nodelist as argument or create it from the network


# 2. Element list creation
# Iterate over all active elements in the network
# Skip sources but include source impedance!
# Warning if not source connected to single impedance!

# 3. Call make_y 









end