abstract type Machine end

function make_power_flow_ac!(machine :: Machine, dict :: Dict{String, Any},
    global_dict :: Dict{String, Any})
    key = string(length(dict["gen"]))
    ((dict["gen"])[string(key)])["pc1"] = 0
    ((dict["gen"])[string(key)])["pc2"] = 0
    ((dict["gen"])[string(key)])["qc1min"] = 0
    ((dict["gen"])[string(key)])["qc1max"] = 0
    ((dict["gen"])[string(key)])["qc2min"] = 0
    ((dict["gen"])[string(key)])["qc2max"] = 0
    ((dict["gen"])[string(key)])["ramp_agc"] = 0
    ((dict["gen"])[string(key)])["ramp_q"] = 0
    ((dict["gen"])[string(key)])["ramp_10"] = 0
    ((dict["gen"])[string(key)])["ramp_30"] = 0
    ((dict["gen"])[string(key)])["apf"] = 0
    ((dict["gen"])[string(key)])["startup"] = 0
    ((dict["gen"])[string(key)])["shutdown"] = 0

    ((dict["gen"])[string(key)])["gen_status"] = 1
    ((dict["gen"])[string(key)])["source_id"] = Any["gen", parse(Int, key)]
    ((dict["gen"])[string(key)])["index"] = parse(Int, key)

    S_base = global_dict["S"] / 1e6
    V_base = global_dict["V"] / 1e3
    ((dict["gen"])[string(key)])["pg"] = machine.P / S_base
    ((dict["gen"])[string(key)])["qg"] = machine.Q / S_base
    ((dict["gen"])[string(key)])["pmin"] = machine.P_min / S_base
    ((dict["gen"])[string(key)])["pmax"] = machine.P_max / S_base
    ((dict["gen"])[string(key)])["qmin"] = machine.Q_min / S_base
    ((dict["gen"])[string(key)])["qmax"] = machine.Q_max / S_base
    ((dict["gen"])[string(key)])["vg"] = machine.V / V_base

    # not using
    ((dict["gen"])[string(key)])["model"] = 1
    ((dict["gen"])[string(key)])["cost"] = 0
    ((dict["gen"])[string(key)])["ncost"] = 0

end

function make_power_flow_dc!(machine :: Machine, dict :: Dict{String, Any},
    global_dict :: Dict{String, Any})
nothing
end

function eval_abcd(machine :: Machine, s :: Complex)
    return eval_y(machine, s)
end

function eval_y(machine :: Machine, s :: Complex)
    Y = eval_parameters(machine, s)
    return Y
end

function make_power_flow!(machine:: Machine, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)

    # Check if AC or DC source (second one not implemented)
    # is_three_phase(elem) ? nothing : error("DC sources are currently not implemented")

    ### MAKE BUSES OUT OF THE NODES
    # Find the nodes not connected to the ground
    ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
    ac_nodes = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(elem.pins)))) #Look in the nodes of this component and convert into tuple
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes, global_dict)
    # Make busses for the non-ground nodes 
    interm_bus = add_interm_bus_ac!(data, global_dict) # No mapping to node, bcs no corresponding node in PowerImpedance

    # Make the generator component for injection
    key = injection_initialization!(data, elem2comp, comp2elem, interm_bus, elem, global_dict)
    key = string(key)

    # Add additional branch & bus for SM transformer (RL-branch)
    
    key_branch = length(data["branch"]) + 1
    key_branch_str = string(key_branch)

    (data["branch"])[key_branch_str] = Dict{String, Any}()
    ((data["branch"])[key_branch_str])["f_bus"] = interm_bus
    ((data["branch"])[key_branch_str])["t_bus"] = ac_bus
    ((data["branch"])[key_branch_str])["source_id"] = Any["branch", key_branch]
    ((data["branch"])[key_branch_str])["index"] = key_branch
    ((data["branch"])[key_branch_str])["rate_a"] = 1
    ((data["branch"])[key_branch_str])["rate_b"] = 1
    ((data["branch"])[key_branch_str])["rate_c"] = 1
    ((data["branch"])[key_branch_str])["br_status"] = 1
    ((data["branch"])[key_branch_str])["angmin"] = ang_min
    ((data["branch"])[key_branch_str])["angmax"] = ang_max
    ((data["branch"])[key_branch_str])["transformer"] = false
    ((data["branch"])[key_branch_str])["tap"] = 1
    ((data["branch"])[key_branch_str])["shift"] = 0
    ((data["branch"])[key_branch_str])["c_rating_a"] = 1

    
    ((data["branch"])[key_branch_str])["br_r"] = machine.rt * (machine.Vᵃᶜ_base^2 / machine.S_base) / global_dict["Z"]
    ((data["branch"])[key_branch_str])["br_x"] = machine.lt * (machine.Vᵃᶜ_base^2 / machine.S_base) / global_dict["Z"]
    ((data["branch"])[key_branch_str])["g_fr"] = 0
    ((data["branch"])[key_branch_str])["b_fr"] = 0
    ((data["branch"])[key_branch_str])["g_to"] = 0
    ((data["branch"])[key_branch_str])["b_to"] = 0

    # Change type of final bus, intermediate bus is PQ-bus
    if isapprox(machine.P_max, machine.P)
        ((data["bus"])[string(interm_bus)]) = set_bus_type((data["bus"])[string(interm_bus)], 1)
    else
        ((data["bus"])[string(interm_bus)]) = set_bus_type((data["bus"])[string(interm_bus)], 2)
    end

    ((data["bus"])[string(interm_bus)])["vm"] = ((data["gen"])[key])["vg"]
    ((data["bus"])[string(interm_bus)])["vmin"] =  0.9*((data["gen"])[key])["vg"]
    ((data["bus"])[string(interm_bus)])["vmax"] =  1.1*((data["gen"])[key])["vg"]
end

