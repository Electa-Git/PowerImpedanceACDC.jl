export power_flow, result, data
function power_flow(net:: Network)
    global ang_min, ang_max, result, nodes2bus, elem2comp, data
    global_dict = PowerModelsACDC.get_pu_bases(1000, net.voltageBase[1]) # 3-PH MVA, LL-RMS, Original setting was 100,320
    global_dict["omega"] = 2π * 50

    ang_min = deg2rad(360)
    ang_max = deg2rad(-360)

    nodes_dict = net.nets
    elem_dict = net.elements

    # PowerModelsACDC network dictoniary
    data = Dict{String, Any}()
    data = data_init(data, global_dict)
   
    ### 2-way dicts so we can have O(1) time complexity (node, elem:PowerImpedance ↔ bus, component:PowerModelsACDC)
    nodes2bus = Dict()
    bus2nodes = Dict() 
    elem2comp = Dict()
    comp2elem = Dict()

    ### Add grounds to the interfaces so we know for following (only do it once)
    ground_nodes = [k for k in keys(net.nets) if startswith(string(k), "gnd")] #TODO: add other ground identifiers (GND, Gnd, Ground, ground)
    push!(nodes2bus, ground_nodes => "gnd")
    push!(bus2nodes, "gnd" => ground_nodes) 
    
    #### 1. Create PowerModelsACDC dictionary and make interface 
    for (elem) in values(elem_dict)
        make_powerflow!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    end

    #### 2. Run PowerModelsACDC power flow
    PowerModelsACDC.process_additional_data!(data)
    ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false)
    result = run_acdcpf(data, ACPPowerModel, ipopt; setting = s)
    println(result["termination_status"])

    #### 3. Update setpoints of active elements
    for (key, element) in net.elements
        #Find the corresponding PowerModels component
        comp_type, key = elem2comp[element.symbol]
        elem_dict = result["solution"][comp_type][string(key)]
        pins = element.pins
        
        if is_converter(element) # In converter bus1 is DC and bus2 is AC
            node1 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "1.")]...)
            node2 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "2.")]...)
            _, dc_bus = nodes2bus[node1]
            _, ac_bus = nodes2bus[node2]
            
            Pdc = elem_dict["pdc"] * global_dict["S"] / 1e6
            Vm = (result["solution"]["bus"][string(ac_bus)]["vm"] * global_dict["V"] / 1e3) * sqrt(2) # Convert the LN-RMS voltage coming from the PF to LN-PK
            θ = result["solution"]["bus"][string(ac_bus)]["va"]
            Vdc = result["solution"]["busdc"][string(dc_bus)]["vm"] * global_dict["V"] / 1e3
            Pac = -elem_dict["pgrid"] * global_dict["S"] / 1e6
            Qac = elem_dict["qgrid"] * global_dict["S"] / 1e6 # Think about this!

            update!(element.element_value, Vm, θ, Pac, Qac, Vdc, Pdc)
            # if isa(element.element_value, MMC)
            #     update_string = "MMC #"
            #     update_mmc(element.element_value, Vm, θ, Pac, Qac, Vdc, Pdc)
            # else
            #     update_string = "TLC #"
            #     update_tlc(element.element_value, Vm, θ, Pac, Qac, Vdc, Pdc)
            # end
            update_string = string(key)
            print(update_string * " Active Power [MW]: ")
            println(Pac)
            print(update_string * " Reactive Power [MVar]: ")
            println(Qac)
            print(update_string * " AC Voltage Magnitude [pu]: ")
            println(result["solution"]["bus"][string(ac_bus)]["vm"])
            print(update_string * " AC Voltage Angle [rad]: ")
            println(θ)
            print(update_string * " DC Voltage [kV]: ")
            println(Vdc)
        elseif is_generator(element) #ac bus is the one with no ground in it's name
            ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
            ac_node = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(element.pins)))) #Look in the nodes of this component and convert into tuple
            bus_type, ac_bus = nodes2bus[ac_node]

            Pgen = elem_dict["pg"] * global_dict["S"] / 1e6 #MW
            Qgen = elem_dict["qg"] * global_dict["S"] / 1e6 #MVAr
            Vm = (result["solution"]["bus"][string(ac_bus)]["vm"] *
                    global_dict["V"] / 1e3) * sqrt(2) # Convert the LN-RMS voltage coming from the PF to LN-PK
            θ = result["solution"]["bus"][string(ac_bus)]["va"]
            update_string = string(key)

            update!(element.element_value, Pgen, Qgen, Vm, θ)
            
            print(update_string * " Active Power [MW]: ")
            println(Pgen)
            print(update_string * " Reactive Power [MVar]: ")
            println(Qgen)
            print(update_string * " AC Voltage Magnitude [pu]: ")
            println(result["solution"]["bus"][string(ac_bus )]["vm"])
            print(update_string * " AC Voltage Angle [rad]: ")
            println(θ)
        end
    end

    return result, data, nodes2bus, elem2comp

end

function set_bus_type(bus_data, type)
    ## This function makes sure we dont overwrite higher bus type (only for AC bus now)
    ## PQ (1) < PV (2) < slack (3)
    current_type = bus_data["bus_type"] #Just to make sure it is int
    if type > current_type
        bus_data["bus_type"] = type
    end
    return bus_data
end

## THis function makes sure we dispatch on the right component
make_powerflow!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict) = make_power_flow!(elem.element_value, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem,global_dict)

function injection_initialization!(data, elem2comp, comp2elem, ac_bus, elem, global_dict)
    ## A lot of initialization for source and machine are the same so combined in here
    

    ### ELEMENT TO COMPONENT
    # Interface
    # Interface element
    key = comp_elem_interface!(data, elem2comp, comp2elem, elem, "gen")
    key = string(key)

    # Network component
    (data["gen"])[key] = Dict{String, Any}()
    ((data["gen"])[key])["mBase"] = global_dict["S"] / 1e6
    ((data["gen"])[key])["gen_bus"] = ac_bus     
    ((data["gen"])[key])["pc1"] = 0
    ((data["gen"])[key])["pc2"] = 0
    ((data["gen"])[key])["qc1min"] = 0
    ((data["gen"])[key])["qc1max"] = 0
    ((data["gen"])[key])["qc2min"] = 0
    ((data["gen"])[key])["qc2max"] = 0
    ((data["gen"])[key])["ramp_agc"] = 0
    ((data["gen"])[key])["ramp_q"] = 0
    ((data["gen"])[key])["ramp_10"] = 0
    ((data["gen"])[key])["ramp_30"] = 0
    ((data["gen"])[key])["apf"] = 0
    ((data["gen"])[key])["startup"] = 0
    ((data["gen"])[key])["shutdown"] = 0

    ((data["gen"])[key])["gen_status"] = 1
    ((data["gen"])[key])["source_id"] = Any["gen", parse(Int, key)]
    ((data["gen"])[key])["index"] = parse(Int, key)

    injecter = elem.element_value
    S_base = global_dict["S"] / 1e6
    V_base = global_dict["V"] / 1e3
    ((data["gen"])[key])["pg"] = injecter.P / S_base
    ((data["gen"])[key])["qg"] = injecter.Q / S_base
    ((data["gen"])[key])["pmin"] = injecter.P_min / S_base
    ((data["gen"])[key])["pmax"] = injecter.P_max / S_base
    ((data["gen"])[key])["qmin"] = injecter.Q_min / S_base
    ((data["gen"])[key])["qmax"] = injecter.Q_max / S_base
    ((data["gen"])[key])["vg"] = injecter.V / V_base

    # not using
    ((data["gen"])[key])["model"] = 1
    ((data["gen"])[key])["cost"] = 0
    ((data["gen"])[key])["ncost"] = 0
    return parse(Int, key)
end
function branch_ac!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    # Add busses for the branch
    pins = elem.pins
    
    # Handle any amount of input and output pins
    node1 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "1.")]...) #(pins[Symbol(1.1)], pins[Symbol(1.2)]) 
    node2 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "2.")]...) #(pins[Symbol(2.1)], pins[Symbol(2.2)])
    bus1 = add_bus_ac!(data, nodes2bus, bus2nodes, node1, global_dict)
    bus2 = add_bus_ac!(data, nodes2bus, bus2nodes, node2, global_dict)

    # Interface element
    key_branch = comp_elem_interface!(data, elem2comp, comp2elem, elem, "branch")

    (data["branch"])[string(key_branch)] = Dict{String, Any}()
    ((data["branch"])[string(key_branch)])["f_bus"] = bus1
    ((data["branch"])[string(key_branch)])["t_bus"] = bus2
    ((data["branch"])[string(key_branch)])["source_id"] = Any["branch", key_branch]
    ((data["branch"])[string(key_branch)])["index"] = key_branch
    ((data["branch"])[string(key_branch)])["rate_a"] = 1
    ((data["branch"])[string(key_branch)])["rate_b"] = 1
    ((data["branch"])[string(key_branch)])["rate_c"] = 1
    ((data["branch"])[string(key_branch)])["br_status"] = 1
    ((data["branch"])[string(key_branch)])["angmin"] = ang_min
    ((data["branch"])[string(key_branch)])["angmax"] = ang_max
    return key_branch
end

function branch_dc!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    # Add busses for the branch
    pins = elem.pins
    
    node1 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "1.")]...) #pins[Symbol(1.1)] 
    node2 = tuple([pins[k] for k in sort(collect(keys(pins))) if startswith(string(k), "2.")]...) #pins[Symbol(2.1)]
    bus1 = add_bus_dc!(data, nodes2bus, bus2nodes, node1, global_dict)
    bus2 = add_bus_dc!(data, nodes2bus, bus2nodes, node2, global_dict)

    # Interface element
    key_branch = comp_elem_interface!(data, elem2comp, comp2elem, elem, "branchdc")

    (data["branchdc"])[string(key_branch)] = Dict{String, Any}()
    ((data["branchdc"])[string(key_branch)])["fbusdc"] = bus1
    ((data["branchdc"])[string(key_branch)])["tbusdc"] = bus2
    ((data["branchdc"])[string(key_branch)])["source_id"] = Any["branchdc", key_branch]
    ((data["branchdc"])[string(key_branch)])["index"] = key_branch
    ((data["branchdc"])[string(key_branch)])["rateA"] = 100
    ((data["branchdc"])[string(key_branch)])["rateB"] = 100
    ((data["branchdc"])[string(key_branch)])["rateC"] = 100
    ((data["branchdc"])[string(key_branch)])["status"] = 1

    #Copied these from element-specific definitions bcs everytime zero
    ((data["branchdc"])[string(key_branch)])["l"] = 0
    ((data["branchdc"])[string(key_branch)])["c"] = 0
    
    return key_branch
end





       

       
function comp_elem_interface!(data, elem2comp, comp2elem, elem, component)
    ## Making component
    key = length(data[component])+1
    push!(comp2elem, (component, key) => elem.symbol)
    push!(elem2comp, elem.symbol => (component, key))
    return key
end

function add_bus_dc!(data, nodes2bus, bus2nodes, node, global_dict)
    # Check if node is already in interface
    if node ∉ values(bus2nodes) 
        
        # Add bus to PowerModelsACDC dict       
        bus = length(data["busdc"]) + 1 # We make new DC bus
        #Update interface
        push!(bus2nodes, ("busdc", bus) => node)
        push!(nodes2bus, node => ("busdc", bus))
        bus = string(bus)
        #UPdata network dict
        (data["busdc"])[bus] = Dict{String, Any}()
        ((data["busdc"])[bus])["busdc_i"] = parse(Int, bus)
        ((data["busdc"])[bus])["source_id"] = Any["busdc", parse(Int, bus)]
        ((data["busdc"])[bus])["grid"] = 1
        ((data["busdc"])[bus])["index"] = parse(Int, bus)
        ((data["busdc"])[bus])["Cdc"] = 0
        ((data["busdc"])[bus])["Vdc"] = 1
        ((data["busdc"])[bus])["Vdcmax"] = 1.1
        ((data["busdc"])[bus])["Vdcmin"] = 0.9
        ((data["busdc"])[bus])["Pdc"] = 0
        ((data["busdc"])[bus])["basekVdc"] = global_dict["V"] / 1e3
        bus = parse(Int, bus)
    else
        #Return bus of this node, nodes correspond to one bus so no risk of same values with different keys
        _, bus = nodes2bus[node]
    end

    return bus
end

function add_bus_ac!(data, nodes2bus, bus2nodes, node, global_dict)

    if node ∉ values(bus2nodes)
        bus = length(data["bus"]) + 1# We make new DC bus, key should be string apparently
        
        #Update interface
        push!(bus2nodes, ("bus", bus) => node)
        push!(nodes2bus, node => ("bus", bus))
        bus = string(bus)

        (data["bus"])[bus] = Dict{String, Any}()
        ((data["bus"])[bus])["source_id"] = Any["bus", parse(Int, bus)]
        ((data["bus"])[bus])["index"] = parse(Int, bus)
        ((data["bus"])[bus])["bus_i"] = parse(Int, bus)
        ((data["bus"])[bus])["zone"] = 1
        ((data["bus"])[bus])["area"] = 1
        ((data["bus"])[bus])["vmin"] = 0.9
        ((data["bus"])[bus])["vmax"] = 1.1
        ((data["bus"])[bus])["vm"] = 1
        ((data["bus"])[bus])["va"] = 0
        ((data["bus"])[bus])["base_kv"] = global_dict["V"] / 1e3
        ((data["bus"])[bus])["bus_type"] = 1 # bus type - depends on components 1 is default PQ
        bus = parse(Int, bus)
    else
        #Return bus of this node, nodes correspond to one bus so no risk of same values with different keys
        _, bus = nodes2bus[node]
    end

    return bus
end
"""
    function add_interm_bus_ac!(data)

Add intermediate AC bus to PowerModels data, which has no corresponding node(s). 
An example is the Synchronous Machine wich has an intermediate bus, connecting the generator and RL-branch (transformer)
"""
function add_interm_bus_ac!(data, global_dict)

    bus = length(data["bus"]) + 1# We make new DC bus, key should be string apparently
    bus = string(bus)

    (data["bus"])[bus] = Dict{String, Any}()
    ((data["bus"])[bus])["source_id"] = Any["bus", parse(Int, bus)]
    ((data["bus"])[bus])["index"] = parse(Int, bus)
    ((data["bus"])[bus])["bus_i"] = parse(Int, bus)
    ((data["bus"])[bus])["zone"] = 1
    ((data["bus"])[bus])["area"] = 1
    ((data["bus"])[bus])["vmin"] = 0.9
    ((data["bus"])[bus])["vmax"] = 1.1
    ((data["bus"])[bus])["vm"] = 1
    ((data["bus"])[bus])["va"] = 0
    ((data["bus"])[bus])["base_kv"] = global_dict["V"] / 1e3
    ((data["bus"])[bus])["bus_type"] = 1 # bus type - depends on components 1 is default PQ
    bus = parse(Int, bus)
    return bus
end

function data_init(data, global_dict)
    data["source_type"] = "matpower"
    data["name"] = "network"
    data["source_version"] = "0.0.0"
    data["per_unit"] = true
    data["dcpol"] = 2 # bipolar converter topologym check in the future
    data["baseMVA"] = global_dict["S"] / 1e6
    data["bus"] = Dict{String, Any}()
    data["busdc"] = Dict{String, Any}()
    data["shunt"] = Dict{String, Any}()     # empty
    data["dcline"] = Dict{String, Any}()    # empty
    data["storage"] = Dict{String, Any}()   # empty
    data["switch"] = Dict{String, Any}()    # empty
    data["load"] = Dict{String, Any}()      # empty
    data["branch"] = Dict{String, Any}()
    data["branchdc"] = Dict{String, Any}()
    data["gen"] = Dict{String, Any}()
    data["convdc"] = Dict{String, Any}()
    return data
end