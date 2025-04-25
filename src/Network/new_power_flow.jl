export new_power_flow
function new_power_flow(net:: Network)
    global global_dict
    global_dict = PowerModelsACDC.get_pu_bases(1000, net.voltageBase[1]) # 3-PH MVA, LL-RMS, Original setting was 100,320
    global_dict["omega"] = 2π * 50

    ang_min = deg2rad(360)
    ang_max = deg2rad(-360)

    nodes_dict = net.nets
    elem_dict = net.elements

    # PowerModelsACDC network dictoniary
    data = Dict{String, Any}()
    data = data_init(data)
   
    ### 2-way dicts so we can have O(1) time complexity (node, elem:PowerImpedance ↔ bus, component:PowerModelsACDC)
    nodes2bus = Dict()
    bus2nodes = Dict() 
    elem2comp = Dict()
    comp2elem = Dict()

    ### Add grounds to the interfaces so we know for following (only do it once)
    ground_nodes = [k for k in keys(net.nets) if startswith(string(k), "gnd")] #TODO: add other ground identifiers (GND, Gnd, Ground, ground)
    push!(nodes2bus, ground_nodes => "gnd")
    push!(bus2nodes, "gnd" => ground_nodes) 
    #### 1. Create interface between PowerImp element and PowerModels dict key
    for (elem) in values(elem_dict)
        make_interface!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    end
    PowerModelsACDC.process_additional_data!(data)

    return data, nodes2bus, bus2nodes, elem2comp, comp2elem

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
make_interface!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict) = _interface!(elem.element_value, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem,global_dict)

function injection_initialization!(data, elem2comp, comp2elem, ac_bus, elem)
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
end
function branch_ac!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)
    
    # Add busses for the branch
    pins = elem.pins
    
    node1 = (pins[Symbol(1.1)], pins[Symbol(1.2)]) 
    node2 = (pins[Symbol(2.1)], pins[Symbol(2.2)])
    bus1 = add_bus_ac!(data, nodes2bus, bus2nodes, node1)
    bus2 = add_bus_ac!(data, nodes2bus, bus2nodes, node2)

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

function branch_dc!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)
    # Add busses for the branch
    pins = elem.pins
    
    node1 = pins[Symbol(1.1)] 
    node2 = pins[Symbol(2.1)]
    bus1 = add_bus_dc!(data, nodes2bus, bus2nodes, node1)
    bus2 = add_bus_dc!(data, nodes2bus, bus2nodes, node2)

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
function  _interface!(t :: Transformer, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    # Initialize an AC branch between both nodes
    key_branch = branch_ac!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)

    # Add transformer data)
    ((data["branch"])[string(key_branch)])["transformer"] = true

    ((data["branch"])[string(key_branch)])["shift"] = 0
    ((data["branch"])[string(key_branch)])["c_rating_a"] = 1

    abcd = eval_abcd(t, global_dict["omega"] * 1im)
    n = 3
    (a, b, c, d) = (abcd[1:n,1:n], abcd[1:n, n+1:end], abcd[n+1:end,1:n], abcd[n+1:end,n+1:end])
    Y = [d*inv(b) c-d*inv(b)*a; -inv(b) a*inv(b)] * global_dict["Z"]

    tap = sqrt(real(Y[n+1,n+1]/Y[1,1]))
    ys = -Y[1,n+1]*tap
    yc = Y[n+1,n+1] - ys

    ((data["branch"])[string(key_branch)])["tap"] = tap
    ((data["branch"])[string(key_branch)])["br_r"] = real(1/ys)
    ((data["branch"])[string(key_branch)])["br_x"] = imag(1/ys)
    ((data["branch"])[string(key_branch)])["g_fr"] = real(yc)/2
    ((data["branch"])[string(key_branch)])["b_fr"] = imag(yc)/2
    ((data["branch"])[string(key_branch)])["g_to"] = real(yc)/2
    ((data["branch"])[string(key_branch)])["b_to"] = imag(yc)/2  
end

function  _interface!(imp:: Impedance, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    if is_three_phase(elem)
        if is_load(elem) #This means it's a ground connected impedance --> shunt impedance
            ### MAKE BUSES OUT OF THE NODES
            # Find the nodes not connected to the ground
            ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
            ac_nodes = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(elem.pins)))) #Look in the nodes of this component and convert into tuple
        
            ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes)
            key = comp_elem_interface!(data, elem2comp, comp2elem, elem, "shunt")

            (data["shunt"])[string(key)] = Dict{String, Any}()
            ((data["shunt"])[string(key)])["source_id"] = Any["bus", ac_bus]
            ((data["shunt"])[string(key)])["index"] = key
            ((data["shunt"])[string(key)])["shunt_bus"]  = ac_bus
            data["shunt"][string(key)]["status"] = 1

            abcd = eval_abcd(imp, global_dict["omega"] * 1im)
            n = 3
            Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]
            data["shunt"][string(key)]["gs"] = real(1/Z)
            data["shunt"][string(key)]["bs"] = imag(1/Z)
        else
            # Initialize an AC branch between both nodes
            key = branch_ac!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)
            ((data["branch"])[string(key)])["transformer"] = false
            ((data["branch"])[string(key)])["tap"] = 1
            ((data["branch"])[string(key)])["shift"] = 0
            ((data["branch"])[string(key)])["c_rating_a"] = 1

            abcd = eval_abcd(imp, global_dict["omega"] * 1im)
            n = 3
            Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"] # Assuming impedance with equal values for all phases 
            ((data["branch"])[string(key)])["br_r"] = real(Z)
            ((data["branch"])[string(key)])["br_x"] = imag(Z)
            ((data["branch"])[string(key)])["g_fr"] = 0
            ((data["branch"])[string(key)])["b_fr"] = 0
            ((data["branch"])[string(key)])["g_to"] = 0
            ((data["branch"])[string(key)])["b_to"] = 0
        end
    else
        ## DC impedance
        key = branch_dc!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)

        abcd = eval_abcd(imp, 1e-6*1im)
        Z = abcd[1,2] / global_dict["Z"]
        ((data["branchdc"])[string(key)])["r"] = real(Z)
    end
    
end

function _interface!(tl :: Transmission_line, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    pins = elem.pins
 
    if is_three_phase(elem)    
        key = branch_ac!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)
        ((data["branch"])[string(key)])["transformer"] = false
        ((data["branch"])[string(key)])["tap"] = 1
        ((data["branch"])[string(key)])["shift"] = 0
        ((data["branch"])[string(key)])["c_rating_a"] = 1

        abcd = eval_abcd(tl, global_dict["omega"] * 1im)
        n = Int(size(abcd, 1)/2)
        Z_ph = (abcd[1:n,n+1:end]) / global_dict["Z"] # Phase domain impedance data
        T_seq = [1 1 1;1 exp(2*pi/3im) exp(4*pi/3im);1 exp(4*pi/3im) exp(2*pi/3im)]/sqrt(3) # Transformation matrix for sequence domain
        Z = (inv(T_seq) * Z_ph * T_seq)[2,2] # Taking the positive sequence impedance
        Y = (abcd[n+1:end,1:n] * inv(abcd[n+1:end,n+1:end]))[1,1] * global_dict["Z"]

        ((data["branch"])[string(key)])["br_r"] = real(Z)
        ((data["branch"])[string(key)])["br_x"] = imag(Z)
        ((data["branch"])[string(key)])["g_fr"] = real(Y)/2
        ((data["branch"])[string(key)])["b_fr"] = imag(Y)/2
        ((data["branch"])[string(key)])["g_to"] = real(Y)/2
        ((data["branch"])[string(key)])["b_to"] = imag(Y)/2

    else
        key = branch_dc!(data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem)
        abcd = eval_abcd(tl, 1e-6*1im)
        n = Int(size(abcd, 1)/2)
        Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]
        Y = (abcd[n+1:end,1:n] * inv(abcd[n+1:end,n+1:end]))[1,1] * global_dict["Z"]
        ((data["branchdc"])[string(key)])["r"] = real(Z)
    end


end

function _interface!(source:: Source, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    # Check if AC or DC source (second one not implemented)
    is_three_phase(elem) ? nothing : error("DC sources are currently not implemented")

    ### MAKE BUSES OUT OF THE NODES
    # Find the nodes not connected to the ground
    ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
    ac_nodes = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(elem.pins)))) #Look in the nodes of this component and convert into tuple
    
    # Make busses for the non-ground nodes 
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes)

    # Make the generator component for injection
    injection_initialization!(data, elem2comp, comp2elem, ac_bus, elem)
    key = string((elem2comp[elem.symbol])[2]) # Of form "gen", 1 so convert 2nd element to string

    # Change bus information
    ((data["bus"])[string(ac_bus)])["vmin"] =  0.9*((data["gen"])[key])["vg"]
    ((data["bus"])[string(ac_bus)])["vmax"] =  1.1*((data["gen"])[key])["vg"]
    ((data["bus"])[string(ac_bus)])["vm"] = ((data["gen"])[key])["vg"]

    ((data["bus"])[string(ac_bus)]) = set_bus_type((data["bus"])[string(ac_bus)], 3)
end
       
function _interface!(machine:: Machine, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)

    # Check if AC or DC source (second one not implemented)
    is_three_phase(elem) ? nothing : error("DC sources are currently not implemented")

    ### MAKE BUSES OUT OF THE NODES
    # Find the nodes not connected to the ground
    ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
    ac_nodes = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(elem.pins)))) #Look in the nodes of this component and convert into tuple
    
    # Make busses for the non-ground nodes 
    interm_bus = add_interm_bus_ac!(data) # No mapping to node, bcs no corresponding node in PowerImpedance

    # Make the generator component for injection
    injection_initialization!(data, elem2comp, comp2elem, interm_bus, elem)
    key = string(elem2comp(elem.symbol))

    # Add additional branch & bus for SM transformer (RL-branch)
    key_branch = String(length(data["branch"])+1)
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes)
    # Interface element
    key_branch = comp_elem_interface!(data, elem2comp, comp2elem, elem, "branch")
    key_branch = string(key_branch)

    (data["branch"])[key_branch] = Dict{String, Any}()
    ((data["branch"])[key_branch])["f_bus"] = interm_bus
    ((data["branch"])[key_branch])["t_bus"] = ac_bus
    ((data["branch"])[key_branch])["source_id"] = Any["branch", parse(Int, key_branch)]
    ((data["branch"])[key_branch])["index"] = parse(Int, key_branch)
    ((data["branch"])[key_branch])["rate_a"] = 1
    ((data["branch"])[key_branch])["rate_b"] = 1
    ((data["branch"])[key_branch])["rate_c"] = 1
    ((data["branch"])[key_branch])["br_status"] = 1
    ((data["branch"])[key_branch])["angmin"] = ang_min
    ((data["branch"])[key_branch])["angmax"] = ang_max
    ((data["branch"])[key_branch])["transformer"] = false
    ((data["branch"])[key_branch])["tap"] = 1
    ((data["branch"])[key_branch])["shift"] = 0
    ((data["branch"])[key_branch])["c_rating_a"] = 1

    
    ((data["branch"])[key_branch])["br_r"] = machine.rt * (machine.Vᵃᶜ_base^2 / machine.S_base) / global_dict["Z"]
    ((data["branch"])[key_branch])["br_x"] = machine.lt * (machine.Vᵃᶜ_base^2 / machine.S_base) / global_dict["Z"]
    ((data["branch"])[key_branch])["g_fr"] = 0
    ((data["branch"])[key_branch])["b_fr"] = 0
    ((data["branch"])[key_branch])["g_to"] = 0
    ((data["branch"])[key_branch])["b_to"] = 0

    # Change type of final bus, intermediate bus is PQ-bus
    if isapprox(source_machine.P_max, source_machine.P)
        ((data["bus"])[string(ac_bus)]) = set_bus_type((data["bus"])[string(ac_bus)], 1)
    else
        ((data["bus"])[string(ac_bus)]) = set_bus_type((data["bus"])[string(ac_bus)], 2)
    end

    ((data["bus"])[string(ac_bus)])["vm"] = ((data["gen"])[key_branch])["vg"]
    ((data["bus"])[string(ac_bus)])["vmin"] =  0.9*((data["gen"])[key_branch])["vg"]
    ((data["bus"])[string(ac_bus)])["vmax"] =  1.1*((data["gen"])[key_branch])["vg"]
end
       
function comp_elem_interface!(data, elem2comp, comp2elem, elem, component)
    ## Making component
    key = length(data[component])+1
    push!(comp2elem, (component, key) => elem.symbol)
    push!(elem2comp, elem.symbol => (component, key))
    return key
end
function _interface!(converter:: Converter, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    pins = elem.pins

    # Busses interface (dc_bus --> 1.1 & ac_bus --> 2.1 and 2.2)
    dc_node = pins[Symbol(1.1)] # This is the node wherefore nodes_dict gives us all the pins, we have to map this to a dc_bus
    ac_nodes = (pins[Symbol(2.1)], pins[Symbol(2.2)]) #Similar AC bus
    dc_bus = add_bus_dc!(data, nodes2bus, bus2nodes, dc_node)
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes)
    
    # Interface element
    key = comp_elem_interface!(data, elem2comp, comp2elem, elem, "convdc")

    (data["convdc"])[string(key)] = Dict{String, Any}()
    ((data["convdc"])[string(key)])["busdc_i"] = string(dc_bus)
    ((data["convdc"])[string(key)])["busac_i"] = string(ac_bus)
    ((data["convdc"])[string(key)])["source_id"] = Any["convdc", key]
    ((data["convdc"])[string(key)])["status"] = 1
    ((data["convdc"])[string(key)])["index"] = key

    ((data["convdc"])[string(key)])["basekVac"] = global_dict["V"] / 1e3

    if in(:vac, keys(converter.controls)) || in(:vac_supp, keys(converter.controls)) 
        ((data["convdc"])[string(key)])["type_ac"] = 2  # PV ac bus
        data["bus"][string(ac_bus)] = set_bus_type(data["bus"][string(ac_bus)], 2)
        # TODO: The line below sometimes gives errors during power flow (NUMERICAL_ERROR)
         # Not entirely sure if this is necessary.
        if in(:vac, keys(converter.controls))
            ((data["convdc"])[string(key)])["Vtar"] = converter.controls[:vac].ref[1] * 1e3 / (global_dict["V"] * sqrt(2))
        else
            ((data["convdc"])[string(key)])["Vtar"] = converter.controls[:vac_supp].ref[1] * 1e3 / (global_dict["V"] * sqrt(2))
        end
        data["bus"][string(ac_bus)]["vm"] = ((data["convdc"])[string(key)])["Vtar"]
    else
        ((data["convdc"])[string(key)])["type_ac"] = 1  # PQ ac bus TODO: Check if this makes sense in the presence of GFM
        ((data["convdc"])[string(key)])["Vtar"] = converter.Vₘ * 1e3 / global_dict["V"]
    end
    if in(:p, keys(converter.controls))
        ((data["convdc"])[string(key)])["type_dc"] = 1  # constant AC active power        
    elseif in(:dc, keys(converter.controls))
        ((data["convdc"])[string(key)])["type_dc"] = 2  # constant DC voltage
    else
        ((data["convdc"])[string(key)])["type_dc"] = 3  # DC voltage droop
    end

    if in(:vac_supp, keys(converter.controls))
        ((data["convdc"])[string(key)])["acq_droop"] = 1
        ((data["convdc"])[string(key)])["kq_droop"] = converter.controls[:vac_supp].Kₚ
    else
        ((data["convdc"])[string(key)])["acq_droop"] = 0
        ((data["convdc"])[string(key)])["kq_droop"] = 0
    end

    # droop control - not implemented
    ((data["convdc"])[string(key)])["droop"] = 0
    ((data["convdc"])[string(key)])["Pdcset"] = converter.P_dc
    ((data["convdc"])[string(key)])["Vdcset"] = converter.Vᵈᶜ * 1e3 / global_dict["V"]
    ((data["convdc"])[string(key)])["dVdcSet"] = 0
    # LCC converter
    ((data["convdc"])[string(key)])["islcc"] = 0

    # without transformer
    ((data["convdc"])[string(key)])["transformer"] = 0
    ((data["convdc"])[string(key)])["rtf"] = 0
    ((data["convdc"])[string(key)])["xtf"] = 0
    ((data["convdc"])[string(key)])["tm"] = 1
    # without filter
    ((data["convdc"])[string(key)])["filter"] = 0
    ((data["convdc"])[string(key)])["bf"] = 0
    # with reactor
    ((data["convdc"])[string(key)])["reactor"] = 1
    ((data["convdc"])[string(key)])["rc"] = (converter.Rᵣ + converter.Rₐᵣₘ / 2) / global_dict["Z"]
    ((data["convdc"])[string(key)])["xc"] = (converter.Lᵣ + converter.Lₐᵣₘ / 2) * global_dict["omega"] / global_dict["Z"]
    converter.ω₀ = global_dict["omega"]

    # default values
    ((data["convdc"])[string(key)])["Vmmax"] = 1.1 * converter.Vₘ * 1e3 / global_dict["V"]
    ((data["convdc"])[string(key)])["Vmmin"] = 0.9 * converter.Vₘ * 1e3 / global_dict["V"]
    ((data["convdc"])[string(key)])["Imax"] = 1.1 * abs(converter.P_max) / converter.Vₘ

    ((data["convdc"])[string(key)])["P_g"] = converter.P
    ((data["convdc"])[string(key)])["Q_g"] = converter.Q
    

    ((data["convdc"])[string(key)])["LossA"] = 0
    ((data["convdc"])[string(key)])["LossB"] = 0
    ((data["convdc"])[string(key)])["LossCrec"] = converter.Rₐᵣₘ / 2
    ((data["convdc"])[string(key)])["LossCinv"] = converter.Rₐᵣₘ / 2

    ((data["convdc"])[string(key)])["Qacmax"] = converter.Q_max
    ((data["convdc"])[string(key)])["Qacmin"] = converter.Q_min
    ((data["convdc"])[string(key)])["Pacmax"] = converter.P_max
    ((data["convdc"])[string(key)])["Pacmin"] = converter.P_min

    
    
    #Voltage limits for bus, these are already set when PV bus
    if (data["bus"][string(ac_bus)]["bus_type"] == 1) #PQ-bus
        data["bus"][string(ac_bus)]["vm"] = ((data["convdc"])[string(key)])["Vtar"]
        ((data["bus"])[string(ac_bus)])["vmin"] = 0.9 * data["bus"][string(ac_bus)]["vm"]
        ((data["bus"])[string(ac_bus)])["vmax"] = 1.1 * data["bus"][string(ac_bus)]["vm"]
    end
    ((data["busdc"])[string(dc_bus)])["Vdc"] = converter.Vᵈᶜ * 1e3 / global_dict["V"]
    ((data["busdc"])[string(dc_bus)])["Vdcmax"] = 1.1 * ((data["busdc"])[string(dc_bus)])["Vdc"]
    ((data["busdc"])[string(dc_bus)])["Vdcmin"] = 0.9 * ((data["busdc"])[string(dc_bus)])["Vdc"]
    
end

function add_bus_dc!(data, nodes2bus, bus2nodes, node)
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

    else
        #Return bus of this node, nodes correspond to one bus so no risk of same values with different keys
        _, bus = nodes2bus[node]
    end

    return bus
end

function add_bus_ac!(data, nodes2bus, bus2nodes, node)

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

    else
        #Return bus of this node, nodes correspond to one bus so no risk of same values with different keys
        _, bus = nodes2bus[node]
    end

    return bus
end

function add_interm_bus_ac!(data)

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

    return bus
end

function data_init(data)
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