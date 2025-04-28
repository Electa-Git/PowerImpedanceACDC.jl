
@with_kw mutable struct Source
    Z :: Union{Float64, Int, Basic} = 0 # source series impedance [Ω]
    V :: Union{Float64, Int} = 0        # DC voltage or voltage magnitude [kV]

    P   :: Union{Float64, Int} = 0      # active power output [MW]
    Q   :: Union{Float64, Int} = 0      # reactive power output [MVAr]
    P_min :: Union{Float64, Int} = 0    # min active power output [MW]
    P_max :: Union{Float64, Int} = 0    # max active power output [MW]
    Q_min :: Union{Float64, Int} = 0    # min reactive power output [MVA]
    Q_max :: Union{Float64, Int} = 0    # max reactive power output [MVA]

    pins :: Int = 1
    ABCD :: Array{Basic} = Basic[]
end

function make_abcd(source :: Source)
    source.ABCD = convert(Array{Basic}, Diagonal([1 for dummy in 1:2source.pins]))
    source.ABCD[1:source.pins, source.pins+1:end] = convert(Array{Basic}, Diagonal([source.Z for dummy in 1:source.pins]))
end

function eval_abcd(source :: Source, s :: Complex)
    abcd = N.(source.ABCD)
    abcd = convert(Array{Float64}, real(abcd)) + 1im*convert(Array{Float64}, imag(abcd))
    abcd = convert(Array{Complex}, abcd)
end

function eval_y(source :: Source, s :: Complex)
    abcd = eval_abcd(source, s)
    if source.Z == 0
        abcd[1:source.pins, source.pins+1:end] = 1e-6 * Diagonal([1 for i in 1:source.pins])
    end
    return abcd_to_y(abcd)
end

function make_power_flow_dc!(source :: Source, dict :: Dict{String, Any},
        global_dict :: Dict{String, Any})
    nothing
end

function make_power_flow_ac!(source :: Source, dict :: Dict{String, Any},
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
    ((dict["gen"])[string(key)])["pg"] = source.P / S_base
    ((dict["gen"])[string(key)])["qg"] = source.Q / S_base
    ((dict["gen"])[string(key)])["pmin"] = source.P_min / S_base
    ((dict["gen"])[string(key)])["pmax"] = source.P_max / S_base
    ((dict["gen"])[string(key)])["qmin"] = source.Q_min / S_base
    ((dict["gen"])[string(key)])["qmax"] = source.Q_max / S_base
    ((dict["gen"])[string(key)])["vg"] = source.V / V_base

    # not using
    ((dict["gen"])[string(key)])["model"] = 1
    ((dict["gen"])[string(key)])["cost"] = 0
    ((dict["gen"])[string(key)])["ncost"] = 0
end
function make_power_flow!(source:: Source, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)
    
    # Check if AC or DC source (second one not implemented)
    is_three_phase(elem) ? nothing : error("DC sources are currently not implemented")

    ### MAKE BUSES OUT OF THE NODES
    # Find the nodes not connected to the ground
    ground_nodes = Set(bus2nodes["gnd"]) #Collect ground nodes and make them a set for faster lookup
    ac_nodes = Tuple(collect(Iterators.filter(x -> !(x in ground_nodes), values(elem.pins)))) #Look in the nodes of this component and convert into tuple
    
    # Make busses for the non-ground nodes 
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes)

    # Make the generator component for injection
    key = injection_initialization!(data, elem2comp, comp2elem, ac_bus, elem)
    key = string(key) # Of form "gen", 1 so convert 2nd element to string

    # Change bus information
    ((data["bus"])[string(ac_bus)])["vmin"] =  0.9*((data["gen"])[key])["vg"]
    ((data["bus"])[string(ac_bus)])["vmax"] =  1.1*((data["gen"])[key])["vg"]
    ((data["bus"])[string(ac_bus)])["vm"] = ((data["gen"])[key])["vg"]

    ((data["bus"])[string(ac_bus)]) = set_bus_type((data["bus"])[string(ac_bus)], 3)
end