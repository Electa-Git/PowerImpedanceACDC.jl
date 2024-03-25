# Main network definition

export Network, add!, connect!, disconnect!, @network,
        composite_element
export power_flow  # for testing

import Base: delete!

const Net = Vector{Tuple{Symbol,Symbol}} # pairs of element designator and pin name

"""
```
struct Network
    elements::OrderedDict{Symbol, Element}
    nets :: Dict{Symbol, Net}
    connections :: Dict{Symbol, Net}
    Network() = new(OrderedDict{Symbol, Element}(), Dict{Symbol, Net}(), Dict{Symbol, Vector{Int}}())
end
```
"""
struct Network
    elements::OrderedDict{Symbol, Element}
    nets :: Dict{Symbol, Net}
    connections :: Dict{Symbol, Net}
    Network() = new(OrderedDict{Symbol, Element}(), Dict{Symbol, Net}(), Dict{Symbol, Vector{Int}}())
end

# adding elements
"""
    add!(n::Network, elem::Element)
Adds the element `elem` to the network `n`, creating and returning a new, unique
reference designator, leaving its pins unconnected. #FP: It will be connected with the pin connection at the end of the code (usually at the end)
"""
function add!(n::Network, elem::Element)
    for (k, v) in n.elements #n.elements =
        if v == elem
            return k
        end
    end
    designator = gensym() #gensym-> generates a symbol that will not conflict with the other variable names
    add!(n, designator, elem) #add the element elem to the network n, with the reference designator: designator
    return designator
end

add!(n::Network, elems::Element...) = ((add!(n, elem) for elem in elems)...,)

"""
    add!(n::Network, designator::Symbol, elem::Element)
Adds the element `elem` to the network `n` with the reference designator
`designator`, leaving its pins unconnected. If the network already contained
an element named `designator`, it is removed first.
"""
function add!(n::Network, designator::Symbol, elem::Element)
    if haskey(n.elements, designator) #haskey -> determine whether a collection (n.elements) has a mapping for a given key (designator)
        delete!(n, designator) #delete the element named designator from the network n (disconnecting all its pins)
    end
    for pin in keys(elem.pins) #keys elem.pins returns an array of keys. Return an iterator over all keys in a dictionary
        add!(n, (designator, pin))
    end
    add!(elem, :symbol, designator)
    n.elements[designator] = elem
end

# adding pins and nets
"""
    add!(n::Network, pin::Tuple{Symbol, Symbol})
Adds the pin `pin` to the nets in `n.nets` with the generic name.
"""
function add!(n::Network, pin::Tuple{Symbol, Symbol})
    for (name, elem_pins) in n.nets
        if (pin in elem_pins)
            return
        end
    end
    n.nets[gensym()] = [pin]
end
add!(n::Network, p::Tuple{Symbol,Any}) = add!(n, (p[1], Symbol(p[2])))

function add!(n::Network, name::Symbol, pins::Union{Tuple{Symbol,Any}}...)
    n.nets[name] = []
    for pin in pins
        append!(n.nets[name], pin)
    end
end

"""
    delete!(n::Network, designator::Symbol)
Deletes the element named `designator` from the network `n` (disconnecting all
its pins).
"""
function delete!(n::Network, designator::Symbol)
    for (sym,net) in n.nets
        filter!(elempin -> elempin[1] != designator, net)
    end
    delete!(n.elements, designator)
end

function netfor!(n::Network, p::Tuple{Symbol,Symbol})
    for (name, net) in n.nets
        p ∈ net && return net
    end
    throw(ArgumentError("Unknown pin $p"))
end
netfor!(n::Network, p::Tuple{Symbol,Any}) = netfor!(n, (p[1], Symbol(p[2])))

function netfor!(n::Network, name::Symbol)
    if !haskey(n.nets, name)
        n.nets[name] = []
    end
    n.nets[name]
end

# get net name
function netname(n::Network, name::Symbol)
    if haskey(n.nets, name)
        return name
    else
        throw(ArgumentError("Unknown net name $name."))
    end
end

function netname(n::Network, pin::Tuple{Symbol,Symbol})
    for (name, pins) in n.nets
        pin ∈ pins && return name
    end
    #throw(ArgumentError("Unknown pin $pin."))
    return Symbol()
end
netname(n::Network, pin::Tuple{Symbol, Any}) = netname(n, Tuple(pin[1], Symbol(pin[2])))

function netname(n::Network, pins::Union{Symbol,Tuple{Symbol, Symbol}}...)
    for (name, net_pins) in n.nets
        all((isa(pin, Symbol) && (pin == name)) || (pin ∈ net_pins) for pin ∈ pins) && return name
    end
    #throw(ArgumentError("Unknown net connected to pins $pins."))
    return Symbol()
end

function netname(n::Network, pins::Array{Tuple{Symbol, Symbol}})
    for (name, net_pins) in n.nets
        all((pin ∈ net_pins) for pin ∈ pins) &&
            length(pins) > 0 && return name
    end
    return Symbol()
end

"""
    connect!(n::Network)
Connects all elements' pins with their node names.
"""
function connect!(n::Network)
    for net in keys(n.nets)
        for pin in n.nets[net]
            n.elements[pin[1]].pins[pin[2]] = net
        end
    end
end

@doc doc"""
    connect!(n::Network, pins::Union{Symbol,Tuple{Symbol,Any}}...)
Connects the given pins (or named nets) to each other in the network `n`. Named
nets are given as `Symbol`s, pins are given as `Tuple{Symbols,Any}`s, where the
first entry is the reference designator of an element in `c`, and the second
entry is the pin name. For convenience, the latter is automatically converted to
a `Symbol` as needed.
# Example
```jldoctest; output = false, setup = :(include("../src/HVDCstability.jl"); using .HVDCstability), filter = r"(HVDCstability\.)?Network\(.*"s
network = Network()
add!(network, :r, impedance(z = 1e3, pins = 1))
add!(network, :src, dc_source(V = 5))
connect!(network, (:src, 2.1), (:r, 2.1), :gnd) # connect to gnd net
network
# output
Network(...)
```
"""
function connect!(n::Network, pins::Union{Symbol,Tuple{Symbol,Any}}...)
    nets = []
    for net in unique([netfor!(n, pin) for pin in pins])
        append!(nets, net)
        delete!(n.nets, netname(n, net))
    end

    # add pins to named net
    if any(isa(pin, Symbol) for pin in pins)
        n.nets[filter(p -> isa(p, Symbol), collect(pins))[]] = nets
    else
        add!(n, nets[1])
        n.nets[netname(n, nets[1])] = nets
    end
end

"""
    disconnect!(n::Network, p::Tuple{Symbol,Symbol})
Disconnects the given pin `p` from anything else in the network `n`. The pin is
given as a `Tuple{Symbols,Any}`, where the first entry is the reference
designator of an element in `n`, and the second entry is the pin name. For
convenience, the latter is automatically converted to a `Symbol` as needed. Note
that if e.g. three pin `p1`, `p2`, and `p3` are connected then
`disconnect!(n, p1)` will disconnect `p1` from `p2` and `p3`, but leave `p2` and
`p3` connected to each other.
"""
function disconnect!(n::Network, pin::Tuple{Symbol,Symbol})
    net = netfor!(n, pin)
    filter!(p -> p != pin, net)

    push!(n.nets, [pin])
end
disconnect!(n::Network, p::Tuple{Symbol,Any}) = disconnect!(n, (p[1], Symbol(p[2])))

"""
    check_lumped_elements(net :: Network)
Checks if pins of all elements are connected.
"""
function check_lumped_elements(net :: Network)
    for (sym, elements) in net.nets
        if occursin("gnd", string(sym))
            continue
        else
            if length(elements) == 1
                (s, p) = elements[1]
                throw(ArgumentError("Element $s has lumped pin $p."))
            end
        end
    end
end

"""
function power_flow(net :: Network)
    Forms the dictionary needed for solving the power flow problem using
    package PowerModelsACDC. After successful power flow solving, it updates
    the operating point of the power converter.
"""
function power_flow(net :: Network)
    global global_dict
    global ang_min, ang_max
    global max_gen
    #TODO: Check if this has to be LN-RMS or LL-RMS. Do the necessary changes internally after validations against PSCAD.
    global_dict = PowerModelsMCDC.get_pu_bases(1000, 380) # 3-PH MVA, LL-RMS, Original setting was 100,320
    global_dict["omega"] = 2π * 50

    ang_min = deg2rad(360)
    ang_max = deg2rad(-360)

    function add_bus_ac(data :: Dict{String, Any})
        key = string(length(data["bus"]) + 1)
        (data["bus"])[key] = Dict{String, Any}()
        ((data["bus"])[key])["source_id"] = Any["bus", parse(Int, key)]
        ((data["bus"])[key])["index"] = parse(Int, key)
        ((data["bus"])[key])["bus_i"] = parse(Int, key)
        ((data["bus"])[key])["zone"] = 1
        ((data["bus"])[key])["area"] = 1
        ((data["bus"])[key])["vmin"] = 0.9
        ((data["bus"])[key])["vmax"] = 1.1
        ((data["bus"])[key])["vm"] = 1
        ((data["bus"])[key])["va"] = 0
        ((data["bus"])[key])["base_kv"] = global_dict["V"] / 1e3
        ((data["bus"])[key])["bus_type"] = 1 # bus type - depends on components
    end
    function add_bus_dc(data :: Dict{String, Any})
        key = string(length(data["busdc"]) + 1)
        (data["busdc"])[key] = Dict{String, Any}()
        ((data["busdc"])[key])["busdc_i"] = parse(Int, key)
        ((data["busdc"])[key])["source_id"] = Any["busdc", parse(Int, key)]
        ((data["busdc"])[key])["grid"] = 1
        ((data["busdc"])[key])["index"] = parse(Int, key)
        ((data["busdc"])[key])["Cdc"] = 0
        ((data["busdc"])[key])["Vdc"] = 1
        ((data["busdc"])[key])["Vdcmax"] = 1.1
        ((data["busdc"])[key])["Vdcmin"] = 0.9
        ((data["busdc"])[key])["Pdc"] = 0
        ((data["busdc"])[key])["basekVdc"] = global_dict["V"] / 1e3
    end

    function make_branch_ac(data :: Dict{String, Any}, element :: Element, dict_ac :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        !isempty(new_i) ? (push!(dict_ac, new_i); add_bus_ac(data); key_i = length(data["bus"])) :
                            key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)
        !isempty(new_o) ? (push!(dict_ac, new_o); add_bus_ac(data); key_o = length(data["bus"])) :
                            key_o = findfirst(p -> in(element.pins[Symbol(2.1)], p), dict_ac)
        key_e = length(data["branch"])+1
        (data["branch"])[string(key_e)] = Dict{String, Any}()
        ((data["branch"])[string(key_e)])["f_bus"] = key_i
        ((data["branch"])[string(key_e)])["t_bus"] = key_o
        ((data["branch"])[string(key_e)])["source_id"] = Any["branch", key_e]
        ((data["branch"])[string(key_e)])["index"] = key_e
        ((data["branch"])[string(key_e)])["rate_a"] = 1
        ((data["branch"])[string(key_e)])["rate_b"] = 1
        ((data["branch"])[string(key_e)])["rate_c"] = 1
        ((data["branch"])[string(key_e)])["br_status"] = 1
        ((data["branch"])[string(key_e)])["angmin"] = ang_min
        ((data["branch"])[string(key_e)])["angmax"] = ang_max

        make_power_flow_ac!(element.element_value, data, global_dict)
    end
    function make_branch_dc(data :: Dict{String, Any}, element :: Element, dict_dc :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        !isempty(new_i) ? (push!(dict_dc, new_i); add_bus_dc(data); key_i = length(data["busdc"])) :
                            key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_dc)
        !isempty(new_o) ? (push!(dict_dc, new_o); add_bus_dc(data); key_o = length(data["busdc"])) :
                            key_o = findfirst(p -> in(element.pins[Symbol(2.1)], p), dict_dc)
        key_e = length(data["branchdc"])+1
        (data["branchdc"])[string(key_e)] = Dict{String, Any}()
        ((data["branchdc"])[string(key_e)])["fbusdc"] = key_i
        ((data["branchdc"])[string(key_e)])["tbusdc"] = key_o
        ((data["branchdc"])[string(key_e)])["source_id"] = Any["branchdc", key_e]
        ((data["branchdc"])[string(key_e)])["index"] = key_e
        ((data["branchdc"])[string(key_e)])["rateA"] = 100
        ((data["branchdc"])[string(key_e)])["rateB"] = 100
        ((data["branchdc"])[string(key_e)])["rateC"] = 100
        ((data["branchdc"])[string(key_e)])["status"] = 1

        make_power_flow_dc!(element.element_value, data, global_dict)
    end

    function make_generator(data :: Dict{String, Any}, element :: Element, dict_dc :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        !isempty(new_i) ? (!any(occursin("gnd", string(x)) for x in new_i) && (push!(dict_ac, new_i); add_bus_ac(data))) : nothing
        key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)
        !isempty(new_o) ? (!any(occursin("gnd", string(x)) for x in new_o) && (push!(dict_ac, new_o); add_bus_ac(data))) : nothing
        key_o = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)

        (key_i === nothing) ? key_b = key_o : key_b = key_i

        key = length(data["gen"])+1
        (data["gen"])[string(key)] = Dict{String, Any}()
        ((data["gen"])[string(key)])["mBase"] = global_dict["S"] / 1e6
        ((data["gen"])[string(key)])["gen_bus"] = key_b     

        make_power_flow_ac!(element.element_value, data, global_dict)

        if isapprox(max_gen, element.element_value.P)
            ((data["bus"])[string(key_b)])["bus_type"] = 3
        else
            ((data["bus"])[string(key_b)])["bus_type"] = 2
        end
        ((data["bus"])[string(key_b)])["bus_type"] = 3
        ((data["bus"])[string(key_b)])["vm"] = ((data["gen"])[string(key)])["vg"]
        ((data["bus"])[string(key_b)])["vmin"] =  0.9*((data["gen"])[string(key)])["vg"]
        ((data["bus"])[string(key_b)])["vmax"] =  1.1*((data["gen"])[string(key)])["vg"]
    end

    function make_syngen(data :: Dict{String, Any}, element :: Element, dict_dc :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        !isempty(new_i) ? (!any(occursin("gnd", string(x)) for x in new_i) && (push!(dict_ac, new_i); add_bus_ac(data))) : nothing
        key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)
        !isempty(new_o) ? (!any(occursin("gnd", string(x)) for x in new_o) && (push!(dict_ac, new_o); add_bus_ac(data))) : nothing
        key_o = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)

        (key_i === nothing) ? key_b = key_o : key_b = key_i

        # !isempty(new_i) ? (push!(dict_ac, new_i); add_bus_ac(data); key_i = length(data["bus"])) :
        #                     key_i_br = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac)
        # !isempty(new_o) ? (push!(dict_ac, new_o); add_bus_ac(data); key_o = length(data["bus"])) :
        #                     key_o_br = findfirst(p -> in(element.pins[Symbol(2.1)], p), dict_ac)
        # Add a new branch to get power flow from the SG to the infinite bus

        add_bus_ac(data)
        key_i = length(data["bus"])
        push!(dict_ac, key_i)

        key_e = length(data["branch"])+1
        (data["branch"])[string(key_e)] = Dict{String, Any}()
        ((data["branch"])[string(key_e)])["f_bus"] = key_i
        ((data["branch"])[string(key_e)])["t_bus"] = key_o
        ((data["branch"])[string(key_e)])["source_id"] = Any["branch", key_e]
        ((data["branch"])[string(key_e)])["index"] = key_e
        ((data["branch"])[string(key_e)])["rate_a"] = 1
        ((data["branch"])[string(key_e)])["rate_b"] = 1
        ((data["branch"])[string(key_e)])["rate_c"] = 1
        ((data["branch"])[string(key_e)])["br_status"] = 1
        ((data["branch"])[string(key_e)])["angmin"] = ang_min
        ((data["branch"])[string(key_e)])["angmax"] = ang_max
        ((data["branch"])[string(key_e)])["transformer"] = false
        ((data["branch"])[string(key_e)])["tap"] = 1
        ((data["branch"])[string(key_e)])["shift"] = 0
        ((data["branch"])[string(key_e)])["c_rating_a"] = 1

        
        ((data["branch"])[string(key_e)])["br_r"] = element.element_value.rt * (element.element_value.Vᵃᶜ_base^2 / element.element_value.S_base) / global_dict["Z"]
        ((data["branch"])[string(key_e)])["br_x"] = element.element_value.lt * (element.element_value.Vᵃᶜ_base^2 / element.element_value.S_base) / global_dict["Z"]
        ((data["branch"])[string(key_e)])["g_fr"] = 0
        ((data["branch"])[string(key_e)])["b_fr"] = 0
        ((data["branch"])[string(key_e)])["g_to"] = 0
        ((data["branch"])[string(key_e)])["b_to"] = 0       

        key_b = key_i
        key = length(data["gen"])+1
        (data["gen"])[string(key)] = Dict{String, Any}()
        ((data["gen"])[string(key)])["mBase"] = global_dict["S"] / 1e6
        ((data["gen"])[string(key)])["gen_bus"] = key_b

        make_power_flow_ac!(element.element_value, data, global_dict)

        # TODO: Apply a reactive power based saturation instead.
        if isapprox(element.element_value.P_max, element.element_value.P)
            ((data["bus"])[string(key_b)])["bus_type"] = 1
        else
            ((data["bus"])[string(key_b)])["bus_type"] = 2
        end
        # ((data["bus"])[string(key_b)])["bus_type"] = 3
        ((data["bus"])[string(key_b)])["vm"] = ((data["gen"])[string(key)])["vg"]
        ((data["bus"])[string(key_b)])["vmin"] =  0.9*((data["gen"])[string(key)])["vg"]
        ((data["bus"])[string(key_b)])["vmax"] =  1.1*((data["gen"])[string(key)])["vg"]
    end

    function make_converter(data :: Dict{String, Any}, element :: Element,
        dict_dc :: Array{Any}, dict_ac :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        !isempty(new_i) ? (push!(dict_dc, new_i); add_bus_dc(data); key_i = length(data["busdc"])) :
                            key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_dc)
        !isempty(new_o) ? (push!(dict_ac, new_o); add_bus_ac(data); key_o = length(data["bus"])) :
                            key_o = findfirst(p -> in(element.pins[Symbol(2.1)], p), dict_ac)

        key = length(data["convdc"])+1
        (data["convdc"])[string(key)] = Dict{String, Any}()
        ((data["convdc"])[string(key)])["busdc_i"] = key_i
        ((data["convdc"])[string(key)])["busac_i"] = key_o

        make_power_flow!(element.element_value, data, global_dict)
    end

    function make_shunt_ac(data :: Dict{String, Any}, element :: Element, dict_ac :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        # merge input and output pins in one bus
        new_i = vcat(new_i, new_o)
        !isempty(new_i) ? (push!(dict_ac, new_i); add_bus_ac(data); key_i = length(data["bus"])) :
                            (key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac))
        key_e = length(data["shunt"])+1
        (data["shunt"])[string(key_e)] = Dict{String, Any}()
        ((data["shunt"])[string(key_e)])["source_id"] = Any["bus", key_i]
        ((data["shunt"])[string(key_e)])["index"] = key_e
        ((data["shunt"])[string(key_e)])["shunt_bus"]  = key_i
        data["shunt"][string(key_e)]["status"] = 1

        make_power_flow_ac!(element.element_value, data, global_dict)
    end

    function make_shunt_ac_impedance(data :: Dict{String, Any}, element :: Element, dict_ac :: Array{Any}, new_i :: Array{Symbol}, new_o :: Array{Symbol})
        # merge input and output pins in one bus
        # TODO: Consider representing as a load instead of a shunt.
        # TODO: Right now, this only works if the first pin of the load is connected to the bus, and the second pin grounded. This can be generalized.
        new_i = vcat(new_i, new_o)
        !isempty(new_i) ? (push!(dict_ac, new_i); add_bus_ac(data); key_i = length(data["bus"])) :
                            (key_i = findfirst(p -> in(element.pins[Symbol(1.1)], p), dict_ac))
        key_e = length(data["shunt"])+1
        (data["shunt"])[string(key_e)] = Dict{String, Any}()
        ((data["shunt"])[string(key_e)])["source_id"] = Any["bus", key_i]
        ((data["shunt"])[string(key_e)])["index"] = key_e
        ((data["shunt"])[string(key_e)])["shunt_bus"]  = key_i
        data["shunt"][string(key_e)]["status"] = 1

        abcd = eval_abcd(element.element_value, global_dict["omega"] * 1im)
        n = 3
        Z = (abcd[1:n,n+1:end])[1,1] / global_dict["Z"]
        data["shunt"][string(key_e)]["gs"] = real(1/Z)
        data["shunt"][string(key_e)]["bs"] = imag(1/Z)

    end

    !any((is_converter(element) || is_generator(element)) for element in values(net.elements)) && return

    dict_ac = Any[]
    dict_dc = Any[]

    data = Dict{String, Any}()
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

    processed = Symbol[]
    max_gen = 0
    for element in values(net.elements)
       if is_source(element)
           max_gen = max(max_gen, element.element_value.P)
       end
    end
    for (symbol, element) in net.elements
        new_i = Symbol[]
        new_o = Symbol[]

        for (key, val) in element.pins
            if !in(val, processed)
                push!(processed, val)
                (parse(Int, string(key)[1]) == 1) ? push!(new_i, val) : push!(new_o, val)
            end
        end
        # TODO: DC sources are not modeled in the power flow. Check with Hakan!
        if is_passive(element)
            if is_three_phase(element)
                if is_load(element) # TODO: Check this generalization in detail.
                    make_shunt_ac_impedance(data, element, dict_ac, new_i, new_o)
                else
                    make_branch_ac(data, element, dict_ac, new_i, new_o)
                end
            else
                make_branch_dc(data, element, dict_dc, new_i, new_o)
            end
        elseif is_source(element) && is_three_phase(element)
            make_generator(data, element, dict_ac, new_i, new_o)
        elseif is_generator(element)
            make_syngen(data, element, dict_ac, new_i, new_o)
        elseif is_converter(element)
            make_converter(data, element, dict_dc, dict_ac, new_i, new_o)
        end
    end

    PowerModelsMCDC.process_additional_data!(data)
    ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false)
    result = run_acdcpf(data, ACPPowerModel, ipopt; setting = s)
    println(result["termination_status"])
    # println(result["solution"]["bus"])
    # println("Power flow solution:")
    # println(result["solution"])
    # println("Calculating Jacobian")
    # data_for_jac = make_basic_network(data)
    # jac = calc_basic_jacobian_matrix(data)
    # writedlm("./files/power_flow_jacobian.csv",  jac, ',')
    id_converter = 1
    for (key, element) in net.elements
        if is_converter(element)
            conv_dict = result["solution"]["convdc"][string(id_converter)]
            I = conv_dict["iconv"] * global_dict["I"]
            Pdc = conv_dict["pdc"] * global_dict["S"] / 1e6
            Vm = (result["solution"]["bus"][string(data["convdc"][string(id_converter)]["busac_i"])]["vm"] *
                    global_dict["V"] / 1e3) * sqrt(2) # Convert the LN-RMS voltage coming from the PF to LN-PK
            θ = result["solution"]["bus"][string(data["convdc"][string(id_converter)]["busac_i"])]["va"]
            Vdc = result["solution"]["busdc"][string(data["convdc"][string(id_converter)]["busdc_i"])]["vm"] *
                    global_dict["V"] / 1e3
            Pac = -conv_dict["pgrid"] * global_dict["S"] / 1e6
            Qac = conv_dict["qgrid"] * global_dict["S"] / 1e6 # Think about this!
            if isa(element.element_value, MMC)
                update_string = "MMC #"
                update_mmc(element.element_value, Vm, θ, Pac, Qac, Vdc, Pdc)
            elseif isa(element.element_value, MMC_BI)
                update_string = "MMC #"
                update_mmc(element.element_value, Vm, θ, Pac, Qac, Vdc/2, -Vdc/2, Pdc)
            else
                update_string = "TLC #"
                update_tlc(element.element_value, Vm, θ, Pac, Qac, Vdc, Pdc)
            end
            print(update_string * string(id_converter) * " Active Power [MW]: ")
            println(Pac)
            print(update_string * string(id_converter) * " Reactive Power [MVar]: ")
            println(Qac)
            print(update_string * string(id_converter) * " AC Voltage Magnitude [kV]: ")
            println(Vm)
            print(update_string * string(id_converter) * " AC Voltage Angle [rad]: ")
            println(θ)
            print(update_string * string(id_converter) * " DC Voltage [kV]: ")
            println(Vdc)
            id_converter += 1
        end
    end
    id_gen = 1
    for (key, element) in net.elements
        # TODO: Generalize this method to find the correct branch.
        # Instead of taking the voltage and power from the generator bus, we need to check the bus of the branch.
        # Here there is a mismatch between id_gen and the id of the SG. We need to find another way to get to the network element that is a SG
        if is_generator(element) || is_source(element)
            if isa(element.element_value, SynchronousMachine)
                gen_dict = result["solution"]["gen"][string(id_gen)]
                # Search based on the active power gives problems for a meshed network 
                # (branch_dict[key]["pf"] == gen_dict["pg"]).
                # Trying now a seach based on the reactive power.
                branch_dict = result["solution"]["branch"]
                gen_branch_id = 0
                for (key_i,value) in branch_dict
                    if isapprox(branch_dict[key_i]["qf"],gen_dict["qg"],atol=1e-3)
                        gen_branch_id = key_i
                    end
                end
                # gen_branch_dict = data["branch"][string(key)]
                # We need to access the other end of the branch
                # gen_bus = data["branch"][string(gen_branch_id)]["t_bus"]
                Pgen = gen_dict["pg"] * global_dict["S"] / 1e6 #MW
                Qgen = gen_dict["qg"] * global_dict["S"] / 1e6 #MVAr
                Vm = (result["solution"]["bus"][string(data["branch"][string(gen_branch_id)]["t_bus"] )]["vm"] *
                        global_dict["V"] / 1e3) * sqrt(2) # Convert the LN-RMS voltage coming from the PF to LN-PK
                θ = result["solution"]["bus"][string(data["branch"][string(gen_branch_id)]["t_bus"] )]["va"]
                print("SG #" * string(id_gen) * " Active Power [MW]: ")
                println(Pgen)
                print("SG #" * string(id_gen) * " Reactive Power [MVar]: ")
                println(Qgen)
                print("SG #" * string(id_gen) * " AC Voltage Magnitude [kV]: ")
                println(Vm)
                print("SG #" * string(id_gen) * " AC Voltage Angle [rad]: ")
                println(θ)
                update_gen(element.element_value, Pgen, Qgen, Vm, θ)
            end
            id_gen += 1
        end
    end    
    return result
end

@doc doc"""
    @network begin #= ... =# end
Provides a simple domain-specific language to decribe networks. The
`begin`/`end` block can hold element definitions of the form
`refdes = elementfunc(params)` and connection specifications of the form
`refdes1[pin1] ⟷ refdes2[pin2]`.
# Example
To create a network with a voltage source connected to a resistor:

```jldoctest; output = false, setup = :(include("../src/HVDCstability.jl"); using .HVDCstability), filter = r"(HVDCstability\.)?Network\(.*"s
@network begin
    src = dc_source(V = 5)
    r = impedance(z = 1000, pins = 1)
    src[1.1] ⟷ r[1.1]
    src[2.1] ⟷ r[2.1]
end
# output
Network(...)
```
Alternatively, connection specifications can be given after an element
specification, separated by commas. In that case, the `refdes` may be omitted,
defaulting to the current element.
# Example
```jldoctest; output = false, setup = :(include("../src/HVDCstability.jl"); using .HVDCstability), filter = r"(HVDCstability\.)?Network\(.*"s
@network begin
    src = dc_source(V = 5)
    r = impedance(z = 1000, pins = 1), src[1.1] ⟷ [1.1], src[2.1] ⟷ [2.1]
end
# output
Network(...)
```
Finally, a connection endpoint may simply be of the form `netname`, to connect
to a named net. (Such named nets are created as needed.)
# Example
```jldoctest; output = false, setup = :(include("../src/HVDCstability.jl"); using .HVDCstability), filter = r"(HVDCstability\.)?Network\(.*"s
@network begin
    src = dc_source(V = 5), [2.1] ⟷ gnd
    r = impedance(z = 1000, pins = 1), [1.1] ⟷ src[1.1], [2.1] ⟷ gnd
end
# output
Network(...)
```
If a net or pin specification is not just a single symbol or number, and has to
be put in quotes (e.g. `"in+"`, `"9V"`)
!!! note
    Instead of `⟷` (`\\longleftrightarrow`), one can also use `==`.
"""
macro network(cdef)
    is_conn_spec(expr::Expr) =
        (expr.head === :call && (expr.args[1] === :(⟷) || expr.args[1] === :(↔) || expr.args[1] === :(==))) ||
        (expr.head === :comparison && all(c -> c === :(==), expr.args[2:2:end]))
    is_conn_spec(::Any) = false

    function elem_spec(expr)
        if !isa(expr, Expr) || expr.head !== :(=)
            error("invalid element specification$locinfo: $(expr)")
        end
        if !isa(expr.args[1], Symbol)
            error("invalid element identifier$locinfo: $(expr.args[1])")
        end
        if isa(expr.args[2], Expr) && expr.args[2].head === :tuple
            if isempty(expr.args[2].args)
                error("invalid element specification$locinfo: $(expr.args[2])")
            end
            elemspec = expr.args[2].args[1]
            conn_exprs = expr.args[2].args[2:end]
        else
            elemspec = expr.args[2]
            conn_exprs = []
        end
        push!(ccode.args, :(add!(network, $(QuoteNode(expr.args[1])), $(esc(elemspec)))))
        for conn_expr in conn_exprs
            if !is_conn_spec(conn_expr)
                error("invalid connection specification$locinfo: $conn_expr")
            end
            push!(ccode.args, Expr(:call, :connect!, :network, extractpins(conn_expr, expr.args[1])...))
        end
    end

    function extractpins(expr::Expr, default_element=nothing)
        if expr.head === :call && (expr.args[1] === :(⟷) || expr.args[1] === :(↔) || expr.args[1] === :(==))
            return vcat((extractpins(a, default_element) for a in expr.args[2:end])...)
        elseif expr.head === :comparison && all(c -> c === :(==), expr.args[2:2:end])
            return vcat((extractpins(a, default_element) for a in expr.args[1:2:end])...)
        elseif expr.head === :ref
            return [:(($(QuoteNode(expr.args[1])), $(QuoteNode(expr.args[2]))))]
        elseif expr.head === :vect && length(expr.args) == 1
            if default_element === nothing
                error("missing element$(locinfo): $expr")
            end
            return [:(($(QuoteNode(default_element)), $(QuoteNode(expr.args[1]))))]
        else
            error("invalid pin specification$(locinfo): $expr")
        end
    end

    function extractpins(netname::Symbol, default_element=nothing)
        return [QuoteNode(netname)]
    end

    extractpins(netname::String, default_element=nothing) =
        extractpins(Symbol(netname), default_element)

    if !isa(cdef, Expr) || cdef.head !== :block
        error("@network must be followed by a begin/end block")
    end
    ccode = Expr(:block)
    push!(ccode.args, :(network = Network()))
    locinfo = ""
    for expr in cdef.args
        if isa(expr, LineNumberNode)
            locinfo = " at $(expr.file):$(expr.line)"
            continue
        end
        if !isa(expr, Expr)
            error("invalid statement in network definition$locinfo: $expr")
        end
        if expr.head === :line
            locinfo = " at $(expr.args[2]):$(expr.args[1])"
        elseif expr.head === :(=)
            elem_spec(expr)
        elseif is_conn_spec(expr)
            push!(ccode.args, Expr(:call, :connect!, :network, extractpins(expr)...))
        else 
            error("invalid statement in network definition$locinfo: $expr")
        end
    end

    # here you can add functions for network before the end of the code
    push!(ccode.args, :(check_lumped_elements(network)))
    push!(ccode.args, :(connect!(network)))
    push!(ccode.args, :(power_flow(network)))
    push!(ccode.args, :(network))
    return ccode
end


@doc doc"""
    function composite_element(subnet::Network, input_pins::Array{Any}, output_pins::Array{Any})
Create a net element from the (sub-)network `net`. The `input_pins` and `output_pin`
define input and output nodes of the element.
# Example
```jldoctest; output = false, setup = :(include("../src/HVDCstability.jl"); using .HVDCstability), filter = r"(HVDCstability\.)?Element\(.*"s
net = @network begin
   r1 = impedance(z = 10e3, pins = 1)
   r2 = impedance(z = 10e3, pins = 1), [1.1] == r1[2.1]
   c = impedance(z = 10e3, pins = 1), [1.1] == r2[1.1], [2.1] == r2[2.1]
   src = dc_source(V = 5), [1.1] == r1[1.1], [2.1] == r2[2.1]
end
composite_element(net, Any[(:r2, Symbol(1.1))], Any[(:r2, Symbol(2.1))])
# output

Element(...)
```
"""
function composite_element(subnet::Network, input_pins::Array{Any}, output_pins::Array{Any})
    element = Element(input_pins = length(input_pins), output_pins = length(output_pins),
            element_value = subnet)
    for i in 1:length(input_pins)
        name = netname(subnet, input_pins[i])
        subnet.connections[Symbol("1.",i)] = subnet.nets[name]
        element.pins[Symbol("1.",i)] = name
    end
    for i in 1:length(output_pins)
        name = netname(subnet, output_pins[i])
        subnet.connections[Symbol("2.",i)] = subnet.nets[name]
        element.pins[Symbol("2.",i)] = name
    end

    return element
end


function eval_abcd(subnet :: Network, s :: Complex)
    start_pins = Symbol[]
    end_pins = Symbol[]
    dict = Dict{Symbol, Array{Union{Symbol,Int}}}(:node_list => Symbol[], :element_list => Symbol[],
        :output_list => Symbol[])
    dict[:element_list] = [elem_symbol for elem_symbol in keys(subnet.elements)]
    for (key, val) in subnet.connections
        if occursin("2.", string(key))
            push!(dict[:output_list], netname(subnet, val))
            push!(end_pins, netname(subnet, val))
        else
            push!(start_pins, netname(subnet, val))
        end

    end
    dict[:node_list] = setdiff([node_symbol for node_symbol in keys(subnet.nets)], vcat(start_pins, end_pins))
    dict[:node_list] = vcat(start_pins, dict[:node_list])

    make_abcd(subnet, dict, start_pins, end_pins, s)
end
