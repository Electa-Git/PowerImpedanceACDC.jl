# Main network definition

export Network, add!, connect!, disconnect!, @network,
        composite_element, eval_abcd
export power_flow, data, result  # for testing

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
    voltageBase :: Array{Float64,1} # Network line-ground RMS voltage base used in the power flow. This is necessary to get a matching power flow result in the presence of voltage controlling converters.
    Network() = new(OrderedDict{Symbol, Element}(), Dict{Symbol, Net}(), Dict{Symbol, Vector{Int}}(), [220/sqrt(3)])
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

# Doctest, should be added again to @doc block
# ```jldoctest; output = false, setup = :(include("../src/PowerImpedanceACDC.jl"); using .PowerImpedanceACDC), filter = r"(PowerImpedanceACDC\.)?Network\(.*"s
# network = Network()
# add!(network, :r, impedance(z = 1e3, pins = 1))
# add!(network, :src, dc_source(V = 5))
# connect!(network, (:src, 2.1), (:r, 2.1), :gnd) # connect to gnd net
# network
# # output
# Network(...)
# ```

@doc doc"""
    connect!(n::Network, pins::Union{Symbol,Tuple{Symbol,Any}}...)
Connects the given pins (or named nets) to each other in the network `n`. Named
nets are given as `Symbol`s, pins are given as `Tuple{Symbols,Any}`s, where the
first entry is the reference designator of an element in `c`, and the second
entry is the pin name. For convenience, the latter is automatically converted to
a `Symbol` as needed.
# Example

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

# Doctest, should be added again to @doc block
# ```jldoctest; output = false, setup = :(include("../src/PowerImpedanceACDC.jl"); using .PowerImpedanceACDC), filter = r"(PowerImpedanceACDC\.)?Network\(.*"s
# @network begin
#     src = dc_source(V = 5)
#     r = impedance(z = 1000, pins = 1)
#     src[1.1] ⟷ r[1.1]
#     src[2.1] ⟷ r[2.1]
# end
# # output
# Network(...)
# ```

# ```jldoctest; output = false, setup = :(include("../src/PowerImpedanceACDC.jl"); ImpedanceACDC), filter = r"(PowerImpedanceACDC\.)?Network\(.*"s
# @network begin
#     src = dc_source(V = 5)
#     r = impedance(z = 1000, pins = 1), src[1.1] ⟷ [1.1], src[2.1] ⟷ [2.1]
# end
# # output
# Network(...)
# ```

# ```jldoctest; output = false, setup = :(include("../src/PowerImpedanceACDC.jl"); ImpedanceACDC), filter = r"(PowerImpedanceACDC\.)?Network\(.*"s
# @network begin
#     src = dc_source(V = 5), [2.1] ⟷ gnd
#     r = impedance(z = 1000, pins = 1), [1.1] ⟷ src[1.1], [2.1] ⟷ gnd
# end
# # output
# Network(...)
# ```

@doc doc"""
    @network begin #= ... =# end
Provides a simple domain-specific language to decribe networks. The
`begin`/`end` block can hold element definitions of the form
`refdes = elementfunc(params)` and connection specifications of the form
`refdes1[pin1] ⟷ refdes2[pin2]`.
# Example
To create a network with a voltage source connected to a resistor:


Alternatively, connection specifications can be given after an element
specification, separated by commas. In that case, the `refdes` may be omitted,
defaulting to the current element.
# Example

Finally, a connection endpoint may simply be of the form `netname`, to connect
to a named net. (Such named nets are created as needed.)
# Example

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
            if expr.args[1] == :voltageBase
                push!(ccode.args, :(network.voltageBase[1] = $(esc(expr.args[2]))))
                return
            end
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

# Doctest, should be added again to @doc
# ```jldoctest; output = false, setup = :(include("../src/PowerImpedanceACDC.jl"); using .PowerImpedanceACDC), filter = r"(PowerImpedanceACDC\.)?Element\(.*"s
# net = @network begin
#    r1 = impedance(z = 10e3, pins = 1)
#    r2 = impedance(z = 10e3, pins = 1), [1.1] == r1[2.1]
#    c = impedance(z = 10e3, pins = 1), [1.1] == r2[1.1], [2.1] == r2[2.1]
#    src = dc_source(V = 5), [1.1] == r1[1.1], [2.1] == r2[2.1]
# end
# composite_element(net, Any[(:r2, Symbol(1.1))], Any[(:r2, Symbol(2.1))])
# # output

# Element(...)
# ```

@doc doc"""
    function composite_element(subnet::Network, input_pins::Array{Any}, output_pins::Array{Any})
Create a net element from the (sub-)network `net`. The `input_pins` and `output_pin`
define input and output nodes of the element.
# Example

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
