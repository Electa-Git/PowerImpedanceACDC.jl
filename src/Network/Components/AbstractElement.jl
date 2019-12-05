__precompile__()

export save_data, plot_data

include("tools/tools.jl")

"""
Struct guarantees representation of the component like a multiport
network using ABCD parameters. It consists of:
- element unique symbol inside constructed network - `symbol`
- dictionary that maps pins inside the network - `pins`
- number of input pins - `input_pins`
- number of output pins - `output_pins`
- component definition - `element_value`
"""
mutable struct Element
  symbol::Symbol
  pins :: Dict{Symbol, Symbol}
  input_pins :: Int
  output_pins :: Int
  element_value :: Any  # component defined type

  transformation :: Bool

  function Element(;args...)
    elem = new()
    for (key, val) in kwargs_pairs(args)
      if (key in propertynames(elem))
        setfield!(elem, key, val)
      else
        throw(ArgumentError("The property name $(key) is not defined."))
      end
    end

    # definition of pins
    if !isdefined(elem, :transformation)
        elem.transformation = false
    elseif (elem.transformation)
        elem.input_pins -= 1
        elem.output_pins -= 1
    end
    if !isdefined(elem, :pins)
      elem.pins = merge(Dict{Symbol, Symbol}(Symbol(string("1.",i)) => Symbol() for i in 1:nip(elem)),
                        Dict{Symbol, Symbol}(Symbol(string("2.",i)) => Symbol() for i in 1:nop(elem)))
    end

    elem
  end
end

for (n,m) in Dict(:nip => :input_pins, :nop => :output_pins)
  @eval ($n)(e::Element) = e.$m
end
np(e::Element) = nip(e) + nop(e) # number pins

function add!(elem::Element, sym::Symbol, value::Any)
  if (sym in propertynames(elem))
    setfield!(elem, sym, value)
  end
end

function get_nodes(element::Element)
    return values(element.pins)
end

function get_nodes(element::Element, pin::Symbol)
    array = Symbol[]
    for (key, val) in element.pins
        (pin != key) && push!(array, val)
        # !occursin(string(pin)[1:2], string(key)) && push!(array, val)
    end
    return array
end

################### ABCD functions ################################
function get_abcd(element::Element, s::Complex)
    abcd = eval_abcd(element.element_value, s)

    if (element.transformation)
        if np(element) == 2
            return transformation_dc(abcd)
        elseif is_three_phase(element)
            ω₀ = 100*π
            abcd₁ = eval_abcd(element.element_value, s + 1im*ω₀)
            abcd₂ = eval_abcd(element.element_value, s - 1im*ω₀)
            return transformation_dq(abcd₁, abcd₂)
        end
    end
    return abcd
end

nip_abcd(e::Element) = size(get_abcd(e, 1im),1)
nop_abcd(e::Element) = size(get_abcd(e, 1im),2)
np_abcd(e::Element) = Int((nip_abcd(e) + nop_abcd(e))/2) # number pins

########################## Y functions #############################
function get_y(element :: Element, s :: Complex)
    y = eval_y(element.element_value, s)
    # if (element.transformation)
    #     return transformation_dc(abcd)
    # end
    return y
end

######################### Element type #############################
function is_passive(element :: Element)
    (isa(element.element_value, MMC) || isa(element.element_value, Source)) && return false
    true
end

function is_source(element :: Element)
    isa(element.element_value, Source)
end

function is_converter(element :: Element)
    isa(element.element_value, MMC)
end

function is_three_phase(element :: Element)
    (np(element) == 6) || (np(element) == 4 && (element.transformation)) && return true
    return false
end

####################### Save/plot ##################################
"""
    function save_data(element :: Element, file_name :: String; omega_range = (-3, 5, 1000),
        scale = :log)
Used for saving data of the specific component defined as `Element`. The function calls component
specific function. The data is saved in csv textual file.

Additional parameter `omega_range` provides possibility to manually add the
frequency scale for saving data. Scale can be given as logarithmic (`:log`) and
linear (`:lin`).
"""
function save_data(element :: Element, file_name :: String; omega_range = (-3, 5, 1000),
    scale = :log)
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    if (scale == :log)
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = [min_ω+i*n for i in 1:Int(n_ω)]
    end

    file_name = string("./files/",  file_name)
    save_data(element.element_value, file_name, omegas)
end

"""
    function plot_data(element :: Element; omega_range = (-3, 5, 1000),
        scale = :log)

It plots the component defined data with the desired frequency range.
It differs from component to component.

Additional parameter `omega_range` provides possibility to manually add the
frequency scale for saving data. Scale can be given as logarithmic (`:log`) and
linear (`:lin`).
"""
function plot_data(element :: Element; omega_range = (-3, 5, 1000),
    scale = :log)
    (min_ω, max_ω, n_ω) = omega_range
    n = (max_ω - min_ω) / n_ω
    if (scale == :log)
        omegas = [exp10(min_ω)*10^(i*n) for i in 1:Int(n_ω)]
    else
        omegas = [min_ω+i*n for i in 1:Int(n_ω)]
    end

    plot_data(element.element_value, omegas)
end

include("../Network.jl")
include("element_types.jl")
