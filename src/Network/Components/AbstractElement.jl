"""
Struct guarantees representation of the component like a multiport
network using ABCD parameters. It consists of:
- element unique symbol inside constructed network - `symbol`
- dictionary that maps pins inside the network (to network nodes) - `pins`
- number of input pins - `input_pins`
- number of output pins - `output_pins`
- component definition - `element_value`
- transformation flag - `transformation`
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
    for (key, val) in pairs(args)
      if (key in propertynames(elem))
        setfield!(elem, key, val)
      else
        throw(ArgumentError("The property name $(key) is not defined."))
      end
    end

    # definition of pins
    if !isdefined(elem, :transformation)
        elem.transformation = false
    elseif (elem.transformation) #TODO: Not generalizable, only makes sense for 3-phase systems
        elem.input_pins -= 1
        elem.output_pins -= 1
    end
    if !isdefined(elem, :pins) # Initialize pins field with empty symbol, if not defined 
      elem.pins = merge(Dict{Symbol, Symbol}(Symbol(string("1.",i)) => Symbol() for i in 1:nip(elem)),
                        Dict{Symbol, Symbol}(Symbol(string("2.",i)) => Symbol() for i in 1:nop(elem)))
    end

    elem
  end
end

for (n,m) in Dict(:nip => :input_pins, :nop => :output_pins)
  @eval ($n)(e::Element) = e.$m # creation of functions nip() and nop(), fetching the input_pins and output_pins parameters within the element structure
end
np(e::Element) = nip(e) + nop(e) # total number of pins

function add!(elem::Element, sym::Symbol, value::Any)
  if (sym in propertynames(elem))
    setfield!(elem, sym, value)
  end
end

function get_nodes(element::Element) # Returns all nodes connected to the element
    return values(element.pins)
end

function get_nodes(element::Element, pin::Symbol) # Returns all nodes connected to the element, except the one specified by pin
    array = Symbol[]
    for (key, val) in element.pins
        (pin != key) && push!(array, val)
        # !occursin(string(pin)[1:2], string(key)) && push!(array, val)
    end
    return array
end

################### ABCD functions ################################
function get_abcd(element::Element, s::Complex)
    
    if (element.transformation)
        if np(element) == 2
            abcd = eval_abcd(element.element_value, s)
            return transformation_dc(abcd)
        elseif is_three_phase(element)
            ω₀ = 100*π
            abcd₁ = eval_abcd(element.element_value, s + 1im*ω₀)
            abcd₂ = eval_abcd(element.element_value, s - 1im*ω₀)
            return transformation_dq(abcd₁, abcd₂)
        end
    else
        abcd = eval_abcd(element.element_value, s)
    end
    return abcd
end

function nip_abcd(e::Element)
    if isa(e.element_value, MMC) || isa(e.element_value, TLC)
        return 3
    else
        return 2nip(e)
    end
    # return 2nip(e)
end

function nop_abcd(e::Element)
    if isa(e.element_value, MMC) || isa(e.element_value, TLC)
        return 3
    else
        return 2nop(e)
    end
    # return 2nop(e)
end
np_abcd(e::Element) = Int((nip_abcd(e) + nop_abcd(e))/2) # number pins

########################## Y functions #############################
# TODO: Crashes when using with MMC, TLC, SynchronousMachine. Still needed?
function get_y(element :: Element, s :: Complex)
    abcd = get_abcd(element, s)
    return abcd_to_y(abcd)
end

######################### Element type #############################
function is_passive(element :: Element)
    (isa(element.element_value, MMC) || isa(element.element_value, TLC) || isa(element.element_value, Source) || isa(element.element_value, SynchronousMachine)) && return false
    true
end

function is_source(element :: Element)
    isa(element.element_value, Source)
end

function is_converter(element :: Element)
    (isa(element.element_value, MMC) || isa(element.element_value, TLC))
end

function is_shunt_reactor(element :: Element)
    isa(element.element_value, Shunt_reactor)
end

function is_generator(element :: Element)
    isa(element.element_value, SynchronousMachine)
end
 

function is_impedance(element :: Element)
    isa(element.element_value, Impedance) && !any(occursin("gnd", string(x)) for x in element.pins)
end

function is_load(element :: Element)
    isa(element.element_value, Impedance) && any(occursin("gnd", string(x)) for x in element.pins)
end

function is_three_phase(element :: Element)
    (np(element) == 6) || (np(element) == 4 && (element.transformation) && !is_converter(element)) && return true
    return false
end



