export mixed_OHL_cable

include("cable.jl")
include("transmission_line.jl")


@with_kw mutable struct Mixed_OHL_cable
    parts :: Dict{Symbol, Union{Cable, Transmission_line}} = Dict{Symbol, Union{Cable, Transmission_line}}()
end

"""
    function mixed_OHL_cable(;args...)
This function creates mixed OHL-cable component. The component consists of dictionary
with overhead line and cable definitions. It creates one component by concatenating
overhead line/cable sections. Overhead line and cable definitions are as already defined.
"""
function mixed_OHL_cable(;args...)
    m = Mixed_OHL_cable()
    pins = 0
    for (key, val) in kwargs_pairs(args)
        m.parts[key] = val.element_value
        pins = nip(val)
    end

    elem = Element(input_pins = pins,
                    output_pins = pins, element_value = m)
end

function create_abcd!(element :: Element, m :: Mixed_OHL_cable)
    n = Int(np(element))
    abcd_tf =  reshape([Basic(string(element.symbol,i)) for i in 1:n^2], n, n)
    return abcd_tf
end

function eval_abcd(m :: Mixed_OHL_cable, s :: Complex)
    abcd = 0
    for (key, val) in m.parts
        if (abcd == 0)
            abcd = eval_abcd(m.parts[key], s)
        else
            abcd = abcd * eval_abcd(m.parts[key], s)
        end
    end
    return abcd
end
