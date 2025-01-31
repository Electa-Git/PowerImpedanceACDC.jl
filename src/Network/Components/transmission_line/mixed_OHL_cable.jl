export mixed_OHL_cable

@with_kw mutable struct Mixed_OHL_cable <: Transmission_line
    parts :: Dict{Symbol, Transmission_line} = Dict{Symbol, Transmission_line}()
end


function mixed_OHL_cable(;args...)
    m = Mixed_OHL_cable()
    pins = 0
    transformation = false
    for (key, val) in pairs(args)
        if (key == :transformation)
            transformation = val
        else
            m.parts[key] = val.element_value
            pins = nip(val)
        end
    end

    elem = Element(input_pins = pins, output_pins = pins, element_value = m,
            transformation = transformation)
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
