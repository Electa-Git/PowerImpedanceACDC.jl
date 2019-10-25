export dc_source

@with_kw mutable struct Source
    value :: Union{Float64, Int} = 0      # impedance value
    pins :: Int = 1
end

"""
    dc_source(voltage, exp)
Creates dc voltage in Volts.

Internal impedance can be added with a command impedance after
the equality sign.

Pins: `1` and `2`

Plus pin is connected to `1` and minus to `2`. To ground the source,
connect the pin to the ground while constructing the network.
"""
function dc_source(; voltage = 0)
    elem = Element(input_pins = 1, output_pins = 1,
        element_value = Source(value = voltage, pins = 1))

    elem
end

function create_abcd!(element :: Element, src :: Source)
    abcd = convert(Array{Basic}, Diagonal([1 for dummy in 1:np(element)]))
end

function eval_abcd(element :: Element, source :: Source, s :: Complex)
    abcd = N.(element.ABCD)
end
