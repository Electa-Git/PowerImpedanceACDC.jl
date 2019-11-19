export dc_source

# include("source.jl")

"""
    dc_source(voltage, exp)
Creates dc voltage in Volts.

Internal impedance can be added with a command impedance after
the equality sign.

Pins: `1.1` and `2.1`

Plus pin is connected to `1.1` and minus to `2.1`. To ground the source,
connect the pin to the ground while constructing the network.
"""
function dc_source(; voltage = 0)
    source = Source(V = voltage, pins = 1)
    source.ABCD = convert(Array{Basic}, Diagonal([1 for dummy in 1:2]))
    elem = Element(input_pins = 1, output_pins = 1,
        element_value = source)
    elem
end
