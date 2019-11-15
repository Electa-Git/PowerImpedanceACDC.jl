export ac_source_3Φ, ac_source

include("source.jl")

"""
    ac_source(;args...)

Creates  n phase ac source with specified amplitude V,
phase shift (in radians) θ, reference active power P, reference reactive power Q,
maximum and minimum generated active and reactive powers.

Pins: `1.1`, `1.2`, `1.3`, and `2.1`, `2.2`, `2.3`

Plus pin is connected to `1.x` and minus to `2.x`. To ground the source,
connect the pin to the ground while constructing the network.
"""
function ac_source(;args...)
    source = Source()
    for (key, val) in kwargs_pairs(args)
        setfield!(source, key, val)
    end
    source.ABCD = convert(Array{Basic}, Diagonal([1 for dummy in 1:6]))
    elem = Element(input_pins = source.pins, output_pins = source.pins, element_value = source)
    elem
end
