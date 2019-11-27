export ac_source

"""
    ac_source(;args...)

Creates  n phase ac source with specified amplitude V,
phase shift (in radians) θ, reference active power P, reference reactive power Q,
maximum and minimum generated active and reactive powers.

Pins: `1.x` and `2.x` for x ∈ {1, ..., pins}

Plus pin is connected to `1.x` and minus to `2.x`. To ground the source,
connect the pin to the ground while constructing the network.
"""
function ac_source(;args...)
    source = Source()
    transformation = false
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(source))
            setfield!(source, key, val)
        elseif (key == :transformation)
            transformation = val
        else
            throw(ArgumentError("Source does not have a property $(key)."))
        end
    end
    make_abcd(source)
    elem = Element(input_pins = source.pins, output_pins = source.pins, element_value = source,
        transformation = transformation)
    elem
end
