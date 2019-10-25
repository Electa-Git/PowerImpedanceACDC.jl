export ac_source_3Φ, ac_source

"""
    ac_source(amplitude, frequency, phase, impedance)

Creates 1 and 3 phase ac source with specified amplitude, frequency and
phase shift (in radians), respectively.

Internal impedance can be added with a command impedance.

Pins: `1.1` and `2.1`

Plus pin is connected to "1.1" and minus to "2.1". To ground the source,
connect the pin to the ground while constructing the network.
"""
function ac_source(; amplitude = 0, frequency = 0, phase = 0, impedance = 0)
    elem = Element(input_pins = 1, output_pins = 1, element_value = Source(pins = 1))
    elem
end

"""
    ac_source_3Φ(;amplitude, frequency, phase, impedance)

Creates  3 phase ac source with specified amplitude, frequency and
phase shift (in radians), respectively.

Pins: `1.1`, `1.2`, `1.3`, and `2.1`, `2.2`, `2.3`

Plus pin is connected to `1.x` and minus to `2.x`. To ground the source,
connect the pin to the ground while constructing the network.
"""
function ac_source_3Φ(; amplitude = 0, frequency = 0, phase = 0, impedance = 0)
    elem = Element(input_pins = 3, output_pins = 3, element_value = Source(pins = 3))
    elem
end
