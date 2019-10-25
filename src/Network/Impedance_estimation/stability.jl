export check_stability

"""
    function check_stability(net :: Network, mmc :: Element, direction :: Symbol = :dc)
"""
function check_stability(net :: Network, mmc :: Element, direction :: Symbol = :dc)
    if !isa(mmc.element_value, MMC)
        throw(ArgumentError("Cannot determine stability of the passive element."))
    end

    node_list = []

    if (direction == :dc)
        mmc.ABCD.direction = 2
        Z_mmc = closing_impedance(get_abcd(mmc), 0)
        Y_mmc = 1/Z_mmc[1]
        for (pin, node_name) in pairs(mmc.pins)
            if (occursin("2.", string(pin)))
                push!(node_list, node_name)
            end
        end
        imp = determine_impedance(net, input_pins = ((mmc.symbol, node_list[1]),),
                                    output_pins = ((mmc.symbol, node_list[2]),))[]
        return [Y_mmc, imp, imp * Y_mmc]
    else
        for (pin, node_name) in pairs(mmc.pins)
            if (occursin("1.", string(pin)))
                push!(node_list, node_name)
            end
        end
        imp = (determine_impedance(net, input_pins = ((mmc.symbol, node_list[1]),),
                                    output_pins = ((mmc.symbol, node_list[2]),)))
    end
end
