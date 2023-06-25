abstract type Machine end

function make_power_flow_ac!(machine :: Machine, dict :: Dict{String, Any},
    global_dict :: Dict{String, Any})
    key = string(length(dict["gen"]))
    ((dict["gen"])[string(key)])["pc1"] = 0
    ((dict["gen"])[string(key)])["pc2"] = 0
    ((dict["gen"])[string(key)])["qc1min"] = 0
    ((dict["gen"])[string(key)])["qc1max"] = 0
    ((dict["gen"])[string(key)])["qc2min"] = 0
    ((dict["gen"])[string(key)])["qc2max"] = 0
    ((dict["gen"])[string(key)])["ramp_agc"] = 0
    ((dict["gen"])[string(key)])["ramp_q"] = 0
    ((dict["gen"])[string(key)])["ramp_10"] = 0
    ((dict["gen"])[string(key)])["ramp_30"] = 0
    ((dict["gen"])[string(key)])["apf"] = 0
    ((dict["gen"])[string(key)])["startup"] = 0
    ((dict["gen"])[string(key)])["shutdown"] = 0

    ((dict["gen"])[string(key)])["gen_status"] = 1
    ((dict["gen"])[string(key)])["source_id"] = Any["gen", parse(Int, key)]
    ((dict["gen"])[string(key)])["index"] = parse(Int, key)

    S_base = global_dict["S"] / 1e6
    V_base = global_dict["V"] / 1e3
    ((dict["gen"])[string(key)])["pg"] = machine.P / S_base
    ((dict["gen"])[string(key)])["qg"] = machine.Q / S_base
    ((dict["gen"])[string(key)])["pmin"] = machine.P_min / S_base
    ((dict["gen"])[string(key)])["pmax"] = machine.P_max / S_base
    ((dict["gen"])[string(key)])["qmin"] = machine.Q_min / S_base
    ((dict["gen"])[string(key)])["qmax"] = machine.Q_max / S_base
    ((dict["gen"])[string(key)])["vg"] = machine.V / V_base

    # not using
    ((dict["gen"])[string(key)])["model"] = 1
    ((dict["gen"])[string(key)])["cost"] = 0
    ((dict["gen"])[string(key)])["ncost"] = 0

end

function make_power_flow_dc!(machine :: Machine, dict :: Dict{String, Any},
    global_dict :: Dict{String, Any})
nothing
end

function eval_abcd(machine :: Machine, s :: Complex)
    return eval_y(machine, s)
end

function eval_y(machine :: Machine, s :: Complex)
    Y = eval_parameters(machine, s)
    return Y
end


