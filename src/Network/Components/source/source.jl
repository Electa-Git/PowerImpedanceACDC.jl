
@with_kw mutable struct Source
    impedance :: Union{Float64, Int, Basic} = 0
    V :: Union{Float64, Int} = 0        # DC voltage or voltage magnitude [kV]

    P   :: Union{Float64, Int} = 0      # active power output [MW]
    Q   :: Union{Float64, Int} = 0      # reactive power output [MVAr]
    P_min :: Union{Float64, Int} = 0    # min active power output
    P_max :: Union{Float64, Int} = 0    # max active power output
    Q_min :: Union{Float64, Int} = 0    # min reactive power output
    Q_max :: Union{Float64, Int} = 0    # max reactive power output

    pins :: Int = 1
    ABCD :: Array{Basic} = Basic[]
end

function eval_abcd(source :: Source, s :: Complex)
    abcd = N.(source.ABCD)
end


function make_power_flow_dc!(source :: Source, dict :: Dict{String, Any},
        global_dict :: Dict{String, Any})
    nothing
end

function make_power_flow_ac!(source :: Source, dict :: Dict{String, Any},
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
    ((dict["gen"])[string(key)])["pg"] = source.P / S_base
    ((dict["gen"])[string(key)])["qg"] = source.Q / S_base
    ((dict["gen"])[string(key)])["pmin"] = source.P_min / S_base
    ((dict["gen"])[string(key)])["pmax"] = source.P_max / S_base
    ((dict["gen"])[string(key)])["qmin"] = source.Q_min / S_base
    ((dict["gen"])[string(key)])["qmax"] = source.Q_max / S_base
    ((dict["gen"])[string(key)])["vg"] = source.V / V_base

    # not using
    ((dict["gen"])[string(key)])["model"] = 1
    ((dict["gen"])[string(key)])["cost"] = 0
    ((dict["gen"])[string(key)])["ncost"] = 0
end
