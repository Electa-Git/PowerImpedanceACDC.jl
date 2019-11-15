abstract type Converter end

function make_power_flow!(converter :: Converter, dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = string(length(dict["convdc"]))
    ((dict["convdc"])[string(key)])["source_id"] = Any["convdc", parse(Int, key)]
    ((dict["convdc"])[string(key)])["status"] = 1
    ((dict["convdc"])[string(key)])["index"] = parse(Int, key)

    ((dict["convdc"])[string(key)])["basekVac"] = global_dict["V"] / 1e3

    ((dict["convdc"])[string(key)])["type_ac"] = 1 # default, PQ ac bus
    if in(:power, converter.controls)
        ((dict["convdc"])[string(key)])["type_dc"] = 1  # constant AC active power
    elseif in(:dc, converter.controls)
        ((dict["convdc"])[string(key)])["type_dc"] = 2  # constant DC voltage
    else
        ((dict["convdc"])[string(key)])["type_dc"] = 3  # DC voltage droop
    end

    # droop control - not implemented
    ((dict["convdc"])[string(key)])["droop"] = 0
    ((dict["convdc"])[string(key)])["Pdcset"] = 0
    ((dict["convdc"])[string(key)])["Vdcset"] = converter.Vᵈᶜ / global_dict["V"]
    ((dict["convdc"])[string(key)])["dVdcSet"] = 0
    # LCC converter
    ((dict["convdc"])[string(key)])["islcc"] = 0

    # without transformer
    ((dict["convdc"])[string(key)])["transformer"] = 0
    ((dict["convdc"])[string(key)])["rtf"] = 0
    ((dict["convdc"])[string(key)])["xtf"] = 0
    ((dict["convdc"])[string(key)])["tm"] = 1
    # without filter
    ((dict["convdc"])[string(key)])["filter"] = 0
    ((dict["convdc"])[string(key)])["bf"] = 0
    # with reactor
    ((dict["convdc"])[string(key)])["reactor"] = 1
    ((dict["convdc"])[string(key)])["rc"] = converter.Rᵣ / global_dict["Z"]
    ((dict["convdc"])[string(key)])["xc"] = converter.Lᵣ * global_dict["omega"] / global_dict["V"]
    converter.ω₀ = global_dict["omega"]

    # default values
    ((dict["convdc"])[string(key)])["Imax"] = 1.1
    ((dict["convdc"])[string(key)])["Vmmax"] = 1.1
    ((dict["convdc"])[string(key)])["Vmmin"] = 0.9

    ((dict["convdc"])[string(key)])["P_g"] = converter.P
    ((dict["convdc"])[string(key)])["Q_g"] = converter.Q
    ((dict["convdc"])[string(key)])["Vtar"] = converter.Vm

    ((dict["convdc"])[string(key)])["LossA"] = 0
    ((dict["convdc"])[string(key)])["LossB"] = 0
    ((dict["convdc"])[string(key)])["LossCrec"] = converter.Rₐᵣₘ / 2
    ((dict["convdc"])[string(key)])["LossCinv"] = converter.Rₐᵣₘ / 2

    ((dict["convdc"])[string(key)])["Qacmax"] = 50
    ((dict["convdc"])[string(key)])["Qacmin"] = -50
    ((dict["convdc"])[string(key)])["Pacmax"] = 100
    ((dict["convdc"])[string(key)])["Pacmin"] = -100

end
