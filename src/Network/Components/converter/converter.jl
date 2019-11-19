abstract type Converter end

function make_power_flow!(converter :: Converter, dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = string(length(dict["convdc"]))
    ((dict["convdc"])[string(key)])["source_id"] = Any["convdc", parse(Int, key)]
    ((dict["convdc"])[string(key)])["status"] = 1
    ((dict["convdc"])[string(key)])["index"] = parse(Int, key)

    ((dict["convdc"])[string(key)])["basekVac"] = global_dict["V"] / 1e3

    ((dict["convdc"])[string(key)])["type_ac"] = 1 # default, PQ ac bus
    if in(:power, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 1  # constant AC active power
        ((dict["convdc"])[string(key)])["type_ac"] = 2  # PV ac bus
    elseif in(:dc, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 2  # constant DC voltage
    else
        ((dict["convdc"])[string(key)])["type_dc"] = 3  # DC voltage droop
    end

    # droop control - not implemented
    ((dict["convdc"])[string(key)])["droop"] = 0
    ((dict["convdc"])[string(key)])["Pdcset"] = 1
    ((dict["convdc"])[string(key)])["Vdcset"] = converter.Vᵈᶜ * 1e3 / global_dict["V"]
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
    ((dict["convdc"])[string(key)])["xc"] = converter.Lᵣ * global_dict["omega"] / global_dict["Z"]
    converter.ω₀ = global_dict["omega"]

    # default values
    ((dict["convdc"])[string(key)])["Vmmax"] = 1.1 * converter.Vₘ * 1e3 / global_dict["V"]
    ((dict["convdc"])[string(key)])["Vmmin"] = 0.9 * converter.Vₘ * 1e3 / global_dict["V"]
    ((dict["convdc"])[string(key)])["Imax"] = 1.1 * converter.P / converter.Vₘ

    ((dict["convdc"])[string(key)])["P_g"] = converter.P
    ((dict["convdc"])[string(key)])["Q_g"] = converter.Q
    ((dict["convdc"])[string(key)])["Vtar"] = converter.Vₘ * 1e3 / global_dict["V"]

    ((dict["convdc"])[string(key)])["LossA"] = 0
    ((dict["convdc"])[string(key)])["LossB"] = 0
    ((dict["convdc"])[string(key)])["LossCrec"] = converter.Rₐᵣₘ / 2
    ((dict["convdc"])[string(key)])["LossCinv"] = converter.Rₐᵣₘ / 2

    ((dict["convdc"])[string(key)])["Qacmax"] = converter.Q_max
    ((dict["convdc"])[string(key)])["Qacmin"] = converter.Q_min
    ((dict["convdc"])[string(key)])["Pacmax"] = converter.P_max
    ((dict["convdc"])[string(key)])["Pacmin"] = converter.P_min

    key_o = ((dict["convdc"])[string(key)])["busac_i"]
    if (dict["bus"][string(key_o)]["bus_type"] == 1)
        dict["bus"][string(key_o)]["vm"] = ((dict["convdc"])[string(key)])["Vtar"]
        ((dict["bus"])[string(key_o)])["vmin"] = 0.9 * dict["bus"][string(key_o)]["vm"]
        ((dict["bus"])[string(key_o)])["vmax"] = 1.1 * dict["bus"][string(key_o)]["vm"]
    end
    key_i = ((dict["convdc"])[string(key)])["busdc_i"]
    ((dict["busdc"])[string(key_i)])["Vdc"] = converter.Vᵈᶜ * 1e3 / global_dict["V"]
    ((dict["busdc"])[string(key_i)])["Vdcmax"] = 1.1 * ((dict["busdc"])[string(key_i)])["Vdc"]
    ((dict["busdc"])[string(key_i)])["Vdcmin"] = 0.9 * ((dict["busdc"])[string(key_i)])["Vdc"]
end
