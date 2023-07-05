abstract type Converter end

function eval_abcd(converter :: Converter, s :: Complex)
    return eval_y(converter, s)
end

function eval_y(converter :: Converter, s :: Complex)
    Y = eval_parameters(converter, s)
    return Y
end

function make_power_flow!(converter :: Converter, dict :: Dict{String, Any},
            global_dict :: Dict{String, Any})
    key = string(length(dict["convdc"]))
    ((dict["convdc"])[string(key)])["source_id"] = Any["convdc", parse(Int, key)]
    ((dict["convdc"])[string(key)])["status"] = 1
    ((dict["convdc"])[string(key)])["index"] = parse(Int, key)

    ((dict["convdc"])[string(key)])["basekVac"] = global_dict["V"] / 1e3

    if in(:vac, keys(converter.controls)) || in(:vac_supp, keys(converter.controls)) 
        ((dict["convdc"])[string(key)])["type_ac"] = 2  # PV ac bus
    else
        ((dict["convdc"])[string(key)])["type_ac"] = 1  # PQ ac bus TODO: Check if this makes sense.
    end
    if in(:p, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 1  # constant AC active power        
    elseif in(:dc, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 2  # constant DC voltage
    else
        ((dict["convdc"])[string(key)])["type_dc"] = 3  # DC voltage droop
    end


    # droop control - not implemented
    ((dict["convdc"])[string(key)])["droop"] = 0
    ((dict["convdc"])[string(key)])["Pdcset"] = converter.P_dc
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
    ((dict["convdc"])[string(key)])["rc"] = (converter.Rᵣ + converter.Rₐᵣₘ / 2) / global_dict["Z"]
    ((dict["convdc"])[string(key)])["xc"] = (converter.Lᵣ + converter.Lₐᵣₘ / 2) * global_dict["omega"] / global_dict["Z"]
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

function save_data(conv :: Converter, file_name :: String, omegas)
    open(string(file_name, "_y.txt"), "w") do f
        for omega in omegas
            Y = eval_parameters(conv, 1im*omega)
            writedlm(f, [omega reshape(Y, 1, length(Y))], ",")
        end
    end
end

function plot_data(conv :: Converter, omegas)
    Y_ac = []
    Y_dc = []
    Y_acdc = []
    for omega in omegas
        Y = (eval_parameters(conv, 1im*omega))
        push!(Y_ac, Y[2:3,2:3])
        push!(Y_dc, Y[1,1])
        push!(Y_acdc, [Y[2:3,1] Y[1,2:3]])
    end

    p = bode(Y_ac, omega = omegas, titles = ["Y_{d,d}" "Y_{d,q}"; "Y_{q,d}" "Y_{q,q}"])
    save_plot(p, "files/ac_side")
    p = bode(Y_dc, omega = omegas, titles = ["Y_{dc}"])
    save_plot(p, "files/dc_side")
    p = bode(Y_acdc, omega = omegas, titles = ["Y_{d,dc}" "Y_{dc,d}"; "Y_{q, dc}" "Y_{dc,q}"])
    save_plot(p, "files/acdc_side")
end
