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
        # TODO: The line below sometimes gives errors during power flow (NUMERICAL_ERROR)
        dict["bus"][string(((dict["convdc"])[string(key)])["busac_i"])]["bus_type"] = 2 # Not entirely sure if this is necessary.
        if in(:vac, keys(converter.controls))
            ((dict["convdc"])[string(key)])["Vtar"] = converter.controls[:vac].ref[1] * 1e3 / (global_dict["V"] * sqrt(2))
        else
            ((dict["convdc"])[string(key)])["Vtar"] = converter.controls[:vac_supp].ref[1] * 1e3 / (global_dict["V"] * sqrt(2))
        end
        dict["bus"][string(((dict["convdc"])[string(key)])["busac_i"])]["vm"] = ((dict["convdc"])[string(key)])["Vtar"]
    else
        ((dict["convdc"])[string(key)])["type_ac"] = 1  # PQ ac bus TODO: Check if this makes sense in the presence of GFM
        ((dict["convdc"])[string(key)])["Vtar"] = converter.Vₘ * 1e3 / global_dict["V"]
    end
    if in(:p, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 1  # constant AC active power        
    elseif in(:dc, keys(converter.controls))
        ((dict["convdc"])[string(key)])["type_dc"] = 2  # constant DC voltage
    else
        ((dict["convdc"])[string(key)])["type_dc"] = 3  # DC voltage droop
    end

    if in(:vac_supp, keys(converter.controls))
        ((dict["convdc"])[string(key)])["acq_droop"] = 1
        ((dict["convdc"])[string(key)])["kq_droop"] = converter.controls[:vac_supp].Kₚ
    else
        ((dict["convdc"])[string(key)])["acq_droop"] = 0
        ((dict["convdc"])[string(key)])["kq_droop"] = 0
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
    ((dict["convdc"])[string(key)])["Imax"] = 1.1 * abs(converter.P) / converter.Vₘ

    ((dict["convdc"])[string(key)])["P_g"] = converter.P
    ((dict["convdc"])[string(key)])["Q_g"] = converter.Q
    

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

function timeDelayPadeMatrices(padeOrderNum,padeOrderDen,t_delay,numberVars)

    # Calculation of the state-space representation of a pade approximation with nth-order numerator and nth-order denominator
    # padeOrderNum = Order of the numerator, padeOrderDen= Order of the denominator, t_delay = Time delay duration, numberVars= Variables to be delayed max.2 !
    #
    size_A=padeOrderDen;
    a_k=factorial(padeOrderNum);
    b_l=((factorial(padeOrderDen)*(-1)^padeOrderNum)*(t_delay^(padeOrderNum-padeOrderDen)))/a_k;
    Ad=zeros(size_A,size_A);
    Bd=zeros(size_A,1);
    Bd[end]=1;
    Cd=zeros(1,size_A);
    Dd=b_l;
    Ad[1:end-1,2:end] = Matrix(1.0I, padeOrderDen-1, padeOrderDen-1);
    for i=0:padeOrderDen-1
        a_i=(t_delay^(i-padeOrderDen)*(factorial(padeOrderNum+padeOrderDen-i)*factorial(padeOrderDen)/(factorial(i)*factorial(padeOrderDen-i))))/a_k;
        b_i=(t_delay^(i-padeOrderDen)*((-1)^i)*(factorial(padeOrderNum+padeOrderDen-i)*factorial(padeOrderNum)/(factorial(i)*factorial(padeOrderNum-i))))/a_k;
        Ad[end,i+1] =-a_i;
        Cd[i+1] = (b_i-(a_i*b_l));
    end
    # Result above gives state space in canonical realization, which is ill-conditioned, especially for higher order pades.
    # In order to improve numerical handling, transform state space into modal realization
    # Alternative # TODO: What about Dd?
    sys = ss(Ad,Bd,Cd,Dd)
    sys_modal = modal_form(sys; C1=true)
    Ad = sys_modal[1].A
    Bd = sys_modal[1].B
    Cd = sys_modal[1].C
    if numberVars == 1
        A_Pade=Ad;
        B_Pade=Bd;
        C_Pade=Cd;
        D_Pade=Dd;
    elseif numberVars == 2
        A_Pade=cat(Ad,Ad;dims=[1,2]);
        B_Pade=cat(Bd,Bd;dims=[1,2]);
        C_Pade=cat(Cd,Cd;dims=[1,2]);
        D_Pade=cat(Dd,Dd;dims=[1,2]);
    end

    return A_Pade,B_Pade,C_Pade,D_Pade
end