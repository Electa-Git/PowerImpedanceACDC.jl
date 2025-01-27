export tlc

@with_kw mutable struct TLC <: Converter
    ω₀ :: Union{Int, Float64} = 100*π

    P :: Union{Int, Float64} = -10              # active power [MW]
    Q :: Union{Int, Float64} = 3                # reactive power [MVA]
    P_dc :: Union{Int, Float64} = 100           # DC power [MW]
    P_min :: Union{Float64, Int} = -100         # min active power output [MW]
    P_max :: Union{Float64, Int} = 100          # max active power output [MW]
    Q_min :: Union{Float64, Int} = -50          # min reactive power output [MVA]
    Q_max :: Union{Float64, Int} = 50           # max reactive power output [MVA]

    θ :: Union{Int, Float64} = 0
    Vₘ :: Union{Int, Float64} = 333             # AC voltage, amplitude [kV]
    Vᵈᶜ :: Union{Int, Float64} = 640            # DC-bus voltage [kV]

    Lₐᵣₘ :: Union{Int, Float64}  = 0        # arm inductance [H]
    Rₐᵣₘ :: Union{Int, Float64}  = 0        # equivalent arm resistance

    Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the phase reactor [H]
    Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the phase reactor [Ω]

    controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
    equilibrium :: Array{Union{Int, Float64}} = [0]
    A :: Array{Complex} = [0]
    B :: Array{Complex} = [0]
    C :: Array{Complex} = [0]
    D :: Array{Complex} = [0]

    timeDelay :: Float64 = 0
    padeOrderNum :: Int = 0
    padeOrderDen :: Int = 0

    vACbase_LL_RMS :: Union{Int, Float64} = 220 # Voltage base in kV
    Sbase :: Union{Int, Float64} = 500 # Power base in MW

    vACbase :: Float64 = 0 # AC voltage base for impedance/admittance calculation
    iACbase :: Float64 = 0 # AC current base for impedance/admittance calculation
end

"""
    function mmc(;args...)
It constructs MMC operating both as a rectifier and an inverter. MMC is constructed as a struct with the
following fields.
```julia
ω₀ :: Union{Int, Float64} = 100*π

    P :: Union{Int, Float64} = -10              # active power [MW]
    Q :: Union{Int, Float64} = 3                # reactive power [MVA]
    P_dc :: Union{Int, Float64} = 100           # DC power [MW]
    P_min :: Union{Float64, Int} = -100         # min active power output [MW]
    P_max :: Union{Float64, Int} = 100          # max active power output [MW]
    Q_min :: Union{Float64, Int} = -50          # min reactive power output [MVA]
    Q_max :: Union{Float64, Int} = 50           # max reactive power output [MVA]

    θ :: Union{Int, Float64} = 0
    Vₘ :: Union{Int, Float64} = 333             # AC voltage, amplitude [kV]
    Vᵈᶜ :: Union{Int, Float64} = 640            # DC-bus voltage [kV]

    Lₐᵣₘ :: Union{Int, Float64}  = 0        # arm inductance [H]
    Rₐᵣₘ :: Union{Int, Float64}  = 0        # equivalent arm resistance

    Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the phase reactor [H]
    Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the phase reactor [Ω]

    controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
    equilibrium :: Array{Union{Int, Float64}} = [0]
    A :: Array{Complex} = [0]
    B :: Array{Complex} = [0]
    C :: Array{Complex} = [0]
    D :: Array{Complex} = [0]

    timeDelay :: Float64 = 0
    padeOrderNum :: Int = 0
    padeOrderDen :: Int = 0

    vACbase_LL_RMS :: Union{Int, Float64} = 220 # Voltage base in kV
    Sbase :: Union{Int, Float64} = 500 # Power base in MW

    vACbase :: Float64 = 0 # AC voltage base for impedance/admittance calculation
    iACbase :: Float64 = 0 # AC current base for impedance/admittance calculation
"""


function tlc(;args...)
    converter = TLC()

    for (key, val) in pairs(args)
        if isa(val, Controller)
            converter.controls[key] = val
        elseif in(key, propertynames(converter))
            setfield!(converter, key, val)
        end
    end

    elem = Element(input_pins = 1, output_pins = 2, element_value = converter)
end

function update_tlc(converter :: TLC, Vm, θ, Pac, Qac, Vdc, Pdc)
    

    wbase = 100*pi
    vAC_base = converter.vACbase_LL_RMS*sqrt(2/3)
    Sbase = converter.Sbase
    iAC_base = 2*Sbase/3/vAC_base
    zAC_base = (3/2)*vAC_base^2/Sbase
    lAC_base = zAC_base/wbase

    converter.vACbase = vAC_base
    converter.iACbase = iAC_base

    Lᵣ = converter.Lᵣ / lAC_base
    Rᵣ = converter.Rᵣ / zAC_base

    ω₀ = converter.ω₀

    Qac *=-1 # Correction for reactive power sign

    converter.Vₘ = Vm
    converter.θ = θ
    converter.Vᵈᶜ = Vdc
    converter.P = Pac
    converter.Q = Qac
    converter.P_dc = Pdc # Has the same sign as Pac

    Vm /= vAC_base
    Vdc /= (2*vAC_base)
    Pac /= Sbase
    Qac /= Sbase
    Pdc /= Sbase
    
    Vᴳd = Vm # Vᴳd = Vm * cos(θ)
    Vᴳq = 0# Vᴳq = -Vm * sin(θ)

    Id = ((Vᴳd * Pac + Vᴳq * Qac) / (Vᴳd^2 + Vᴳq^2)) 
    Iq = ((Vᴳq * Pac - Vᴳd * Qac) / (Vᴳd^2 + Vᴳq^2)) 

    # setup control parameters and equations
    init_x = zeros(2, 1)

    # TODO: These have to be updated to be compliant with the PU model!
    for (key, val) in (converter.controls)
        # fix coefficients
        if (val.Kₚ == 0) && (val.Kᵢ == 0)
            if (key == :occ)                            # pole placement
                val.Kᵢ = Lₑ * val.bandwidth^2
                val.Kₚ = 2 * val.ζ * val.bandwidth * Lₑ - Rₑ
            end
        end

        # fix reference values
        if (key == :occ)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Id Iq]
            end
            init_x[1] = val.ref[1]
            init_x[2] = val.ref[2]
        elseif (key == :p)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Pac]
            end
        elseif (key == :q)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Qac]
            end
        elseif (key == :dc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Vdc]
            end
        elseif (key == :vac) || (key == :vac_supp)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Vm]
            end
        end
    end

    index = 2
    index_PLL_angle = 0

    init_x[1] = Pac
    init_x[2] = Qac

    exp = Expr(:block)
    # Represent measurements here before possible voltage filtering
    if in(:pll, keys(converter.controls))
        index_PLL_angle = index + 3
        if in(:v_meas_filt, keys(converter.controls))
            index_PLL_angle +=2
        end        
        push!(exp.args, :(
            T_θ = [cos(x[$index_PLL_angle]) -sin(x[$index_PLL_angle]); sin(x[$index_PLL_angle]) cos(x[$index_PLL_angle])];
            I_θ = [cos(x[$index_PLL_angle]) sin(x[$index_PLL_angle]); -sin(x[$index_PLL_angle]) cos(x[$index_PLL_angle])];          
        ))
    else
        push!(exp.args, :(
            T_θ = [1 0; 0 1];
            I_θ = [1 0; 0 1];
        ))
    end
    push!(exp.args, :(
        Vdc = inputs[1];
        (Vᴳd, Vᴳq) = T_θ * [inputs[2]; inputs[3]];
    ))

    # add voltage measurement filter
    if in(:v_meas_filt, keys(converter.controls))
        push!(exp.args, :(
            Vᴳd_f = x[$index+1];
            Vᴳq_f = x[$index+2];
            F[$index+1] = $(converter.controls[:v_meas_filt].ω_f) * (Vᴳd - Vᴳd_f);
            F[$index+2] = $(converter.controls[:v_meas_filt].ω_f) * (Vᴳq - Vᴳq_f);
            ))
            indexVᴳdf = index + 1
            indexVᴳqf = index + 2
        index +=2
    else
        push!(exp.args, :(
            (Vᴳd_f, Vᴳq_f) = (Vᴳd, Vᴳq);
            ))
    end
    

    # add PLL
    if in(:pll, keys(converter.controls))
        if (converter.controls[:pll].ω_f != 0) # A PLL filter is implemented
            push!(exp.args, :(
                vₚₗₗ = x[$index+1];
                F[$index+1] = $(converter.controls[:pll].ω_f) * (Vᴳq_f - vₚₗₗ);
            ))
            index +=1
        else
            push!(exp.args, :(
                vₚₗₗ = Vᴳq
            ))
        end
        push!(exp.args, :(
            F[$index+1] = -vₚₗₗ*$(converter.controls[:pll].Kᵢ);
            Δω = $(converter.controls[:pll].Kₚ) * (-vₚₗₗ) + x[$index+1];
            ω = $(converter.ω₀)/$wbase + Δω;
            F[$index+2] = $wbase*Δω;
        ))
        index += 2
    else
        push!(exp.args, :(
            Δω = 0;
        ))
    end
    push!(exp.args, :((i_d_pcc_c, i_q_pcc_c) = T_θ * [x[1]; x[2]];))

    # add current measurement filter
    if in(:i_meas_filt, keys(converter.controls))
        push!(exp.args, :(
            (i_d_pcc_f, i_q_pcc_f) = x[$index+1:$index+2];
            F[$index+1] = $(converter.controls[:i_meas_filt].ω_f) * (i_d_pcc_c - i_d_pcc_f);
            F[$index+2] = $(converter.controls[:i_meas_filt].ω_f) * (i_q_pcc_c - i_q_pcc_f);
            ))
        index +=2
    else
        push!(exp.args, :(
            (i_d_pcc_f, i_q_pcc_f) = (i_d_pcc_c, i_q_pcc_c);
            ))
    end

    # DC voltage control not yet implemented!!
    # TODO: Think about DC voltage control and generalize the case for the absence of power controllers
    if in(:p, keys(converter.controls))
        # add frequency support
        if in(:f_supp, keys(converter.controls))
            push!(exp.args, :(
                F[$index+1] = $(converter.controls[:f_supp].ω_f) *(-$(converter.controls[:f_supp].Kₚ)*Δω - x[$index+1]);
                p_ref = $(converter.controls[:p].ref[1]) + x[$index+1]))
            index +=1
        else
            push!(exp.args, :(
                p_ref = $(converter.controls[:p].ref[1])))
        end
        # active power control
        push!(exp.args, :(
            P_ac = (Vᴳd_f * i_d_pcc_f + Vᴳq_f * i_q_pcc_f);
            # iΔd_ref = (Kp_Pac * (Pac_ref - Pac) + Ki_Pac * xiPac);
            id_ref = ($(converter.controls[:p].Kₚ) * (p_ref - P_ac) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:p].Kᵢ) *(p_ref - P_ac)))
        index += 1
    end
    if in(:q, keys(converter.controls))
        # add voltage support
        if in(:vac_supp, keys(converter.controls))
            push!(exp.args, :(
                Vᴳ_mag = sqrt(Vᴳd_f^2 + Vᴳq_f^2);
                Δq_unf = $(converter.controls[:vac_supp].Kₚ)*($(converter.controls[:vac_supp].ref[1])-Vᴳ_mag);
                F[$index+1] = $(converter.controls[:vac_supp].ω_f) *(Δq_unf - x[$index+1]);
                q_ref = $(converter.controls[:q].ref[1]) + x[$index+1]))
            index +=1
        else
            push!(exp.args, :(
                q_ref = $(converter.controls[:q].ref[1])))
        end
        # reactive power control
        push!(exp.args, :(
            Q_ac =  (-Vᴳq_f * i_d_pcc_f + Vᴳd_f * i_q_pcc_f);
            iq_ref = ($(converter.controls[:q].Kₚ) * (q_ref - Q_ac) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:q].Kᵢ) *(q_ref - Q_ac)))
        index += 1
    end
    # add control equations
    for (key, val) in (converter.controls)                
        if (key == :occ)
            # output current control
            push!(exp.args, :(
                        F[$index+1] = $(converter.controls[:occ].Kᵢ) * (id_ref - i_d_pcc_f);
                        F[$index+2] = $(converter.controls[:occ].Kᵢ) * (iq_ref - i_q_pcc_f);

                        md_c = 2 * ( x[$index+1] +
                                    $(converter.controls[:occ].Kₚ) * (id_ref - i_d_pcc_f) + $Lᵣ * (1 + Δω) * i_q_pcc_f + Vᴳd_f) / Vdc;
                        mq_c = 2 * ( x[$index+2] +
                                    $(converter.controls[:occ].Kₚ) * (iq_ref - i_q_pcc_f) - $Lᵣ * (1 + Δω) * i_d_pcc_f + Vᴳq_f) / Vdc;
                        (md, mq) = I_θ * [md_c; mq_c]))
            index += 2
        end
    end

    if !in(:occ, keys(converter.controls))
        push!(exp.args,
                    :(md = 0;
                      mq = 0;))
    end

    # TODO: Not yet finalized
    # add time delays here, if there are controllers implemented
    if (converter.timeDelay != 0.0) && (in(:occ, keys(converter.controls)))
        push!(exp.args,
            :(
            T_ab_dq=0.5*[1 im;-im 1];# from alpha-beta to dq
            T_dq_ab=0.5*[1 -im;im 1];#from dq to alpha-beta
            ))
        if in(:occ, keys(converter.controls))
            push!(exp.args,
            :(
                timeDelayIn = [vMΔ_ref;vMΔ_ref];
                statesDelay = x[$index + 1 : $index + 2*$converter.padeOrderDen]; 
                timeDelayOut = timeDelayPadeMatrices($converter.padeOrderNum,$converter.padeOrderDen,$converter.timeDelay,length(timeDelayIn));
                
                A_delay = timeDelayOut[1];
                B_delay = timeDelayOut[2];
                C_delay = timeDelayOut[3];
                D_delay = timeDelayOut[4];
                F[$index + 1 : $index + 2*$converter.padeOrderDen] = A_delay*statesDelay + B_delay*timeDelayIn;
                # timeDelayOut = C_delay*statesDelay + D_delay*timeDelayIn;
                # Implement phase shifts by transforming the dq voltage references to alpha-beta
                vM_ab_ref = (cos($converter.ω₀*$converter.timeDelay)-sin($converter.ω₀*$converter.timeDelay)*im)*(T_dq_ab*(C_delay*statesDelay + D_delay*timeDelayIn));
                vM_dq_ref = real(T_ab_dq * conj(vM_ab_ref) + conj(T_ab_dq) * vMΔ_ab_ref);
                vMd_ref = vM_dq_ref[1];
                vMq_ref = vM_dq_ref[2];
            ))
            index += 2*converter.padeOrderDen
        end
        
    end

    # add state variables
    push!(exp.args,
    :(
        (vMd, vMq) = 0.5 * Vdc * [md; mq];
        
        # dw neglected here
        F[1] = (vMd - inputs[2] - $Rᵣ*x[1] - $Lᵣ*x[2])/$Lᵣ;             
        F[2] = (vMq - inputs[3] - $Rᵣ*x[2] + $Lᵣ*x[1])/$Lᵣ;       
        F[1:2] *= $wbase;

        ))

    exp_equilibrium = copy(exp)
    # add outputs (DC current and dq currents)
    push!(exp.args,
    :(
        F[$index+1] = (vMd * x[1] + vMq * x[2]) / Vdc ; # Power balance
        F[$index+2] = x[1] ; 
        F[$index+3] = x[2] ; 
        ))

    function f!(expr, F, x, inputs) # F derivative of state variable x state variable vector, inputs input vqlue expr equation of mmc
       f = eval(:((F,x,inputs) -> $expr))
       return Base.invokelatest(f, F,x,inputs)
    end

    vector_inputs = [Vdc, Vᴳd, Vᴳq]
    init_x = [init_x; zeros(index-2,1)]
    if in(:v_meas_filt, keys(converter.controls))
        init_x[indexVᴳdf] =  Vᴳd
        init_x[indexVᴳqf] =  1e-3 # Initialize to a small non-zero value to avoid Inf or NaN problems with nlsolve
    end

    g!(F,x) = f!(exp_equilibrium, F, x, vector_inputs)
    k = nlsolve(g!, init_x, autodiff = :forward, iterations = 100, ftol = 1e-6, xtol = 1e-3, method = :trust_region)
    if converged(k)
        println("TLC steady-state solution found!")
    end
    converter.equilibrium = k.zero

    number_output = 3
    number_input =3

    h(F,x) = f!(exp, F, x[1:end-number_input], x[end-number_input+1:end])
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index + number_output); h(F, x); return F) # 
    jac = zeros(index + number_output , index + number_input)
    ForwardDiff.jacobian!(jac, ha, [k.zero' vector_inputs'])
    converter.A = jac[1:index, 1:index] # index indicates the number of state variables
    converter.B = jac[1:index, index+1:end]
    converter.C = jac[index+1:end, 1:index]
    converter.D = jac[index+1:end, index+1:end]

end

function eval_parameters(converter :: TLC, s :: Complex)
    # numerical
    I = Matrix{Complex}(Diagonal([1 for dummy in 1:size(converter.A,1)]))
    # Y = (converter.C*inv(s*I-converter.A))*converter.B + converter.D # This matrix is in pu
    Y = converter.C * ((s*I-converter.A) \ converter.B) + converter.D # This matrix is in pu
    Y *= converter.iACbase / converter.vACbase
    return Y
end

