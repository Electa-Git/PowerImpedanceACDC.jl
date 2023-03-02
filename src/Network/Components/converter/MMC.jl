export mmc

include("converter.jl")
include("controller.jl")

@with_kw mutable struct MMC <: Converter
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

    Lₐᵣₘ :: Union{Int, Float64}  = 50e-3        # arm inductance [H]
    Rₐᵣₘ :: Union{Int, Float64}  = 1.07         # equivalent arm resistance
    Cₐᵣₘ :: Union{Int, Float64}  = 10e-3        # capacitance per submodule [F]
    N :: Int = 400                              # number of submodules per arm

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
end

"""
    function mmc(;args...)
It constructs MMC operating both as a rectifier and an inverter. MMC is constructed as a struct with the
following fields.
```julia
ω₀ :: Union{Int, Float64} = 100*π

P :: Union{Int, Float64} = -10              # active power [MW]
Q :: Union{Int, Float64} = 3                # reactive power [MVA]
P_dc :: Union{Int, Float64} = 100           # DC power [kW]
P_min :: Union{Float64, Int} = -100         # min active power output [MW]
P_max :: Union{Float64, Int} = 100          # max active power output [MW]
Q_min :: Union{Float64, Int} = -50          # min reactive power output [MVA]
Q_max :: Union{Float64, Int} = 50           # max reactive power output [MVA]

θ :: Union{Int, Float64} = 0
Vₘ :: Union{Int, Float64} = 333             # AC voltage [kV]
Vᵈᶜ :: Union{Int, Float64} = 640            # DC-bus voltage [kV]

Lₐᵣₘ :: Union{Int, Float64}  = 50e-3        # arm inductance [H]
Rₐᵣₘ :: Union{Int, Float64}  = 1.07         # equivalent arm resistance
Cₐᵣₘ :: Union{Int, Float64}  = 10e-3        # capacitance per submodule [F]
N :: Int = 401                              # number of submodules per arm

Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the phase reactor [H]
Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the phase reactor [Ω]

controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
equilibrium :: Array{Union{Int, Float64}} = [0]
A :: Array{Complex} = [0]
B :: Array{Complex} = [0]
C :: Array{Complex} = [0]
D :: Array{Complex} = [0]
```

The constructed MMC has 2 pins on the AC side: `2.1`, `2.2`, and 1 pin on its
DC-side: `1.1`.
"""
function mmc(;args...)
    converter = MMC()

    for (key, val) in kwargs_pairs(args)
        if isa(val, Controller)
            converter.controls[key] = val
        elseif in(key, propertynames(converter))
            setfield!(converter, key, val)
        end
    end

    elem = Element(input_pins = 1, output_pins = 2, element_value = converter)
end

function update_mmc(converter :: MMC, Vm, θ, Pac, Qac, Vdc, Pdc)
    Lₑ = converter.Lₐᵣₘ / 2 + converter.Lᵣ
    Rₑ = converter.Rₐᵣₘ / 2 + converter.Rᵣ
    N = converter.N
    Lₐᵣₘ = converter.Lₐᵣₘ
    Rₐᵣₘ = converter.Rₐᵣₘ
    Cₐᵣₘ = converter.Cₐᵣₘ
    Cₑ = 6*Cₐᵣₘ / N
    ω₀ = converter.ω₀

    converter.Vₘ = Vm
    converter.θ = θ
    converter.Vᵈᶜ = Vdc
    converter.P = Pac
    converter.Q = Qac
    converter.P_dc = Pdc

    Vm *= 1e3
    Vdc *= 1e3
    Pac *= 1e6
    Qac *= 1e6
    Pdc *= 1e6

    Vᴳd = Vm * cos(θ)
    Vᴳq = -Vm * sin(θ)
    Id = 2/3*(Vᴳd * Pac + Vᴳq * Qac) / (Vᴳd^2 + Vᴳq^2)
    Iq = 2/3*(Vᴳq * Pac - Vᴳd * Qac) / (Vᴳd^2 + Vᴳq^2)

    # println(Pac)
    # println(Qac) 
    # setup control parameters and equations
    init_x = zeros(12, 1)

    for (key, val) in (converter.controls)
        # fix coefficients
        if (val.Kₚ == 0) && (val.Kᵢ == 0)
            if (key == :occ)                            # pole placement
                val.Kᵢ = Lₑ * val.bandwidth^2
                val.Kₚ = 2 * val.ζ * val.bandwidth * Lₑ - Rₑ

            elseif (key == :ccc) || (key == :zcc)       # pole placement
                val.Kᵢ = converter.Lₐᵣₘ * val.bandwidth^2
                val.Kₚ = 2 * val.ζ * val.bandwidth * converter.Lₐᵣₘ - converter.Rₐᵣₘ
            end
        end

        # fix reference values
        if (key == :occ)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Id Iq]
            end
            init_x[1] = val.ref[1]
            init_x[2] = val.ref[2]
        elseif (key == :energy)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 3 * (Cₐᵣₘ * Vdc^2)/ N
            end
        elseif (key == :zcc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 3*Vᴳd*Id/6/Vdc
            end
        elseif (key == :ccc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [0 0]
            end
            init_x[3] = val.ref[1]
            init_x[4] = val.ref[2]
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
        end
    end
    init_x[5] = Pdc/3/Vdc
    init_x[12] = Vdc

    vdc_position = 12

    exp = Expr(:block)

    # add PLL
    if in(:pll, keys(converter.controls))
        push!(exp.args, :(# θ = x[14]
                        T_θ = [cos(x[14]) -sin(x[14]); sin(x[14]) cos(x[14])];
                        I_θ = [cos(x[14]) sin(x[14]); -sin(x[14]) cos(x[14])];
                        T_2θ = [cos(-2x[14]) -sin(-2x[14]); sin(-2x[14]) cos(-2x[14])];
                        I_2θ = [cos(-2x[14]) sin(-2x[14]); -sin(-2x[14]) cos(-2x[14])];
                        (Vᴳd, Vᴳq) = T_θ * [inputs[2]; inputs[3]];
                        F[13] = -Vᴳq;
                        Δω = $(converter.controls[:pll].Kₚ) * (-Vᴳq) + $(converter.controls[:pll].Kᵢ) * x[13];
                        ω = $(converter.ω₀) + Δω;
                        F[14] = Δω;
                        ))
        index = 14
    else
        push!(exp.args, :(
                        T_θ = [1 0; 0 1];
                        I_θ = [1 0; 0 1];
                        T_2θ = [1 0; 0 1];
                        I_2θ = [1 0; 0 1];
                        Vᴳd = inputs[2];
                        Vᴳq = inputs[3];
                        ω = $(converter.ω₀)))
        index = 12
    end
    push!(exp.args, :((iΔd, iΔq) = T_θ * [x[1]; x[2]];
                      (iΣd, iΣq) = T_2θ * [x[3]; x[4]];))
    if in(:dc, keys(converter.controls))
        push!(exp.args, :(Vdc = x[$index+1]; Idc = inputs[1];
                    F[$index+1] = (Idc - 3*x[5]) / $Cₑ;
                    F[$index+2] = $(converter.controls[:dc].ref[1]) - Vdc;

                    iΔd_ref = -($(converter.controls[:dc].Kₚ) * ($(converter.controls[:dc].ref[1]) - Vdc) +
                                $(converter.controls[:dc].Kᵢ) * x[$index+2])))
        vdc_position = index + 1
        index += 2
    else
        push!(exp.args, :(Vdc = inputs[1]))
    end

    if in(:p, keys(converter.controls))
        # active power control
        push!(exp.args, :(
            # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
            P_ac = 1.5 * (Vᴳd * iΔd + Vᴳq * iΔq);
            # iΔd_ref = (Kp_Pac * (Pac_ref - Pac) + Ki_Pac * xiPac);
            iΔd_ref = ($(converter.controls[:p].Kₚ) * ($(converter.controls[:p].ref[1]) - P_ac) +
                        $(converter.controls[:p].Kᵢ) * x[$index+1]);
            F[$index+1] = $(converter.controls[:p].ref[1]) - P_ac))
        index += 1
    elseif !in(:dc, keys(converter.controls))
        push!(exp.args, :(
            iΔd_ref = $(converter.controls[:occ].ref[1])))
    end
    if in(:q, keys(converter.controls))
        # reactive power control
        push!(exp.args, :(

            # Q_ac = (3/2)*(Vgq*iΔd - Vgd*iΔq)
            Q_ac = 1.5 * (Vᴳq * iΔd - Vᴳd * iΔq);
            # iΔq_ref = -(Kp_Qac * (Qac_ref - Qac) + Ki_Qac * xiQac);
            iΔq_ref = -($(converter.controls[:q].Kₚ) * ($(converter.controls[:q].ref[1]) - Q_ac) +
                        $(converter.controls[:q].Kᵢ) * x[$index+1]);
            F[$index+1] = $(converter.controls[:q].ref[1]) - Q_ac))
        index += 1
    else
        push!(exp.args, :(
            iΔq_ref = $(converter.controls[:occ].ref[2])))
    end
    # add control equations
    for (key, val) in (converter.controls)
        if (key == :ccc)
            # circulating current control
            push!(exp.args, :(
                        (iΣd_ref, iΣq_ref) = T_θ * [$(converter.controls[:ccc].ref[1]); $(converter.controls[:ccc].ref[2])];
                        F[$index+1] = iΣd_ref - iΣd;
                        F[$index+2] = iΣq_ref - iΣq;
                        # vMΣd_ref = - Ki_Σ * xiΣd - Kp_Σ * (iΣd_ref -  iΣd) + 2*w*Larm*iΣq
                        vMΣd_ref_c = (-$(converter.controls[:ccc].Kᵢ) * x[$index+1] -
                                $(converter.controls[:ccc].Kₚ) * (iΣd_ref - iΣd) + 2 * ω * $Lₐᵣₘ * iΣq);
                        # vMΣq_ref = - Ki_Σ * xiΣq - Kp_Σ * (iΣq_ref -  iΣq) - 2*w*Larm*iΣd
                        vMΣq_ref_c = (-$(converter.controls[:ccc].Kᵢ) * x[$index+2] -
                                $(converter.controls[:ccc].Kₚ) * (iΣq_ref - iΣq) - 2 * ω * $Lₐᵣₘ * iΣd);
                        (vMΣd_ref, vMΣq_ref) = I_2θ * [vMΣd_ref_c; vMΣq_ref_c]))
            index += 2
        elseif (key == :zcc) && !in(:energy, keys(converter.controls))
            # zero current control
            push!(exp.args,
                        :(F[$index+1] = $(converter.controls[:zcc].ref[1]) - x[5];
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (Vdc/2 - $(converter.controls[:zcc].Kₚ) *
                            ($(converter.controls[:zcc].ref[1]) - x[5]) - $(converter.controls[:zcc].Kᵢ) * x[$index+1])))
            index += 1
        elseif (key == :energy)
            # zero energy control
            push!(exp.args, :(
                        # wΣz = 3*Carm * (vCΔd^2 + vCΔq^2 + vCΔZd^2 + vCΔZq^2 + vCΣd^2 + vCΣq^2 + 2*vCΣz^2)/(2*N)
                        wΣz = 3 * $Cₐᵣₘ / 2 / $N * (x[6]^2 + x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2 + x[11]^2 + 2x[12]^2);
                        F[$index+1] = $(converter.controls[:energy].ref[1]) - wΣz;
                        # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
                        P_ac = 1.5 * (Vᴳd * iΔd + Vᴳq * iΔq);
                        #iΣz_ref = (Kp_wΣ * (wΣz_ref - wΣz) + Ki_wΣ * xwΣz + Pac) / 3 / Vdc,
                        iΣz_ref = ($(converter.controls[:energy].Kₚ) * ($(converter.controls[:energy].ref[1]) - wΣz) +
                            $(converter.controls[:energy].Kᵢ) * x[$index+1] + P_ac) / 3 / Vdc;
                        F[$index+2] = iΣz_ref - x[5];
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (Vdc/2 - $(converter.controls[:zcc].Kₚ) *
                            (iΣz_ref - x[5]) -  $(converter.controls[:zcc].Kᵢ) * x[$index+2])))
            index += 2            
        elseif (key == :occ)
            # output current control
            push!(exp.args, :(
                        (iΔd_ref, iΔq_ref) = T_θ * [iΔd_ref; iΔq_ref];
                        F[$index+1] = iΔd_ref - iΔd;
                        F[$index+2] = iΔq_ref - iΔq;
                        # vMΔd_ref = Ki_Δ * xiΔd + Kp_Δ * (iΔd_ref -  iΔd) + w*L_eqac*iΔq + Vᴳd
                        vMΔd_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+1] +
                                    $(converter.controls[:occ].Kₚ) * (iΔd_ref - iΔd) + ω * $Lₑ * iΔq + Vᴳd);
                        # vMΔq_ref = Ki_Δ * xiΔq + Kp_Δ * (iΔq_ref -  iΔq) - w*Leqac*iΔd + Vᴳq
                        vMΔq_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+2] +
                                    $(converter.controls[:occ].Kₚ) * (iΔq_ref - iΔq) - ω * $Lₑ * iΔd + Vᴳq);
                        (vMΔd_ref, vMΔq_ref) = I_θ * [vMΔd_ref_c; vMΔq_ref_c]))
            index += 2
        end
    end
    push!(exp.args,
                :(vMΔZd_ref = 0;
                  vMΔZq_ref = 0;))
    if !in(:occ, keys(converter.controls))
        push!(exp.args,
                    :(vMΔd_ref = 0;
                      vMΔq_ref = 0;))
    end
    if !in(:ccc, keys(converter.controls))
        push!(exp.args,
                    :(vMΣd_ref = 0;
                      vMΣq_ref = 0;))
    end
    if !in(:zcc, keys(converter.controls))
        push!(exp.args, :(vMΣz_ref = Vdc/2))
    end

    # add time delays here, if there are controllers implemented
    if (converter.timeDelay != 0.0) && (in(:occ, keys(converter.controls)) || in(:ccc, keys(converter.controls)) || in(:zcc, keys(converter.controls)))
        push!(exp.args,
            :(
            T_ab_dq=0.5*[1 im;-im 1];# from alpha-beta to dq
            T_dq_ab=0.5*[1 -im;im 1];#from dq to alpha-beta
            ))
        if in(:occ, keys(converter.controls))
            push!(exp.args,
            :(
                timeDelayIn = [vMΔd_ref;vMΔq_ref];
                statesDelay = x[$index + 1 : $index + 2*$converter.padeOrderDen]; 
                timeDelayOut = timeDelayPadeMatrices($converter.padeOrderNum,$converter.padeOrderDen,$converter.timeDelay,length(timeDelayIn));
                
                A_delay = timeDelayOut[1];
                B_delay = timeDelayOut[2];
                C_delay = timeDelayOut[3];
                D_delay = timeDelayOut[4];
                F[$index + 1 : $index + 2*$converter.padeOrderDen] = A_delay*statesDelay + B_delay*timeDelayIn;
                # timeDelayOut = C_delay*statesDelay + D_delay*timeDelayIn;
                # Implement phase shifts by transforming the dq voltage references to alpha-beta
                vMΔ_ab_ref = (cos($converter.ω₀*$converter.timeDelay)-sin($converter.ω₀*$converter.timeDelay)*im)*(T_dq_ab*(C_delay*statesDelay + D_delay*timeDelayIn));
                vMΔ_dq_ref = real(T_ab_dq * conj(vMΔ_ab_ref) + conj(T_ab_dq) * vMΔ_ab_ref);
                vMΔd_ref = vMΔ_dq_ref[1];
                vMΔq_ref = vMΔ_dq_ref[2];
            ))
            index += 2*converter.padeOrderDen
        end
        if in(:ccc, keys(converter.controls))
            push!(exp.args,
            :(
                timeDelayIn = [vMΣd_ref;vMΣq_ref];
                statesDelay = x[$index + 1 : $index + 2*$converter.padeOrderDen]; 
                timeDelayOut = timeDelayPadeMatrices($converter.padeOrderNum,$converter.padeOrderDen,$converter.timeDelay,length(timeDelayIn));
                
                A_delay = timeDelayOut[1];
                B_delay = timeDelayOut[2];
                C_delay = timeDelayOut[3];
                D_delay = timeDelayOut[4];
                F[$index + 1 : $index + 2*$converter.padeOrderDen] = A_delay*statesDelay + B_delay*timeDelayIn;
                # timeDelayOut = C_delay*statesDelay + D_delay*timeDelayIn;
                # Phase shifts
                vMΣ_ab_ref = (cos(-2*$converter.ω₀*$converter.timeDelay)-sin(-2*$converter.ω₀*$converter.timeDelay)*im)*(T_dq_ab*(C_delay*statesDelay + D_delay*timeDelayIn));
                vMΣ_dq_ref = real(T_ab_dq * conj(vMΣ_ab_ref) + conj(T_ab_dq) * vMΣ_ab_ref);
                vMΣd_ref = vMΣ_dq_ref[1];
                vMΣq_ref = vMΣ_dq_ref[2];
            ))
            index += 2*converter.padeOrderDen
        end
        if in(:zcc, keys(converter.controls))
            push!(exp.args,
            :(
                timeDelayIn = vMΣz_ref;
                statesDelay = x[$index + 1 : $index + $converter.padeOrderDen]; 
                timeDelayOut = timeDelayPadeMatrices($converter.padeOrderNum,$converter.padeOrderDen,$converter.timeDelay,1);
                
                A_delay = timeDelayOut[1];
                B_delay = timeDelayOut[2];
                C_delay = timeDelayOut[3];
                D_delay = timeDelayOut[4];

                F[$index + 1 : $index + $converter.padeOrderDen] = A_delay*statesDelay + B_delay*timeDelayIn;
                timeDelayOut = dot(C_delay,statesDelay) + D_delay*timeDelayIn; # Get rid of 1-element array.
                vMΣz_ref = timeDelayOut;
            ))
            index += converter.padeOrderDen
        end
    end
    # add state variables
    # x = [iΔd, iΔq, iΣd, iΣq, iΣz, vCΔd, vCΔq, vCΔZd, vCΔZq, vCΣd, vCΣq, vCΣz] = [x[1], x[2], ...]
    # add corresponding differential equations [diΔd_dt, diΔq_dt, ...] = [F[1], F[2], ...]
    # m = [mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz], vM = [vMΔd, vMΔq, vMΔZd, vMΔZq,vMΣd, vMΣq, vMΣz]
    # vM_ref = [vMΔd_ref, vMΔq_ref, vMΔZd_ref, vMΔZq_ref, vMΣd_ref, vMΣq_ref, vMΣz_ref] = s
    push!(exp.args,
    :(
    # [mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz] = 2/Vᵈᶜ * [-vMΔd_ref, -vMΔq_ref, -vMΔZd_ref, -vMΔZq_ref, vMΣd_ref, vMΣq_ref, vMΣz_ref]
    (mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz) = 2/Vdc * [-vMΔd_ref; -vMΔq_ref; -vMΔZd_ref; -vMΔZq_ref; vMΣd_ref; vMΣq_ref; vMΣz_ref];

    vMΔd =(mΔq*x[11])/4 - (mΔd*x[12])/2 - (mΔd*x[10])/4 - (mΔZd*x[10])/4 + (mΔZq*x[11])/4 - (mΣd*x[6])/4 - (mΣz*x[6])/2 +
            (mΣq*x[7])/4 - (mΣd*x[8])/4 + (mΣq*x[9])/4;
    vMΔq =(mΔd*x[11])/4 + (mΔq*x[10])/4 - (mΔq*x[12])/2 - (mΔZd*x[11])/4 - (mΔZq*x[10])/4 + (mΣd*x[7])/4 + (mΣq*x[6])/4 -
            (mΣz*x[7])/2 - (mΣd*x[9])/4 - (mΣq*x[8])/4;
    vMΔZd =- (mΔd*x[10])/4 - (mΔq*x[11])/4 - (mΔZd*x[12])/2 - (mΣd*x[6])/4 - (mΣq*x[7])/4 - (mΣz*x[8])/2;
    vMΔZq =(mΔd*x[11])/4 - (mΔq*x[10])/4 - (mΔZq*x[12])/2 - (mΣd*x[7])/4 + (mΣq*x[6])/4 - (mΣz*x[9])/2;

    # vMΣd = (mΔd*vCΔd)/4 - (mΔq*vCΔq)/4 + (mΔd*vCΔZd)/4 + (mΔZd*vCΔd)/4 + (mΔq*vCΔZq)/4 + (mΔZq*vCΔq)/4 + (mΣd*vCΣz)/2 + (mΣz*vCΣd)/2;
    vMΣd =mΔd*x[6]/4 - mΔq*x[7]/4 + mΔd*x[8]/4 + mΔZd*x[6]/4 + mΔq*x[9]/4 + mΔZq*x[7]/4 + mΣd*x[12]/2 + mΣz*x[10]/2;
    # vMΣq = (mΔq*vCΔZd)/4 - (mΔq*vCΔd)/4 - (mΔd*vCΔZq)/4 - (mΔd*vCΔq)/4 + (mΔZd*vCΔq)/4 - (mΔZq*vCΔd)/4 + (mΣq*vCΣz)/2 + (mΣz*vCΣq)/2;
    vMΣq = mΔq*x[8]/4 - mΔq*x[6]/4 - mΔd*x[9]/4 - mΔd*x[7]/4 + mΔZd*x[7]/4 - mΔZq*x[6]/4 + mΣq*x[12]/2 + mΣz*x[11]/2;
    # vMΣz = (mΔd*vCΔd)/4 + (mΔq*vCΔq)/4 + (mΔZd*vCΔZd)/4 + (mΔZq*vCΔZq)/4 + (mΣd*vCΣd)/4 + (mΣq*vCΣq)/4 + (mΣz*vCΣz)/2;
    vMΣz = mΔd*x[6]/4 + mΔq*x[7]/4 + mΔZd*x[8]/4 + mΔZq*x[9]/4 + mΣd*x[10]/4 + mΣq*x[11]/4 + mΣz*x[12]/2;

    F[1] = -(inputs[2] - vMΔd + $Rₑ*x[1] + $Lₑ*$ω₀*x[2])/$Lₑ;                 # diΔd_dt =-(Vgd - vMΔd + Rₑ*iΔd + Lₑ*iΔq*w)/Lₑ
    F[2] = -(inputs[3] - vMΔq + $Rₑ*x[2] - $Lₑ*$ω₀*x[1])/$Lₑ;                 # diΔq_dt =-(Vgq - vMΔq + Rₑ*iΔq - Lₑ*iΔd*w)/Lₑ
    F[3] = -(vMΣd + $Rₐᵣₘ*x[3] - 2*$Lₐᵣₘ*$ω₀*x[4])/$Lₐᵣₘ;                                  # diΣd_dt =-(vMΣd + Rₐᵣₘ*iΣd - 2*Lₐᵣₘ*iΣq*w)/Lₐᵣₘ
    F[4] = -(vMΣq + $Rₐᵣₘ*x[4] + 2*$Lₐᵣₘ*$ω₀*x[3])/$Lₐᵣₘ;                                  # diΣq_dt =-(vMΣq + Rₐᵣₘ*iΣq + 2*Lₐᵣₘ*iΣd*w)/Lₐᵣₘ
    F[5] = -(vMΣz - Vdc/2 + $Rₐᵣₘ*x[5])/$Lₐᵣₘ;                                     # diΣz_dt =-(vMΣz - Vᵈᶜ/2 + Rₐᵣₘ*iΣz)/Lₐᵣₘ
    # dvCΔd_dt =(N*(iΣz*mΔd - (iΔq*mΣq)/4 + iΣd*(mΔd/2 + mΔZd/2) - iΣq*(mΔq/2 + mΔZq/2) + iΔd*(mΣd/4 + mΣz/2) - (2*Cₐᵣₘ*vCΔq*w)/N))/(2*Cₐᵣₘ)
    F[6] = ($N*(x[5]*mΔd - x[2]*mΣq/4 + x[3]*(mΔd/2 + mΔZd/2) - x[4]*(mΔq/2 + mΔZq/2) + x[1]*(mΣd/4 + mΣz/2) - 2*$Cₐᵣₘ*x[7]*$ω₀/$N))/2/$Cₐᵣₘ;
    # dvCΔq_dt =-(N*((iΔd*mΣq)/4 - iΣz*mΔq + iΣq*(mΔd/2 - mΔZd/2) + iΣd*(mΔq/2 - mΔZq/2) + iΔq*(mΣd/4 - mΣz/2) - (2*Cₐᵣₘ*vCΔd*w)/N))/(2*Cₐᵣₘ)
    F[7] = -($N*((x[1]*mΣq)/4 - x[5]*mΔq + x[4]*(mΔd/2 - mΔZd/2) + x[3]*(mΔq/2 - mΔZq/2) + x[2]*(mΣd/4 - mΣz/2) - 2*$Cₐᵣₘ*x[6]*$ω₀/$N))/2/$Cₐᵣₘ;
    # dvCΔZd_dt =(N*(iΔd*mΣd + 2*iΣd*mΔd + iΔq*mΣq + 2*iΣq*mΔq + 4*iΣz*mΔZd))/(8*Cₐᵣₘ) - 3*vCΔZq*w
    F[8] = ($N*(x[1]*mΣd + 2*x[3]*mΔd + x[2]*mΣq + 2*x[4]*mΔq + 4*x[5]*mΔZd))/(8*$Cₐᵣₘ) - 3*x[9]*$ω₀;
    # dvCΔZq_dt =3*vCΔZd*w + (N*(iΔq*mΣd - iΔd*mΣq + 2*iΣd*mΔq - 2*iΣq*mΔd + 4*iΣz*mΔZq))/(8*Cₐᵣₘ)
    F[9] = 3*x[8]*$ω₀ + ($N*(x[2]*mΣd - x[1]*mΣq + 2*x[3]*mΔq - 2*x[4]*mΔd + 4*x[5]*mΔZq))/(8*$Cₐᵣₘ);
    # dvCΣd_dt =(N*(iΣd*mΣz + iΣz*mΣd + iΔd*(mΔd/4 + mΔZd/4) - iΔq*(mΔq/4 - mΔZq/4) + (4*Cₐᵣₘ*vCΣq*w)/N))/(2*Cₐᵣₘ)
    F[10] = ($N*(x[3]*mΣz + x[5]*mΣd + x[1]*(mΔd/4 + mΔZd/4) - x[2]*(mΔq/4 - mΔZq/4) + 4*$Cₐᵣₘ*x[11]*$ω₀/$N))/(2*$Cₐᵣₘ);
    # dvCΣq_dt =-(N*(iΔq*(mΔd/4 - mΔZd/4) - iΣz*mΣq - iΣq*mΣz + iΔd*(mΔq/4 + mΔZq/4) + (4*Cₐᵣₘ*vCΣd*w)/N))/(2*Cₐᵣₘ)
    F[11] = -($N*(x[2]*(mΔd/4 - mΔZd/4) - x[5]*mΣq - x[4]*mΣz + x[1]*(mΔq/4 + mΔZq/4) + 4*$Cₐᵣₘ*x[10]*$ω₀/$N))/(2*$Cₐᵣₘ);
    # dvCΣz_dt =(N*(iΔd*mΔd + iΔq*mΔq + 2*iΣd*mΣd + 2*iΣq*mΣq + 4*iΣz*mΣz))/(8*Cₐᵣₘ)
    F[12] = ($N*(x[1]*mΔd + x[2]*mΔq + 2*x[3]*mΣd + 2*x[4]*mΣq + 4*x[5]*mΣz))/(8*$Cₐᵣₘ)))

    function f!(expr, F, x, inputs) # F derivative of state variable x state variable vector, inputs input vqlue expr equation of mmc
       f = eval(:((F,x,inputs) -> $expr))
       return Base.invokelatest(f, F,x,inputs)
    end

    if in(:dc, keys(converter.controls))
        vector_inputs = [Pdc/Vdc, Vᴳd, Vᴳq]
        init_x = [init_x; zeros(index-12,1)]
        init_x[vdc_position] = Vdc
    else
        vector_inputs = [Vdc, Vᴳd, Vᴳq]
        init_x = [init_x; zeros(index-12,1)]
    end

    g!(F,x) = f!(exp, F, x, vector_inputs)
    k = nlsolve(g!, init_x, autodiff = :forward, iterations = 100, ftol = 1e-6, xtol = 1e-3, method = :newton)
    converter.equilibrium = k.zero

    h(F,x) = f!(exp, F, x[1:end-3], x[end-2:end])
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index+3); h(F, x); return F)
    A = zeros(index+3,index+3)
    ForwardDiff.jacobian!(A, ha, [k.zero' vector_inputs'])
    converter.A = A[1:end-3, 1:end-3]
    converter.B = A[1:end-3, end-2:end]

    converter.C = zeros(3, size(converter.A,1))
    converter.C[2,1] = 1
    converter.C[3,2] = 1
    !in(:dc, keys(converter.controls)) ? converter.C[1,5] = 3 : converter.C[1, vdc_position] = 1
    converter.D = zeros(3,3)
end

function eval_parameters(converter :: MMC, s :: Complex)
    # numerical
    I = Matrix{Complex}(Diagonal([1 for dummy in 1:size(converter.A,1)]))
    Y = (converter.C*inv(s*I-converter.A))*converter.B + converter.D

    if in(:dc, keys(converter.controls))
        (m11, m12, m21, m22) = (Y[1,1], Y[1,2:3], Y[2:3,1], Y[2:3,2:3])
        m11 = 1/m11
        Y = [m11 -transpose(m12)*m11; m21*m11 m22-m21*m11*transpose(m12)]
    end

    return Y
end

function timeDelayPadeMatrices(padeOrderNum,padeOrderDen,t_delay,numberVars)

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
    # Convert from controllable canonical form to diagonal form - This does not work as such a state-space representation is complex-valued
    d_not_sorted=eigvals(Ad);
    v_not_sorted=eigvecs(Ad);
    d=zeros(ComplexF64,padeOrderDen,1);
    v=zeros(ComplexF64,padeOrderDen,padeOrderDen);
    matrix_ind=1;
    # Sort the eigenvalues and eigenvectors
    for i=1:padeOrderDen
        if sign(imag(d_not_sorted[i])==0) #Real eigenvalue, copy as is
            d[matrix_ind]=d_not_sorted[i];
            v[:,matrix_ind]=v_not_sorted[:,i];
            matrix_ind+=1;
        else # Complex eigenvalue
            if matrix_ind < padeOrderDen
                if sign(imag(d_not_sorted[i]))==-1
                    d[matrix_ind]    = real(d_not_sorted[i])     -imag(d_not_sorted[i])*im;
                    d[matrix_ind+1]  = real(d_not_sorted[i])     +imag(d_not_sorted[i])*im;
                    v[:,matrix_ind]  = real(v_not_sorted[:,i])   - imag(v_not_sorted[:,i])*im;
                    v[:,matrix_ind+1]= real(v_not_sorted[:,i])   + imag(v_not_sorted[:,i])*im;
                else
                    d[matrix_ind]    = real(d_not_sorted[i])     +imag(d_not_sorted[i])*im;
                    d[matrix_ind+1]  = real(d_not_sorted[i])     -imag(d_not_sorted[i])*im;
                    v[:,matrix_ind]  = real(v_not_sorted[:,i])   + imag(v_not_sorted[:,i])*im;
                    v[:,matrix_ind+1]= real(v_not_sorted[:,i])   - imag(v_not_sorted[:,i])*im;
                end
                matrix_ind+=2;
            end
        end
    end
    matrix_ind=1;
    T_inv=zeros(padeOrderDen,padeOrderDen);
    for i=1:padeOrderDen
        if imag(d[i])==0 #Real eigenvalue
            T_inv[:,matrix_ind]=v[:,i];
            matrix_ind+=1;
        else # Complex eigenvalue
            if matrix_ind < padeOrderDen
                T_inv[:,matrix_ind]  =real(v[:,i]);
                T_inv[:,matrix_ind+1]=imag(v[:,i]);
                matrix_ind+=2;
            end
        end
    end
    T=inv(T_inv);
    Ad=T*Ad*T_inv;
    Bd=T*Bd;
    Cd=Cd*T_inv;
    # Concatenate the Pade matrices: nDelta_d, nDelta_q, nSigma_d, nSigma_q, nSigma_z
    # A_Pade=cat(Ad,Ad,Ad,Ad,Ad;dims=[1,2]);
    # B_Pade=cat(Bd,Bd,Bd,Bd,Bd;dims=[1,2]);
    # C_Pade=cat(Cd,Cd,Cd,Cd,Cd;dims=[1,2]);
    # D_Pade=cat(Dd,Dd,Dd,Dd,Dd;dims=[1,2]);
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
