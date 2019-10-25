export mmc

include("controller.jl")

@with_kw mutable struct MMC
    f₀ :: Union{Int, Float64}  = 50             # nominal frequency
    ω₀ :: Union{Int, Float64} = 100*π

    S_base :: Union{Int, Float64} = 1000e6      # power rating of the system
    V_base :: Union{Int, Float64} = 333e3       # base voltage rating
    I_base :: Union{Int, Float64} = 0           # base current rating
    Vᵈᶜ :: Union{Int, Float64} = 640e3          # DC-bus voltage
    Vᴳd :: Union{Int, Float64} = 0
    Vᴳq :: Union{Int, Float64} = 0

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
end

"""
    function mmc(;args...)
It constructs MMC operating both as a rectifier and an inverter. MMC is constructed as a struct with the
following fields.
```
f₀ :: Union{Int, Float64}  = 50             # nominal frequency
ω₀ :: Union{Int, Float64} = 100*π

S_base :: Union{Int, Float64} = 1000e6      # power rating of the system
V_base :: Union{Int, Float64} = 333e3       # base voltage rating
I_base :: Union{Int, Float64} = 0           # base current rating
Vₘ :: Union{Int, Float64} = 0               # AC-grid voltage in abc frame [V]
Vᵈᶜ :: Union{Int, Float64} = 640e3          # DC-bus voltage
Vᴳd :: Union{Int, Float64} = 0
Vᴳq :: Union{Int, Float64} = 0

Lₐᵣₘ :: Union{Int, Float64}  = 50e-3        # arm inductance [H]
Rₐᵣₘ :: Union{Int, Float64}  = 1.07         # equivalent arm resistance
Cₐᵣₘ :: Union{Int, Float64}  = 10e-3        # capacitance per submodule [F]
N :: Int = 400                              # number of submodules per arm

Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the phase reactor [H]
Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the phase reactor [Ω]

controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
equilibrium :: Array{Union{Int, Float64}} = [0]
A :: Array{Union{Int, Float64}} = [0]
B :: Array{Union{Int, Float64}} = [0]
C :: Array{Union{Int, Float64}} = [0]
D :: Array{Union{Int, Float64}} = [0]

```

The constructed MMC has 3 pins on the AC side: `2.1`, `2.2`, and 2 pins on its
DC-side: `1.1` and `1.2`.

The component is described using two ABCD parameters' matrix sized 4×4.
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


    converter.Vᴳq = 0

    Lₑ = converter.Lₐᵣₘ / 2 + converter.Lᵣ
    Rₑ = converter.Rₐᵣₘ / 2 + converter.Rᵣ
    N = converter.N
    Lₐᵣₘ = converter.Lₐᵣₘ
    Rₐᵣₘ = converter.Rₐᵣₘ
    Cₐᵣₘ = converter.Cₐᵣₘ
    ω₀ = converter.ω₀

    if (converter.V_base != 0)
        converter.Vᴳd = converter.V_base * sqrt(2/3)
        converter.I_base = 2converter.S_base / 3 / converter.V_base / sqrt(2/3)
    end

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
                val.ref = [-converter.I_base 0]
            end
            init_x[1] = val.ref[1]
            init_x[2] = val.ref[2]
        elseif (key == :energy)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 3 * (Cₐᵣₘ * converter.Vᵈᶜ^2)/ N
            end
        elseif (key == :zcc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = - 3*converter.Vᴳd*converter.I_base/6/converter.Vᵈᶜ
            end
            init_x[5] = val.ref[1]
        elseif (key == :ccc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [0 0]
            end
            init_x[3] = val.ref[1]
            init_x[4] = val.ref[2]
        elseif (key == :power)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [-converter.S_base 0]
            end
        end
    end
    init_x[12] = converter.Vᵈᶜ

    exp = Expr(:block)

    # add PLL
    if in(:pll, keys(converter.controls))
        push!(exp.args, :(# θ = x[14]
                        T_θ = [cos(x[14]) -sin(x[14]); sin(x[14]) cos(x[14])];
                        I_θ = [cos(x[14]) sin(x[14]); -sin(x[14]) cos(x[14])];
                        T_2θ = [cos(-2x[14]) -sin(-2x[14]); sin(-2x[14]) cos(-2x[14])];
                        I_2θ = [cos(-2x[14]) sin(-2x[14]); -sin(-2x[14]) cos(-2x[14])];
                        (Vᴳd, Vᴳq) = T_θ * [inputs[1]; inputs[2]];
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
                        Vᴳd = inputs[1];
                        Vᴳq = inputs[2];
                        ω = $(converter.ω₀)))
        index = 12
    end

    # add control equations
    for (key, val) in (converter.controls)
        if (key == :occ) && !in(:power, keys(converter.controls))
            # output current control
            push!(exp.args, :((iΔd, iΔq) = T_θ * [x[1]; x[2]];
                        (iΔd_ref, iΔq_ref) = T_θ * [$(converter.controls[:occ].ref[1]); $(converter.controls[:occ].ref[2])];
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
        elseif (key == :ccc)
            # circulating current control
            push!(exp.args, :((iΣd, iΣq) = T_2θ * [x[3]; x[4]];
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
                        :(F[$index+1] = $(converter.controls[:zcc].ref[1]) - x[5],
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (inputs[3]/2 - $(converter.controls[:zcc].Kₚ) *
                            ($(converter.controls[:zcc].ref[1]) - x[5]) - $(converter.controls[:zcc].Kᵢ) * x[$index+1])))
            index += 1
        elseif (key == :energy)
            # zero energy control
            push!(exp.args,
                        :((iΔd, iΔq) = T_θ * [x[1]; x[2]];
                        # wΣz = 3*Carm * (vCΔd^2 + vCΔq^2 + vCΔZd^2 + vCΔZq^2 + vCΣd^2 + vCΣq^2 + 2*vCΣz^2)/(2*N)
                        wΣz = 3 * $Cₐᵣₘ / 2 / $N * (x[6]^2 + x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2 + x[11]^2 + 2x[12]^2);
                        F[$index+1] = $(converter.controls[:energy].ref[1]) - wΣz;
                        # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
                        P_ac = 1.5 * (Vᴳd * iΔd + Vᴳq * iΔq);
                        #iΣz_ref = (Kp_wΣ * (wΣz_ref - wΣz) + Ki_wΣ * xwΣz + Pac) / 3 / Vdc,
                        iΣz_ref = ($(converter.controls[:energy].Kₚ) * ($(converter.controls[:energy].ref[1]) - wΣz) +
                            $(converter.controls[:energy].Kᵢ) * x[$index+1] + P_ac) / 3 / inputs[3];
                        F[$index+2] = iΣz_ref - x[5];
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (inputs[3]/2 - $(converter.controls[:zcc].Kₚ) *
                            (iΣz_ref - x[5]) -  $(converter.controls[:zcc].Kᵢ) * x[$index+2])))
            index += 2
        elseif (key == :power)
            # active and reactive power control
            push!(exp.args, :((iΔd, iΔq) = T_θ * [x[1]; x[2]];
                        # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
                        P_ac = 1.5 * (Vᴳd * iΔd + Vᴳq * iΔq);
                        # Q_ac = (3/2)*(Vgq*iΔd - Vgd*iΔq)
                        Q_ac = 1.5 * (Vᴳq * iΔd - Vᴳd * iΔq);
                        # iΔd_ref = (Kp_Pac * (Pac_ref - Pac) + Ki_Pac * xiPac);
                        iΔd_ref = ($(converter.controls[:power].Kₚ) * ($(converter.controls[:power].ref[1]) - P_ac) +
                                    $(converter.controls[:power].Kᵢ) * x[$index+1]);
                        # iΔq_ref = -(Kp_Qac * (Qac_ref - Qac) + Ki_Qac * xiQac);
                        iΔq_ref = -($(converter.controls[:power].Kₚ) * ($(converter.controls[:power].ref[2]) - Q_ac) +
                                    $(converter.controls[:power].Kᵢ) * x[$index+2]);
                        F[$index+1] = $(converter.controls[:power].ref[1]) - P_ac;
                        F[$index+2] = $(converter.controls[:power].ref[2]) - Q_ac;

                        # OCC
                        F[$index+3] = iΔd_ref - iΔd;
                        F[$index+4] = iΔq_ref - iΔq;
                        # vMΔd_ref = Ki_Δ * xiΔd + Kp_Δ * (iΔd_ref -  iΔd) + w*L_eqac*iΔq + Vᴳd
                        vMΔd_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+3] +
                                    $(converter.controls[:occ].Kₚ) * (iΔd_ref - iΔd) + ω * $Lₑ * iΔq + Vᴳd);
                        # vMΔq_ref = Ki_Δ * xiΔq + Kp_Δ * (iΔq_ref -  iΔq) - w*Leqac*iΔd + Vᴳq
                        vMΔq_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+4] +
                                    $(converter.controls[:occ].Kₚ) * (iΔq_ref - iΔq) - ω * $Lₑ * iΔd + Vᴳq);
                        (vMΔd_ref, vMΔq_ref) = I_θ * [vMΔd_ref_c; vMΔq_ref_c]))
            index += 4
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
        push!(exp.args, :(vMΣz_ref = inputs[3]))
    end

    # add state variables
    # x = [iΔd, iΔq, iΣd, iΣq, iΣz, vCΔd, vCΔq, vCΔZd, vCΔZq, vCΣd, vCΣq, vCΣz] = [x[1], x[2], ...]
    # add corresponding differential equations [diΔd_dt, diΔq_dt, ...] = [F[1], F[2], ...]
    # m = [mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz], vM = [vMΔd, vMΔq, vMΔZd, vMΔZq,vMΣd, vMΣq, vMΣz]
    # vM_ref = [vMΔd_ref, vMΔq_ref, vMΔZd_ref, vMΔZq_ref, vMΣd_ref, vMΣq_ref, vMΣz_ref] = s
    push!(exp.args,
    :(
    # [mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz] = 2/Vᵈᶜ * [-vMΔd_ref, -vMΔq_ref, -vMΔZd_ref, -vMΔZq_ref, vMΣd_ref, vMΣq_ref, vMΣz_ref]
    (mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz) = 2/inputs[3] * [-vMΔd_ref; -vMΔq_ref; -vMΔZd_ref; -vMΔZq_ref; vMΣd_ref; vMΣq_ref; vMΣz_ref];

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

    F[1] = -(inputs[1] - vMΔd + $Rₑ*x[1] + $Lₑ*$ω₀*x[2])/$Lₑ;                 # diΔd_dt =-(Vgd - vMΔd + Rₑ*iΔd + Lₑ*iΔq*w)/Lₑ
    F[2] = -(inputs[2] - vMΔq + $Rₑ*x[2] - $Lₑ*$ω₀*x[1])/$Lₑ;                 # diΔq_dt =-(Vgq - vMΔq + Rₑ*iΔq - Lₑ*iΔd*w)/Lₑ
    F[3] = -(vMΣd + $Rₐᵣₘ*x[3] - 2*$Lₐᵣₘ*$ω₀*x[4])/$Lₐᵣₘ;                                  # diΣd_dt =-(vMΣd + Rₐᵣₘ*iΣd - 2*Lₐᵣₘ*iΣq*w)/Lₐᵣₘ
    F[4] = -(vMΣq + $Rₐᵣₘ*x[4] + 2*$Lₐᵣₘ*$ω₀*x[3])/$Lₐᵣₘ;                                  # diΣq_dt =-(vMΣq + Rₐᵣₘ*iΣq + 2*Lₐᵣₘ*iΣd*w)/Lₐᵣₘ
    F[5] = -(vMΣz - inputs[3]/2 + $Rₐᵣₘ*x[5])/$Lₐᵣₘ;                                     # diΣz_dt =-(vMΣz - Vᵈᶜ/2 + Rₐᵣₘ*iΣz)/Lₐᵣₘ
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

    function f!(expr, F, x, inputs)
       f = eval(:((F,x,inputs) -> $expr))
       return Base.invokelatest(f, F,x,inputs)
    end

    vector_inputs = [converter.Vᴳd, converter.Vᴳq, converter.Vᵈᶜ]
    init_x = [init_x; zeros(index-12,1)]
    g!(F,x) = f!(exp, F, x, vector_inputs)

    k = nlsolve(g!, init_x, autodiff = :forward, iterations = 10000, ftol = 1e-6, method = :newton)
    converter.equilibrium = k.zero[1:12]

    h(F,x) = f!(exp, F, x[1:end-3], x[end-2:end])
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index+3); h(F, x); return F)
    A = zeros(index+3,index+3)
    ForwardDiff.jacobian!(A, ha, [k.zero' vector_inputs'])
    converter.A = A[1:12, 1:12]
    converter.B = A[1:12, end-2:end]

    converter.C = zeros(3, 12)
    converter.C[1,1] = 1
    converter.C[2,2] = 1
    converter.C[3,5] = 1
    converter.D = zeros(3,3)

    elem = Element(input_pins = 2, output_pins = 2, element_value = converter)
end

function create_abcd!(element :: Element, mmc :: MMC)
    n = 4
    symbol = "mmc"
    abcd_tf =  reshape([Basic(string(symbol,i)) for i in 1:n^2], n, n)
    return abcd_tf
end

function eval_parameters(converter :: MMC, s :: Complex)
    # numerical
    I = Matrix{Complex}(Diagonal([1 for dummy in 1:size(converter.A,1)]))
    Y = (converter.C*inv(s*I-converter.A))*converter.B + converter.D
    return Y
end

function eval_abcd(element :: Element, mmc :: MMC, complex_s :: Complex)
    n = 4
    Y = eval_parameters(mmc,complex_s)

    abcd = [-Y[1,1]/2/Y[1,3]  -Y[1,2]/2/Y[1,3]   1/2/Y[1,3]   0;
            Y[2,1]/2/Y[2,3]   Y[2,2]/2/Y[2,3]    0            -1/2/Y[2,3];
            (3Y[3,1]-3Y[1,1]*Y[3,3]/Y[1,3]) (3Y[3,2]-3Y[3,3]*Y[1,2]/Y[1,3]) 3Y[3,3]/Y[1,3] 0;
            -(3Y[3,1]-3Y[2,1]*Y[3,3]/Y[2,3]) -(3Y[3,2]-3Y[3,3]*Y[2,2]/Y[2,3]) 0 -3Y[3,3]/Y[2,3]]

    # println(complex_s)
    # println(Y)
    # println(abcd)
    return abcd
end


function save_data(mmc :: MMC, file_name :: String, omegas)
    open(string(file_name, "_y.txt"), "w") do f
        for omega in omegas
            Y = eval_parameters(mmc, 1im*omega)
            writedlm(f, [omega reshape(Y, 1, length(Y))], ",")
        end
    end
end

function plot_data(mmc :: MMC, omegas)
    Y_ac = []
    Y_dc = []
    Y_acdc = []
    for omega in omegas
        Y = (eval_parameters(mmc, 1im*omega))
        push!(Y_ac, Y[1:2,1:2])
        push!(Y_dc, Y[3,3])
        push!(Y_acdc, [Y[1:2,3] Y[3,1:2]])
    end

    bode(Y_ac, omega = omegas, titles = ["Y_{dd}" "Y_{dq}"; "Y_{qd}" "Y_{qq}"])
    save_plot("files/ac_side")
    bode(Y_dc, omega = omegas, titles = ["Y_{zz}"])
    save_plot("files/dc_side")
    bode(Y_acdc, omega = omegas, titles = ["Y_{dz}" "Y_{zd}"; "Y_{qz}" "Y_{zq}"])
    save_plot("files/acdc_side")
end
