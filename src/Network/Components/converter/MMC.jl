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
    Rₐᵣₘ :: Union{Int, Float64}  = 1.07         # equivalent arm resistance[Ω]
    Cₐᵣₘ :: Union{Int, Float64}  = 10e-3        # capacitance per submodule [F]
    N :: Int = 400                              # number of submodules per arm

    Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the converter transformer at the converter side [H]
    Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the converter transformer at the converter side [Ω]


    gfm :: Bool = false                         # If variable is false, GFL is assumed
    controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
    equilibrium :: Array{Union{Int, Float64}} = [0]
    # TODO: The state-space matrices were previously of type Complex. See if this will create any incompatibility issues.
    A :: Array{Float64} = [0]
    B :: Array{Float64} = [0]
    C :: Array{Float64} = [0]
    D :: Array{Float64} = [0]

    timeDelay :: Float64 = 0
    padeOrderNum :: Int = 0
    padeOrderDen :: Int = 0

    vACbase_LL_RMS :: Union{Int, Float64} = 380 # Voltage base in kV
    Sbase :: Union{Int, Float64} = 1000 # Power base in MW
    vDCbase :: Union{Int, Float64} = 640 # DC voltage base in kV

    turnsRatio :: Union{Int, Float64} = 1 # Turns ratio of the converter transformer, converter side/AC side

    vACbase :: Float64 = 0 # AC voltage base for impedance/admittance calculation
    iACbase :: Float64 = 0 # AC current base for impedance/admittance calculation
    iDCbase :: Float64 = 0 # DC current base for impedance/admittance calculation
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

    # TODO: Error handling for further unreasonable controller implementations, such as....
#=     
    if haskey(converter.controls,....) && haskey(converter.controls,....)

        println("Two power controllers are defined! Program will crash!")
        return

    end 
 =#


    elem = Element(input_pins = 1, output_pins = 2, element_value = converter)
end

function update_mmc(converter :: MMC, Vm, θ, Pac, Qac, Vdc, Pdc)
    

    wbase = 100*pi
    vAC_base = converter.vACbase_LL_RMS*sqrt(2/3)
    vDC_base = converter.vDCbase
    Sbase = converter.Sbase
    iAC_base = 2*Sbase/3/vAC_base
    iDC_base = Sbase/vDC_base
    zAC_base = (3/2)*vAC_base^2/Sbase
    zDC_base = vDC_base/iDC_base
    lAC_base = zAC_base/wbase
    lDC_base = zDC_base/wbase
    cbase = 1/wbase/zDC_base

    converter.vACbase = vAC_base
    converter.iACbase = iAC_base
    converter.vDCbase = vDC_base
    converter.iDCbase = iDC_base

    Lₑ = (converter.Lₐᵣₘ / 2 + converter.Lᵣ) / lAC_base
    Rₑ = (converter.Rₐᵣₘ / 2 + converter.Rᵣ) / zAC_base
    N = converter.N
    Lₐᵣₘ = converter.Lₐᵣₘ / lDC_base
    Rₐᵣₘ = converter.Rₐᵣₘ / zDC_base
    Cₐᵣₘ = converter.Cₐᵣₘ / cbase
    Cₑ = 1e-6/ cbase
    ω₀ = converter.ω₀

    baseConv1 = vAC_base/vDC_base;# AC to DC voltage
    baseConv2 = vDC_base/vAC_base;# DC to AC voltage
    baseConv3 = iAC_base/iDC_base;# AC to DC current

    converter.Vₘ = Vm
    converter.θ = θ
    converter.Vᵈᶜ = Vdc
    converter.P = Pac
    converter.Q = Qac
    converter.P_dc = Pdc # Has the same sign as Pac

    Vm /= vAC_base
    Vdc /= vDC_base
    Pac /= Sbase
    Qac /= Sbase
    Pdc /= Sbase
    
    Vᴳd = Vm * cos(θ)
    Vᴳq = -Vm * sin(θ)
    # Vᴳd = Vm
    # Vᴳq = 0

    Id = ((Vᴳd*converter.turnsRatio * Pac + Vᴳq*converter.turnsRatio * Qac) / ((Vᴳd*converter.turnsRatio)^2 + (Vᴳq*converter.turnsRatio)^2)) 
    Iq = ((Vᴳq*converter.turnsRatio * Pac - Vᴳd*converter.turnsRatio * Qac) / ((Vᴳd*converter.turnsRatio)^2 + (Vᴳq*converter.turnsRatio)^2)) 

    # setup control parameters and equations
    init_x = zeros(12, 1) #pu

    # TODO: These have to be updated to be compliant with the PU model!
    for (key, val) in (converter.controls)
        

        # fix coefficients. Only valid for PI controllers and specific plants
        if isa(val, PI_control)
            
            if (val.Kₚ == 0) && (val.Kᵢ == 0)
                if (key == :occ)                            # pole placement
                    val.Kᵢ = Lₑ * val.bandwidth^2
                    val.Kₚ = 2 * val.ζ * val.bandwidth * Lₑ - Rₑ
    
                elseif (key == :ccc) || (key == :zcc)       # pole placement
                    val.Kᵢ = converter.Lₐᵣₘ * val.bandwidth^2
                    val.Kₚ = 2 * val.ζ * val.bandwidth * converter.Lₐᵣₘ - converter.Rₐᵣₘ
                end
            end 


        end

        # TODO: Decide on default reference for the grid forming controls and fill here
        if (key == :occ)
            if (length(val.ref) == 1) && (val.ref[1] == 0) # no reference existent 
                val.ref = [Id Iq]
            end
            init_x[1] = val.ref[1]
            init_x[2] = val.ref[2]
        elseif (key == :energy)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 1 #3 * (Cₐᵣₘ * Vdc^2)/ N #per unit
            end
        elseif (key == :zcc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 3*Vᴳd*Id/6/Vdc #per-unit 
            end
        elseif (key == :ccc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [0 0]
            end
            init_x[3] = val.ref[1]
            init_x[4] = val.ref[2]
        elseif (key == :p)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Pac*Sbase]                    # Comes from Powerflow definition P, conversion to SI. Will be converted back, later in the code 
            end
        elseif (key == :q)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Qac*Sbase]                    # Comes from Powerflow definition Q, conversion to SI. Will be converted back, later in the code 
            end
        elseif (key == :dc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Vdc]                          # Comes from Powerflow definition Vdc, in pu.
            end
        elseif (key == :vac) || (key == :vac_supp)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                # TODO: Decide what the default reference has to be
                # val.ref = [Vm*converter.turnsRatio]
                # val.ref = [1.0]
            end
        end
    end
    Idc_in = Pdc/Vdc
    init_x[5] = Pdc/3/Vdc
    init_x[12] = Vdc

    vdc_position = 12
    
    exp = Expr(:block) # Start construction of the state-space equations

    # TODO: Add code before to define transformation matrices generally, like in TLC.jl. Determine the statenumber based on if statements
    index = 12;
    if in(:pll, keys(converter.controls))
        

        if (((converter.controls[:pll].n_f)) == 1 ) && (converter.gfm) #1st order butterworth only Francesco's model so far 

            push!(exp.args, :(
            Δθ_pll=x[14];
            Δθᵥ = x[19];
            T_θ_pll = [cos(Δθ_pll) -sin(Δθ_pll); sin(Δθ_pll) cos(Δθ_pll)];
            T_θ = [cos(Δθᵥ) -sin(Δθᵥ); sin(Δθᵥ) cos(Δθᵥ)];
            I_θ = [cos(Δθᵥ) sin(Δθᵥ); -sin(Δθᵥ) cos(Δθᵥ)];
            T_2θ = [cos(-2Δθᵥ) -sin(-2Δθᵥ); sin(-2Δθᵥ) cos(-2Δθᵥ)];
            I_2θ = [cos(-2Δθᵥ) sin(-2Δθᵥ); -sin(-2Δθᵥ) cos(-2Δθᵥ)];
            (Vᴳd_pll, Vᴳq_pll) = T_θ_pll * [inputs[2] * $converter.turnsRatio; inputs[3] * $converter.turnsRatio];  # Vgd: Input 1 and Vgq: Input 2
            F[13] = x[15]*$(converter.controls[:pll].Kᵢ);
            Δω = ($(converter.controls[:pll].Kₚ) * x[15]) + x[13]; #Delta omega_pll [pu]
            ω = $(converter.ω₀)/$wbase + Δω;
            F[14] = $wbase*Δω; 
            F[15] = -(Vᴳq_pll+ x[15])*$(converter.controls[:pll].ω_f);
            (Vᴳd, Vᴳq) = T_θ * [inputs[2] * $converter.turnsRatio; inputs[3] * $converter.turnsRatio]; #Vd_grid input 2 and Vq_grid input 3 both expressed in the grid frame and at grid side
            ))
            init_x = [init_x;zeros(index-length(init_x))];
            init_x = [init_x;0;θ;0];
            index += 3;
            
        
        
        # TODO: Implement PLL with arbitary butterworth filtering
        else
            push!(exp.args, :(# θ = x[14]
            T_θ = [cos(x[14]) -sin(x[14]); sin(x[14]) cos(x[14])];
            I_θ = [cos(x[14]) sin(x[14]); -sin(x[14]) cos(x[14])];
            T_2θ = [cos(-2x[14]) -sin(-2x[14]); sin(-2x[14]) cos(-2x[14])];
            I_2θ = [cos(-2x[14]) sin(-2x[14]); -sin(-2x[14]) cos(-2x[14])];
            (Vᴳd, Vᴳq) = T_θ * [inputs[2] * $converter.turnsRatio; inputs[3] * $converter.turnsRatio];  # Vgd: Input 1 and Vgq: Input 2
            F[13] = -Vᴳq*$(converter.controls[:pll].Kᵢ);
            Δω = $(converter.controls[:pll].Kₚ) * (-Vᴳq) + x[13];
            ω = $(converter.ω₀)/$wbase + Δω;
            F[14] = $wbase*Δω;

            ))

            index += 2;
        end

                        
        
    else
        push!(exp.args, :(
                        T_θ = [1 0; 0 1];
                        I_θ = [1 0; 0 1];
                        T_2θ = [1 0; 0 1];
                        I_2θ = [1 0; 0 1];
                        Vᴳd = inputs[2] * $converter.turnsRatio;
                        Vᴳq = inputs[3] * $converter.turnsRatio;
                        ω = $(converter.ω₀)))
        
    end

    push!(exp.args, :((iΔd, iΔq) = T_θ * [x[1]; x[2]]; # Currents in grid dq frame defined: x1 and x2, see circuit equations far below 
    (iΣd, iΣq) = T_2θ * [x[3]; x[4]];
    Vdc = inputs[1])) #Vdc voltage input 1 


    if in(:p, keys(converter.controls)) # P control
    
        converter.controls[:p].ref[1] /= Sbase # conversion back to pu again
        push!(exp.args, :(
                P_ac = (Vᴳd * iΔd + Vᴳq * iΔq);))

       if isa(converter.controls[:p], VSE) && (converter.gfm) #VSE with grid-forming converter

            # TODO: Implement arbitary order of butterworth filtering
            if ((converter.controls[:p].n_f)) == 2 # 2nd order butterworth P filtering like Francesco's model
            push!(exp.args, :(
                F[$index+1]=x[$index+2];  #eta_1
                F[$index+2]=-($(converter.controls[:p].ω_f))^2 *x[$index+1]-sqrt(2)*$(converter.controls[:p].ω_f)*x[$index+2]+($(converter.controls[:p].ω_f))^2 *P_ac; #eta_2
                P_ac_f=x[$index+1]))
                
            init_x = [init_x;zeros(index-length(init_x))];
            init_x = [init_x;Pac;0];
            index += 2
            
            else #No P filtering
            

            end
        # Swing equation
        push!(exp.args, :(
            ωᵥ= x[$index+1]; # Actually delta Omega_VSM, relative angle [pu]. With initial conditions it becomes absolute: ωᵥ = Δωᵥ + ωᵥ(0)
            F[$index+1] =($(converter.controls[:p].ref[1]) - P_ac_f - $(converter.controls[:p].K_d)*(ωᵥ-(Δω+1)) -  $(converter.controls[:p].K_ω)*(ωᵥ - $(converter.controls[:p].ref_ω)))/(2*$(converter.controls[:p].H)) ;
            F[$index+2] =$wbase*(x[$index+1]-1); 
            ))
        init_x = [init_x;zeros(index-length(init_x))];
        init_x = [init_x;1;θ]; # Setting Delta Theta_VSM could be a problem for larger angles?
        index += 2

       elseif isa(converter.controls[:p], VSE) && !(converter.gfm) #VSE with grid-following converter



       elseif isa(converter.controls[:p], PI_control) && !(converter.gfm) #P control with grid-following converter
        
            # active power control
            # TODO: Implement arbitary order of butterworth filtering
            push!(exp.args, :(
                # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
                #P_ac = (Vᴳd * iΔd + Vᴳq * iΔq);
                # iΔd_ref = (Kp_Pac * (Pac_ref - Pac) + Ki_Pac * xiPac);
                iΔd_ref = ($(converter.controls[:p].Kₚ) * ($(converter.controls[:p].ref[1]) - P_ac) +
                            x[$index+1]);
                F[$index+1] = $(converter.controls[:p].Kᵢ) *($(converter.controls[:p].ref[1]) - P_ac)))
            index += 1

       end


       

    elseif in(:dc, keys(converter.controls)) # DC voltage control
        # With a DC voltage filter
        # push!(exp.args, :(
        #         F[$index+1] = $(converter.controls[:dc].ω_f) *(Vdc - x[$index+1]);
        #         F[$index+2] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - x[$index+1]);
        #             iΔd_ref = -($(converter.controls[:dc].Kₚ) * ($(converter.controls[:dc].ref[1]) - x[$index+1]) +
        #                          x[$index+2])))
        # vdc_position = index + 1
        # epsilon_vdc_index = index + 2
        # index += 2

        push!(exp.args, :(
                F[$index+1] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - Vdc);
                    iΔd_ref = -($(converter.controls[:dc].Kₚ) * ($(converter.controls[:dc].ref[1]) - Vdc) +
                                 x[$index+1])))
        vdc_position = index + 1
        epsilon_vdc_index = index + 1
        index += 1
        # Additional state variable to represent the DC voltage
        # push!(exp.args, :(
        #         Vdc = x[$index+1]; Idc = inputs[1];
        #         F[$index+1] = $wbase * (Idc - 3*x[5]) / $Cₑ;
        #         F[$index+2] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - Vdc);
        #             iΔd_ref = -($(converter.controls[:dc].Kₚ) * ($(converter.controls[:dc].ref[1]) - Vdc) +
        #                          x[$index+2])))
        # vdc_position = index + 1
        # index += 2
    else
        push!(exp.args, :(
            iΔd_ref = $Pac))
    end
    if in(:q, keys(converter.controls))
        converter.controls[:q].ref[1] /= -Sbase # The minus sign corrects for the Q convention used in the model.
        
        # Voltage droop control
        if in(:vac_supp, keys(converter.controls))
            push!(exp.args, :(
                Vᴳ_mag = sqrt(Vᴳd^2+Vᴳq^2);
                Δq_unf = $(converter.controls[:vac_supp].Kₚ)*($(converter.controls[:vac_supp].ref[1])-Vᴳ_mag);
                F[$index+1] = $(converter.controls[:vac_supp].ω_f) *(Δq_unf - x[$index+1]);
                q_ref = $(converter.controls[:q].ref[1]) - x[$index+1]))
            index +=1
        else
            push!(exp.args, :(
                q_ref = $(converter.controls[:q].ref[1])))
        end
        push!(exp.args, :(
            Q_ac =  (Vᴳq * iΔd - Vᴳd * iΔq);))
        # Reactive power control
        if (converter.gfm) # Q control for GFM

            # TODO: Implement arbitary order of butterworth filtering
            if ((converter.controls[:q].n_f)) == 2 # Q Filtering second order butterworth
                push!(exp.args, :(
                    F[$index+1]=x[$index+2];  #eta_1
                    F[$index+2]=-($(converter.controls[:q].ω_f))^2 *x[$index+1]-sqrt(2)*$(converter.controls[:q].ω_f)*x[$index+2]+($(converter.controls[:q].ω_f))^2 *Q_ac; #eta_2
                    Q_ac_f=x[$index+1]))
                init_x = [init_x;zeros(index-length(init_x))];
                init_x = [init_x;Qac;0];
                index += 2
                
            else #TODO: Implement without filtering
                #Q_ac_f = Q_ac
            end
            push!(exp.args, :(
                F[$index+1]=$(converter.controls[:q].Kᵢ)*(q_ref - Q_ac_f);
                Vⱽd_ref=-($(converter.controls[:q].Kₚ)*(q_ref - Q_ac_f) + x[$index+1]);
            ))
            index += 1
            
            





        else # TODO: Implement arbitary order of butterworth filtering
            push!(exp.args, :(

            # Q_ac = (3/2)*(Vgq*iΔd - Vgd*iΔq)
            # iΔq_ref = -(Kp_Qac * (Qac_ref - Qac) + Ki_Qac * xiQac);
            iΔq_ref = -($(converter.controls[:q].Kₚ) * (q_ref - Q_ac) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:q].Kᵢ) *(q_ref - Q_ac)))
            index += 1

        end
        
    
    elseif in(:vac, keys(converter.controls))
        push!(exp.args, :(
            Vᴳ_mag = sqrt(Vᴳd^2+Vᴳq^2);
            iΔq_ref = ($(converter.controls[:vac].Kₚ) * ($(converter.controls[:vac].ref[1]) - Vᴳ_mag) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:vac].Kᵢ) *($(converter.controls[:vac].ref[1]) - Vᴳ_mag)
        ))
        index +=1
        # epsilon_vac_index = index + 1
    else 
        push!(exp.args, :(
            iΔq_ref = $Qac))
    end

    # Virtual impedance

    if in(:VI, keys(converter.controls))
        # The implementation is equivalent with the one before, since results are matching
        # Hence the state-space realization, which is different now impacts the results!
        # TODO: 
        if ((converter.controls[:VI].n_f)) >=1  # Voltage filtering with butterworth filter

            Abutt, Bbutt, Cbutt, Dbutt =  butterworthMatrices(converter.controls[:VI].n_f, converter.controls[:VI].ω_f, 2);
            #Voltage filtering of Vᴳd
            push!(exp.args, :(
                voltagesIn = [Vᴳd;Vᴳq];
                statesButt= x[$index + 1 : $index + 2*$(converter.controls[:VI].n_f)]; 
                F[$index + 1 : $index + 2*$(converter.controls[:VI].n_f)] = $Abutt*statesButt + $Bbutt*voltagesIn;
                voltagesOut=$Cbutt*statesButt+$Dbutt*voltagesIn;
                Vᴳd_f=voltagesOut[1];
                Vᴳq_f=voltagesOut[2];
                ))

            init_x = [init_x;zeros(index-length(init_x))];
            #init_x = [init_x; (Vm * cos(θ))/Cbutt[1];zeros(converter.controls[:VI].n_f-1);(-Vm * sin(θ))/Cbutt[1];zeros(converter.controls[:VI].n_f-1)];
            init_x = [init_x; 2*zeros(converter.controls[:VI].n_f)];
            index += 2*(converter.controls[:VI].n_f) 

        else  # No voltage filtering

            push!(exp.args, :(
                Vᴳd_f=Vᴳd;
                Vᴳq_f=Vᴳq;
            ))

        end

        
        if isa(converter.controls[:VI], CCQSEM)
            
            push!(exp.args, :(
                
                iΔd_ref=($(converter.controls[:VI].Rᵥ)*(($(converter.controls[:VI].ref_vd)+Vⱽd_ref)-Vᴳd_f) + ωᵥ*$(converter.controls[:VI].Lᵥ)*(Vᴳq_f-$(converter.controls[:VI].ref_vq)))/(($(converter.controls[:VI].Rᵥ))^2+ωᵥ^2*($(converter.controls[:VI].Lᵥ))^2);
                iΔq_ref=($(converter.controls[:VI].Rᵥ)*($(converter.controls[:VI].ref_vq)-Vᴳq_f) + ωᵥ*$(converter.controls[:VI].Lᵥ)*(-Vᴳd_f+($(converter.controls[:VI].ref_vd)+Vⱽd_ref)))/(($(converter.controls[:VI].Rᵥ))^2+ωᵥ^2*($(converter.controls[:VI].Lᵥ))^2)

                ))
        
        else

        end

    



    end

    # add control equations
    for (key, val) in (converter.controls)
        if (key == :ccc)
            # circulating current control
            push!(exp.args, :(
                        (iΣd_ref, iΣq_ref) = [$(converter.controls[:ccc].ref[1]); $(converter.controls[:ccc].ref[2])];
                        F[$index+1] = $(converter.controls[:ccc].Kᵢ) * (iΣd_ref - iΣd);
                        F[$index+2] = $(converter.controls[:ccc].Kᵢ) * (iΣq_ref - iΣq);
                        # vMΣd_ref = - Ki_Σ * xiΣd - Kp_Σ * (iΣd_ref -  iΣd) + 2*w*Larm*iΣq
                        # vMΣq_ref = - Ki_Σ * xiΣq - Kp_Σ * (iΣq_ref -  iΣq) - 2*w*Larm*iΣd
                        # Assuming constant w
                        vMΣd_ref_c = (- x[$index+1] -
                                $(converter.controls[:ccc].Kₚ) * (iΣd_ref - iΣd) + 2 * $Lₐᵣₘ * iΣq);
                        vMΣq_ref_c = (- x[$index+2] -
                                $(converter.controls[:ccc].Kₚ) * (iΣq_ref - iΣq) - 2 * $Lₐᵣₘ * iΣd);
                        (vMΣd_ref, vMΣq_ref) = I_2θ * [vMΣd_ref_c; vMΣq_ref_c]))
            index += 2
        elseif (key == :zcc) && !in(:energy, keys(converter.controls)) # TODO: Doesnt make sense to have a ZCC with fixed current reference in my opinion ?!
            # zero current control
            push!(exp.args,
                        :(F[$index+1] = $(converter.controls[:zcc].Kᵢ) *($(converter.controls[:zcc].ref[1]) - x[5]);
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (Vdc/2 - $(converter.controls[:zcc].Kₚ) *
                            ($(converter.controls[:zcc].ref[1]) - x[5]) -  x[$index+1])))
            index += 1
        elseif (key == :energy)
            # zero energy control
            push!(exp.args, :(
                        # wΣz = 3*Carm * (vCΔd^2 + vCΔq^2 + vCΔZd^2 + vCΔZq^2 + vCΣd^2 + vCΣq^2 + 2*vCΣz^2)/(2*N)
                        wΣz = (x[6]^2 + x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2 + x[11]^2 + 2x[12]^2)/2;
                        F[$index+1] = $(converter.controls[:energy].Kᵢ) * ($(converter.controls[:energy].ref[1]) - wΣz);
                        # Pac = (3/2)*(Vgd*iΔd + Vgq*iΔq)
                        P_ac = (Vᴳd * iΔd + Vᴳq * iΔq);
                        #iΣz_ref = (Kp_wΣ * (wΣz_ref - wΣz) + Ki_wΣ * xwΣz + Pac) / 3 / Vdc,
                        iΣz_ref = ($(converter.controls[:energy].Kₚ) * ($(converter.controls[:energy].ref[1]) - wΣz) +
                             x[$index+1] + P_ac_f) / 3 / Vdc;
                        F[$index+2] = $(converter.controls[:zcc].Kᵢ) *(iΣz_ref - x[5]);
                        # vMΣz_ref = Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz,
                        vMΣz_ref = (Vdc/2 - $(converter.controls[:zcc].Kₚ) *
                            (iΣz_ref - x[5]) -   x[$index+2])))
            index += 2            
        elseif (key == :occ)
            # output current control


            if ((converter.controls[:occ].n_f)) >=1  # Filtering of the voltage in the voltage feedforward with nth-order butterworth filter with gain 1 

                Abutt_fc, Bbutt_fc, Cbutt_fc, Dbutt_fc =  butterworthMatrices(converter.controls[:occ].n_f, converter.controls[:occ].ω_f, 2);
                #Voltage filtering of Vᴳd
                push!(exp.args, :(
                    voltagesIn = [Vᴳd;Vᴳq];
                    statesButt_fc= x[$index + 1 : $index + 2*$(converter.controls[:occ].n_f)]; 
                    F[$index + 1 : $index + 2*$(converter.controls[:occ].n_f)] = $Abutt_fc*statesButt_fc + $Bbutt_fc*voltagesIn;
                    voltagesOut_fc=$Cbutt_fc*statesButt_fc+$Dbutt_fc*voltagesIn;
                    Vᴳd_fc=voltagesOut_fc[1];
                    Vᴳq_fc=voltagesOut_fc[2];
                    ))
                    init_x = [init_x;zeros(index-length(init_x))];
                    init_x = [init_x; 2*zeros(converter.controls[:occ].n_f)];
                    index += 2*(converter.controls[:occ].n_f) 

            else # No filtering of the voltage in the voltage feedforward

                push!(exp.args, :(
                    Vᴳd_fc=1*Vᴳd;
                    Vᴳq_fc=1*Vᴳq;
                ))


            end
                push!(exp.args, :(
                    # (iΔd_ref, iΔq_ref) = T_θ * [iΔd_ref; iΔq_ref];
                    F[$index+1] = $(converter.controls[:occ].Kᵢ) * (iΔd_ref - iΔd);
                    F[$index+2] = $(converter.controls[:occ].Kᵢ) * (iΔq_ref - iΔq);
                    # vMΔd_ref = Ki_Δ * xiΔd + Kp_Δ * (iΔd_ref -  iΔd) + w*L_eqac*iΔq + Vᴳd
                    # vMΔq_ref = Ki_Δ * xiΔq + Kp_Δ * (iΔq_ref -  iΔq) - w*Leqac*iΔd + Vᴳq
                    
                    # Original case, w provided by the PLL
                    # vMΔd_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+1] +
                    #             $(converter.controls[:occ].Kₚ) * (iΔd_ref - iΔd) + ω * $Lₑ * iΔq + Vᴳd);
                    # vMΔq_ref_c = ($(converter.controls[:occ].Kᵢ) * x[$index+2] +
                    #             $(converter.controls[:occ].Kₚ) * (iΔq_ref - iΔq) - ω * $Lₑ * iΔd + Vᴳq);
                    # Assuming constant w. Having (1+deltaw) instead of 1 results in mismatches above 100 Hz in y_dq and y_qq.
                    vMΔd_ref_c = ( x[$index+1] +
                                $(converter.controls[:occ].Kₚ) * (iΔd_ref - iΔd) + $Lₑ * iΔq + 1*Vᴳd_fc);
                    vMΔq_ref_c = ( x[$index+2] +
                                $(converter.controls[:occ].Kₚ) * (iΔq_ref - iΔq) - $Lₑ * iΔd + 1*Vᴳq_fc);
                    (vMΔd_ref, vMΔq_ref) = I_θ * [vMΔd_ref_c; vMΔq_ref_c]))  # Transformation from converter frame to grid dq frame 
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
    (mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz) = 2/Vdc * [-vMΔd_ref * $baseConv1; -vMΔq_ref * $baseConv1; -vMΔZd_ref * $baseConv1; -vMΔZq_ref * $baseConv1; vMΣd_ref; vMΣq_ref; vMΣz_ref];

    vMΔd = $baseConv2 * ((mΔq*x[11])/4 - (mΔd*x[12])/2 - (mΔd*x[10])/4 - (mΔZd*x[10])/4 + (mΔZq*x[11])/4 - (mΣd*x[6])/4 - (mΣz*x[6])/2 +
            (mΣq*x[7])/4 - (mΣd*x[8])/4 + (mΣq*x[9])/4);
    vMΔq = $baseConv2 * ((mΔd*x[11])/4 + (mΔq*x[10])/4 - (mΔq*x[12])/2 - (mΔZd*x[11])/4 - (mΔZq*x[10])/4 + (mΣd*x[7])/4 + (mΣq*x[6])/4 -
            (mΣz*x[7])/2 - (mΣd*x[9])/4 - (mΣq*x[8])/4);

    # vMΣd = (mΔd*vCΔd)/4 - (mΔq*vCΔq)/4 + (mΔd*vCΔZd)/4 + (mΔZd*vCΔd)/4 + (mΔq*vCΔZq)/4 + (mΔZq*vCΔq)/4 + (mΣd*vCΣz)/2 + (mΣz*vCΣd)/2;
    vMΣd =mΔd*x[6]/4 - mΔq*x[7]/4 + mΔd*x[8]/4 + mΔZd*x[6]/4 + mΔq*x[9]/4 + mΔZq*x[7]/4 + mΣd*x[12]/2 + mΣz*x[10]/2;
    # vMΣq = (mΔq*vCΔZd)/4 - (mΔq*vCΔd)/4 - (mΔd*vCΔZq)/4 - (mΔd*vCΔq)/4 + (mΔZd*vCΔq)/4 - (mΔZq*vCΔd)/4 + (mΣq*vCΣz)/2 + (mΣz*vCΣq)/2;
    vMΣq = mΔq*x[8]/4 - mΔq*x[6]/4 - mΔd*x[9]/4 - mΔd*x[7]/4 + mΔZd*x[7]/4 - mΔZq*x[6]/4 + mΣq*x[12]/2 + mΣz*x[11]/2;
    # vMΣz = (mΔd*vCΔd)/4 + (mΔq*vCΔq)/4 + (mΔZd*vCΔZd)/4 + (mΔZq*vCΔZq)/4 + (mΣd*vCΣd)/4 + (mΣq*vCΣq)/4 + (mΣz*vCΣz)/2;
    vMΣz = mΔd*x[6]/4 + mΔq*x[7]/4 + mΔZd*x[8]/4 + mΔZq*x[9]/4 + mΣd*x[10]/4 + mΣq*x[11]/4 + mΣz*x[12]/2;

    F[1] = -(inputs[2] * $converter.turnsRatio - vMΔd + $Rₑ*x[1] + $Lₑ*x[2])/$Lₑ;                 # diΔd_dt =-(Vgd - vMΔd + Rₑ*iΔd + Lₑ*iΔq*w)/Lₑ, grid frame, pu
    F[2] = -(inputs[3] * $converter.turnsRatio - vMΔq + $Rₑ*x[2] - $Lₑ*x[1])/$Lₑ;                 # diΔq_dt =-(Vgq - vMΔq + Rₑ*iΔq - Lₑ*iΔd*w)/Lₑ, grid frame, pu
    F[3] = -(vMΣd + $Rₐᵣₘ*x[3] - 2*$Lₐᵣₘ*x[4])/$Lₐᵣₘ;                                  # diΣd_dt =-(vMΣd + Rₐᵣₘ*iΣd - 2*Lₐᵣₘ*iΣq*w)/Lₐᵣₘ, grid 2w frame
    F[4] = -(vMΣq + $Rₐᵣₘ*x[4] + 2*$Lₐᵣₘ*x[3])/$Lₐᵣₘ;                                  # diΣq_dt =-(vMΣq + Rₐᵣₘ*iΣq + 2*Lₐᵣₘ*iΣd*w)/Lₐᵣₘ,  grid 2w frame
    F[5] = -(vMΣz - Vdc/2 + $Rₐᵣₘ*x[5])/$Lₐᵣₘ;                                     # diΣz_dt =-(vMΣz - Vᵈᶜ/2 + Rₐᵣₘ*iΣz)/Lₐᵣₘ
    # dvCΔd_dt =(N*(iΣz*mΔd - (iΔq*mΣq)/4 + iΣd*(mΔd/2 + mΔZd/2) - iΣq*(mΔq/2 + mΔZq/2) + iΔd*(mΣd/4 + mΣz/2) - (2*Cₐᵣₘ*vCΔq*w)/N))/(2*Cₐᵣₘ)
    F[6] = ($N*(x[5]*mΔd - x[2]*$baseConv3*mΣq/4 + x[3]*(mΔd/2 + mΔZd/2) - x[4]*(mΔq/2 + mΔZq/2) + x[1]*$baseConv3*(mΣd/4 + mΣz/2) - 2*$Cₐᵣₘ*x[7]/$N))/2/$Cₐᵣₘ;
    # dvCΔq_dt =-(N*((iΔd*mΣq)/4 - iΣz*mΔq + iΣq*(mΔd/2 - mΔZd/2) + iΣd*(mΔq/2 - mΔZq/2) + iΔq*(mΣd/4 - mΣz/2) - (2*Cₐᵣₘ*vCΔd*w)/N))/(2*Cₐᵣₘ)
    F[7] = -($N*((x[1]*$baseConv3*mΣq)/4 - x[5]*mΔq + x[4]*(mΔd/2 - mΔZd/2) + x[3]*(mΔq/2 - mΔZq/2) + x[2]*$baseConv3*(mΣd/4 - mΣz/2) - 2*$Cₐᵣₘ*x[6]/$N))/2/$Cₐᵣₘ;
    # dvCΔZd_dt =(N*(iΔd*mΣd + 2*iΣd*mΔd + iΔq*mΣq + 2*iΣq*mΔq + 4*iΣz*mΔZd))/(8*Cₐᵣₘ) - 3*vCΔZq*w
    F[8] = ($N*(x[1]*$baseConv3*mΣd + 2*x[3]*mΔd + x[2]*$baseConv3*mΣq + 2*x[4]*mΔq + 4*x[5]*mΔZd))/(8*$Cₐᵣₘ) - 3*x[9];
    # dvCΔZq_dt =3*vCΔZd*w + (N*(iΔq*mΣd - iΔd*mΣq + 2*iΣd*mΔq - 2*iΣq*mΔd + 4*iΣz*mΔZq))/(8*Cₐᵣₘ)
    F[9] = 3*x[8] + ($N*(x[2]*$baseConv3*mΣd - x[1]*$baseConv3*mΣq + 2*x[3]*mΔq - 2*x[4]*mΔd + 4*x[5]*mΔZq))/(8*$Cₐᵣₘ);
    # dvCΣd_dt =(N*(iΣd*mΣz + iΣz*mΣd + iΔd*(mΔd/4 + mΔZd/4) - iΔq*(mΔq/4 - mΔZq/4) + (4*Cₐᵣₘ*vCΣq*w)/N))/(2*Cₐᵣₘ)
    F[10] = ($N*(x[3]*mΣz + x[5]*mΣd + x[1]*$baseConv3*(mΔd/4 + mΔZd/4) - x[2]*$baseConv3*(mΔq/4 - mΔZq/4) + 4*$Cₐᵣₘ*x[11]/$N))/(2*$Cₐᵣₘ);
    # dvCΣq_dt =-(N*(iΔq*(mΔd/4 - mΔZd/4) - iΣz*mΣq - iΣq*mΣz + iΔd*(mΔq/4 + mΔZq/4) + (4*Cₐᵣₘ*vCΣd*w)/N))/(2*Cₐᵣₘ)
    F[11] = -($N*(x[2]*$baseConv3*(mΔd/4 - mΔZd/4) - x[5]*mΣq - x[4]*mΣz + x[1]*$baseConv3*(mΔq/4 + mΔZq/4) + 4*$Cₐᵣₘ*x[10]/$N))/(2*$Cₐᵣₘ);
    # dvCΣz_dt =(N*(iΔd*mΔd + iΔq*mΔq + 2*iΣd*mΣd + 2*iΣq*mΣq + 4*iΣz*mΣz))/(8*Cₐᵣₘ)
    F[12] = ($N*(x[1]*$baseConv3*mΔd + x[2]*$baseConv3*mΔq + 2*x[3]*mΣd + 2*x[4]*mΣq + 4*x[5]*mΣz))/(8*$Cₐᵣₘ);
    F[1:12] *= $wbase))

    function f!(expr, F, x, inputs) # Creating callable function that evaluates the function "F" with respect to "x" and "inputs". F(x,inputs)= x^dot
       f = eval(:((F,x,inputs) -> $expr))
       return Base.invokelatest(f, F,x,inputs)
    end

    # if in(:dc, keys(converter.controls))
    #     vector_inputs = [Pdc/Vdc, Vᴳd, Vᴳq]
    #     init_x = [init_x; zeros(index-12,1)]
    #     init_x[vdc_position] = Vdc
    # else
    #     vector_inputs = [Vdc, Vᴳd, Vᴳq]
    #     init_x = [init_x; zeros(index-12,1)]
    # end

    vector_inputs = [Vdc, Vᴳd, Vᴳq] # Inputs to the system for the steady state solution and linearization
    init_x=[init_x;zeros(index-length(init_x))]; ################################ATTENTION!!!!!#####################

    # If there is a dc voltage controller, add an additional equation to represent the dc voltage, only for the steady-state solution
    exp_steadyState = copy(exp) # copy the state-space formulation f 
    if in(:dc, keys(converter.controls))  
        init_x =[init_x;Vdc]
        push!(exp_steadyState.args,
        :(
            F[$index+1] = $wbase * ($Idc_in - 3*x[5]) / $Cₑ;
            F[$epsilon_vdc_index] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - x[end]);
        ))
    end
    # if in(:vac, keys(converter.controls))
    #     init_x =[init_x;Vm;0]
    #     push!(exp_steadyState.args,
    #     :(
    #         F[$index+1] = $wbase * ($Id - x[1]) / 1e-9;
    #         F[$index+2] = $wbase * ($Iq - x[2]) / 1e-9;
    #         Vᴳ_mag_SS = sqrt(x[end-1]^2 + x[end]^2);
    #         F[$epsilon_vac_index] = $(converter.controls[:vac].Kᵢ) *($(converter.controls[:vac].ref[1]) - Vᴳ_mag_SS);
    #     ))
    # end

    g!(F,x) = f!(exp_steadyState, F, x, vector_inputs) # g is the state-space formulation used to obtain the steady-state operation point, copy from f, see some lines above
    # TODO: Check if it makes sense to use newton or trust_region
    # Newton seems to give really bad estimates for control-related states, while trust_region cannot even get id and iq correct.
    # k = nlsolve(g!, init_x, autodiff = :forward, iterations = 100, method = :newton)
    k = nlsolve(g!, init_x, autodiff = :forward, iterations = 100, ftol = 1e-6, xtol = 1e-3, method = :trust_region) # solve for steady state with initial operating point "init_x" 
    if converged(k)
        println("MMC steady-state solution found!")
    end
    # if in(:dc, keys(converter.controls))
    #     converter.equilibrium = converter.equilibrium[1:end-1]
    # end
    # if in(:vac, keys(converter.controls))
    #     converter.equilibrium = converter.equilibrium[1:end-2]
    # end
    if in(:dc, keys(converter.controls))
        converter.equilibrium = k.zero[1:end-1]
    else
        converter.equilibrium = k.zero
    end
    

    h(F,x) = f!(exp, F, x[1:end-3], x[end-2:end]) # See it as a function call
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index+3); h(F, x); return F)
    A = zeros(index+3,index+3)
    ForwardDiff.jacobian!(A, ha, [converter.equilibrium' vector_inputs'])
    converter.A = real(A[1:end-3, 1:end-3])
    converter.B = real(A[1:end-3, end-2:end])

    converter.C = real(zeros(3, size(converter.A,1)))
    converter.C[2,1] = 1  #iΔd in grid frame 
    converter.C[3,2] = 1  #iΔq in grid frame 
    # !in(:dc, keys(converter.controls)) ? converter.C[1,5] = 3 : converter.C[1, vdc_position] = 1
    converter.C[1,5] = 3  # iDC voltage ??
    converter.D = real(zeros(3,3))

    # converter.B[:,2:3] *= converter.turnsRatio
    # converter.D[:,2:3] *= converter.turnsRatio
end

function eval_parameters(converter :: MMC, s :: Complex)
    # numerical
    I = Matrix{Complex}(Diagonal([1 for dummy in 1:size(converter.A,1)]))
    # Y = (converter.C*inv(s*I-converter.A))*converter.D + converter.D # This matrix is in pu
    Y = converter.C * ((s*I-converter.A) \ converter.B) + converter.D # This matrix is in pu

    # if in(:dc, keys(converter.controls))
    #     (m11, m12, m21, m22) = (Y[1,1], Y[1,2:3], Y[2:3,1], Y[2:3,2:3])
    #     m11 = 1/m11
    #     Y = [m11 -transpose(m12)*m11; m21*m11 m22-m21*m11*transpose(m12)]
    # end

    Y[1,:] *= converter.iDCbase
    Y[:,1] /= converter.vDCbase
    Y[2:3,:] *= converter.iACbase # Base current of the converter side 
    # # # Multiplication with the AC voltage base converts the pu admittance to SI.
    # # # The double division with the turns ratio is actually a multiplication,
    # # # and is needed to bring the grid-side voltage to the converter side.
    Y[:,2:3] /= (converter.vACbase / converter.turnsRatio) # Base voltage at the grid side 
    

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
    # d_not_sorted=eigvals(Ad);
    # v_not_sorted=eigvecs(Ad);
    # d=zeros(ComplexF64,padeOrderDen,1);
    # v=zeros(ComplexF64,padeOrderDen,padeOrderDen);
    # matrix_ind=1;
    # # Sort the eigenvalues and eigenvectors
    # for i=1:padeOrderDen
    #     if sign(imag(d_not_sorted[i])==0) #Real eigenvalue, copy as is
    #         d[matrix_ind]=d_not_sorted[i];
    #         v[:,matrix_ind]=v_not_sorted[:,i];
    #         matrix_ind+=1;
    #     else # Complex eigenvalue
    #         if matrix_ind < padeOrderDen
    #             if sign(imag(d_not_sorted[i]))==-1
    #                 d[matrix_ind]    = real(d_not_sorted[i])     -imag(d_not_sorted[i])*im;
    #                 d[matrix_ind+1]  = real(d_not_sorted[i])     +imag(d_not_sorted[i])*im;
    #                 v[:,matrix_ind]  = real(v_not_sorted[:,i])   - imag(v_not_sorted[:,i])*im;
    #                 v[:,matrix_ind+1]= real(v_not_sorted[:,i])   + imag(v_not_sorted[:,i])*im;
    #             else
    #                 d[matrix_ind]    = real(d_not_sorted[i])     +imag(d_not_sorted[i])*im;
    #                 d[matrix_ind+1]  = real(d_not_sorted[i])     -imag(d_not_sorted[i])*im;
    #                 v[:,matrix_ind]  = real(v_not_sorted[:,i])   + imag(v_not_sorted[:,i])*im;
    #                 v[:,matrix_ind+1]= real(v_not_sorted[:,i])   - imag(v_not_sorted[:,i])*im;
    #             end
    #             matrix_ind+=2;
    #         end
    #     end
    # end
    # matrix_ind=1;
    # T_inv=zeros(padeOrderDen,padeOrderDen);
    # for i=1:padeOrderDen
    #     if imag(d[i])==0 #Real eigenvalue
    #         T_inv[:,matrix_ind]=v[:,i];
    #         matrix_ind+=1;
    #     else # Complex eigenvalue
    #         if matrix_ind < padeOrderDen
    #             T_inv[:,matrix_ind]  =real(v[:,i]);
    #             T_inv[:,matrix_ind+1]=imag(v[:,i]);
    #             matrix_ind+=2;
    #         end
    #     end
    # end
    #Original implementation, resulting in a SingularException for Pade orders larger than 3.
    # T=inv(T_inv);
    # Ad=T*Ad*T_inv;
    # Bd=T*Bd;
    # Cd=Cd*T_inv;
    # Alternative
    sys = ss(Ad,Bd,Cd,Dd)
    sys_modal = modal_form(sys; C1=true)
    Ad = sys_modal[1].A
    Bd = sys_modal[1].B
    Cd = sys_modal[1].C
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

function butterworthMatrices(buttOrder,ω_c,numberVars)

    # Calculation of the state-space representation of a n-order butterworth filter with a gain of 1.
    # buttOrder = Order of butterworth filter, ω_c= Cutoff frequency of the filter in [rad/s]
    # 
    # TODO: Reference to equations
   
    size_A=buttOrder;
    Ab=zeros(size_A,size_A);
    Bb=zeros(size_A,1);
    Bb[end]=1;
    Cb=zeros(1,size_A);
    Cb[1]=1;
    Db=0;
    Ab[1:end-1,2:end] = Matrix(1.0I, buttOrder-1, buttOrder-1);
    
    γ=pi/(2*buttOrder)
    
    # Calculation of the matrix entries in A 
    # Calculation of the coefficients of the denominator polynominal aₙ*sⁿ+...+a₀
    for i=0:buttOrder-1
    
        
        if i==0
    
            a_i = 1 
            Ab[end, i+1] = -a_i;
        
        else 
    
            a_i = 1; 
            for μ=1:i
    
                a_i=a_i*cos((μ-1)γ)/(sin(μ*γ));
            
            end
            Ab[end, i+1] = -a_i * (1/ω_c)^(i)
    
        end
    
    
    end
    
    # Convert from aₙ*sⁿ+...+a₀ to sⁿ+...+a₀ by dividing numerator and denominator by 1/aₙ
    Ab[end, 1:end]=Ab[end, 1:end]*(ω_c)^buttOrder;
    Cb[1]=Cb[1]*(ω_c)^buttOrder;
    

    # Controllable canonical form can create problems while solving for MMC steady-state, espec. when high order pade approximation is used. 
    # Probably related to high conditioning number of controllable canonical form.
    # Transformation to modal form, which results in lower conditioning number.
    sys = ss(Ab,Bb,Cb,Db)
    sys_modal = modal_form(sys;C1 = true)
    Ab= sys_modal[1].A
    Bb = sys_modal[1].B
    Cb = sys_modal[1].C
    Db = sys_modal[1].D


    # Adjust matrices for multiple,independent inputs, so far only up to 2 possible
    if numberVars == 1 #One input, one output
        A_butt=Ab;
        B_butt=Bb;
        C_butt=Cb;
        D_butt=Db;
    elseif numberVars == 2 #Two inputs, two outputs 
        A_butt=cat(Ab,Ab;dims=[1,2]);
        B_butt=cat(Bb,Bb;dims=[1,2]);
        C_butt=cat(Cb,Cb;dims=[1,2]);
        D_butt=cat(Db,Db;dims=[1,2]);
    end
    
    return A_butt,B_butt,C_butt,D_butt
   
end 




