export mmc

@with_kw mutable struct MMC <: Converter
    ω₀ :: Union{Int, Float64} = 100*π           # Base angular frequncy [rad/s]

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
    N :: Int = 400                              # number of submodules per arm [-]

    Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the converter transformer at the converter side [H]
    Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the converter transformer at the converter side [Ω]


    gfm :: Bool = false                         # Grid-following or grid-froming control. True= GFM, False=GFL
    controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
    equilibrium :: Array{Union{Int, Float64}} = [0]
    
    A :: Array{Float64} = [0]
    B :: Array{Float64} = [0]
    C :: Array{Float64} = [0]
    D :: Array{Float64} = [0]

    timeDelay :: Float64 = 0                    # Time delay [s]
    padeOrderNum :: Int = 0                     # Order of the numerator polynominal of pade approximation [-]
    padeOrderDen :: Int = 0                     # Order of the denominator polynominal of pade approximation [-]

    vACbase_LL_RMS :: Union{Int, Float64} = 380 # Voltage base [kV]
    Sbase :: Union{Int, Float64} = 1000         # Power base [MW]
    vDCbase :: Union{Int, Float64} = 640        # DC voltage base [kV]

    turnsRatio :: Union{Int, Float64} = 1       # Turns ratio of the converter transformer, converter side/AC side [-]

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
Rₐᵣₘ :: Union{Int, Float64}  = 1.07         # equivalent arm resistance [Ω]
Cₐᵣₘ :: Union{Int, Float64}  = 10e-3        # capacitance per submodule [F]
N :: Int = 401                              # number of submodules per arm [-]

Lᵣ :: Union{Int, Float64}  = 60e-3          # inductance of the phase reactor [H]
Rᵣ :: Union{Int, Float64}  = 0.535          # resistance of the phase reactor [Ω]

controls :: OrderedDict{Symbol, Controller} = OrderedDict{Symbol, Controller}()
equilibrium :: Array{Union{Int, Float64}} = [0]
A :: Array{Complex} = [0]
B :: Array{Complex} = [0]
C :: Array{Complex} = [0]
D :: Array{Complex} = [0]
........................

```

The constructed MMC has 2 pins on the AC side: `2.1`, `2.2`, and 1 pin on its
DC-side: `1.1`.
"""
function mmc(;args...) #Constructor 
    converter = MMC()

    for (key, val) in pairs(args)
        if isa(val, Controller)
            converter.controls[key] = val
        elseif in(key, propertynames(converter))
            setfield!(converter, key, val)
        end
    end
    elem = Element(input_pins = 1, output_pins = 2, element_value = converter)
end

function update!(converter :: MMC, Vm, θ, Pac, Qac, Vdc, Pdc) #Function to calculate state space and impedance of MMC with respect to power flow solution
    

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

    Qac *=-1 # Correction for reactive power sign

    converter.Vₘ = Vm
    converter.θ = θ
    converter.Vᵈᶜ = Vdc
    converter.P = Pac
    converter.Q = Qac   
    converter.P_dc = Pdc # Has the same sign as Pac

    Vm /= vAC_base       # Grid side voltage (peak,phase) perunitized by converter-side base voltage (peak,phase)
    Vdc /= vDC_base
    Pac /= Sbase
    Qac /= Sbase
    Pdc /= Sbase
    
    Vᴳd = Vm * cos(θ)   
    Vᴳq = -Vm * sin(θ)  
 
    #TODO: Equations correct, see Q compare to TLC
    Id = ((Vᴳd*converter.turnsRatio * Pac - Vᴳq*converter.turnsRatio * Qac) / ((Vᴳd*converter.turnsRatio)^2 + (Vᴳq*converter.turnsRatio)^2)) 
    Iq = ((Vᴳq*converter.turnsRatio * Pac + Vᴳd*converter.turnsRatio * Qac) / ((Vᴳd*converter.turnsRatio)^2 + (Vᴳq*converter.turnsRatio)^2)) 

###############################Manipulating controller references and defining steady-state values for the MMC model###################
    init_x = zeros(12, 1) # Initial steady-state values for MMC model

    
    for (key, val) in (converter.controls)
                
        if (key == :occ)
            if (length(val.ref) == 1) && (val.ref[1] == 0) # no reference existent 
                val.ref = [Id Iq]
            end
            init_x[1] = val.ref[1]
            init_x[2] = val.ref[2]
        elseif (key == :energy)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref[1] = 1 
            end
        elseif (key == :ccc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [0 0]
            end
            init_x[3] = val.ref[1]
            init_x[4] = val.ref[2]
        elseif (key == :p)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Pac*Sbase]   # Comes from Powerflow definition P, conversion to SI. Will be converted back, later in the code 
            end
        elseif (key == :q)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Qac*Sbase]   # Comes from Powerflow definition Q, conversion to SI. Will be converted back, later in the code 
            end
        elseif (key == :dc)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                val.ref = [Vdc]         # Comes from Powerflow definition Vdc, in pu.
            end
        elseif (key == :vac) || (key == :vac_supp)
            if (length(val.ref) == 1) && (val.ref[1] == 0)
                
                val.ref = [sqrt(2)*converter.Vₘ]


            end
        end
    end
    Idc_in = Pdc/Vdc
    init_x[5] = Pdc/3/Vdc
    init_x[12] = Vdc

    index = 12
    exp = Expr(:block) # Start construction of the state-space equations




#################################### Define transformation matrices######################################################################
    # Determine position of state variable of VSE Δθᵥ
    # Determine position of state variable of PLL
    index_theta_VSM = 14
    index_pll=14

    if in(:pll, keys(converter.controls))
           
        index_pll+=converter.controls[:pll].n_f

    else

        
    end

    if converter.gfm

        index_theta_VSM=index_pll+2 #PLL states
        index_theta_VSM+=converter.controls[:p].n_f #P filter states
    else


    end
    
    # Define transformation angle and matrix with or without PLL
    if in(:pll, keys(converter.controls)) # PLL

        push!(exp.args, :(
            Δθ_pll=x[$index_pll];
            T_θ_pll = [cos(Δθ_pll) -sin(Δθ_pll); sin(Δθ_pll) cos(Δθ_pll)];
        ))

    # No PLL implementation
    else
    
        push!(exp.args, :(
        
            Δθ_pll=θ;
        
        ))

    end 
   
    
    if converter.gfm #Grid-forming
        push!(exp.args, :(
            Δθᵥ=x[$index_theta_VSM];
            T_θ = [cos(Δθᵥ) -sin(Δθᵥ); sin(Δθᵥ) cos(Δθᵥ)];
            I_θ = [cos(Δθᵥ) sin(Δθᵥ); -sin(Δθᵥ) cos(Δθᵥ)];
            T_2θ = [cos(-2Δθᵥ) -sin(-2Δθᵥ); sin(-2Δθᵥ) cos(-2Δθᵥ)];
            I_2θ = [cos(-2Δθᵥ) sin(-2Δθᵥ); -sin(-2Δθᵥ) cos(-2Δθᵥ)];
        ))

    else #Grid-following
        push!(exp.args, :(
            T_θ = [cos(Δθ_pll) -sin(Δθ_pll); sin(Δθ_pll) cos(Δθ_pll)];
            I_θ = [cos(Δθ_pll) sin(Δθ_pll); -sin(Δθ_pll) cos(Δθ_pll)];
            T_2θ = [cos(-2Δθ_pll) -sin(-2Δθ_pll); sin(-2Δθ_pll) cos(-2Δθ_pll)];
            I_2θ = [cos(-2Δθ_pll) sin(-2Δθ_pll); -sin(-2Δθ_pll) cos(-2Δθ_pll)];
            ))

    end

    # Transformation of required variables for control purposes
    push!(exp.args, :(
        (iΔd, iΔq) = T_θ * [x[1]; x[2]]; # Currents in grid dq frame defined: x1 and x2, see circuit equations far below 
        (iΣd, iΣq) = T_2θ * [x[3]; x[4]];
        (Vᴳd, Vᴳq) = T_θ * [inputs[2] * $converter.turnsRatio; inputs[3] * $converter.turnsRatio]; #Vd_grid input 2 and Vq_grid input 3 both expressed in the grid frame and at grid side
        Vdc = inputs[1];)) #Vdc voltage input 1  
    

 ##################################################PLL###################################################################################
    if in(:pll, keys(converter.controls))
        
        push!(exp.args, :(
            (Vᴳd_pll, Vᴳq_pll) = T_θ_pll * [inputs[2] * $converter.turnsRatio; inputs[3] * $converter.turnsRatio];  # Vgd: Input 1 and Vgq: Input 2
            ))

        if (((converter.controls[:pll].n_f)) >= 1 ) #Filtering of Vᴳq_pll

            Abutt_pll, Bbutt_pll, Cbutt_pll, Dbutt_pll =  butterworthMatrices(converter.controls[:pll].n_f, converter.controls[:pll].ω_f, 1);
            push!(exp.args, :(
                
                statesButt_pll= x[$index + 1 : $index + 1*$(converter.controls[:pll].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:pll].n_f)] = $Abutt_pll*statesButt_pll + $Bbutt_pll*-Vᴳq_pll;
                Vᴳq_pll_f=dot($Cbutt_pll,statesButt_pll)+$Dbutt_pll*-Vᴳq_pll;# Get rid of 1-element array
            ))

            # init_x = [init_x;zeros(index-length(init_x))];
            # init_x = [init_x; 1*zeros(converter.controls[:pll].n_f)];
            index += 1*(converter.controls[:pll].n_f)

        else #No filtering of Vᴳq_pll
            
            push!(exp.args, :(
                
                Vᴳq_pll_f=-1*Vᴳq_pll;
                
                ))

        end

        #PLL equations
        push!(exp.args, :(
            
            F[$index+1] = Vᴳq_pll_f*$(converter.controls[:pll].Kᵢ);
            Δω = $(converter.controls[:pll].Kₚ) * (Vᴳq_pll_f) + x[$index+1]; #Delta omega_pll [pu]
            F[$index+2] = $wbase*Δω;

            ))
            # init_x = [init_x;zeros(index-length(init_x))];
            # init_x = [init_x;0;θ];
            index += 2;

    end

############################################ P control###################################################################################
    if in(:p, keys(converter.controls)) 
    
        converter.controls[:p].ref[1] /= Sbase # conversion back to pu again
        push!(exp.args, :(
                P_ac = (Vᴳd * iΔd + Vᴳq * iΔq);))

        if ((converter.controls[:p].n_f)) >= 1 # Filtering of p_ac

            Abutt_p, Bbutt_p, Cbutt_p, Dbutt_p =  butterworthMatrices(converter.controls[:p].n_f, converter.controls[:p].ω_f, 1);
            push!(exp.args, :(
                
                statesButt_p= x[$index + 1 : $index + 1*$(converter.controls[:p].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:p].n_f)] = $Abutt_p*statesButt_p + $Bbutt_p*P_ac;
                P_ac_f=dot($Cbutt_p,statesButt_p)+$Dbutt_p*P_ac;
        
            ))
            # init_x = [init_x;zeros(index-length(init_x))];
            # init_x = [init_x; 1*zeros(converter.controls[:p].n_f)];
            index += 1*(converter.controls[:p].n_f)

        else # No filtering of P_ac
            push!(exp.args, :(

            P_ac_f=P_ac;

            ))


        end

       # P controller equations
       if isa(converter.controls[:p], VSE) && (converter.gfm) #VSE with grid-forming converter


            if in(:pll, keys(converter.controls))
            
                # Swing equation with PLL
                push!(exp.args, :(
                    ωᵥ= x[$index+1]; # Actually delta Omega_VSM, relative angle [pu]. With initial conditions it becomes absolute: ωᵥ = Δωᵥ + ωᵥ(0)
                    F[$index+1] =($(converter.controls[:p].ref[1]) - P_ac_f - $(converter.controls[:p].K_d)*(ωᵥ-(Δω+1)) -  $(converter.controls[:p].K_ω)*(ωᵥ - $(converter.controls[:p].ref_ω)))/(2*$(converter.controls[:p].H)) ;
                    F[$index+2] =$wbase*(x[$index+1]-1); 
                    ))

            else
                # Swing equation without PLL no explicit damping term
                push!(exp.args, :(
                    ωᵥ= x[$index+1]; # Actually delta Omega_VSM, relative angle [pu]. With initial conditions it becomes absolute: ωᵥ = Δωᵥ + ωᵥ(0)
                    F[$index+1] =($(converter.controls[:p].ref[1]) - P_ac_f -  $(converter.controls[:p].K_ω)*(ωᵥ - $(converter.controls[:p].ref_ω)))/(2*$(converter.controls[:p].H)) ;
                    F[$index+2] =$wbase*(x[$index+1]-1); 
                    ))


            end
            init_x = [init_x;zeros(index-length(init_x))];
            init_x = [init_x;1;θ]; # Approximating Delta Theta_VSM with θ
            index += 2

       #elseif isa(converter.controls[:p], VSE) && !(converter.gfm) #VSE with grid-following converter



       elseif isa(converter.controls[:p], PI_control) && !(converter.gfm) #Typical PI control with grid-following converter
        
            push!(exp.args, :(

                iΔd_ref = ($(converter.controls[:p].Kₚ) * ($(converter.controls[:p].ref[1]) - P_ac_f) +
                            x[$index+1]);
                F[$index+1] = $(converter.controls[:p].Kᵢ) *($(converter.controls[:p].ref[1]) - P_ac_f);
                
                ))
            
            index += 1

       end

     
#########################################DC voltage control#############################################################################

    elseif in(:dc, keys(converter.controls)) # DC voltage control

        if ((converter.controls[:dc].n_f)) >= 1 # Filtering of Vdc

            Abutt_vdc, Bbutt_vdc, Cbutt_vdc, Dbutt_vdc =  butterworthMatrices(converter.controls[:dc].n_f, converter.controls[:dc].ω_f, 1);
            push!(exp.args, :(
                
                statesButt_vdc= x[$index + 1 : $index + 1*$(converter.controls[:vdc].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:vdc].n_f)] = $Abutt_vdc*statesButt_vdc + $Bbutt_vdc*Vdc;
                Vdc_f=dot($Cbutt_vdc,statesButt_vdc)+$Dbutt_vdc*Vdc;
        
            ))
            # init_x = [init_x;zeros(index-length(init_x))];
            # init_x = [init_x; 1*zeros(converter.controls[:dc].n_f)];
            index += 1*(converter.controls[:dc].n_f)


        else # No filtering of Vdc
            push!(exp.args, :(

            Vdc_f= Vdc;

            ))
            

        end


        # DC voltage controller equations
        push!(exp.args, :(
                F[$index+1] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - Vdc_f);
                    iΔd_ref = -($(converter.controls[:dc].Kₚ) * ($(converter.controls[:dc].ref[1]) - Vdc_f) +
                                 x[$index+1]);))
        epsilon_vdc_index = index + 1
        index += 1    
    
    else # No explicit control of P for grid-following
        push!(exp.args, :(
            iΔd_ref = $(converter.controls[:p].ref[1])/Vᴳd;))
    end

################################################Q control###############################################################################
    if in(:q, keys(converter.controls))
        converter.controls[:q].ref[1] /= Sbase # The minus sign corrects for the Q convention used in the model.
        
        
        push!(exp.args, :(
            Q_ac =  (-Vᴳq * iΔd + Vᴳd * iΔq);))

        if ((converter.controls[:q].n_f)) >= 1 # Filtering of q_ac

            Abutt_q, Bbutt_q, Cbutt_q, Dbutt_q =  butterworthMatrices(converter.controls[:q].n_f, converter.controls[:q].ω_f, 1);
            push!(exp.args, :(
                
                statesButt_q= x[$index + 1 : $index + 1*$(converter.controls[:q].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:q].n_f)] = $Abutt_q*statesButt_q + $Bbutt_q*Q_ac;
                Q_ac_f=dot($Cbutt_q,statesButt_q)+$Dbutt_q*Q_ac;
        
            ))
            # init_x = [init_x;zeros(index-length(init_x))];
            # init_x = [init_x; 1*zeros(converter.controls[:q].n_f)];
            index += 1*(converter.controls[:q].n_f)

        else # No filtering of Q_ac
            push!(exp.args, :(

            Q_ac_f=Q_ac;

            ))


        end


        if in(:vac_supp, keys(converter.controls)) # Voltage droop control
            converter.controls[:vac_supp].ref[1] /= (vAC_base / converter.turnsRatio)
            
            push!(exp.args, :(
                Vᴳ_mag = sqrt(Vᴳd^2+Vᴳq^2);))
                
            if ((converter.controls[:vac_supp].n_f)) >= 1 # Filtering of Vac

                Abutt_vac, Bbutt_vac, Cbutt_vac, Dbutt_vac =  butterworthMatrices(converter.controls[:vac_supp].n_f, converter.controls[:vac_supp].ω_f, 1);
                push!(exp.args, :(
                    
                    statesButt_vac= x[$index + 1 : $index + 1*$(converter.controls[:vac_supp].n_f)]; 
                    F[$index + 1 : $index + 1*$(converter.controls[:vac_supp].n_f)] = $Abutt_vac*statesButt_vac + $Bbutt_vac*Vᴳ_mag;
                    Vᴳ_mag_f=dot($Cbutt_vac,statesButt_vac)+$Dbutt_vac*Vᴳ_mag;
            
                ))
                index += 1*(converter.controls[:vac_supp].n_f)

            else# No Filtering of Vac
                push!(exp.args, :(

                Vᴳ_mag_f=Vᴳ_mag;
    
                ))

            end
                push!(exp.args, :(
                
                q_ref = $(converter.controls[:vac_supp].Kₚ)*($(converter.controls[:vac_supp].ref[1])-Vᴳ_mag_f) +$(converter.controls[:q].ref[1]);
                
                ))
           

        else
            push!(exp.args, :(
                q_ref = $(converter.controls[:q].ref[1])))
        end

        
        if (converter.gfm) # Q control for GFM

            push!(exp.args, :(
                F[$index+1]=$(converter.controls[:q].Kᵢ)*(q_ref - Q_ac_f);
                Vⱽd_ref=($(converter.controls[:q].Kₚ)*(q_ref - Q_ac_f) + x[$index+1]);
            ))
            index += 1
            
            
        else 
            push!(exp.args, :(

            iΔq_ref = ($(converter.controls[:q].Kₚ) * (q_ref - Q_ac_f) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:q].Kᵢ) *(q_ref - Q_ac_f);))
            index += 1

        end
        

#TODO: Implement filter for voltage filtering here, discuss with Ozgur...🤐. I think the minus sign is missing here!
    elseif in(:vac, keys(converter.controls))
        
        converter.controls[:vac].ref[1] /= (vAC_base / converter.turnsRatio)
        
        push!(exp.args, :(
            Vᴳ_mag = sqrt(Vᴳd^2+Vᴳq^2);))

        if ((converter.controls[:vac].n_f)) >= 1 # Filtering of Vac

            Abutt_vac, Bbutt_vac, Cbutt_vac, Dbutt_vac =  butterworthMatrices(converter.controls[:vac].n_f, converter.controls[:vac].ω_f, 1);
            push!(exp.args, :(
                
                statesButt_vac= x[$index + 1 : $index + 1*$(converter.controls[:vac].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:vac].n_f)] = $Abutt_vac*statesButt_vac + $Bbutt_vac*Vᴳ_mag;
                Vᴳ_mag_f=dot($Cbutt_vac,statesButt_vac)+$Dbutt_vac*Vᴳ_mag;
        
            ))
            index += 1*(converter.controls[:vac].n_f)

        else# No Filtering of Vac
            push!(exp.args, :(

            Vᴳ_mag_f=Vᴳ_mag;

            ))

        end

        push!(exp.args, :(
            iΔq_ref = ($(converter.controls[:vac].Kₚ) * ($(converter.controls[:vac].ref[1]) -  Vᴳ_mag_f) +
                         x[$index+1]);
            F[$index+1] = $(converter.controls[:vac].Kᵢ) *($(converter.controls[:vac].ref[1]) -  Vᴳ_mag_f)
        ))

        index +=1
        # epsilon_vac_index = index + 1

    else # No explicit control of Q for grid-following

        push!(exp.args, :(
            iΔq_ref = $(converter.controls[:q].ref[1])/Vᴳd;))
    end
###############################################################Virtual impedance########################################################
    
    if in(:VI, keys(converter.controls))

        if ((converter.controls[:VI].n_f)) >=1  # Voltage filtering

            Abutt, Bbutt, Cbutt, Dbutt =  butterworthMatrices(converter.controls[:VI].n_f, converter.controls[:VI].ω_f, 2);
            push!(exp.args, :(
                voltagesIn = [Vᴳd;Vᴳq];
                statesButt= x[$index + 1 : $index + 2*$(converter.controls[:VI].n_f)]; 
                F[$index + 1 : $index + 2*$(converter.controls[:VI].n_f)] = $Abutt*statesButt + $Bbutt*voltagesIn;
                voltagesOut=$Cbutt*statesButt+$Dbutt*voltagesIn;
                Vᴳd_f=voltagesOut[1];
                Vᴳq_f=voltagesOut[2];
                ))
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
###############################################################Circulating current controller###########################################
    
    for (key, val) in (converter.controls)
        if (key == :ccc)
            # circulating current control
            push!(exp.args, :(
                        (iΣd_ref, iΣq_ref) = [$(converter.controls[:ccc].ref[1]); $(converter.controls[:ccc].ref[2])];
                        F[$index+1] = $(converter.controls[:ccc].Kᵢ) * (iΣd_ref - iΣd);
                        F[$index+2] = $(converter.controls[:ccc].Kᵢ) * (iΣq_ref - iΣq);
                        # vMΣd_ref = 2/Vdc*(- Ki_Σ * xiΣd - Kp_Σ * (iΣd_ref -  iΣd) + 2*Larm*iΣq)
                        # vMΣq_ref = 2/Vdc*(- Ki_Σ * xiΣq - Kp_Σ * (iΣq_ref -  iΣq) - 2*Larm*iΣd)
                        # Assuming constant w
                        vMΣd_ref_c =2/Vdc* (- x[$index+1] -
                                $(converter.controls[:ccc].Kₚ) * (iΣd_ref - iΣd) + 2 * $Lₐᵣₘ * iΣq);
                        vMΣq_ref_c = 2/Vdc*(- x[$index+2] -
                                $(converter.controls[:ccc].Kₚ) * (iΣq_ref - iΣq) - 2 * $Lₐᵣₘ * iΣd);
                        (vMΣd_ref, vMΣq_ref) = I_2θ * [vMΣd_ref_c; vMΣq_ref_c]))
            index += 2

 ######################################################Total energy controller##########################################################
        elseif (key == :energy)
            # zero energy control
            push!(exp.args, :(
                        # wΣz = (vCΔd^2 + vCΔq^2 + vCΔZd^2 + vCΔZq^2 + vCΣd^2 + vCΣq^2 + 2*vCΣz^2)/(2)
                        wΣz = (x[6]^2 + x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2 + x[11]^2 + 2x[12]^2)/2;
                        F[$index+1] = $(converter.controls[:energy].Kᵢ) * ($(converter.controls[:energy].ref[1]) - wΣz);
                        #iΣz_ref = (Kp_wΣ * (wΣz_ref - wΣz) + Ki_wΣ * xwΣz + Pac_f) / 3 / Vdc,
                        iΣz_ref = ($(converter.controls[:energy].Kₚ) * ($(converter.controls[:energy].ref[1]) - wΣz) +
                             x[$index+1] + P_ac_f) / 3 / Vdc;
                        F[$index+2] = $(converter.controls[:zcc].Kᵢ) *(iΣz_ref - x[5]);
                        # vMΣz_ref = 2/Vdc*(Vdc/2 - Kp_Σz*(iΣz_ref - iΣz) - Ki_Σz * xiΣz),
                        vMΣz_ref = 2/Vdc*(Vdc/2 - $(converter.controls[:zcc].Kₚ) *
                            (iΣz_ref - x[5]) -   x[$index+2])))
            index += 2            
 
 ###################################################################Output current controller###########################################
        elseif (key == :occ)
            # output current control
            if ((converter.controls[:occ].n_f)) >=1  # Filtering of the voltage in the voltage feedforward with nth-order butterworth filter with gain 1 

                Abutt_fc, Bbutt_fc, Cbutt_fc, Dbutt_fc =  butterworthMatrices(converter.controls[:occ].n_f, converter.controls[:occ].ω_f, 2);
                push!(exp.args, :(
                    voltagesIn = [Vᴳd;Vᴳq];
                    statesButt_fc= x[$index + 1 : $index + 2*$(converter.controls[:occ].n_f)]; 
                    F[$index + 1 : $index + 2*$(converter.controls[:occ].n_f)] = $Abutt_fc*statesButt_fc + $Bbutt_fc*voltagesIn;
                    voltagesOut_fc=$Cbutt_fc*statesButt_fc+$Dbutt_fc*voltagesIn;
                    Vᴳd_fc=voltagesOut_fc[1];
                    Vᴳq_fc=voltagesOut_fc[2];
                    ))
                    index += 2*(converter.controls[:occ].n_f) 

            else # No filtering of the voltage in the voltage feedforward

                push!(exp.args, :(
                    Vᴳd_fc=1*Vᴳd;
                    Vᴳq_fc=1*Vᴳq;
                ))


            end
            push!(exp.args, :(
                F[$index+1] = $(converter.controls[:occ].Kᵢ) * (iΔd_ref - iΔd);
                F[$index+2] = $(converter.controls[:occ].Kᵢ) * (iΔq_ref - iΔq);
                # vMΔd_ref = 2/Vdc*(Ki_Δ * xiΔd + Kp_Δ * (iΔd_ref -  iΔd) + Leqac*iΔq + Vᴳd)
                # vMΔq_ref = 2/Vdc*(Ki_Δ * xiΔq + Kp_Δ * (iΔq_ref -  iΔq) - Leqac*iΔd + Vᴳq)
                vMΔd_ref_c = 2/Vdc*( x[$index+1] +
                            $(converter.controls[:occ].Kₚ) * (iΔd_ref - iΔd) + $Lₑ * iΔq + 1*Vᴳd_fc);
                vMΔq_ref_c = 2/Vdc*( x[$index+2] +
                            $(converter.controls[:occ].Kₚ) * (iΔq_ref - iΔq) - $Lₑ * iΔd + 1*Vᴳq_fc);
                (vMΔd_ref, vMΔq_ref) = I_θ * [vMΔd_ref_c; vMΔq_ref_c]))  # Transformation from converter dq frame to grid dq frame 
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
        push!(exp.args, :(vMΣz_ref = 1))
    end
###########################################################Time delays##################################################################
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
###########################################################MMC State variables##########################################################
    # add state variables
    # x = [iΔd, iΔq, iΣd, iΣq, iΣz, vCΔd, vCΔq, vCΔZd, vCΔZq, vCΣd, vCΣq, vCΣz] = [x[1], x[2], ...]
    # add corresponding differential equations [diΔd_dt, diΔq_dt, ...] = [F[1], F[2], ...]
    # m = [mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz], vM = [vMΔd, vMΔq, vMΔZd, vMΔZq,vMΣd, vMΣq, vMΣz]
    # vM_ref = [vMΔd_ref, vMΔq_ref, vMΔZd_ref, vMΔZq_ref, vMΣd_ref, vMΣq_ref, vMΣz_ref] = s
    push!(exp.args,
    :(
    (mΔd, mΔq, mΔZd, mΔZq, mΣd, mΣq, mΣz) = 1 * [-vMΔd_ref * $baseConv1; -vMΔq_ref * $baseConv1; -vMΔZd_ref * $baseConv1; -vMΔZq_ref * $baseConv1; vMΣd_ref; vMΣq_ref; vMΣz_ref];

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
    # diΔd_dt =-(Vgd - vMΔd + Rₑ*iΔd + Lₑ*iΔq*w)/Lₑ, grid frame
    F[1] = -(inputs[2] * $converter.turnsRatio - vMΔd + $Rₑ*x[1] + $Lₑ*x[2])/$Lₑ;                 
    # diΔq_dt =-(Vgq - vMΔq + Rₑ*iΔq - Lₑ*iΔd*w)/Lₑ, grid frame
    F[2] = -(inputs[3] * $converter.turnsRatio - vMΔq + $Rₑ*x[2] - $Lₑ*x[1])/$Lₑ;                 
    # diΣd_dt =-(vMΣd + Rₐᵣₘ*iΣd - 2*Lₐᵣₘ*iΣq*w)/Lₐᵣₘ, grid 2w frame
    F[3] = -(vMΣd + $Rₐᵣₘ*x[3] - 2*$Lₐᵣₘ*x[4])/$Lₐᵣₘ;                                 
    # diΣq_dt =-(vMΣq + Rₐᵣₘ*iΣq + 2*Lₐᵣₘ*iΣd*w)/Lₐᵣₘ,  grid 2w frame
    F[4] = -(vMΣq + $Rₐᵣₘ*x[4] + 2*$Lₐᵣₘ*x[3])/$Lₐᵣₘ;                                  
    # diΣz_dt =-(vMΣz - Vᵈᶜ/2 + Rₐᵣₘ*iΣz)/Lₐᵣₘ
    F[5] = -(vMΣz - Vdc/2 + $Rₐᵣₘ*x[5])/$Lₐᵣₘ;                                     
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

#######################################################Creating function call and inputs#################################################
    function f!(expr, F, x, inputs) # Creating callable function that takes "x" and ""inputs" as arguments and evaluates function "F" within  the expression "expr"
       f = eval(:((F,x,inputs) -> $expr)) # F(x,inputs)= x^dot
       return Base.invokelatest(f, F,x,inputs) 
    end

    vector_inputs = [Vdc, Vᴳd, Vᴳq] # Inputs to the system for the steady-state solution and linearization. Comes from powerflow results
    init_x=[init_x;zeros(index-length(init_x))]; ################################ATTENTION!!!!!#####################

    # If there is a dc voltage controller, add an additional equation to represent the dc voltage, only for the steady-state solution
    exp_steadyState = copy(exp) # copy the state-space formulation f, only used to obtain steady state
    if in(:dc, keys(converter.controls))  
        init_x =[init_x;Vdc]
        push!(exp_steadyState.args,
        :(
            F[$index+1] = $wbase * ($Idc_in - 3*x[5]) / $Cₑ;
            F[$epsilon_vdc_index] = $(converter.controls[:dc].Kᵢ) * ($(converter.controls[:dc].ref[1]) - x[end]);
        ))
    end

##################################################Steady state solution###############################################################
    # TODO: Define termination conditions
    g!(du,u,p,t) = f!(exp_steadyState, du, u, vector_inputs) # g is the state-space formulation used to obtain the steady-state operation point, copy from f, see some lines above
    println("Starting to solve for Steady-State Solution!")
    prob = SteadyStateProblem(g!, init_x)
    sol=solve(prob,SSRootfind(TrustRegion()),maxiters=20,abstol = 1e-8,reltol = 1e-8)
    
    # Command to show solver results
    # sol.trace
    if SciMLBase.successful_retcode(sol)
        println("MMC steady-state solution found!")
    else
        println("MMC steady-state solution not found!")
    end

    # Delete solution for additional equation in case of DC voltage control
    if in(:dc, keys(converter.controls))
        converter.equilibrium = sol.u[1:end-1] 
    else
        converter.equilibrium = sol.u
    end
    
########################################################Linearization###################################################################
    h(F,x) = f!(exp, F, x[1:end-3], x[end-2:end]) # Using callable function from above 
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index+3); h(F, x); return F)
    A = zeros(index+3,index+3) # Create empty matrix that shall contain the Jacobian of F with respect to x and inputs
    ForwardDiff.jacobian!(A, ha, [converter.equilibrium' vector_inputs']) # Calculating the Jacobian with respect to x and inputs
    converter.A = real(A[1:end-3, 1:end-3])
    converter.B = real(A[1:end-3, end-2:end])

    converter.C = real(zeros(3, size(converter.A,1)))
    converter.C[2,1] = 1  # iΔd in grid frame, generator sign convention exiting converter
    converter.C[3,2] = 1  # iΔq in grid frame, generator sign convention exiting converter
    converter.C[1,5] = 3  # iDC = 3*iΣz, load sign convention entering converter
    converter.D = real(zeros(3,3))

end
###############################################################Calculate admittance######################################################
function eval_parameters(converter :: MMC, s :: Complex)
    # numerical
    I = Matrix{Complex}(Diagonal([1 for dummy in 1:size(converter.A,1)]))
    Y = converter.C * ((s*I-converter.A) \ converter.B) + converter.D # C*((sI-A)^-1)*B+D. This matrix is in pu
    
    #Conversion of admittance from pu to SI
    Y[1,:] *= converter.iDCbase
    Y[:,1] /= converter.vDCbase
    Y[2:3,:] *= converter.iACbase # Base current of the converter side 
    # # # Multiplication with the AC voltage base converts the pu admittance to SI.
    # # # The double division with the turns ratio is actually a multiplication,
    # # # and is needed to bring the grid-side voltage to the converter side.
    Y[:,2:3] /= (converter.vACbase / converter.turnsRatio) # Base voltage at the grid side 
    
    # Structure of Y
    #         Ydcdc   Ydcd   Ydcq
    #    Y=   Yddc    Ydd    Ydq
    #         Yqdc    Yqd    Yqq

    return Y
end







