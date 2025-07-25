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

    Lᵣ :: Union{Int, Float64}  = 60e-3         # inductance of the converter transformer at the converter side [H]
    Rᵣ :: Union{Int, Float64}  = 0.535         # resistance of the converter transformer at the converter side [H]

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
    vDCbase :: Union{Int, Float64} = 640        # DC voltage base [kV]

    vACbase :: Float64 = 0 # AC voltage base for impedance/admittance calculation
    iACbase :: Float64 = 0 # AC current base for impedance/admittance calculation

    debug = nothing
end

"""
    function tlc(;args...)
It constructs tlc operating both as a rectifier and an inverter. TLC is constructed as a struct with the
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
            if key ∈ (:vac_supp, :v_ac)
                error("AC voltage support is not yet implemented")
            end
            converter.controls[key] = val
        elseif in(key, propertynames(converter))
            setfield!(converter, key, val)
        end
    end

    elem = Element(input_pins = 2, output_pins = 2, element_value = converter)
end

function update!(converter :: TLC, Pac, Qac, Vm, θ)
    

    wbase = 100*pi
    vAC_base = converter.vACbase_LL_RMS*sqrt(2/3)
    Sbase = converter.Sbase
    vDC_base = converter.vDCbase
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
    # converter.Vᵈᶜ = Vdc
    Vdc = converter.Vᵈᶜ
    converter.P = Pac
    converter.Q = Qac
    # converter.P_dc = Pdc # Has the same sign as Pac
    Pdc = converter.P_dc

    Vm /= vAC_base
    Vdc /= vDC_base
    Pac /= Sbase
    Qac /= Sbase
    Pdc /= Sbase
    
    Vᴳd = Vm * cos(θ)
    Vᴳq = -Vm * sin(θ)

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
        index_PLL_angle = index + 1*(converter.controls[:pll].n_f) +2
        if in(:v_meas_filt, keys(converter.controls))
            index_PLL_angle +=  2*(converter.controls[:v_meas_filt].n_f)  
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
        Abutt, Bbutt, Cbutt, Dbutt =  butterworthMatrices(converter.controls[:v_meas_filt].n_f, converter.controls[:v_meas_filt].ω_f, 2);
        push!(exp.args, :(
            voltagesIn = [Vᴳd;Vᴳq];
            statesButt= x[$index + 1 : $index + 2*$(converter.controls[:v_meas_filt].n_f)]; 
            F[$index + 1 : $index + 2*$(converter.controls[:v_meas_filt].n_f)] = $Abutt*statesButt + $Bbutt*voltagesIn;
            voltagesOut=$Cbutt*statesButt+$Dbutt*voltagesIn;
            Vᴳd_f=voltagesOut[1];
            Vᴳq_f=voltagesOut[2];
            ))
        index += 2*(converter.controls[:v_meas_filt].n_f) 
        init_x = [init_x;Vᴳd;Vᴳq] #Initalize to avoid steady-state solver problems
    else
        push!(exp.args, :(
            (Vᴳd_f, Vᴳq_f) = (Vᴳd, Vᴳq);
            ))
    end
    

    # add PLL
    if in(:pll, keys(converter.controls))
        if (converter.controls[:pll].ω_f != 0) # A PLL filter is implemented
            Abutt_pll, Bbutt_pll, Cbutt_pll, Dbutt_pll =  butterworthMatrices(converter.controls[:pll].n_f, converter.controls[:pll].ω_f, 1);
            push!(exp.args, :(               
                statesButt_pll= x[$index + 1 : $index + 1*$(converter.controls[:pll].n_f)]; 
                F[$index + 1 : $index + 1*$(converter.controls[:pll].n_f)] = $Abutt_pll*statesButt_pll + $Bbutt_pll*Vᴳq_f;
                vₚₗₗ=dot($Cbutt_pll,statesButt_pll)+$Dbutt_pll*Vᴳq_f;# Get rid of 1-element array
            ))

           
            index += 1*(converter.controls[:pll].n_f)
            
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
        Abutt_i, Bbutt_i, Cbutt_i, Dbutt_i =  butterworthMatrices(converter.controls[:i_meas_filt].n_f, converter.controls[:i_meas_filt].ω_f, 2);
        push!(exp.args, :(
            currentsIn = [i_d_pcc_c;i_q_pcc_c];
            statesButt_i= x[$index + 1 : $index + 2*$(converter.controls[:i_meas_filt].n_f)]; 
            F[$index + 1 : $index + 2*$(converter.controls[:i_meas_filt].n_f)] = $Abutt_i*statesButt_i + $Bbutt_i*currentsIn;
            currentsOut=$Cbutt_i*statesButt_i+$Dbutt_i*currentsIn;
            i_d_pcc_f=currentsOut[1];
            i_q_pcc_f=currentsOut[2];
            ))
        # init_x = [init_x;zeros(index-length(init_x))]
        index += 2*(converter.controls[:i_meas_filt].n_f)
        # init_x = [init_x;Id;Iq]

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
            converter.controls[:vac_supp].ref[1] /= vAC_base #Per unitize voltage reference
            push!(exp.args, :(
                Vᴳ_mag = sqrt(Vᴳd_f^2 + Vᴳq_f^2);
                Δq_unf = $(converter.controls[:vac_supp].Kₚ)*($(converter.controls[:vac_supp].ref[1])-Vᴳ_mag);
                F[$index+1] = $(converter.controls[:vac_supp].ω_f) *(Δq_unf - x[$index+1]);
                q_ref = $(converter.controls[:q].ref[1]) + x[$index+1]))
            # init_x = [init_x;zeros(index-length(init_x))] #Initalize states before voltage support to zero
            index +=1
            # init_x = [init_x;converter.controls[:vac_supp].ref[1]]
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
         #Small value to converge
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
                timeDelayIn = [md;mq];
                statesDelay = x[$index + 1 : $index + 2*$converter.padeOrderDen]; 
                timeDelayOut = timeDelayPadeMatrices($converter.padeOrderNum,$converter.padeOrderDen,$converter.timeDelay,length(timeDelayIn));
                
                A_delay = timeDelayOut[1];
                B_delay = timeDelayOut[2];
                C_delay = timeDelayOut[3];
                D_delay = timeDelayOut[4];
                F[$index + 1 : $index + 2*$converter.padeOrderDen] = A_delay*statesDelay + B_delay*timeDelayIn;
                # timeDelayOut = C_delay*statesDelay + D_delay*timeDelayIn;
                # Implement phase shifts by transforming the dq voltage references to alpha-beta
                m_ab_ref = (cos($converter.ω₀*$converter.timeDelay)-sin($converter.ω₀*$converter.timeDelay)*im)*(T_dq_ab*(C_delay*statesDelay + D_delay*timeDelayIn));
                m_dq_ref = real(T_ab_dq * conj(m_ab_ref) + conj(T_ab_dq) * m_ab_ref);
                md = m_dq_ref[1];
                mq = m_dq_ref[2];
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

    function g(u)
        du = similar(u)
        f!(exp_equilibrium, du, u, vector_inputs)
        return du
    end

    vector_inputs = [Vdc, Vᴳd, Vᴳq]
    init_x = [init_x; zeros(index-length(init_x))]
    # if in(:v_meas_filt, keys(converter.controls))
    #     init_x[indexVᴳdf] =  Vᴳd
    #     init_x[indexVᴳqf] =  1e-3 # Initialize to a small non-zero value to avoid Inf or NaN problems with nlsolve
    # end
    ##################################################Steady state solution###############################################################
    
    g!(du,u,p,t) = f!(exp_equilibrium, du, u, vector_inputs) # g is the state-space formulation used to obtain the steady-state operation point, copy from f, see some lines above
 
    println("Starting to solve for Steady-State Solution!")
    prob = SteadyStateProblem(g!, init_x)
    sol=solve(prob,SSRootfind(TrustRegion()),maxiters=20,abstol = 1e-8,reltol = 1e-8)
    
    steady_state_jacobian = ForwardDiff.jacobian(g, init_x)
    converter.debug = [steady_state_jacobian, sol]

    # Command to show solver results
    # println(sol.trace)
    # converter.debug = sol
    if SciMLBase.successful_retcode(sol)
        println("TLC steady-state solution found!")
    else
        println("TLC steady-state solution not found!")
    end
    
    converter.equilibrium = sol.u

    number_output = 3
    number_input =3

    h(F,x) = f!(exp, F, x[1:end-number_input], x[end-number_input+1:end])
    ha = x -> (F = fill(zero(promote_type(eltype(x), Float64)), index + number_output); h(F, x); return F) # 
    jac = zeros(index + number_output , index + number_input)
    ForwardDiff.jacobian!(jac, ha, [converter.equilibrium' vector_inputs'])
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

function make_power_flow!(tlc:: TLC, data, nodes2bus, bus2nodes, elem2comp, comp2elem, elem, global_dict)

    # Check if AC or DC source (second one not implemented)
    # is_three_phase(elem) ? nothing : error("DC sources are currently not implemented")

    ### MAKE BUSES OUT OF THE NODES
    # Find the nodes not connected to the ground
    ac_nodes = make_non_ground_node(elem, bus2nodes) 
    ac_bus = add_bus_ac!(data, nodes2bus, bus2nodes, ac_nodes, global_dict)
    # Make busses for the non-ground nodes 
    interm_bus = add_interm_bus_ac!(data, global_dict) # No mapping to node, bcs no corresponding node in PowerImpedance

    # Make the load component for injection
    key_load = comp_elem_interface!(data, elem2comp, comp2elem, elem, "load")

    (data["load"])[string(key_load)] = Dict{String, Any}()
    ((data["load"])[string(key_load)])["status"] = 1
    data["load"][string(key_load)]["load_bus"] = interm_bus
    data["load"][string(key_load)]["pd"] = -tlc.P / (global_dict["S"] / 1e6)
    data["load"][string(key_load)]["qd"] = -tlc.Q / (global_dict["S"] / 1e6)

    # Add additional branch & bus for SM transformer (RL-branch)
    
    key_branch = length(data["branch"]) + 1
    key_branch_str = string(key_branch)

    (data["branch"])[key_branch_str] = Dict{String, Any}()
    ((data["branch"])[key_branch_str])["f_bus"] = interm_bus
    ((data["branch"])[key_branch_str])["t_bus"] = ac_bus
    ((data["branch"])[key_branch_str])["source_id"] = Any["branch", key_branch]
    ((data["branch"])[key_branch_str])["index"] = key_branch
    ((data["branch"])[key_branch_str])["rate_a"] = 1
    ((data["branch"])[key_branch_str])["rate_b"] = 1
    ((data["branch"])[key_branch_str])["rate_c"] = 1
    ((data["branch"])[key_branch_str])["br_status"] = 1
    ((data["branch"])[key_branch_str])["angmin"] = ang_min
    ((data["branch"])[key_branch_str])["angmax"] = ang_max
    ((data["branch"])[key_branch_str])["transformer"] = false
    ((data["branch"])[key_branch_str])["tap"] = 1
    ((data["branch"])[key_branch_str])["shift"] = 0
    ((data["branch"])[key_branch_str])["c_rating_a"] = 1

    
    ((data["branch"])[key_branch_str])["br_r"] = tlc.Rᵣ / global_dict["Z"] # Convert SI to pu
    ((data["branch"])[key_branch_str])["br_x"] = tlc.Lᵣ*global_dict["omega"] / global_dict["Z"]
    ((data["branch"])[key_branch_str])["g_fr"] = 0
    ((data["branch"])[key_branch_str])["b_fr"] = 0
    ((data["branch"])[key_branch_str])["g_to"] = 0
    ((data["branch"])[key_branch_str])["b_to"] = 0

    # Only PQ-bus implemented
    
    ((data["bus"])[string(interm_bus)]) = set_bus_type((data["bus"])[string(interm_bus)], 1)
    

    # ((data["bus"])[string(interm_bus)])["vm"] = ((data["gen"])[key])["vg"]
    # ((data["bus"])[string(interm_bus)])["vmin"] =  0.9*((data["gen"])[key])["vg"]
    # ((data["bus"])[string(interm_bus)])["vmax"] =  1.1*((data["gen"])[key])["vg"]
end
