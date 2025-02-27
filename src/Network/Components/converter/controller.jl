export Controller, PI_control, VSE, CCQSEM

abstract type Controller end


@with_kw mutable struct PI_control <: Controller
    bandwidth :: Union{Int, Float64}  = 0
    ζ :: Union{Int, Float64}  = 0
    Kₚ :: Union{Int, Float64}  = 0                     # proportional gain [pu]
    Kᵢ :: Union{Int, Float64}  = 0                     # integral gain [pu/s]
    ref :: Array{Union{Int, Float64}}  = [0]           # reference value
    n_f :: Int = 1                                     # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                   # Filter cutoff frequency [rad/s]
end

@with_kw mutable struct VSE <: Controller
    H :: Union{Int, Float64}  = 0                       # Virtual Inertia [s]
    K_d :: Union{Int, Float64}  = 0                     # Damping coefficient [-]
    K_ω    :: Union{Int, Float64}  = 0                  # Droop coefficient [-]
    ref :: Array{Union{Int, Float64}}  = [0]            # Active power reference [MW]
    ref_ω :: Union{Int, Float64}  = 0                   # Angular frequency reference [pu]
    n_f :: Int = 0                                      # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                    # Butterworth Filter cutoff frequency [rad/s]
end

@with_kw mutable struct CCQSEM <: Controller                 
    Rᵥ :: Union{Int, Float64}  = 0                       # Virtual Resistance [pu]
    Lᵥ :: Union{Int, Float64}  = 0                       # Virtual inductance [pu]
    ref_vd :: Union{Int, Float64}  = 0                   # Vd voltage reference [pu]
    ref_vq :: Union{Int, Float64}  = 0                   # Vq voltage reference [pu]
    n_f :: Int  = 0                                      # Butterworth Filter order [-]
    ω_f  :: Union{Int, Float64}  = 0                     # Butterworth Filter cutoff frequency [rad/s]


end

#############################################Additional functions used within the model#################################################


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

function butterworthMatrices(buttOrder,ω_c,numberVars)

    # Calculation of the state-space representation of a n-order butterworth filter with a gain of 1.
    # buttOrder = Order of butterworth filter, ω_c= Cutoff frequency of the filter in [rad/s], numberVars= Variables to be filtered max.2 !
   
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
    

    # Controllable canonical form can create problems while solving for MMC steady-state, espec. when high-order pade delay approximation is used. 
    # Related to high conditioning number of controllable canonical form.
    # Transformation to modal form, which results in lower conditioning number.
    # sys = ss(Ab,Bb,Cb,Db)
    # sys_modal = modal_form(sys;C1 = true)
    # Ab= sys_modal[1].A
    # Bb = sys_modal[1].B
    # Cb = sys_modal[1].C
    # Db = sys_modal[1].D


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