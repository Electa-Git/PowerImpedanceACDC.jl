export autotransformer
s = symbols(:s)
#Data from Deliverable 1
@with_kw mutable struct Autotransformer
    ABCD :: Union{Array{Any}, Array{Basic}}  = []
    Xₕₗ :: Union{Int, Float64,Complex{Float64}} = 0.15      # leakage reactance HL [p.u.]
    Xₕₜ :: Union{Int, Float64,Complex{Float64}} = 0.02      # leakage reactance HT [p.u.]
    Xₜₗ :: Union{Int, Float64,Complex{Float64}} = 0.02      # leakage reactance TL [p.u.]
    #Vₕᵥ :: Union{Int, Float64} = 0    # Voltage HV winding [kV?]
    #Vₗᵥ :: Union{Int, Float64} = 0    #Voltage LV winding[kV?]
    #Vₜᵥ :: Union{Int, Float 64} = 0   #Voltage tertiary winding[kV?]
    S :: Union{Int, Float64} = 0      # nominal power [MVA]
    f :: Union{Int, Float64} = 50      # base frequency [Hz]
    ω :: Union{Int, Float64} = 2*π*f   # rated frequency in [rad/s]
    pins :: Int = 3       #why 9? 3 HV, 3LV, 3Tertiary???
end
"""
function autotransformer( )
    Implements autotransformer from the zero, positive and negative sequence.
"""
function autotransformer(;args...)
    transformation = false
    Xₕₗ :: Union{Int, Float64,Complex{Float64}} = 0.15      # leakage reactance HL [p.u.]
    Xₕₜ :: Union{Int, Float64,Complex{Float64}} = 0.02      # leakage reactance HT [p.u.]
    Xₜₗ :: Union{Int, Float64,Complex{Float64}} = 0.02      # leakage reactance TL [p.u.]
    #Vₕᵥ :: Union{Int, Float64} = 0    # Voltage HV winding [kV?]
    #Vₗᵥ :: Union{Int, Float64} = 0    #Voltage LV winding[kV?]
    #Vₜᵥ :: Union{Int, Float 64} = 0   #Voltage tertiary winding[kV?]
    S :: Union{Int, Float64} = 400e6      # nominal power [MVA]
    f :: Union{Int, Float64} = 50      # base frequency [Hz]
    ω :: Union{Int, Float64} = 2*π*f   # rated frequency in [rad/s]
    pins :: Int = 3       #why 9? 3 HV, 3LV, 3Tertiary???
    at = Autotransformer()
    at.Xₕₗ = Xₕₗ;
    at.Xₕₜ = Xₕₜ;
    at.Xₜₗ = Xₜₗ;
    at.S = S;
    at.f = f;
    at.ω = ω;
    at.pins = pins;
    element = Element(element_value = at, input_pins = at.pins, output_pins = at.pins, transformation = transformation)
    return element
end
function eval_abcd(at :: Autotransformer, s :: Complex)
    Xₕₗ = at.Xₕₗ;
    Xₕₜ = at.Xₕₜ;
    Xₜₗ = at.Xₜₗ;
    S = at.S;
    f = at.f;
    ω = at.ω;
    pins = at.pins;

    #Sequence transformation operator
    a = exp(1im*2*π/3)
    #Transformation matrix
    dqtrfbase = [1  1   1   ;
                1  a^2   a ;
                1  a   a^2 ]

    #Creation of a 12*12 sequence transformation matrix
    dqtrf = zeros(Complex{Float64},12,12)
    for i in 1:4
        if i==1
            k=0;
        elseif i==2
            k=3;
        elseif i==3
            k=6;
        elseif i==4
            k=9;
        end
        for j in 1:3
            for z in 1:3
                dqtrf[k+j,k+z]=dqtrfbase[j,z]
            end
        end

    end

    #inverse of the 12*12 sequence transformation matrix
    dqtrfinv = inv(dqtrf)

    ZHX = s*(at.Xₕₗ)/at.ω
    ZHY = s*(at.Xₕₜ)/at.ω
    ZXY = s*(at.Xₜₗ)/at.ω

    #pag 155 book Protective Relaying Principles and Applications & Willem's matlab code
    ZY = 1/2*(ZXY+ZHY-ZHX)
    ZX = 1/2*(ZXY+ZHX-ZHY)
    ZH = 1/2*(ZHY+ZHX-ZXY)

    ZX0 = 1/2*(ZHX-ZHY+ZXY)
    ZH0 = 1/2*(ZHX+ZHY-ZXY)
    ZY0 = 1/2*(-ZHX+ZHY+ZXY)
    #Finding the corresponding admittance values
    YH = 1/ZH;
    YX = 1/ZX;
    YY = 1/ZY;
    YX0 = 1/ZX0;
    YH0 = 1/ZH0;
    YY0 = 1/ZY0;

    Yseq = [YH  0   0  0  0  0   0  0  0   -YH       0          0;
            0  YH   0  0  0  0   0  0  0    0       -YH         0;
            0   0  YH0 0  0  0   0  0  0    0        0        -YH0;
            0   0   0  YX 0  0   0  0  0   -YX       0          0;
            0   0   0  0  YX 0   0  0  0    0       -YX         0;
            0   0   0  0  0 YX0  0  0  0    0        0        -YX0;
            0   0   0  0  0  0   YY 0  0   -YY       0          0;
            0   0   0  0  0  0   0  YY 0    0       -YY         0;
            0   0   0  0  0  0   0  0  0    0        0          0;
            -YH 0   0 -YX 0  0  -YY 0  0 YH+YX+YY    0          0;
            0  -YH  0  0 -YX 0   0 -YY 0    0     YH+YX+YY      0;
            0   0 -YH0 0  0 -YX0 0  0  0    0        0     YH0+YX0+YY0]

    Yphase = dqtrf*Yseq*dqtrfinv

    #Kron reduction, it brings the Yseq[12x12] matrix into Y[9x9] matrix
    #Not used the Aleksandra's kron function. Checked the Ymatrix with the one from Willem's code -> same result
    red = 6;
    Yaa = Yphase[1:red, 1:red]
    Ye1 = Yphase[1:red, red+1:12]
    Ye2 = Yphase[red+1:12, 1:red]
    Yee = Yphase[red+1:12, red+1:12]
    Y= Yaa - Ye1*inv(Yee)*Ye2

    #Passing from y to ABCD formulation
    Ypp = Yps = Ysp = Yss = zeros(Complex{Float64},3,3) #definition of the matrices needed for ABCD transformation

    m = Int16(sqrt(length(Y))/2)

    Ypp = Y[1:m,1:m]
    Yps = Y[1:m,m+1:2m]
    Ysp = Y[m+1:2m,1:m]
    Yss = Y[m+1:2m,m+1:2m]

    #Formulas from section 3.3 simulator_tutorial
    A = -inv(Ysp)*Yss
    B = -inv(Ysp)
    C = Yps - Ypp*inv(Ysp)*Yss
    D = -inv(Ysp)*Ypp
    at.ABCD= [A B; C D] #ABCD autotransformer formulation -> verified with the MATLAB fnc ytoabcd -> same result
    return ABCD
end



function eval_y(at :: Autotransformer, s :: Complex)
    return abcd_to_y(eval_abcd(t, s))
end
#function eval_abcd(at :: Autotransformer, s :: Complex)
    #Y = eval_parameters(at,s)
    #abcd = y_to_abcd(Y)
#end

#function eval_y(at :: Autotransformer, s :: Complex)
    #Y = eval_parameters(at,s)
#end
