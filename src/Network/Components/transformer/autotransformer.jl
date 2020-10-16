export autotransformer

@with_kw mutable struct Autotransformer
    ABCD :: Array{Basic} = []          # ABCD value

    pins :: Int = 3                     # marks single or three phase


    Xₕₗ :: Union{Int, Float64} = 0      # leakage reactance HL [p.u.]
    Xₕₜ :: Union{Int, Float64} = 0      # leakage reactance HT [p.u.]
    Xₜₗ :: Union{Int, Float64} = 0      # leakage reactance TL [p.u.]

    S :: Union{Int, Float64} = 0      # nominal power [MVA]
    f :: Union{Int, Float64} = 50      # base frequency [Hz]
    ω :: Union{Int, Float64} = 2*π*50   # rated frequency in [rad/s]

end

"""
function autotransformer(;args...)
    Implements autotransformer from the zero, positive and negative sequence.
"""
function autotransformer(;args...)
    at = Autotransformer()
    transformation = false
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(at))
            setfield!(at, key, val)
        elseif (key == :transformation)
            transformation = val
        else
            throw(ArgumentError("Autotransformer does not have a property $(key)."))
        end
    end

# abcd parameters calculation
s = symbols(:s)

a = exp(1im*2*π/3)
dqtrfbase = [1  1   1   ;
             1  a^2   a ;
             1  a   a^2 ]
dqtrf =convert(Array{Basic}, zeros(12,12))
for i in 1:4, j in 1:4
    if i==j
    dqtrf[1+3(i-1):3i , 1+3(j-1):3j ] = dqtrfbase
    end
end
dqtrfinv = inv(dqtrf)


    ZHX = 1im*(at.Xₕₗ)
    ZHY = 1im*(at.Xₕₜ)
    ZXY = 1im*(at.Xₜₗ)

    ZY = 1/2*(ZXY+ZHY-ZHX)
    ZX = 1/2*(ZXY+ZHX-ZHY)
    ZH = 1/2*(ZHY+ZHX-ZXY)

    ZX0 = 1/2*(ZHX-ZHY+ZXY)
    ZH0 = 1/2*(ZHX+ZHY-ZXY)
    ZT0 = 1/2*(-ZHX+ZHY+ZXY)

    YH = 1/ZH;
    YX = 1/ZX;
    YY = 1/ZY;
    YX0 = 1/ZX0;
    YH0 = 1/ZH0;
    YT0 = 1/ZT0;
# calculating Y for posetive, negative and zero sequence of 3 windings and Neutral point


    Yseq = [YH  0  0   0  0  0   0  0  0   -YH    0       0;
            0  YH  0   0  0  0   0  0  0   0     -YH      0;
            0   0 YH0  0  0  0   0  0  0   0      0    -YH0;
            0   0  0   YX 0  0   0  0  0   -YX    0       0;
            0   0  0   0  YX 0   0  0  0   0     -YX      0;
            0   0  0   0  0  YX0 0  0  0   0      0    -YX0;
            0   0  0   0  0  0   YY 0  0   -YY    0       0;
            0   0  0   0  0  0   0  YY 0   0     -YY      0;
            0   0  0   0  0  0   0  0  0   0      0       0;
          -YH   0  0  -YX 0  0  -YY 0  0 YH+YX+YY 0       0;
            0 -YH  0   0 -YX 0   0 -YY 0   0   YH+YX+YY   0;
            0   0 -YH0 0  0 -YX0 0  0  0   0      0   YH0+YX0+YT0]
# calculating 3 phase (abc) from sequence Y
    Yphase     = dqtrf*Yseq*dqtrfinv;
    YTemp           = Yphase;
    iB=[i for i in 1:6]
    iT=[i for i in 1:6]
    n_iB =length(iB);
    n_iT =length(iT);

    YTemp[:,1:n_iT] = YTemp[:,iT];          #move columns of test buses to first columns
    YTemp[:,iT]     =Yphase[:,1:n_iT];          #move first columns to columns of test buses

    YTemp2          = YTemp;
    YTemp[1:n_iB,:] = YTemp[iB,:];          #move rows of injection buses to first rows
    YTemp[iB,:]     = YTemp2[1:n_iB,:];    #move first rows to rows of injection buses

    #perform actual elimination
    #See also https://wiki.openelectrical.org/index.php?title=Kron_Reduction
    Yaa             = YTemp[1:n_iT,1:n_iT];
    Ye1             = YTemp[1:n_iT,n_iT+1:end];
    Ye2             = YTemp[n_iT+1:end,1:n_iT];
    Yee             = YTemp[n_iT+1:end,n_iT+1:end];
    Y_red           = Yaa-(Ye1*inv(Yee)*Ye2)


    m = 3; #2N/2 dimension

    # Allocate memory for the ABCD-parameters
    #abcd_params = zeros(Complex{Float64}, 6,6);

            # Get the Y-parameters
           y11 = Y_red[1:m,1:m];
           y12 = Y_red[1:m,(m+1):(2*m)];
           y21 = Y_red[(m+1):(2*m),1:m];
           y22 = Y_red[(m+1):(2*m),(m+1):(2*m)];

            # Calculate the ABCD-parameters

         A= -y21\y22;
         B=-inv(y21);
         C= y12-y11*inv(y21)*y22;
         D=-y11*inv(y21);

         ABCD=[A B; C D]
         at.ABCD=real(ABCD)+imag(ABCD)*s/at.ω

elem = Element(input_pins = at.pins, output_pins = at.pins, element_value = at,
             transformation = transformation)
end

function eval_abcd(at :: Autotransformer, s :: Complex)
    value = N.(subs.(at.ABCD, symbols(:s), s))
    return convert(Array{Complex}, value)
end

function eval_y(at :: Autotransformer, s :: Complex)
    return abcd_to_y(eval_abcd(at, s))
end
