@with_kw mutable struct Autotransformer
    Zₕₗ :: Union{Int, Float64} = 0      # leakage reactance HL [p.u.]
    Zₕₜ :: Union{Int, Float64} = 0      # leakage reactance HT [p.u.]
    Zₜₗ :: Union{Int, Float64} = 0      # leakage reactance TL [p.u.]

    S :: Union{Int, Float64} = 0      # nominal power [MVA]
    f :: Union{Int, Float64} = 50      # base frequency [Hz]

    pins :: Int = 9
end

function autotransformer(;args...)
    t = Autotransformer()
    for (key, val) in kwargs_pairs(args)
        if in(key, propertynames(t))
            setfield!(t, key, val)
        end
    end

    elem = Element(input_pins = t.pins, output_pins = t.pins, element_value = t)
end

function eval_parameters(t :: Autotransformer, s :: Complex)
    a = exp(1im*2*π/3)
    dqtrfbase = [1  1   1   ;
                 1  a^2   a ;
                 1  a   a^2 ]
    dqtf = zeros(12,12)
    for i in 1:4
        dqtf[1+3(i-1):3i] = dqtrfbase
    end
    dqtrfinv = inv(dqtrf)

    ZHX = real(t.Zₕₗ) + imag(Z.ₕₗ) * s / 1im / 2 / π / t.f
    ZHY = real(t.Zₕₜ) + imag(Z.ₕₜ) * s / 1im / 2 / π / t.f
    ZXY = real(t.Zₜₗ) + imag(Z.ₜₗ) * s / 1im / 2 / π / t.f

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

    Yseq = [YH0 0   0  0   0   0   0   0   0   -YH0        0   0;
            0  YH   0  0   0   0   0   0   0   0         -YH   0;
            0   0  YH  0   0   0   0   0   0   0           0   -YH;
            0   0   0  YX0 0   0   0   0   0   -YX0        0   0;
            0   0   0  0  YX   0   0   0   0   0         -YX   0;
            0   0   0  0   0  YX   0   0   0   0           0 -YX;
            0   0   0  0   0   0   0   0   0   0           0   0;
            0   0   0  0   0   0   0  YY   0   0         -YY   0;
            0   0   0  0   0   0   0   0  YY   0          0  -YY ;
            -YH 0   0  -YX 0   0   0   0   0   YH0+YX0+YT0 0    0;
            0 -YH   0  0 -YX   0   0 -YY   0   0   YH+YX+YY    0;
            0   0 -YH0 0   0 -YX   0   0 -YY   0          0  YH+YX+YY]

    Yseq = dqtrf*Yseq*dqtrfinv
    no_eliminate = [i for i in 1:9]
    Y = kron(Yseq, no_eliminate)
end

function eval_abcd(t :: Autotransformer, s :: Complex)
    Y = eval_parameters(t,s)
    abcd = y_to_abcd(Y)
end

function eval_y(t :: Autotransformer, s :: Complex)
    Y = eval_parameters(t,s)
end
