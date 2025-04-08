const Parameters_types = Union{Array{Basic}, Array{Complex},
        Array{Int}, Array{Float64}, Array{ComplexF64}}

function connect_series!(a::Parameters_types, b::Parameters_types)
    return a*b
end

function connect_parallel!(ABCD₁::Parameters_types, ABCD₂::Parameters_types)
    n = Int(size(ABCD₁, 1)/2)
    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    if (n == 1)
        a = (b₂[1] * a₁[1] + b₁[1] * a₂[1]) / (b₁[1] + b₂[1])
        b = b₁[1] * b₂[1] / (b₁[1] + b₂[1])
        c = c₁[1] + c₂[1] + (d₂[1] - d₁[1]) * (a₁[1] - a₂[1]) / (b₁[1] + b₂[1])
        d = d₁[1] + (d₂[1] - d₁[1]) * b₁[1] / (b₁[1] + b₂[1])
    else
        I = convert(Array{Basic}, Diagonal([1 for dummy in 1:n]))
        if all(b₁[i] == 0 for i in 1:n)
            a = a₂
            b = zeros(Basic, n, n)
            c = c₁ + c₂ + (d₂ - d₁) * (b₂ \ (a₁ - a₂))
            d = d₁
        elseif all(b₂[i] == 0 for i in 1:n)
            a = a₁
            b = zeros(Basic, n, n)
            c = c₁ + c₂ + (d₂ - d₁) * (b₁ \ (a₁ - a₂))
            d = d₁
        else
            a = inv(inv(b₁) + inv(b₂))  * (inv(b₁) * a₁ + inv(b₂) * a₂)
            b = inv(inv(b₁) + inv(b₂))
            c = c₁ + c₂ + (d₂ - d₁) * inv(b₁ + b₂) * (a₁ - a₂)
            d = d₁ + (d₂ - d₁) * inv(b₁ + b₂) * b₁
        end
    end

    ABCD = vcat(hcat(a,b), hcat(c,d))
    return ABCD
end

function closing_impedance(ABCD :: Array{Complex}, Zₜ :: Union{Array{Complex}, Int, Float64, Complex}, direction = :output)
    n = Int(size(ABCD, 1)/2)
    m = Int(size(ABCD, 2)/2)
    (a, b, c, d) = (ABCD[1:n,1:m], ABCD[1:n,m+1:end], ABCD[n+1:end,1:m], ABCD[n+1:end, m+1:end])

    Zₑ = 0
    if (length(Zₜ) == 1)
        if (direction == :output)
            Zₑ = (a .* Zₜ + b) ./ (c .* Zₜ + d)
        else
            Zₑ = (d .* Zₜ - b) ./ (c .* Zₜ - a)
        end
    else
        I = convert(Array{Complex}, Diagonal([1 for dummy in 1:n]))
        if (direction == :output)
            Zₑ = (a * Zₜ + b) * pinv(c * Zₜ + d)
        else
            Zₑ = pinv(Zₜ * c - a) * (Zₜ * d - b)
        end
    end
    return Zₑ
end

# making 2×2 matrix (modal domain) from 4×4 matrix (phase domain)
function transformation_dc(ABCD :: Parameters_types)
    n = Int(size(ABCD, 1)/2)
    (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])

    ABCD = [(a[1,1]+a[2,2]-a[1,2]-a[2,1])/2 (b[1,1]+b[2,2]-b[1,2]-b[2,1])
            (c[1,1]+c[2,2]-c[1,2]-c[2,1])/4 (d[1,1]+d[2,2]-d[1,2]-d[2,1])/2]
end
# making 4×4 matrix (dq domain) from 6x6 matrix (phase domain)
function transformation_dq(ABCD₁, ABCD₂)
    n = Int(size(ABCD₁, 1)/2)
    (a₁, b₁, c₁, d₁) = (ABCD₁[1:n,1:n], ABCD₁[1:n,n+1:end], ABCD₁[n+1:end,1:n], ABCD₁[n+1:end, n+1:end])
    (a₂, b₂, c₂, d₂) = (ABCD₂[1:n,1:n], ABCD₂[1:n,n+1:end], ABCD₂[n+1:end,1:n], ABCD₂[n+1:end, n+1:end])

    # ϕ = 2π/3
    # e = exp(1im*ϕ)
    # a = [1 e e^2;
    #     1im 1im*e 1im*e^2;
    #     0 0 0]

    # a_dq = (1/3 * (a * a₁ + conj(a) * a₂) * transpose(real(a)))[1:2,1:2]
    # b_dq = (1/3 * (a * b₁ + conj(a) * b₂) * transpose(real(a)))[1:2,1:2]
    # c_dq = (1/3 * (a * c₁ + conj(a) * c₂) * transpose(real(a)))[1:2,1:2]
    # d_dq = (1/3 * (a * d₁ + conj(a) * d₂) * transpose(real(a)))[1:2,1:2]
    # SOURCE: Impedance transformation from dq to alpha beta and positive negative sequence frame by Eros Avdiaj, Philippe De Rua
    T = 0.5 * [1 -1im;-1im -1]
    CK = (2/3)*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2]
    CKinv = [1 0;-1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]

    a_dq = T * CK * a₂ * CKinv * conj(T) + conj(T) * CK * a₁ * CKinv * T
    b_dq = T * CK * b₂ * CKinv * conj(T) + conj(T) * CK * b₁ * CKinv * T
    c_dq = T * CK * c₂ * CKinv * conj(T) + conj(T) * CK * c₁ * CKinv * T
    d_dq = T * CK * d₂ * CKinv * conj(T) + conj(T) * CK * d₁ * CKinv * T

    abcd = [a_dq b_dq; c_dq d_dq]
end

function y_to_abcd(Y :: Parameters_types)
    n = Int(size(Y,1)/2)
    (Yee, Yei, Yie, Yii) = (Y[1:n,1:n], Y[1:n,n+1:end], Y[n+1:end,1:n], Y[n+1:end,n+1:end])
    Yie = inv(Yie)
    abcd = [-Yie*Yii -Yie; Yei-Yee*Yie*Yii -Yie*Yee]
end

function abcd_to_y(ABCD :: Parameters_types)
    n = Int(size(ABCD,1)/2)
    (a, b, c, d) = (ABCD[1:n,1:n], ABCD[1:n,n+1:end], ABCD[n+1:end,1:n], ABCD[n+1:end, n+1:end])
    b = inv(b)
    Y = [d*b c-d*b*a; -b b*a]
end
