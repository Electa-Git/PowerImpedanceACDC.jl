function kron(matrix :: Array{Complex}, no_eliminate :: Array{Int})
    n = size(matrix, 1)
    eliminate = setdiff(1:n, no_eliminate)
    matrix = matrix[[no_eliminate; eliminate],:]
    matrix = matrix[:,[no_eliminate; eliminate]]

    n_noElim = length(no_eliminate)
    matrix = matrix[1:n_noElim,1:n_noElim] - matrix[1:n_noElim,1+n_noElim:end]*
            inv(matrix[1+n_noElim:end,1+n_noElim:end])*matrix[1+n_noElim:end,1:n_noElim]
end

function kron_abcd(matrix :: Array{Complex}, no_eliminate :: Array{Int})
    n = Int(size(matrix, 1)/2)
    eliminate = setdiff(1:n, no_eliminate)

    (a,b,c,d) = (matrix[1:n,1:n], matrix[1:n,n+1:end], matrix[n+1:end,1:n], matrix[n+1:end, n+1:end])
    (a11, a12, a21, a22) = (a[no_eliminate,no_eliminate], a[no_eliminate, eliminate], a[eliminate, no_eliminate], a[eliminate, eliminate])
    (b11, b12, b21, b22) = (b[no_eliminate,no_eliminate], b[no_eliminate, eliminate], b[eliminate, no_eliminate], b[eliminate, eliminate])
    (c11, c12, c21, c22) = (c[no_eliminate,no_eliminate], c[no_eliminate, eliminate], c[eliminate, no_eliminate], c[eliminate, eliminate])
    (d11, d12, d21, d22) = (d[no_eliminate,no_eliminate], d[no_eliminate, eliminate], d[eliminate, no_eliminate], d[eliminate, eliminate])

    e = inv(b22)
    a = a11 - b12*e*a21
    b = b11 - b12*e*b21
    c = c11 - d12*e*a21
    d = d11 - d12*e*b21

    return [a b; c d]
end
