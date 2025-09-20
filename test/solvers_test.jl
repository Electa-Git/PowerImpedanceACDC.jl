
# Example 7.8 from Grainger
Y= [-im*16.75 im*11.75 im*2.50 im*2.50;
    im*11.75 -im*19.25 im*2.5 im*5.00;
    im*2.50 im*2.50 -im*5.80 im*0.00;
    im*2.50 im*5.00 im*0.00 -im*8.30]

PowerImpedanceACDC.kron(Y, [1,3,4])

# Test set with validation values from Grainger, slightly different due to rounding errors
Y_validation=[ -9.57792*im  4.02597*im   5.55195*im;
               4.02597*im  -5.47532*im   0.649351*im;
               5.55195*im  0.649351*im  -7.0013*im]

@test PowerImpedanceACDC.kron(Y, [1,3,4]) ≈ Y_validation atol=1e-3

