# This is an PowerImpedanceACDC element in order to model a transmission line based on 
# line parameters (Z and Y) or ABCD matrix.
# It is a blackbox model, meaning that the internal structure of the line is not modeled.
# The line is represented by its imported frequency-dependent parameters, e.g. line costants (Z,Y)
# or the ABCD two-port parameters. The parameters are in total length.
#
export blackbox_line

@with_kw mutable struct Blackbox_line <: Transmission_line 

    n :: Int = 3 # Number of phases
    data_type :: Symbol = :ABCD # :ABCD or :ZY # Type of data provided by the user
    path :: String = "" # Path to the file containing the data
    Y :: Array{Basic} = [] # Imported shunt admittance matrix [Ohms]
    Z :: Array{Basic} = [] # Imported series impedance matrix [Siemens]
    ABCD :: Array{Basic} = [] # Imported ABCD matrix
    f:: Array{Basic} = [] # Frequency vector for the imported data [Hz]

end


function blackbox_line(;args...)
    
    bbl = Blackbox_line()
    transformation = false #variable transformation is defined false (park transformation not operated)
    for (key, val) in pairs(args)
        if in(key, propertynames(bbl))
            setfield!(bbl, key, val) 
        elseif (key == :transformation)
            transformation = val
        else
            throw(ArgumentError("Unknown power line property name.")) #If no one of the value specified above -> display an error
        end
    end
    # Importing happening here, read data from file and store it in the struct
    

    # f =
    # Y=[zeros(ComplexF64,nodes, nodes) for _ in 1:f]
    # Z=[zeros(ComplexF64,nodes, nodes) for _ in 1:f]
    
    # First read in frequency vector, then read in the matrices
    #bbl.f =f

    if bbl.data_type == :ABCD
    
        #bbl.ABCD =
    
    elseif bbl.data_type == :ZY

        #bbl.Y =
        #bbl.Z =
    else 

        throw(ArgumentError("Unknown data type. Use :ABCD or :ZY"))

    end 


    elem = Element(input_pins = bbl.n, output_pins = bbl.n, element_value = bbl,
        transformation = transformation)
    
end


# Never called externally, only used inside the element!
function eval_parameters(bbl :: Blackbox_line, s :: Complex) #function eval_parameters with in input two variables: blackbox_line c and s variable


    # Take the imported Y and Z values 
    if s/(2pi*1im) ∈ freq

        Z=bbl.Z[findfirst(==(s/(2pi*1im)), bbl.f)] #Find the index of the frequency point in the frequency vector f that is equal to s/(2pi*1im) and use this index to get the corresponding Z matrix from the list of Z matrices
        Y=bbl.Y[findfirst(==(s/(2pi*1im)), bbl.f)] 
    else # Inter- or extrapolation needed
        #TODO: Implement in a more efficient way, as this creates a new interpolation object every time the function is called.
        Z_interp = linear_interpolation(bbl.f, bbl.Z,extrapolation_bc=Line())
        Y_interp = linear_interpolation(bbl.f, bbl.Y,extrapolation_bc=Line())
        Z=Z_interp(real(s/(2pi*1im))) #
        Y=Y_interp(real(s/(2pi*1im))) # 

    end

    return (Z,Y) #This function returns the value of series impedance Z and shunt admittance Y.

end


# Only function that is called externally, the interface with the rest of the code!
# Two cases: Either the user provides Z and Y, or the user provides the ABCD matrix directly.
function eval_abcd(bbl :: Blackbox_line, s :: Complex)
    
    
    if data_type == :ZY
    # Case 1: User provides Z and Y directly, calculation of the ABCD matrix as in a normal transmission line!
    (Z, Y) = eval_parameters(bbl, s) #function eval_parameters receives in input variables c and s and gives in output the series impedance matrix Z and the shunt admittance matrix Y
    γ = sqrt(convert(Array{Complex}, Z*Y)) # conversion of arrays product Z*Y and then Γ=sqrt(Z*Y) -> Simulator tutorial pag 19, bottom page (Calculations same as transmission line case)
    # γ = sqrt(Z*Y) #TODO: This line is added to replace the one on the top, which was giving problems for single DC cables (no transformation). Need to see if this will cause any issue with other components, but so far so good.
    Yc = Z \ γ #Always ottom page 19 simulator_tutorial
    n = Int(size(Yc,1)) #Saves in n the number of rows present in the vector Yc
    abcd = zeros(Complex, 2n, 2n) #Creation of the ABCD matrix-> zeros matrix with dimension 2n*2n

    abcd[1:n, 1:n] = cosh(γ*1) #Eq 33 pag 19 simulator_tutorial. Matrix A[n*n]=cosh(Γl)
    abcd[1:n,n+1:end] = Yc \ sinh(γ*1) #Eq 33 pag 19 simulator_tutorial. Matrix B[n*n]=Yc^-1sinh(Γl)
    abcd[n+1:end,1:n] = Yc * sinh(γ*1) #Eq 33 pag 19 simulator_tutorial. Matrix C[n*n]=Ycsinh(Γl)
    abcd[n+1:end, n+1:end] = cosh(γ*1) #Eq 33 pag 19 simulator_tutorial. Matrix D[n*n]=cosh(Γl)

    # Case 2: User provides the ABCD matrix directly, no calculation needed!

    else

    # Take the imported ABCD matrix 
        if s/(2pi*1im) ∈ freq

            Z=bbl.ABCD[findfirst(==(s/(2pi*1im)), bbl.f)] #Find the index of the frequency point in the frequency vector f that is equal to s/(2pi*1im) and use this index to get the corresponding ABDC matrix from the list of ABCD matrices

        else # Inter- or extrapolation needed
            #TODO: Implement in a more efficient way, as this creates a new interpolation object every time the function is called.
            ABCD_interp = linear_interpolation(bbl.f, bbl.ABCD,extrapolation_bc=Line())
            abcd=ABCD_interp(real(s/(2pi*1im))) #

        end



    end

    return abcd

end