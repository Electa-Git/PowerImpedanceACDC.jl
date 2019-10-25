
@testset "determine_impedance" begin
    # @testset "parallel impedance" begin
    #     include("impedance_tests/test_1.jl")
    #     @test isequal(imp[], 1.5s+2)
    # end
    #
    # @testset "cable model" begin
    #     include("impedance_tests/cable_test.jl")
    #     c3, c4 = symbols("c3,c4")
    #     @test isequal(imp[], c3/c4)
    # end
    #
    # @testset "transmission line" begin
    #     include("impedance_tests/tl_thomas.jl")
    #     ABCD = reshape([Sym(string("l",i)) for i in 1:16],4,4)
    #     (a,b,c,d) = (ABCD[1:2,1:2],ABCD[1:2,3:4],ABCD[3:4,1:2], ABCD[3:4,3:4])
    #     @test isequal(simplify.(imp), simplify.(b*inv(d)))
    # end

end
