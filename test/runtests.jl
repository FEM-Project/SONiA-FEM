include("..//src//FemElements.jl")
using .FemElements
using Test

# using SONiA.FemElements

@testset "FemElements.jl" begin
    # Tests for EM1.
    @test EM1_N(-1.0) == [1.0, 0.0]
    @test EM1_N(1.0) == [0.0, 1.0]
    @test EM1_dN(0.0) == [-0.5, 0.5]

    # Tests for EM2.
    @test EM2_N(-1.0) == [1.0, 0.0, 0.0]
    @test EM2_N(0.0) == [0.0, 1.0, 0.0]
    @test EM2_N(1.0) == [0.0, 0.0, 1.0]
    @test EM2_dN(0.0) == [-0.5, 0.0, 0.5]

    # Tests for TL1.
    @test TL1_N(0.0, 0.0) == [1.0, 0.0, 0.0]
    @test TL1_N(1.0, 0.0) == [0.0, 1.0, 0.0]
    @test TL1_N(0.0, 1.0) == [0.0, 0.0, 1.0]
    @test TL1_dN(1/3, 1/3) == ([-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0])

    # Tests for QL1.
    @test QL1_N(-1.0, -1.0) == [1.0, 0.0, 0.0, 0.0]
    @test QL1_N(1.0, -1.0) == [0.0, 1.0, 0.0, 0.0]
    @test QL1_N(1.0, 1.0) == [0.0, 0.0, 1.0, 0.0]
    @test QL1_N(-1.0, 1.0) == [0.0, 0.0, 0.0, 1.0]
    @test QL1_dN(0.0, 0.0) == ([-0.25, 0.25, 0.25, -0.25], [-0.25, -0.25, 0.25, 0.25])
end
