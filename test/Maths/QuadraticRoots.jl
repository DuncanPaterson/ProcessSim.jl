@testset "QuadraticRoots" begin
    @test ProcessSim.QuadraticRoots(-10.0, 25.0)[1] ≈ 5.0
    @test ProcessSim.QuadraticRoots(-10.0, 25.0)[2] ≈ 5.0
    @test ProcessSim.QuadraticRoots(-6.0, 8.0)[1] ≈ 2.0
    @test ProcessSim.QuadraticRoots(-6.0, 8.0)[2] ≈ 4.0
end