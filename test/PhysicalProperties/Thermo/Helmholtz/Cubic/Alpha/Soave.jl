
@testset "AlphaSoave" begin
    @test ProcessSim.SoaveAlpham(0.0) ≈ 0.48
    @test ProcessSim.PRAlpham(0.0) ≈ 0.17464

    tc = [120.0, 150.0]
    c = [ProcessSim.SoaveAlpham(0.1) ProcessSim.SoaveAlpham(0.2)]
    params = ProcessSim.SoaveAlphaPolynomial(c, tc)
    
    alpha = ProcessSim.Alpha(120.0, params)
    @test alpha[1] ≈ 1.0
    alpha = ProcessSim.Alpha(150.0, params)
    @test alpha[2] ≈ 1.0
end