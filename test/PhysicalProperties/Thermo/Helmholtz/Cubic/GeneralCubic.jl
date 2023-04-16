
struct MockD end

struct MockB end

function ProcessSim.CubicD(n, T, Model::MockD)
    return 1.0
end

function ProcessSim.CubicB(n, Model::MockB)
    return 1.0
end


@testset "General Cubic" begin
    params = ProcessSim.GeneralCubic(1.0, 0.0, MockD(), MockB())

    T = 100.0
    V = 2.0
    n = [1.0]
    @test ProcessSim.Helmholtz(T, V, n, params) ≈ -8.314462618 * T * log(0.5) - log(1.5) 
end

@testset "General cubic with vdw1f known results" begin
    tc = [190.564, 305.32] # methane, ethane
    pc = [4.599e6, 4.872e6] 
    ω = [0.0115, 0.0995]
    kij = [0.0 0.0; 0.0 0.0]
    alpha_params = ProcessSim.SoaveAlpham(ω)
    BFromCrit(Ω_b, T_c, P_c)
end