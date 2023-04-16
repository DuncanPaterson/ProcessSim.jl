function Volume(T, P, n, Model::HelmModel, Vini) where{HelmModel <: HelmholtzModel}
    function HelmFromV(V)
        return Helmholtz(T, V, n, Model.helmholtz_model)
    end

    nRT = sum(n) * R * T
    result = DiffResults.HessianResult([Vini])

    function fJ!(F, J, V)
        if !(J === nothing)
            result = ForwardDiff.hessian!(result, HelmFromV, V)
            Pcalc = nRT / V[1] - DiffResults.gradient(result)
            if !(F === nothing)
                F[1] = P - Pcalc
            end
            J[1,1] = -nRT / V[1]^2 - DiffResults.Hessian(result)
        elseif !(F===nothing)
            result = ForwardDiff.hessian!(result, HelmFromV, V)
            Pcalc = nRT / V[1] - DiffResults.gradient(result)
            F[1] = P - Pcalc
        end
    end

    solver_results = nlsolve(only_fj!(fJ), [Vini])
    if !(converged(solver_results))
        throw(DomainError(Vini, "Initial value for V did not generate a converged result"))
    end
    return solver_results.zero[1]
end

function Volume(T, P, n, Model::GeneralCubic, Vini)
    sumn = sum(n)
    D = CubicD(T, n, Model.D_model)
    B = CubicB(n, Model.B_model)

    a = D / (sumn * sumn)
    b = B / sumn

    #P = RT / (v - b) - a / ((v+δ₁b)(v + δ₂b))
    #P(v - b)(v+δ₁b)(v + δ₂b) - RT (v+δ₁b)(v + δ₂b) + a(v - b) = 0
    #P(v-b)(v² + vδ₁b + vδ₂b + δ₁δ₂b²) - RTv² - RTvδ₁b - RTvδ₂b - RTδ₁δ₂b² + av - ab = 0
    #v³P + v²δ₁bP + v²δ₂bP + vδ₁δ₂b²P - v²bP - vδ₁b²P - vδ₂b²P - δ₁δ₂b³P - RTv² - vRTδ₁b - vRTδ₂b - RTδ₁δ₂b² + av - ab = 0
    #v³P + v²(δ₁bP + δ₂bP - bP - RT) + v(δ₁δ₂b²P - δ₁b²P - δ₂b²P + a) + (-δ₁δ₂b³P - RTδ₁δ₂b² - ab) = 0
    #v³ + v²(δ₁b + δ₂b - b - RT/P) + v(δ₁δ₂b² - δ₁b² - δ₂b² + a/P) + (-δ₁δ₂b³ - RTδ₁δ₂b²/P - ab/P) = 0

    δ₁ = Model.δ₁
    δ₂ = Model.δ₂
    c1 = b * (δ₁ + δ₂ - 1) - RT/P
    c2 = b² * (δ₁δ₂ - δ₁ - δ₂) + a/P
    c3 = -δ₁δ₂b³ - RTδ₁δ₂b²/P - ab/P

    return sumn .* CubicRoots(c1, c2, c3), D, B
end