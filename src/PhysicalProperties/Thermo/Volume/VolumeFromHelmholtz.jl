struct VolumeFromHelmholtzModel{HelmholtzModel} <: VolumeModel
    helmholtz_model::HelmholtzModel
end

function Volume(T, P, n, Model::HelmholtzModel, Vini)
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
        else if !(F===nothing)
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