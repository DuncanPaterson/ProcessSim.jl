

struct GibbsFromHelmholtzModel{HelmholtzModel} <: GibbsModel
    helmholtz_model::HelmholtzModel
end

function Gibbs(T, P, n, Model::GibbsFromHelmholtzModel, Vini = 0.0)

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

    if (Vini != 0.0)
        nlsolve(only_fj!(fJ), [Vini])
    else
        V = [nRT / P]
        nlsolve()
    end


end