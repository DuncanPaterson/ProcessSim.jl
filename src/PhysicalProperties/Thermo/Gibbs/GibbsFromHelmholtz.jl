
GFromArCubic(T, sumn, B, V, D, δ₁, δ₂, P) = Ar_Cubic(T, sumn, B, V, D, δ₁, δ₂) + P * V - sumn * R * T * (1 + log(P * V / (sumn * R * T)))

function Gibbs(T, P, n, Model::GeneralCubic, Vini = nothing)
    v_min, v_mid, v_max, D, B = Volume(T, P, n, Model)
    sumn = sum(n)
    δ₁ = Model.δ₁
    δ₂ = Model.δ₂
    if !(isnan(v_max))
        if (Vini === nothing)
            return min(GFromArCubic(T, sumn, B, v_min, D, δ₁, δ₂, P), GFromArCubic(T, sumn, B, v_max, D, δ₁, δ₂, P))
        end
        dv1 = abs(v_min / V - 1)
        dv2 = abs(v_max / V - 1)
        if (dv1 < dv2)
            return GFromArCubic(T, sumn, B, v_min, D, δ₁, δ₂, P)
        else
            return GFromArCubic(T, sumn, B, v_max, D, δ₁, δ₂, P)
        end
    end
    return GFromArCubic(T, sumn, B, v_min, D, δ₁, δ₂, P)
end

function Gibbs(T, P, n, Model::HelmModel, Vini = nothing) where{HelmModel <: HelmholtzModel}
end