
abstract type CubicModel <: HelmholtzModel end

struct GeneralCubic{TReal, DModel, BModel} <: CubicModel
    δ₁::TReal
    δ₂::TReal
    D_model::DModel
    B_model::BModel
end

Ar_Cubic(T, sumn, B, V, D, δ₁, δ₂) = -R * T * sumn * log(1.0 - B / V) - D / (B * (δ₁ -δ₂)) * log((1.0 + δ₁ * B / V) / (1.0 + δ₂ * B / V))

function Helmholtz(T, V, n, Model::GeneralCubic)
    D = CubicD(T, n, Model.D_model)
    B = CubicB(n, Model.B_model)
    δ₁ = Model.δ₁
    δ₂ =  Model.δ₂
    sumn = sum(n)
    return Ar_Cubic(T, sumn, B, V, D, δ₁, δ₂)
end