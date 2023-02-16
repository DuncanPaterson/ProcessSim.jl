
abstract type Cubic_D end

include("vdw1f.jl")

function AFromCrit(Ω_a, T_c, P_c)
    return Ω_a * (R * T_c)^2 / P_c
end