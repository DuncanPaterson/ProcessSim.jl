
abstract type Cubic_B end

include("Linear.jl")

function BFromCrit(Ω_b, T_c, P_c)
    return Ω_b * T_c / P_c * R
end