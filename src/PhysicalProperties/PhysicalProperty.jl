import PhysicalConstants.CODATA2018: N_A, k_B
using Unitful

const R = N_A * k_B / 1u"J" * 1u"K" * 1u"mol" |> NoUnits

include("Thermo/Thermo.jl")