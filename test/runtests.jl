using Test
using SpecialFunctions
using PyPlot
using FastGaussQuadrature

include("onTheFlyGH.jl")

@testset "Gauss-Hermite on the fly" begin
    compare()
end

include("GaussLaguerrePlots.jl")
