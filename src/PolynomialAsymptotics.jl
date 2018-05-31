module PolynomialAsymptotics

using FastGaussQuadrature
using SpecialFunctions
using Roots

import FastGaussQuadrature: besselroots

include("laguerre.jl")

include("gauss/gausslaguerre.jl")
include("gauss/gaussfreud.jl")
include("gauss/gaussjacobi.jl")


export asy_gausslaguerre
export asy_gaussfreud
export asy_gaussjacobi

end # module
