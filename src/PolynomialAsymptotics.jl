module PolynomialAsymptotics

using FastGaussQuadrature

include("laguerre.jl")

include("gauss/gausslaguerre.jl")
include("gauss/gaussfreud.jl")
include("gauss/gaussjacobi.jl")


export asy_gausslaguerre
export asy_gaussfreud
export asy_gaussjacobi

end # module
