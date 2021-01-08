__precompile__()

module OrbitTomographyAddOns

# Import standard packages
using JLD2
using FileIO
using Interpolations
using ProgressMeter
using Distributions

# Import special packages
using Equilibrium
using GuidingCenterOrbits
using OrbitTomography

greet() = print("Hello World!")

include("add_ons.jl")
export add_noise, get_orbel_volume

end # module
