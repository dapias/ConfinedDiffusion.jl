__precompile__(true)

module ConfinedDiffusion

using ForwardDiff
import ForwardDiff.derivative
using Roots
using HDF5

export diffusionsine, rms, savedata

include("cdtypes.jl")
include("cdmethods.jl")
include("diffusion.jl")
include("savedata.jl")

end
