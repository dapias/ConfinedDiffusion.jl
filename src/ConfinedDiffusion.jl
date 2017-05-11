__precompile__(true)

module ConfinedDiffusion

using HDF5

export diffusion, rms, savedata, Parameters

include("cdtypes.jl")
include("diffusion.jl")
include("savedata.jl")

end
