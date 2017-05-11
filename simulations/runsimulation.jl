include("../src/ConfinedDiffusion.jl")
include("parameters.jl")

using ConfinedDiffusion

#filename = "0$(y0)lambda$(lambda)sigma$(sigma)"
filename = "Dx=$(Dx)Dy=$(Dy)-1"
pos = diffusion(nparticles, nsteps, nsampling, dt, Dx, Dy, lambda, shape, L, x0, y0)

t,D =  rms(pos, nsampling, dt)

parameters = Parameters(nsteps, nsampling, nparticles, Dx, Dy, dt, lambda, sigma, L, x0, y0)

savedata(filename, parameters, pos, t, D)
