include("../src/ConfinedDiffusion.jl")
include("parameters.jl")

using ConfinedDiffusion

filename = "prueba"
#filename = "straight100"

pos = diffusion(nparticles, nsteps, nsampling, dt, Dx, Dy, lambda, shape, L)

t,D =  rms(pos, nsampling, dt)

parameters = Parameters(nsteps, nsampling, nparticles, Dx, Dy, dt, lambda, sigma, L)

savedata(filename, parameters, pos, t, D)
