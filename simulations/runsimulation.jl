include("../src/ConfinedDiffusion.jl")
include("parameters.jl")

using ConfinedDiffusion

filename = "pruebasigma$(sigma)lambda$(lambda)"

pos = diffusion(nparticles, nsteps, nsampling, dt, Dx, Dy, lambda, sigma, shape)

t,D =  rms(pos, nsampling, dt)

parameters = Parameters(nsteps, nsampling, nparticles, Dx, Dy, dt, lambda, sigma)

savedata(filename, parameters, pos, t, D)
