include("./particleandboundaries.jl")

p = Particle([0.,0.], [0.,0.])
nsteps = 10000
Dx = 1.
Dy = 1.
l = 0.8
s = 2.0


positions = zeros(nsteps,2)
dt = 0.005

positions[1,:] = p.r

for i in 2:nsteps
    p.rprevious = positions[i-1,:]
    p.r += sqrt(2*dt)*[Dx*randn(), Dy*randn()]
    boundary1(p)
    boundary2(p, lambda = l, sigma = s)
    positions[i,:] = p.r 
end


