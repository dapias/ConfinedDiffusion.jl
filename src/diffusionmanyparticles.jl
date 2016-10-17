include("./particleandboundaries.jl")

nsteps = 100000
Dx = 1.
Dy = 1.
l = 0.8
s = 3.0
nparticles = 100


positions = zeros(nsteps,2, nparticles)
dt = 0.005

for j in 1:nparticles
    p = Particle([0.,0.], [0.,0.])
    positions[1,:,j] = p.r

    for i in 2:nsteps
        p.rprevious =   positions[i-1,:,j]
        p.r += sqrt(2*dt)*[Dx*randn(), Dy*randn()]
      #  boundary1(p)
      #  boundary2(p, lambda = l, sigma = s)
        positions[i,:,j] = p.r 
    end
    println("Particle $j done")
end


