function singletrajectory(nsteps::Int64, nsampling::Int64, dt::Float64, Dx::Float64, Dy::Float64, sigma::Float64, b::Boundary)

    positions = zeros(nsteps,2) 
    temporary = zeros(nsampling,2)

    p = Particle([0.,0.], [0.,0.])
    positions[1,:] = p.r

    for i in 2:nsteps
        for k in 1:nsampling
            if k ==1
                p.rprevious =   positions[i-1,:]
                temporary[k,:] = p.rprevious
            else
                p.rprevious =   temporary[k-1,:]
            end
            p.r += sqrt(2*dt)*[sqrt(Dx)*randn(), sqrt(Dy)*randn()]
            boundary(p, b, sigma)
            temporary[k,:] = p.r
        end
        positions[i,:] = p.r 
    end
    xpositions = positions[:,1]
end


function diffusionsine(nparticles::Int64, nsteps::Int64, nsampling::Int64, dt::Float64, Dx::Float64, Dy::Float64, lambda::Float64, sigma::Float64)

    shape(s) = x->sin(s*x)
    b = Boundary(shape(sigma), lambda)

    xpositions = zeros(nsteps, nparticles)
    xpositions[:,1] = singletrajectory(nsteps, nsampling, dt, Dx, Dy, sigma, b)
    println("Particle 1 done")

    try
        for j in 2:nparticles
            xpositions[:,j] = singleparticle(nsteps, nsampling, dt, Dx, Dy, sigma, b)
            println("Particle $j done")
        end
        return xpositions
    catch y
        if isa(y, InterruptException)
            return xpositions
        end
    end
end


function rms(positions::Array{Float64,2}, nsampling::Int64, dt::Float64)
    xsq = positions.^2
    nsteps = length(positions[:,1])
    nparticles = length(positions[1,:])

    t = 0.0:nsampling*dt:nsampling*dt*(nsteps-1)
    t = collect(t)

    xsquare = [mean(xsq[i,:]) for i in 1:nsteps]
#    xmean = [mean(positions[i,:]) for i in 1:nsteps]

    D = (xsquare)./(2*t)

    return t, D

end


# nsteps = 1000
# nsampling = 1000
# nparticles = 70000
# Dx = Dy = 1.
# dt = 5.e-3
# lambda = 0.1
# sigma = 9.0

# pos = diffusionsine(nparticles, nsteps, nsampling, dt, Dx, Dy, lambda, sigma)
# t,D =  rms(pos, nsampling, dt)
