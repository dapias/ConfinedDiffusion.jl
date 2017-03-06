type Particle{T}
    r::Array{T,1}
    rprevious::Array{T,1}
end

type Boundary
    shape::Function
    lambda::Float64
        
    function Boundary(shape, lambda)
        s(x) = 1. + lambda*shape(x)
        new(s,lambda)
    end
end

type Parameters
    nsteps::Int64
    nsampling::Int64
    nparticles::Int64
    Dx::Float64
    Dy::Float64
    dt::Float64
    lambda::Float64
    sigma::Float64
end

