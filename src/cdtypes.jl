type Particle{T}
    r::Array{T,1}
    rprevious::Array{T,1}
end

type Boundary
    shape::Function
    lambda::Float64
    L::Float64

    #We describe a periodic non-homogeneous channel of the form L+lambda*periodic_function
    function Boundary(shape, lambda, L)
        s(x) = L + lambda*shape(x)
        new(s,lambda, L)
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
    L::Float64
end

