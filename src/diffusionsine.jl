#include("./particleandboundaries.jl")

using ForwardDiff
import ForwardDiff.derivative
using Roots
#using Objects

type Particle{T}
    r::Array{T,1}
    rprevious::Array{T,1}
end

type Boundary
    shape::Function
    tangent::Function
    normal::Function
    lambda::Float64
    
    function normal(b::Function, x)
        f = y -> derivative(b,y)
        normalization = (1 + f(x)^2)^(1/2)
        [f(x)/normalization, -1/normalization]
    end
    
    function tangent(b::Function, x)
        f = y -> derivative(b,y)
        normalization = (1 + f(x)^2)^(1/2)
        [1/normalization, f(x)/normalization]
    end
    
    function Boundary(shape, lambda)
        s(x) = 1. + lambda*shape(x)
        n(x) = normal(s, x)
        t(x) = tangent(s, x)
        
        new(s, t, n, lambda)
    end
    
end

function cellchange(p :: Particle, sigma::Float64)

    function cell(r::Array{Float64,1})
        if -pi/(2*sigma) < r[1] < 3*pi/(2*sigma)
            cellnumber = 0
        elseif r[1] > 3*pi/(2*sigma)
            a = r[1] - 3*pi/(2*sigma)
            cellnumber = div(a,(2*pi/sigma)) + 1  #2pi/sigma is the period of the wall
        else 
            a = r[1] + pi/(2*sigma)
            cellnumber = div(a, (2*pi/sigma)) - 1
        end
        cellnumber
    end

    prevcell = cell(p.rprevious)
    currentcell = cell(p.r)

    return abs(currentcell - prevcell)

end

function incell(p::Particle, b::Boundary)
    if p.r[2] < b.shape(p.r[1])
        return true
    else
        return false
    end
end


function boundarysine(p :: Particle, b::Boundary, sigma::Float64)
    if  incell(p, b)
        if (p.rprevious[2] > 1.-b.lambda) && (p.r[2] > 1.- b.lambda)  #If it crosses the wall with a tunneling-like effect
            if cellchange(p, sigma) >= 1
                p.r = p.rprevious  #Rejection method
            end
        end
    else
        #        reflectingboundaries!(p, b)
        p.r = p.rprevious
    end
end    


function boundary1(p :: Particle, b::Boundary)
    if p.r[2] < 0.0
        p.r[2] *= -1.0
    end

    if !incell(p,b)
        p.r = p.rprevious
    end


end

function diffusionsine(nparticles::Int64, nsteps::Int64, nsampling::Int64, dt::Float64, Dx::Float64, Dy::Float64, lambda::Float64, sigma::Float64)

    shape(s) = x->sin(s*x)
    b = Boundary(shape(sigma), lambda)

    xpositions = zeros(nsteps, nparticles)
    xpositions[:,1] = singleparticle(nsteps, nsampling, dt, Dx, Dy, sigma, b)
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

function singleparticle(nsteps::Int64, nsampling::Int64, dt::Float64, Dx::Float64, Dy::Float64, sigma::Float64, b::Boundary)

    positions = zeros(nsteps,2) #I am just interested in x-positions
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
            boundarysine(p, b, sigma)
            boundary1(p, b)
            temporary[k,:] = p.r
        end
        positions[i,:] = p.r 
    end
    xpositions = positions[:,1]
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

function straightlineparticle(p::Particle)
    f(x) = (p.r[2] - p.rprevious[2])/(p.r[1] - p.rprevious[1])*(x - p.r[1]) + p.r[2]
end

function reflectingboundaries!(p::Particle, b::Boundary)
    intersection(x) = b.shape(x) - straightlineparticle(p)(x)
    xint = fzero(intersection,[p.r[1], p.rprevious[1]])
    yint = b.shape(xint)
    rint = p.r - [xint, yint]
    tanvector = dot(rint, b.tangent(xint))*b.tangent(xint)
    norvector = -dot(rint, b.normal(xint))*b.normal(xint) ##Minus to change the sign of the normal vector
#    p.rprevious = [xint, yint]  #This was thinking when multiple reflections were allowed
    p.r = [xint, yint] + tanvector + norvector
end

nsteps = 1000
nsampling = 1000
nparticles = 70000
Dx = Dy = 1.
dt = 5.e-3
lambda = 0.1
sigma = 1.0

pos = diffusionsine(nparticles, nsteps, nsampling, dt, Dx, Dy, lambda, sigma)
t,D =  rms(pos, nsampling, dt)
