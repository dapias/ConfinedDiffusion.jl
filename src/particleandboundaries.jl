module Objects
export Particle, boundary1, boundary2

using NLsolve

type Particle{T}
    r::Array{T,1}
    rprevious::Array{T,1}
end

function boundary1(p :: Particle)
    if p.r[2] < 0.0
        p.r[2] *= -1.0
    end
end

function boundary2(p :: Particle; lambda = 0.8, sigma = 0.1)
    function gaussian(x)
        1.0 + lambda*exp(-x^2/(2*sigma^2))
    end
    
    function slope(particle::Particle)
        (particle.r[2] - particle.rprevious[2])/(particle.r[1] - particle.rprevious[1])
    end

    
    function f!(x, fvec; m = 0.1) 
        fvec[1] = 1. + lambda*exp(-x[1]^2/(2*sigma^2)) - m*(x[1] - p.rprevious[1]) - p.rprevious[2] 
    end

    function normalslope(x)
        sigma^2/(lambda*x)*exp(x^2/(2*sigma^2))
    end

    function normal(x)
        deltax = 1.0
        xpoints = collect(linspace(x-deltax,x+deltax,10))
        ypoints = normalslope(x)*(xpoints - x) + gaussian(x)
        return xpoints, ypoints
    end

    function unitnormalvector(x)
        xn,yn = normal(x)
        g =  gaussian(x)
        if yn[end] < g
            n = norm([xn[end] - x, yn[end] - g] )
            xunit, yunit = x + (xn[end] - x)/n  , g + (yn[end] - g)/n         else
            n = norm([xn[1] - x, yn[1] - g] )
            xunit, yunit = (xn[1] -x)/n +x , (yn[1] -g)/n + g
        end
        [xunit, yunit]
    end

    function newvector(particle::Particle, intersection::Vector{Float64}, n::Vector{Float64})
        particle.r - 2*(dot((particle.r - intersection), n))*n
    end 

    if p.r[2] > gaussian(p.r[1])
        if abs(p.r[1]) .> 5*sigma
            y = 1.0 - (p.r[2] - 1.0)
            p.r[2] = y
        else
            slo = slope(p)
            zeroes(x, fvec) = f!(x, fvec, m = slo)
            sol = nlsolve(zeroes, [p.r[1]])
            xsol = sol.zero[1]
            ysol = gaussian(xsol)
            unitario = unitnormalvector(xsol)
            p.r= newvector(p,[xsol, ysol], [unitario[1] - xsol, unitario[2] - ysol])
        end
    end
end

end
