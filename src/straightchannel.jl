include("./particleandboundaries.jl")

importall Objects

function straightboundaries(p :: Particle)
    if p.r[2] < 0.0
        p.r[2] *= -1.0
  
    elseif p.r[2] > 1.0
        y = 1.0 - (p.r[2] - 1.0)
        p.r[2] = y
    end
end


