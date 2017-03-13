"""Function that computes the difference of the cell number associated to the particle by comparing its current x-position with the previous one"""
function cellchange(p :: Particle, sigma::Float64)
    """It returns a number of cell based on a x-position"""
    function cell(r::Array{Float64,1})
        if -pi/(2*sigma) < r[1] < 3*pi/(2*sigma)
            cellnumber = 0.0
        elseif r[1] > 3*pi/(2*sigma)
            a = r[1] - 3*pi/(2*sigma)
            cellnumber = div(a,(2*pi/sigma)) + 1.0  #2pi/sigma is the period of the boundary
        else 
            a = r[1] + pi/(2*sigma)
            cellnumber = div(a, (2*pi/sigma)) - 1.0
        end
        cellnumber
    end

    prevcell = cell(p.rprevious)
    currentcell = cell(p.r)

    return abs(currentcell - prevcell)

end

function incell(p::Particle, b::Boundary)
    if 0.0 < p.r[2] < b.shape(p.r[1])
        return true
    else
        return false
    end
end

function boundary(p :: Particle, b::Boundary)
    if  !incell(p, b)
        p.r = p.rprevious
    end
end    

# function boundary(p :: Particle, b::Boundary, sigma::Float64, L::Float64)
#     if  incell(p, b)
#         if (p.rprevious[2] > L - b.lambda) && (p.r[2] > L - b.lambda)  #If it crosses the wall with a tunneling-like effect
#             if cellchange(p, sigma) >= 1
#                 p.r = p.rprevious  #Rejection method
#             end
#         end
#     else
#         p.r = p.rprevious
#     end
# end    

