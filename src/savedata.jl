function writeattributes(file, p::Parameters)
    attrs(file)["nsampling"] = p.nsampling
    attrs(file)["nparticles"] = p.nparticles
    attrs(file)["nsteps"] = p.nsteps
    attrs(file)["Dx"] = p.Dx
    attrs(file)["Dy"] = p.Dy
    attrs(file)["dt"] = p.dt
    attrs(file)["lambda"] = p.lambda
    attrs(file)["sigma"] = p.sigma
    attrs(file)["L"] = p.L
    attrs(file)["x0"] = p.x0
    attrs(file)["y0"] = p.y0
end

function savedata(filename::String, p::Parameters,positions::Array, t::Array, D::Array)
    try
        mkdir("../data/")
    end
    
    file = h5open("../data/$(filename).hdf5", "w")
    writeattributes(file,p)

    file["positions"] = positions
    d = hcat(t,D)
    file["D"] = d

    attrs(file)["Deff"] = mean(D[end-100:end])

    close(file)

    println("Filename ../data/$(filename).hdf5 created")

   

end

