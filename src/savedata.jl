function writeattributes(file, p::Parameters)
    attrs(file)["nsampling"] = p.nsampling
    attrs(file)["nparticles"] = p.nparticles
    attrs(file)["nsteps"] = p.nsteps
    attrs(file)["Dx"] = p.Dx
    attrs(file)["Dy"] = p.Dy
    attrs(file)["dt"] = p.dt
    attrs(file)["lambda"] = p.lambda
    attrs(file)["sigma"] = p.sigma
end

function savedata(filename::String, p::Parameters,positions::Array, t::Array, D::Array)
    try
        mkdir("../data/")
    end
    
    file = h5open("../data/$(filename).hdf5", "w")
    writeattributes(file,p)

    file["positions"] = positions
    deff = hcat(t,D)
    file["Deff"] = deff

    close(file)

    println("Filename ../data/$(filename).hdf5 created")

end

