include("../src/diffusionstraight.jl")

using HDF5

xarray = positions[:,1,:]
xsquare = [mean(xarray[i,:].^2) for i in 1:length(xarray[:,1])]
xmean =  [mean(xarray[i,:]) for i in 1:length(xarray[:,1])]

file = h5open("../data/$(nparticles)particlesstraight.hdf5", "w")

attrs(file)["dt"] = dt
attrs(file)["Dx"] = Dx
attrs(file)["Dy"] = Dy
attrs(file)["nsampling"] = nsampling
attrs(file)["nsteps"] = nsteps
attrs(file)["nparticles"] = nparticles

file["positions"] = positions
file["xsquare"] = xsquare
file["xmean"] = xmean

close(file)

println("File $(nparticles)particlesstraight.hdf5 succesfully generated. See file in ../data/")
