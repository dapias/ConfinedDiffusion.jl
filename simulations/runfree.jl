include("../src/diffusionfree.jl")

using HDF5

xarray = positions[:,1,:]

xsquare = [mean(xarray[i,:].^2) for i in 1:length(xarray[:,1])]
xmean =  [mean(xarray[i,:]) for i in 1:length(xarray[:,1])]

file = h5open("../data/$(nparticles)particlesfree.hdf5", "w")

attrs(file)["dt"] = dt
attrs(file)["Dx"] = Dx
attrs(file)["Dy"] = Dy
attrs(file)["nparticles"] = nparticles
attrs(file)["nsampling"] = nsampling
attrs(file)["nsteps"] = nsteps

file["positions"] = positions
file["xsquare"] = xsquare
file["xmean"] = xmean

close(file)

println("File $(nparticles)particlesfree.hdf5 succesfully generated. See file in ../data/")
