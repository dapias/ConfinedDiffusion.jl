using HDF5
using PyPlot

s1 = h5open("../data/1000particless=1.0.hdf5","r")
s2 = h5open("../data/1000particless=2.0.hdf5","r")
s3 = h5open("../data/1000particless=3.0.hdf5","r")
free = h5open("../data/1000free.hdf5","r")
straight = h5open("../data/1000particlesstraight.hdf5","r")

xs1 = read(s1, "xsquare")
xs2 = read(s2, "xsquare")
xs3 = read(s3, "xsquare")

xm1 = read(s1, "xmean")
xm2 = read(s2, "xmean")
xm3 = read(s3, "xmean")

xsfree = read(free, "xsquare")
xmfree = read(free, "xmean")

xsstraight = read(straight, "xsquare")
xmstraight = read(straight, "xmean")

dt  = read(attrs(free)["dt"]);
nsteps = read(attrs(free)["nsteps"])
nsampling = read(attrs(free)["nsampling"])
t = [dt*nsampling*i for i in 1:nsteps]

fig, ax = plt[:subplots]()

a, = ax[:plot](t, xs1 - xm1.^2, label = L"\sigma = 1.0")
b, = ax[:plot](t, xs2 - xm2.^2,  label = L"\sigma = 2.0")
c, = ax[:plot](t, xs3 - xm3.^2,  label = L"\sigma = 3.0")
#d, = ax[:plot](t, xsfree - xmfree.^2,  label = L"\sigma = free")
e, = ax[:plot](t, xsstraight - xmstraight.^2,  label = L"\sigma = straight")

ax[:set_xlabel](L"t")
ax[:set_ylabel](L"\langle x^2 \rangle_t - \langle x \rangle^2_t")
handles, labels = ax[:get_legend_handles_labels]()
#fig[:legend]((a,b,c,d,e),labels,"upper left")
fig[:legend]((a,b,c,e),labels,"upper left")

plt[:savefig]("../plots/meansquare.png")

