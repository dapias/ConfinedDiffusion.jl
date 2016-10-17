using HDF5
using PyPlot

s1 = h5open("../data/100particless=1.0.hdf5","r")
s2 = h5open("../data/100particless=2.0.hdf5","r")
s3 = h5open("../data/100particless=3.0.hdf5","r")
free = h5open("../data/100particlesfree.hdf5","r")

xs1 = read(s1, "xsquare")
xs2 = read(s2, "xsquare")
xs3 = read(s3, "xsquare")

xm1 = read(s1, "xmean")
xm2 = read(s2, "xmean")
xm3 = read(s3, "xmean")

xsfree = read(free, "xsquare")
xmfree = read(free, "xmean")

fig, ax = plt[:subplots]()

a, = ax[:plot](xs1 - xm1.^2, label = L"\sigma = 1.0")
b, = ax[:plot](xs2 - xm2.^2,  label = L"\sigma = 2.0")
c, = ax[:plot](xs3 - xm3.^2,  label = L"\sigma = 3.0")
d, = ax[:plot](xsfree - xmfree.^2,  label = L"\sigma = free")

ax[:set_xlabel](L"n_{steps}")
ax[:set_ylabel](L"\langle x^2 \rangle_t - \langle x \rangle^2_t")
handles, labels = ax[:get_legend_handles_labels]()
fig[:legend]((a,b,c,d),labels,"upper left")

plt[:savefig]("../plots/meansquare.png")

