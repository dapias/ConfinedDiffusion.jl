#Plot a single trajectory
include("./diffusionsingletrajectory.jl")
function gaussian(x)
    l*exp(-x^2/(2*s^2))+1.0
end

using PyPlot
N = 10000
x_inf = 0.0
x = linspace(-8.0, 8.0, N)
y1 = ones(N)*x_inf;
y2 = [gaussian(i) for i in x]
plot(x,y1,"g",linewidth=2,x,y2,"g", linewidth=2)
plot(positions[:,1], positions[:,2],".-")
plt[:ylim](x_inf -0.05, maximum(y2) +0.05)

