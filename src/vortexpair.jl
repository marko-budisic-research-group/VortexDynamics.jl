using Random

using NeuralODETDA
using ComponentArrays

# 1. Define initial configuration of vortices.
nVortices = 10
u0 = ComponentArray(
    x = [0.0; 0.0], 
    y = 0.75 .+ [-1.0; 1.0] ./ 5,
    )
γ = [1.0; -0.5]

# 2. Simulate the differential equation
using OrdinaryDiffEq

t0 = 0.
tf = 20.
Δ = 1e-1;

diffeq = ODEProblem( vortex_biot_savart!, u0, (t0,tf),γ )
sol = solve( diffeq, Tsit5(); saveat=Δ )

# 3. Evaluate the velocity field on a coarse grid, and vorticity field on a fine grid 
subselect = 1:100
us = @view sol.u[subselect]
ts = @view sol.t[subselect]

# extract horizontal and vertical positions of vortices for each snapshot and merge (hcat) them into a matrix
ux = hcat(getindex.(us,:x)...)
uy = hcat(getindex.(us,:y)...)

# get velocity and vorticity fields
grid_x = LinRange(-1.2,1.2,101)
grid_y = LinRange(-1.2,1.2,101)
vxs, vys, Ωs = getfields( us, ts; px=grid_x, py=grid_y, γ, ν=1e-3 );

using NPZ
npzwrite("vortex-pair-tracks.npz", 
    ux = ux,
    uy = uy,
    t = ts,
    gamma = γ)


npzwrite("vortex-pair-fields.npz", 
    grid_x = collect(grid_x),
    grid_y = collect(grid_y),
    t = ts,
    gamma = γ,
    Vx = cat(vxs...; dims=3),
    Vy = cat(vys...; dims=3),
    Omega = cat(Ωs...; dims=3))


using Makie
ax = Makie.heatmap(grid_x, grid_y, Ωs[1], colormap=:balance )
for k = 1:length(γ)
    lines!( getindex.( getindex.(us, :x), k ),
            getindex.( getindex.(us, :y), k ) )
end
ax
