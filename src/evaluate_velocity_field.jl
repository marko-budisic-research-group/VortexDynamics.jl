using VortexDynamics
using CairoMakie

idx = 110

px = LinRange(-.2,.2,41)
py = LinRange(-.2,.2,41)
grid = [ ComponentVector(x=x, y=y) for x in px, y in py]

px_fine = LinRange(-.2,.2,101)
py_fine = LinRange(-.2,.2,101)
grid_fine = [ ComponentVector(x=x, y=y) for x in px_fine, y in py_fine]

# define functions of a single parameter that evaluate the velocity field
v = velocity( γ, 0, sol.u[idx] )
ω = vorticity_burgers(γ, 0, sol.u[idx]; ν=1e-3)


Ω = ω.(grid_fine)
v = v.(grid) 
s = hypot.( getindex.(v,:x), getindex.(v, :y) )

ax = Makie.heatmap(px_fine, py_fine, Ω, colormap=:balance )
Makie.arrows!(ax.axis,px, py, 
    getindex.(v,:x) .* Δ, 
    getindex.(v,:y) .* Δ;
    transparency=true  )
ax

function 