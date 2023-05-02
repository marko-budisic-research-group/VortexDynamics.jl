using NeuralODETDA
using Makie

idx = 110

px = LinRange(-.2,.2,31)
py = LinRange(-.2,.2,31)
grid = [ ComponentVector(x=x, y=y) for x in px, y in py]


vf(p) = vortex_biot_savart( p, γ, 0, sol.u[idx] )
ω = vorticity(sol.u[idx], γ, 0; ν=1e-3)


Ω = ω.(grid)
v = vf.(grid) 
s = hypot.( getindex.(v,:x), getindex.(v, :y) )

ax = Makie.heatmap(px, py, Ω, colormap=:balance )
Makie.arrows!(ax.axis,px, py, 
    getindex.(v,:x) .* Δ, 
    getindex.(v,:y) .* Δ;
    transparency=true  )
ax