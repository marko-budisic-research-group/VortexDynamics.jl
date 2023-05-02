using NeuralODETDA
using Makie

idx = 200

px = LinRange(-.2,.2,31)
py = LinRange(-.2,.2,31)

D = vortex_biot_savart( (px, py), γ, 0, sol.u[idx]; core=0.005 )

shapeme(v) = reshape( v, (length(px), length(py)) )

ax = Makie.heatmap(px, py, shapeme( hypot.( D.x, D.y ) ) )
Makie.arrows!(ax.axis,px, py, shapeme(D.x * Δ), shapeme(D.y * Δ);transparency=true)
ax