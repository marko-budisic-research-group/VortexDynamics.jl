using NeuralODETDA
using Makie

idx = 200

px = LinRange(-.2,.2,31)
py = LinRange(-.2,.2,31)

p = [ ComponentVector(x=x, y=y) for x in u[1], y in u[2]]
p = ComponentArray( 
    x= vec( map( pp -> pp.x, p) ),  
    y= vec( map( pp -> pp.y, p) ) )


D2 = vortex_biot_savart( (px, py), γ, 0, sol.u[idx]; core=0.005 )

arrows(px, py', D.x .* Δ./100, D.y .* Δ./100)