using Random

using VortexDynamics
using ComponentArrays

# 1. Define initial configuration of vortices.
nVortices = 10
a = round.(rand(Float64, 3), digits=1)
b = round.(rand(Float64, 3), digits=1)
u0 = ComponentArray(
    x = [a[1]; a[2]; a[3]],
    y = [b[1]; b[2]; b[3]] ./ 10,
    )
γ = [1.0; 1.1; 1.0]

# 2. Simulate the differential equation
using OrdinaryDiffEq

t0 = 0.
tf = 5.
Δ = 1e-3;

diffeq = ODEProblem( vortex_biot_savart!, u0, (t0,tf),γ )
sol = solve( diffeq, Tsit5(); saveat=Δ )


# store the data on the drive
ux = hcat( (s.x for s in sol.u)... )
uy = hcat( (s.y for s in sol.u)... )

using NPZ
using Dates

name = "results/vortex-" * string(today())*"-"*randstring(2) *".npz"
npzwrite( name, 
    ux = ux,
    uy = uy,
    t = sol.t,
    gamma = γ)


using Plots
plot()
for k = 1:length(γ)
    plot!( ux[k, :], uy[k,:] )
end
current()

