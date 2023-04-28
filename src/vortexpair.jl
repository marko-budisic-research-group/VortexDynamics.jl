using Random

using NeuralODETDA
using ComponentArrays

# 1. Define initial configuration of vortices.
nVortices = 10
u0 = ComponentArray(
    x = [0.0; 0.0], 
    y = [-1.0; 1.0] ./ 10,
    )
γ = [1.0; 1.1]

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
npzwrite("vortex-unipole.npz", 
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

