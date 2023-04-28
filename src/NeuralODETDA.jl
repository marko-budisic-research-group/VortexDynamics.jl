module NeuralODETDA

using Random

include("VortexODE.jl")

export vortex_biot_savart
export vortex_biot_savart!
# 1. Define initial configuration of vortices.
# 2. Simulate the differential equation
# 3. Store the raw trajectory data.
# 4. Store the (mollified) vorticity.
# 5. Compute the PH.
# 6. Store the PH.



end # module NeuralODETDA
