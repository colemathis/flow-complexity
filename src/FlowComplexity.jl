module FlowComplexity

#==============================================================================#
# IMPORTS
#==============================================================================#

# general & pipeline imports
include("ArgParse.jl")
include("pipeline/Params.jl")
include("Helper.jl")
include("pipeline/IO.jl")
include("pipeline/Analysis.jl")
include("pipeline/SLURM.jl")

# core imports
include("core/Chemostats.jl")
include("core/Ensemble.jl")
include("core/Simulation.jl")
include("core/EvolveStochastic.jl")
include("core/EvolveTauleaping.jl")

# (order matters!)

end

#==============================================================================#
# END OF FILE
#==============================================================================#