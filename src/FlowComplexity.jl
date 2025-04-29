module FlowComplexity

#==============================================================================#
# IMPORTS
#==============================================================================#

include("ArgParse.jl")
include("pipeline/Params.jl")
include("pipeline/SLURM.jl")
include("core/Chemostats.jl")
include("core/Ensemble.jl")
include("core/Simulation.jl")
include("core/EvolveStochastic.jl")
include("core/EvolveTauleaping.jl")
include("pipeline/Extract.jl")
include("pipeline/Analysis.jl")
include("Helper.jl")
# (order matters!)

parse_arguments()

end

#==============================================================================#
# END OF FILE
#==============================================================================#