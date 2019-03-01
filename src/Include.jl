# - LOAD SYSTEM PACKAGES ------------------------------------------- #
using MAT
using RowEchelon
using Sobol
using LinearAlgebra
using SparseArrays
using Distributions
using JSON
using GLPK
using Gurobi
using DelimitedFiles
using Optim
using Logging
using ProgressMeter
# ------------------------------------------------------------------ #

# constants -
const path_to_package = dirname(pathof(@__MODULE__))
const glucose_uptake_index = 202   # HEX1

# local code -
include("./flux/Flux.jl")
include("./flux/Data.jl")
include("./flux/Sample.jl")
include("./flux/Utility.jl")
include("./flux/Types.jl")
include("./flux/Export.jl")
include("./flux/Rules.jl")
include("./flux/Solve.jl")
include("./flux/MOMA.jl")
