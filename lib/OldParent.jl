#__precompile__()
#=
module EAGOParametricInterval

using IntervalArithmetic
using EAGODomainReduction

type Param_Bisect_Opts
  DAGflag::Bool
  LPflag::Bool
  kmax_main::Int64
  kmax_cntr::Int64
  style::String
  display::String
  ptol::Float64
  etol::Float64
  rtol::Float64
  DAGpass::Int64
  p_rel_bisect::Bool
  DAGh
  DAGg
  DAGgL
  DAGgU
  DAGsym
end
Param_Bisect_Opts() = Param_Bisect_Opts(false, #DAGflag::Bool
                                        false, #LPflag::Bool
                                        1E2, #kmax_main
                                        1E2, #kmax_cntr
                                        "KrawczykCW", #style
                                        "Full", #display
                                        1E-3, #ptol
                                        1E-3, #etol
                                        1E-6, #rtol
                                        5, #DAGpass
                                        false, # p_rel_bisect
                                        [],
                                        [],
                                        [],
                                        [],
                                        [])

include("src/Parametric_Utility.jl")
include("src/Parametric_Test.jl")
include("src/Parametric_Contractor_STD.jl")
#include("src/Parametric_Contractor_Inplace.jl")
include("src/Parametric_Bisection.jl")
include("src/Parametric_Main.jl")

#function __init__()
#end

export Generalized_Param_Bisection, Param_Bisect_Opts,
       GenerateH, setprec, PI_NewtonGS, PI_KrawczykCW
end
=#
