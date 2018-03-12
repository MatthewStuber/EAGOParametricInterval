module EAGOParametricInterval

#using ForwardDiff
#using ReverseDiff
using Calculus
using ValidatedNumerics
#using EAGODAGContractor

type Param_Bisect_Opts
  DAGflag::Bool
  LPflag::Bool
  kmax_main
  kmax_cntr
  style
  display
  ptol
  etol
  rtol
  DAGpass
  p_rel_bisect
  DAGh
  DAGg
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
                                        [])

include("lib/Parametric_Utility.jl")
include("lib/Parametric_Test.jl")
include("lib/Parametric_Contractor.jl")
include("lib/Parametric_Bisection.jl")
include("lib/Parametric_Main.jl")

#include("lib/utility.jl")
#include("lib/param_iter.jl")
#include("lib/param_test.jl")
#include("lib/param_bisect.jl")

# export Param_Intv_Contactor, Param_Intv_NewtonGS, Param_Intv_Newton,
#        Param_Intv_KrawczykCW, Param_Intv_Krawczyk, GenerateJacobianX,
#        Preconditioner, Exclusion_Test, Inclusion_Test,
#        Param_Intv_ContactorTest, GenerateH, Generalized_Param_Bisection

export Generalized_Param_Bisection, Param_Bisect_Opts, GenerateJacobianX,
       GenerateH, setprec, Miranda, MirandaExc, partialIncTop, partialIncBot,
       Strict_XinY, isEqual, extDivide, extProcess, MC_NewtonGS, MC_KrawczykCW

end
