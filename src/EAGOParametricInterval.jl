module EAGOParametricInterval

using ForwardDiff
using ReverseDiff
using Calculus
using ValidatedNumerics

include("lib/utility.jl")
include("lib/param_iter.jl")
include("lib/param_test.jl")
#include("lib/param_bisect.jl") to add later

 export Param_Intv_Contactor, Param_Intv_NewtonGS, Param_Intv_Newton,
        Param_Intv_KrawczykCW, Param_Intv_Krawczyk, GenerateJacobianX,
        Preconditioner, Exclusion_Test, Inclusion_Test,
        Param_Intv_ContactorTest, GenerateH, kdel, Expr_SubSymbols

end
