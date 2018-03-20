#!/usr/bin/env julia

using EAGOParametricInterval

# write your own tests here
println("Testing Parametric Contractor...")
t = @elapsed include("ParametricContractor.jl")
println("done (took $t seconds).")

println("Testing Sparse Preconditioner...")
t = @elapsed include("SparseCntr.jl")
println("done (took $t seconds).")
