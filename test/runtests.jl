#!/usr/bin/env julia

using EAGOParametricInterval

# write your own tests here
println("Testing ParametricContractor...")
t = @elapsed include("ParametricContractor.jl")
println("done (took $t seconds).")
