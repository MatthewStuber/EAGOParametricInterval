# EAGOParametricInterval
A library for bounding functions via parametric interval methods

## Installation

```julia
julia> Pkg.add("EAGOParametricInterval")
```

## Capabilities

**EAGOParametricInterval.jl** provides methods for performing parametric interval calculations such as (Parametric Interval Newton/Krawczyk) as well as a series of tests to verify the (non)existence of unique enclosed functions. This routine
are used extensively in the `EAGO.jl`(https://github.com/MatthewStuber/EAGO.jl) package solver

Currently, it exports four types of contractor functions for use. This contractors
include embedded test for the guaranteed existence or non-existence. Please see the example.jl file for usage cases. 

## Future Work

* A parametric bisection routine will be updated that can divide the `(X,P)` space
into a a series of boxes that all contain unique branches of the implicit function
`p->y(p)`.
* Minor modifications to the contractors are planned to improve computational
performance. Namely, row handling for the sparse LU factorization for interval
midpoint computation.

## Related Packages
- **EAGO.jl**(https://github.com/MatthewStuber/EAGO.jl): A package containing global and robust solvers based mainly on McCormick relaxations.
This package supports a JuMP and MathProgBase interface.
- **IntervalRootFinding.jl**(https://github.com/JuliaIntervals/IntervalRootFinding.jl): Provides root finding routines using Interval Newton
and Krawczyk methods but don't include parametric method, methods with embedded test,
or handling of large sparse matrices for preconditioner calculation.

## References
- E. R. Hansen and G. W. Walster. Global Optimization Using Interval Analysis. Marcel Dekker, New York, second edition, 2004.
- R. Krawczyk. Newton-algorithmen zur bestimmung con nullstellen mit fehler-schranken. Computing, 4:187–201, 1969.
- R. Krawczyk. Interval iterations for including a set of solutions. Computing, 32:13–31, 1984.
- C. Miranda. Un’osservatione su un teorema di brower. Boll. Un. Mat. Ital., 3:5–7, 1940.
- A. Neumaier. Interval Methods for Systems of Equations. Cambridge University Press, Cambridge, 1990.
- R. E. Moore. A test for existence of solutions to nonlinear systems. SIAM Journal on Numerical Analysis, 14(4):611–615, 1977.
