module ParamChk_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGOParametricInterval

# Tests strict inclusion procedure for interval vectors
Y = [Interval(0,5),Interval(0,5),Interval(0,5)]
X1 = [Interval(1,2),Interval(1,2),Interval(1,2)]
X2 = [Interval(1,2),Interval(-10,10),Interval(1,2)]
X3 = [Interval(1,2),Interval(1,2),Interval(0,5)]
flag1 = EAGOParametricInterval.Strict_XinY(X1,Y)
flag2 = EAGOParametricInterval.Strict_XinY(X2,Y)
flag3 = EAGOParametricInterval.Strict_XinY(X3,Y)
@test flag1 == true
@test flag2 == false
@test flag3 == false

# Checks extended division routine
B = Interval(1,2)
C = Interval(3,4)
A1 = Interval(0)
A2 = Interval(0,3)
A3 = Interval(-2,0)
A4 = Interval(-3,2)
ind1,B1,C1 = EAGOParametricInterval.extDivide(A1,B,C)
ind2,B2,C2 = EAGOParametricInterval.extDivide(A2,B,C)
ind3,B3,C3 = EAGOParametricInterval.extDivide(A3,B,C)
ind4,B4,C4 = EAGOParametricInterval.extDivide(A4,B,C)
@test ind1 == 0
@test ind2 == 1
@test ind3 == 2
@test ind4 == 3
@test B1 == Interval(-Inf,Inf)
@test 0.33333 - 1E-4 <= B2.lo <= 0.33333 + 1E-4
@test B2.hi == Inf
@test B3 == Interval(-Inf,-0.5)
@test B4.lo == -Inf
@test -0.33333 - 1E-4 <= B4.hi <= -0.33333 + 1E-4
@test C1 == Interval(-Inf,Inf)
@test C2 == Interval(Inf,Inf)
@test C3 == Interval(-Inf,-Inf)
@test C4 == Interval(0.5,Inf)

#=
N =
X =
Mii =
S1 =
S2 =
B =
rtol = 1E-4
indx1,box1,box2 = extProcess(N,X,Mii,S1,S2,B,rtol)
=#

end
