module ParametricContractor

using Compat
using Compat.Test
using IntervalArithmetic
using EAGOParametricInterval

# sets options for contractor
opt1 = Any[100    #Int64: Number of iterations
       1.0E-30 #Float64: Tolerance for equality of
       1.0E-6 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]
# Test Problem #1
P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2*z[1]+p[1]) (2*z[2]);
              1              p[2]]
Eflag = false
Iflag = false
eDflag = false
newtonGS1 = PI_NewtonGS(Z1,P1,hj1,h1,opt1,Eflag,Iflag,eDflag)
krawczykCW1 = PI_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

@test -1E-4 <= newtonGS1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= newtonGS1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= newtonGS1[1][2].hi - 0.208486 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].lo + 1.04243 <= 1E-4
@test -1E-4 <= krawczykCW1[1][1].hi + 0.492759 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].lo - 0.0473789 <= 1E-4
@test -1E-4 <= krawczykCW1[1][2].hi - 0.208486 <= 1E-4

end
