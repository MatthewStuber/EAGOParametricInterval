#workspace()
using BenchmarkTools
using IntervalArithmetic
using EAGOParametricInterval

# sets options for contractor
opt1 = Any[10    #Int64: Number of iterations
       1.0E-30 #Float64: Tolerance for equality of
       1.0E-6 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]
Eflag = false
Iflag = false
eDflag = false

# Test Problem #1
P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4;
           z[1]+p[2]*z[2]]
hj1(z,p) = [(2*z[1]+p[1]) (2*z[2]);
              1              p[2]]

newtonGS1 = PI_NewtonGS(Z1,P1,hj1,h1,opt1,Eflag,Iflag,eDflag)
krawczykCW1 = PI_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

function h1!(hh,z,p)
       hh[1] = z[1]^2+z[2]^2+p[1]*z[1]+4
       hh[2] = z[1]+p[2]*z[2]
end
function hj1!(hh,z,p)
       hh[1,1] = 2*z[1]+p[1]
       hh[1,2] = 2*z[2]
       hh[2,1] = 1
       hh[2,2] = p[2]
end

P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
InnewtonGS1 = PIn_NewtonGS(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag,eDflag)
InkrawczykCW1 = PIn_KrawczykCW(Z1,P1,hj1!,h1!,opt1,Eflag,Iflag)
