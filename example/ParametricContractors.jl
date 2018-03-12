workspace()
using IntervalArithmetic
using EAGOParametricInterval

# sets options for contractor
opt1 = [100    #Int64: Number of iterations
       1.0E-6 #Float64: Tolerance for equality of
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
newtonGS1 = MC_NewtonGS(Z1,P1,hj1,h1,opt1,Eflag,Iflag,eDflag)
krawczykCW1 = MC_KrawczykCW(Z1,P1,hj1,h1,opt1,Eflag,Iflag)

# Test Problem #2
#=
Z2 = [Interval(-5.0,5.0),Interval(-5.0,5.0),Interval(-5.0,5.0)]
P2 = [Interval(0.6020, 0.7358),Interval(1.2110, 1.4801),Interval(3.6, 4.4)]
a = [37.3692 18.5805 6.25]
c = [0.602 1.211 3.6]
#f(z,p) = sum([()^2] for j=1:3)
h2(z,p) = [1.00E-9*(exp(38*z[1])-1) + p[1]*z[1] - 1.6722*z[2] + 0.6689*z[3] - 8.0267
          1.98E-9*(exp(38*z[2])-1) + 0.6622*z[1] + p[2]*z[2] + 0.6622*z[3] + 4.0535
          1.00E-9*(exp(38*z[3])-1) + z[1]-z[2] + p[3]*z[3] - 6.0]
hj2(z,p) = [(1.00E-9*(38*exp(z[1]))+p[1]) (-1.6722) 0.6689;
            0.6622 (1.98E-9*(38*exp(38*z[2]))+p[2]) 0.6622;
            1 (-1) (1.00E-9*(38*exp(38*z[3]))+p[3])]
H2 = h(Z,P)
HJ2 = hj(Z,P)
newtonGS2 = MC_NewtonGS(Z2,P2,hj2,h2,opt,Eflag,Iflag,eDflag)
krawczykCW2 = MC_KrawczykCW(Z2,P2,hj2,h2,opt,Eflag,Iflag)
=#
