workspace()
#using EAGOSemiInfinite
#using EAGOGlobalSolver
using IntervalArithmetic
using EAGOBranchBound
#using EAGOSmoothMcCormickGrad
#using StaticArrays


function FDI_invApprox(x::Interval)
  Ap = 0.00005
  g = 9.8
  low = -(x^2)/(2*g*Ap^2)
  high = (x^2)/(2*g*Ap^2)
  if x.lo>0
    return high
  elseif x.hi<0
    return low
  else
    return Interval(0,max(low.hi,high.hi))
  end
end

function FDI_invApprox(x::Float64)
  Ap = 0.00005
  g = 9.8
  if x>=0
    return (x^2)/(2*g*Ap^2)
  else x<0
    return -(x^2)/(2*g*Ap^2)
  end
end

function xbnds(u,p)

         # parameters
         Ap = 0.00005
         g = 9.8
         z = []
         # bound z via McCormick operators
         z1 = (((u[1]/p[3])/Ap)^2)/(2.0*g)
         z3 = (1.0/(2.0*g))*((u[1]+u[2])/(p[4]*Ap))^2
         z2 = (-((u[2]-p[4]*Ap*sqrt(2*g*z3))/(p[5])/Ap)^2)/(2.0*g)
         z4 = (((u[1]*p[1]/p[3])/Ap)^2)/(2.0*g)
         z6 = (1.0/(2.0*g))*((p[1]*u[1]+u[2])/(p[4]*Ap))^2
         z5 = (-((u[2]-p[4]*Ap*sqrt(2*g*z6))/(p[5])/Ap)^2)/(2.0*g)
         z7 = (((u[1]/p[3])/Ap)^2)/(2.0*g)
         z9 = (1.0/(2.0*g))*((u[1]+u[2])/(p[4]*(Ap+3.14159*p[2]^2)))^2
         z8 = FDI_invApprox((u[2]-p[4]*(Ap+3.14159*p[2]^2)*sqrt(2.0*g*z9))/(p[5]))
         push!(z,z1,z2,z3,z4,z5,z6,z7,z8,z9)

         # rotate coordinates back to X & intersects

         x = []
         x1 = min(max(z1-z2+z3,0.0),0.75)
         x2 = min(max(z3,0.0),0.75)
         x3 = min(max(z3-z2,0.0),0.75)
         x4 = min(max(z4-z5+z6,0.0),0.75)
         x5 = min(max(z6,0.0),0.75)
         x6 = min(max(z6-z5,0.0),0.75)
         x7 = min(max(z7-z8+z9,0.0),0.75)
         x8 = min(max(z9,0.0),0.75)
         x9 = min(max(z9-z8,0.0),0.75)
         push!(x,x1,x2,x3,x4,x5,x6,x7,x8,x9)

         return x
end

# sets up SIP constraint with imbedded x(p) calculation
function SIP_g(un,p)
  u = un[1:2]
  xbnd = xbnds(u,p)
  r_xbnds = [zero(xbnd[1]) for i=1:3,j=1:3]
  r_xbnds[:,1] = xbnd[1:3]
  r_xbnds[:,2] = xbnd[4:6]
  r_xbnds[:,3] = xbnd[7:9]
  out = un[3]
  for i=1:3
    for j=1:2
      for k=(j+1):3
        out = out - (r_xbnds[i,j]-r_xbnds[i,k])^2
      end
    end
  end
  return out
end

# sets up SIP objective
function SIP_f(un)
  return -un[3]
end
# sets up SIP
ex_SIP_X = [(1E-5)..(1E-4),(1E-5)..(1E-4),(-10.0)..(10.0)]
ex_SIP_P = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]

#Intv = SIP_g(ex_SIP_X,ex_SIP_P)
#Trial1 = [1E-5,1E-5,0.0]
#Intv1 = SIP_g(Trial1,ex_SIP_P)
#Trial2 = [1E-5,1E-4,0.0]
#Intv2 = SIP_g(Trial2,ex_SIP_P)
#Trial3 = [1E-4,1E-5,0.0]
#Intv3 = SIP_g(Trial3,ex_SIP_P)
#Trial4 = [1E-4,1E-4,0.0]
#Intv4 = SIP_g(Trial4,ex_SIP_P)
#Trial5 = [(1E-5+1E-4)/2,(1E-5+1E-4)/2,0.0]
#Intv5 = SIP_g(Trial5,ex_SIP_P)
#Trial6 = [(1E-5+1E-4)/3,(1E-5+1E-4)/3,0.0]
#Intv6 = SIP_g(Trial6,ex_SIP_P)
#Trial7 = [2*(1E-5+1E-4)/3,2*(1E-5+1E-4)/3,0.0]
#Intv7 = SIP_g(Trial7,ex_SIP_P)
#Trial8 = [(1E-5+1E-4)/2,(1E-5+1E-4)/3,0.0]
#Intv8 = SIP_g(Trial8,ex_SIP_P)
#Trial9 = [(1E-5+1E-4)/2,2*(1E-5+1E-4)/3,0.0]
#Intv9 = SIP_g(Trial9,ex_SIP_P)
#Trial10 = [(1E-5+1E-4)/3,(1E-5+1E-4)/2,0.0]
#Intv10 = SIP_g(Trial10,ex_SIP_P)
Trial11 = [2*(1E-5+1E-4)/3,(1E-5+1E-4)/2,0.0]
Intv11 = SIP_g(Trial11,ex_SIP_P)

ex_SIP_P1a = [0.66,(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P2a = [(0.54)..(0.66),0.005,(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P3a = [(0.54)..(0.66),(0.0005)..(0.005),1.15,(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P4a = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.95,(0.85)..(1.15)]
ex_SIP_P5a = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),1.15]
Intv11a1 = SIP_g(Trial11,ex_SIP_P1a)
Intv11b1 = SIP_g(Trial11,ex_SIP_P2a)
Intv11c1 = SIP_g(Trial11,ex_SIP_P3a)
Intv11d1 = SIP_g(Trial11,ex_SIP_P4a)
Intv11e1 = SIP_g(Trial11,ex_SIP_P5a)

ex_SIP_P1b = [0.54,(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P2b = [(0.54)..(0.66),0.0005,(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P3b = [(0.54)..(0.66),(0.0005)..(0.005),0.85,(0.65)..(0.95),(0.85)..(1.15)]
ex_SIP_P4b = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.65,(0.85)..(1.15)]
ex_SIP_P5b = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),0.85]
Intv11a2 = SIP_g(Trial11,ex_SIP_P1b)
Intv11b2 = SIP_g(Trial11,ex_SIP_P2b)
Intv11c2 = SIP_g(Trial11,ex_SIP_P3b)
Intv11d2 = SIP_g(Trial11,ex_SIP_P4b)
Intv11e2 = SIP_g(Trial11,ex_SIP_P5b)
ex_SIP_P4ca = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.70,(0.85)..(1.15)]
ex_SIP_P4c = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.75,(0.85)..(1.15)]
ex_SIP_P4da = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.80,(0.85)..(1.15)]
ex_SIP_P4d = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),0.85,(0.85)..(1.15)]
Intv11d2a = SIP_g(Trial11,ex_SIP_P4ca)
Intv11d3 = SIP_g(Trial11,ex_SIP_P4c)
Intv11d3a = SIP_g(Trial11,ex_SIP_P4da)
Intv11d4 = SIP_g(Trial11,ex_SIP_P4d)


disc1 = [0.6,0.0005,1.0,0.65,1.00]
ex_SIP_Xt = [(1E-5)..(1E-4),(1E-5)..(1E-4),(0.0)..(0.0)]
Xtest = SIP_g(ex_SIP_Xt,disc1)
m = BnBModel([(1E-5)..(1E-4),(1E-5)..(1E-4),(-10.0)..(10.0)])
s = BnBSolver()
set_to_default!(s)
s.BnB_atol = 1E-2

# solve optimization problem via interval extension
function test_lower(X,k,pos,opt_inner,UBD)
  FInt = -X[3]
  disc1 = [0.6,0.005,1.0,0.65,1.00]
  GInt = SIP_g(X,disc1)
  feas = true
  for i=1:length(GInt)
    if (GInt[i].lo>0.0)
      feas = false
      break
    end
  end
  val = FInt.lo
  pnt = mid.(X)
  return pnt, val, feas, X, [FInt]
end

function test_upper(X,k,pos,opt_inner,temp)
  feas = true
  val = temp[1].hi
  pnt = mid.(X)
  return pnt, val, feas, X
end

s.Lower_Prob = test_lower
s.Upper_Prob = test_upper
solveBnB!(s,m)

mi = BnBModel([(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)])
si = BnBSolver()
set_to_default!(si)
si.BnB_atol = 5E-2

function inner_lower(X,k,pos,opt_inner,UBD)
  disc1 = [5.5E-5,5.5E-5,4.99]
  FInt = -SIP_g(disc1,X)
  feas = true
  val = FInt.lo
  pnt = mid.(X)
  return pnt, val, feas, X, [FInt]
end

function inner_upper(X,k,pos,opt_inner,temp)
  feas = true
  val = temp[1].hi
  pnt = mid.(X)
  return pnt, val, feas, X
end

si.Lower_Prob = inner_lower
si.Upper_Prob = inner_upper
#solveBnB!(si,mi)

m1 = BnBModel([(1E-5)..(1E-4),(1E-5)..(1E-4),(-10.0)..(10.0)])
s1 = BnBSolver()
set_to_default!(s1)
s1.BnB_atol = 1E-2

# solve optimization problem via interval extension
function test_lower(X,k,pos,opt_inner,UBD)
  FInt = -X[3]
  disc1 = [0.6,0.005,1.0,0.65,1.00]
  GInt = SIP_g(X,disc1) +0.5
  feas = true
  for i=1:length(GInt)
    if (GInt[i].lo>0.0)
      feas = false
      break
    end
  end
  val = FInt.lo
  pnt = mid.(X)
  return pnt, val, feas, X, [FInt]
end

function test_upper(X,k,pos,opt_inner,temp)
  feas = true
  val = temp[1].hi
  pnt = mid.(X)
  return pnt, val, feas, X
end

s1.Lower_Prob = test_lower
s1.Upper_Prob = test_upper
solveBnB!(s1,m1)

#=
SIPopt = SIP_opts()
set_to_default!(SIPopt)
SIPopt.Verbosity = "Full" #"Normal"

# solves sample LBP
vals = 1
count = 1
p_sto = zeros(vals^5,5)
for i=1:vals
  for j=1:vals
    for k=1:vals
      for p=1:vals
        for q=1:vals
          t1 = ex_SIP_P[1].lo + (ex_SIP_P[1].hi - ex_SIP_P[1].lo)*(i-1)/(vals-1)
          t2 = ex_SIP_P[2].lo + (ex_SIP_P[2].hi - ex_SIP_P[2].lo)*(j-1)/(vals-1)
          t3 = ex_SIP_P[3].lo + (ex_SIP_P[3].hi - ex_SIP_P[3].lo)*(k-1)/(vals-1)
          t4 = ex_SIP_P[4].lo + (ex_SIP_P[4].hi - ex_SIP_P[4].lo)*(p-1)/(vals-1)
          t5 = ex_SIP_P[5].lo + (ex_SIP_P[5].hi - ex_SIP_P[5].lo)*(q-1)/(vals-1)
          p_sto[count,:] = [t1,t2,t3,t4,t5]
          count += 1
        end
      end
    end
  end
end

bnd_z2 = ex_SIP_X[2]-ex_SIP_P[4]*0.00005*sqrt((ex_SIP_P[1]*ex_SIP_X[1]+ex_SIP_X[2])/(ex_SIP_P[5]*(ex_SIP_P[4]*0.00005)^2))
bnd_z5 = ex_SIP_X[2]-ex_SIP_P[4]*0.00005*sqrt((ex_SIP_P[1]*ex_SIP_X[1]+ex_SIP_X[2])/(ex_SIP_P[5]*(ex_SIP_P[4]*0.00005)^2))
bnd_z8 = (ex_SIP_X[2]-ex_SIP_P[4]*(0.00005+3.14159*ex_SIP_P[2]^2)*sqrt((((ex_SIP_X[1]+ex_SIP_X[2])/(ex_SIP_P[4]*(0.00005+3.14159*ex_SIP_P[2]^2)))^2)))/(ex_SIP_P[5])

SIPopt.tol = 1E-1
SIPopt.eps_g0 = 1.0
SIPopt.P_LBD = [p_sto[i,:] for i=1:vals^5]
SIPopt.P_UBD = [p_sto[i,:] for i=1:vals^5]
SIPopt.LBP_Opt.DAG_depth = -1
SIPopt.LBP_Opt.f = SIP_f
SIPopt.LBP_Opt.g = SIP_g
SIPopt.LBP_Opt.X = ex_SIP_X
SIPopt.LBP_Opt.g = x -> BndProb_reform(x,[],SIP_g,[],0.0)
#SIPopt.LBP_Opt.LBD_func_relax = "Diff2-MV-OFF"
#SIPopt.LBP_Opt.LBD_problem_relax = "NLP2"
#SIPopt.LBP_Opt.LBD_problem_solver = "Ipopt"
SIPopt.LBP_Opt.LBD_func_relax = "Interval"
SIPopt.LBP_Opt.LBD_problem_relax = "Interval"
SIPopt.LBP_Opt.LBD_problem_solver = "Interval"
SIPopt.LBP_Opt.UBD_func_relax = "Interval"
SIPopt.LBP_Opt.UBD_problem_relax = "Interval"
SIPopt.LBP_Opt.UBD_problem_solver = "Interval"
SIPopt.LBP_Opt.BnB_options.Verbosity = "Normal" #"Full"
SIPopt.LBP_Opt.BnB_options.BnB_tol = 1E-2
#LBDg, xbar, feas, tLBP1,tLBP2 = EAGO_Global_Explicit(SIPopt.LBP_Opt)
#time_lower = @time EAGO_Global_Explicit(SIPopt.LBP_Opt)
#time_lower = @time EAGO_Global_Explicit(SIPopt.LBP_Opt)

xbar = mid.(ex_SIP_X)
SIPopt.LLP_Opt.DAG_depth = -1
SIPopt.LLP_Opt.f = p -> -SIP_g(xbar,p)
SIPopt.LLP_Opt.X = ex_SIP_P
#SIPopt.LLP_Opt.LBD_func_relax = "Diff2-MV-OFF"
#SIPopt.LLP_Opt.LBD_problem_relax = "NLP2"
#SIPopt.LLP_Opt.LBD_problem_solver = "Ipopt"
SIPopt.LLP_Opt.LBD_func_relax = "Interval"
SIPopt.LLP_Opt.LBD_problem_relax = "Interval"
SIPopt.LLP_Opt.LBD_problem_solver = "Interval"
SIPopt.LLP_Opt.UBD_func_relax = "Interval"
SIPopt.LLP_Opt.UBD_problem_relax = "Interval"
SIPopt.LLP_Opt.UBD_problem_solver = "Interval"
SIPopt.LLP_Opt.BnB_options.Verbosity = "Normal"
SIPopt.LLP_Opt.BnB_options.BnB_tol = 1E-2
#INNg1, pbar,feas,tLLP1,tLLP2 = EAGO_Global_Explicit(SIPopt.LLP_Opt)
#time_inner = @time EAGO_Global_Explicit(SIPopt.LLP_Opt)
#time_inner = @time EAGO_Global_Explicit(SIPopt.LLP_Opt)

SIPopt.UBP_Opt.DAG_depth = -1
SIPopt.UBP_Opt.f = SIP_f
SIPopt.UBP_Opt.g = SIP_g
SIPopt.UBP_Opt.X = ex_SIP_X
SIPopt.UBP_Opt.g = x -> BndProb_reform(x,[],SIP_g,[],1.0)
#SIPopt.UBP_Opt.LBD_func_relax = "Diff2-MV-OFF"
#SIPopt.UBP_Opt.LBD_problem_relax = "NLP2"
#SIPopt.UBP_Opt.LBD_problem_solver = "Ipopt"
SIPopt.UBP_Opt.LBD_func_relax = "Interval"
SIPopt.UBP_Opt.LBD_problem_relax = "Interval"
SIPopt.UBP_Opt.LBD_problem_solver = "Interval"
SIPopt.UBP_Opt.UBD_func_relax = "Interval"
SIPopt.UBP_Opt.UBD_problem_relax = "Interval"
SIPopt.UBP_Opt.UBD_problem_solver = "Interval"
SIPopt.UBP_Opt.BnB_options.Verbosity = "Normal"
SIPopt.UBP_Opt.BnB_options.BnB_tol = 1E-2
#UBD_temp, xbar, feas,tUBP1,tUBP2 = EAGO_Global_Explicit(SIPopt.UBP_Opt)
#time_upper = @time EAGO_Global_Explicit(SIPopt.UBP_Opt)
#time_upper = @time EAGO_Global_Explicit(SIPopt.UBP_Opt)
#pbar = [0.66, 0.00275, 1.13595, 0.95, 1.15]
#ubar = [5.5e-05, 0.000100001] × [5.5e-05, 0.000100001] × [9.9781, 9.97815]

xbnd_range1 = xbnds(ex_SIP_X[1:2],ex_SIP_P)
SIP_range = SIP_g(ex_SIP_X,ex_SIP_P)
ex_out = Explicit_SIP_Solve(SIP_f,[],SIP_g,ex_SIP_X,ex_SIP_P,SIPopt)
=#
