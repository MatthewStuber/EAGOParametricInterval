# To Do List:
# Add Forward/Reverse Diff Iterative Method
# Finish Documentation

################# Basic Krawczyk Iterations #####################
"""Param_Intv_KrawczykCW(z,M,x,p,h,Y) defines a single iteration of
a standard calculation of the Interval Krawczyk method
"""
function Param_Intv_Krawczyk(z,X,P,h,hj,Y)
  I = eye(length(X),length(X))
  K = z-Y*h(z,P)+(I-Y*hj(X,P))*(X-z)
  K_out = [Interval(0.0) for i=1:length(X)]
  for i = 1:length(X)
    K_out[i] = K[i] ∩ X[i]
  end
  return K_out
end

"""Param_Intv_KrawczykCW(z,M,x,p,h,Y;jac="ForwardDiff") defines a single
iteration of a standard calculation of the Interval Krawczyk method
"""
function Param_Intv_Krawczyk(z,X,P,h,Y;jac="ForwardDiff")
  tempf_jxfd(x,p) = JacobianXFD(h,x,p)
  Param_Intv_Krawczyk(z,X,P,h,tempf_jxfd,Y)
end

"""Param_Intv_KrawczykCW(z,M,x,p,h,Y;jac="ReverseDiff") defines a single
iteration of a standard calculation of the Interval Krawczyk method using
ReverseDiff.jl to compute the Jacobian
"""
function Param_Intv_Krawczyk(z,X,P,h,Y;jac="ReverseDiff")
  tempf_jxrd(x,p) = JacobianXRD(h,x,p)
  Param_Intv_Krawczyk(z,X,P,h,tempf_jxrd,Y)
end

################# Componentwise Krawczyk Iterations #####################
"""Param_Intv_KrawczykCW(z,M,x,p,h,Y) defines a single iteration of
a component-wise calculation of the Interval Krawczyk method
"""
function Param_Intv_KrawczykCW(z,X,P,h,hj,Y)
  K_int = [Interval(0.0) for i=1:length(X)]
  X_out = copy(K_int)
  Bk = Y*h(z,P)
  Ak = eye(length(X),length(X))-Y*hj(X,P)
  for i=1:length(X)
    for j=1:length(X)
      if (j<i)
        K_int[i] += Ak[i,j]*(X_out[j]-z[j])
      else (j>=i)
        K_int[i] += Ak[i,j]*(X[j]-z[j])
      end
    end
    K_int[i] = z[i] - Bk[i] + K_int[i]
    X_out[i] = K_int[i] ∩ X[i]
  end
  return X_out
end

"""Param_Intv_KrawczykCW(z,X,P,h,Y;jac="ForwardDiff") defines a single iteration
of a component-wise calculation of the Interval Krawczyk method using
ForwardDiff.jl to compute the Jacobian
"""
function Param_Intv_KrawczykCW(z,X,P,h,Y;jac="ForwardDiff")
  hjxfd(x,p) = JacobianXFD(h,x,p)
  Param_Intv_KrawczykCW(z,X,P,h,hjxfd,Y)
end
"""Param_Intv_KrawczykCW(z,X,P,h,Y;jac="ForwardDiff") defines a single iteration
of a component-wise calculation of the Interval Krawczyk method using
ReverseDiff.jl to compute the Jacobian
"""
function Param_Intv_KrawczykCW(z,X,P,h,Y;jac="ReverseDiff")
  hjxrd(x,p) = JacobianXRD(h,x,p)
  Param_Intv_KrawczykCW(z,X,P,h,hjxrd,Y)
end

################# Basic Newton Iteration #####################
"""Param_Intv_Newton(z,M,x,p,h,Y) defines a single iteration of
of the standard Interval Newton method
"""
function Param_Intv_Newton(z,X,P,h,hj,Y)
  N_int = [Interval(0.0) for i=1:length(X)]
  X_out = copy(N_int)
  Bk = Y*h(z,P)
  Ak = Y*hj(X,P)
  for i=1:length(X)
    for j=1:length(X)
      if (j<i)||(j>i)
        N_int[i] += Ak[i,j]*(X[j]-z[j])
      end
    end
    if in(0,Ak[i,i])
      Xout1 = copy(Xout)
      Xout2 = copy(Xout)
      Xout1[i] = X[i] ∩ (z[i] - (Bk[i] + N_int[i])/Interval(Ak[i,i].lo,0))
      Xout2[i] = X[i] ∩ (z[i] - (Bk[i] + N_int[i])/Interval(0,Ak[i,i].hi))
      return Xout1,Xout2
    else
      N_int[i] = z[i] - (Bk[i] + N_int[i])/Ak[i,i]
    end
    X_out[i] = N_int[i] ∩ X[i]
  end
  return X_out
end

function Param_Intv_Newton(z,X,P,h,Y;jac="ReverseDiff")
  hjxrd(x,p) = JacobianXRD(h,x,p)
  PParam_Intv_Newton(z,X,P,h,hjxrd,Y)
end

function Param_Intv_Newton(z,X,P,h,Y;jac="ForwardDiff")
  hjxfd(x,p) = JacobianXFD(h,x,p)
  Param_Intv_Newton(z,X,P,h,hjxfd,Y)
end

################# Gauss-Siedel Newton Iteration #####################
"""Param_Intv_NewtonGS(z,M,x,p,h,Y) defines a single iteration of
of the Gauss-Seidel Interval Newton method
"""
function Param_Intv_NewtonGS(z,X,P,h,hj,Y)
  N_int = [Interval(0.0) for i=1:length(X)]
  X_out = copy(N_int)
  Bk = Y*h(z,P)
  Ak = Y*hj(X,P)
  for i=1:length(X)
    for j=1:length(X)
      if (j<i)
        N_int[i] += Ak[i,j]*(X_out[j]-z[j])
      elseif (j>i)
        N_int[i] += Ak[i,j]*(X[j]-z[j])
      end
    end
    if in(0,Ak[i,i])
      Xout1 = copy(X_out)
      Xout2 = copy(X_out)
      Xout1[i] = X[i] ∩ (z[i] - (Bk[i] + N_int[i])/Interval(Ak[i,i].lo,0))
      Xout2[i] = X[i] ∩ (z[i] - (Bk[i] + N_int[i])/Interval(0,Ak[i,i].hi))
      return Xout1,Xout2
    else
      N_int[i] = z[i] - (Bk[i] + N_int[i])/Ak[i,i]
    end
    X_out[i] = N_int[i] ∩ X[i]
  end
  return X_out
end

function Param_Intv_NewtonGS(z,X,P,h,Y;jac="ReverseDiff")
  hjxrd(x,p) = JacobianXRD(h,x,p)
  PParam_Intv_NewtonGS(z,X,P,h,hjxrd,Y)
end

function Param_Intv_NewtonGS(z,X,P,h,Y;jac="ForwardDiff")
  hjxfd(x,p) = JacobianXFD(h,x,p)
  Param_Intv_NewtonGS(z,X,P,h,hjxfd,Y)
end

############### Defines Interval Contactor Program ###################
function Param_Intv_Contactor(k,X,P,h,hj,style)
  Xnext = copy(X)
  for i = 1:k
    # calculates Preconditioner
    z = mid.(Xnext)
    Y = Preconditioner(hj,X,P,jac="User")
    # advances iteration
    if (style=="Krawczyk")
      Xnext = Param_Intv_KrawczykCW(z,X,P,h,hj,Y)
    elseif (style=="KrawczykCW")
      Xnext = Param_Intv_KrawczykCW(z,X,P,h,hj,Y)
    elseif (style=="Newton")
      Xnext = Param_Intv_Newton(z,X,P,h,hj,Y)
    elseif (style=="NewtonGS")
      Xnext = Param_Intv_NewtonGS(z,X,P,h,hj,Y)
    else
      error("Unsupported style of contractor")
    end
  end
  return Xnext
end
#=
function Param_Intv_Contactor(k,X,P,h,hj,style;jac="ReverseDiff")
  hjxrd(x,p) = JacobianXRD(h,x,p)
  Param_Intv_Contactor(k,X,P,h,hjxrd,style)
end

function Param_Intv_Contactor(k,X,P,h,hj,style;jac="ForwardDiff")
  hjxfd(x,p) = JacobianXFD(h,x,p)
  Param_Intv_Contactor(k,X,P,h,hjxfd,style)
end
=#
