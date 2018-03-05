"""
--------------------------------------------------------------------------------
Function: PI_NewtonGS
--------------------------------------------------------------------------------
Description:
Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
Gauss-Siedel Newton.
--------------------------------------------------------------------------------
Inputs:
X0:        IntervalBox{N,Float64} - Bounds for dependent variables
P:         IntervalBox{N,Float64} - Bounds for independent variables
hj:        function - Jacobian of h(x,p) with respect to x
h:         function - Equations defining potential implicit function
opt:       Array{Float64,1} - Containing the following options
                              opt[1] - Int64: Number of iterations
                              opt[2] - Float64: Tolerance for equality of
                                                intervals
                              opt[3] - Float64: Interval to add to M[i,i] when
                                                processing extended interval
                                                division.
Eflag:     Bool - True if no implicit function exists in box
Iflag:     Bool - True if unique implicit function exists in box
eDflag:    Bool - True if extended division occured in iteration
--------------------------------------------------------------------------------
Returns:
A tuple (X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh) where
X:                IntervalBox{N,Float64} - Output interval box
Xtemp:            IntervalBox{N,Float64} - Extra box formed from extended division
Eflag:            Bool - True if no implicit function exists in box
Iflag:            Bool - True if unique implicit function exists in box
inclusionLow:     Array{Bool} - Each component is true if corresponding variable
                                contracted from lower bound
inclusionHigh:    Array{Bool} - Each component is true if corresponding variable
                                contracted from upper bound
--------------------------------------------------------------------------------
"""
function PI_NewtonGS(X0,P,hj,h,opt,Eflag,Iflag,eDflag)
  # unpacks option file
  kmax = opt[1]
  etol = opt[2]
  rtol = opt[3]

  # Initializes variables
  nx = length(X0)
  S1::Interval = Interval(0.0)
  S2::Interval = Interval(0.0)
  inclusion = [false for i=1:nx]
  inclusionHigh = [false for i=1:nx]
  inclusionLow = [false for i=1:nx]
  eD = 0
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false

  # Initializes storage parameters
  X = copy(X0)
  Xi = copy(X0)
  N = copy(X)
  x_mid = mid.(X)
  k = 1
  H = h(x_mid,P)
  J = hj(X,P)
  Y = eye(nx,nx)
  try
    Y = Preconditioner(hj,X,P,jac="User")
  catch
    Y = eye(nx,nx)
  end
  B = Y*H
  M = Y*J

  # Runs one iteration of parameteric interval NewtonGS
  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    if M[i,i].lo*M[i,i].hi > 0.0
      N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
    else
      Ntemp = N
      eD,N[i],Ntemp[i] = extProcess(N[i],X[i],M[i,i],S1,S2,B[i],rtol)
      if eD == 1
        eDflag = true
        Xtemp = X
        Xtemp[i] = Ntemp[i] ∩ X[i]
        X[i] = N[i] ∩ X[i]
        return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
      end
    end
    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  # checks if all components are included
  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end

  # Runs more iterations until limit hit or intervals are approximately equal
  while ((k<kmax) && isEqual(X,Xi,etol) == false)

    Xi = X0
    x_mid = mid.(X)
    H = h(x_mid,P)
    J = hj(X,P)
    Y = eye(nx,nx)
    try
      Y = Preconditioner(hj,X,P,jac="User")
    catch
      Y = eye(nx,nx)
    end
    B = Y*H
    M = Y*J

    for i=1:nx
      S1 = Interval(0.0)
      S2 = Interval(0.0)
      for j=1:nx
        if (j<i)
          S1 += M[i,j]*(X[j]-x_mid[j])
        elseif (j>i)
          S2 += M[i,j]*(X[j]-x_mid[j])
        end
      end
      if M[i,i].lo*M[i,i].hi > 0.0
        N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
      else
        Ntemp = N
        eD,N[i],Ntemp[i] = extProcess(N[i],X[i],M[i,i],S1,S2,B[i],rtol)
        if eD == 1
          eDflag = true
          Xtemp = X
          Xtemp[i] = Ntemp[i] ∩ X[i]
          X[i] = N[i] ∩ X[i]
          return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
        end
      end
      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1

    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  # sets exclusion flag
  if exclusion
    Eflag = true
  end

  # returns outputs
  Xtemp = copy(X)
  return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
end

"""
--------------------------------------------------------------------------------
Function: PI_KrawczykCW
--------------------------------------------------------------------------------
Description:
Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
componentwise Krawczyk.
--------------------------------------------------------------------------------
Inputs:
X0:        IntervalBox{N,Float64} - Bounds for dependent variables
P:         IntervalBox{N,Float64} - Bounds for independent variables
hj:        function - Jacobian of h(x,p) with respect to x
h:         function - Equations defining potential implicit function
opt:       Array{Float64,1} - Containing the following options
                              opt[1] - Int64: Number of iterations
                              opt[2] - Float64: Tolerance for equality of
                                                intervals
Eflag:     Bool - True if no implicit function exists in box
Iflag:     Bool - True if unique implicit function exists in box
--------------------------------------------------------------------------------
Returns:
A tuple (X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh) where
X:                IntervalBox{N,Float64} - Output interval box
Eflag:            Bool - True if no implicit function exists in box
Iflag:            Bool - True if unique implicit function exists in box
inclusionLow:     Array{Bool} - Each component is true if corresponding variable
                                contracted from lower bound
inclusionHigh:    Array{Bool} - Each component is true if corresponding variable
                                contracted from upper bound
--------------------------------------------------------------------------------
"""
function PI_KrawczykCW(X0,P,hj,h,opt,Eflag,Iflag)
  # unpacks option file
  kmax = opt[1]
  etol = opt[2]

  # Initializes variables
  nx = length(X0)
  S1::Interval = Interval(0.0)
  S2::Interval = Interval(0.0)
  inclusion = [false for i=1:nx]
  inclusionHigh = [false for i=1:nx]
  inclusionLow = [false for i=1:nx]
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false
  X = copy(X0)
  Xi = copy(X)
  N = copy(X)
  x_mid = mid.(X)
  k = 1

  H = h(x_mid,P)
  J = hj(X,P)
  Y = eye(nx,nx)
  try
    Y = Preconditioner(hj,X,P,jac="User")
  catch
    Y = Float64.(eye(nx,nx))
  end

  B = Y*H
  M = eye(nx)-Y*J

  # One iteration of parametric interval Krawczyk
  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>=i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    N[i] = x_mid[i] - B[i] + S1 + S2
    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  # Inclusion tests
  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end

  # Run more iterations of parametric interval Krawczyk
  while ((k<kmax) && isEqual(X,Xi,etol) == false)
    Xi = X0
    x_mid = mid.(X)
    H = h(x_mid,P)
    J = hj(X,P)
    try
      Y = Preconditioner(hj,X,P,jac="User")
    catch
      Y = eye(nx,nx)
    end
    B = Y*H
    M = eye(nx)-Y*J

    for i=1:nx
      S1 = Interval(0.0)
      S2 = Interval(0.0)
      for j=1:nx
        if (j<i)
          S1 += M[i,j]*(X[j]-x_mid[j])
        elseif (j>=i)
          S2 += M[i,j]*(X[j]-x_mid[j])
        end
      end
      N[i] = x_mid[i] - B[i]+S1+S2
      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1
    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  # sets exclusion flag
  if exclusion
    Eflag = true
  end

  # Returns outputs
  return X,Eflag,Iflag,inclusionLow,inclusionHigh
end
