# Newton Gauss Siedel with Parametric Tests
function MC_NewtonGS(X0,P,hj,h,opt,Eflag,Iflag,eDflag)
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

  X = copy(X0)
  Xi = copy(X0)
  N = copy(X)
  #println("X:  ", X)
  x_mid = mid.(X)
  k = 1
  #println("x_mid:  ", x_mid)
  #println("P:  ", P)
  #println("typeof(x_mid):  ", typeof(x_mid))
  #println("typeof(P):  ", typeof(P))
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
  #println("J: ", J)
  #println("Y: ", Y)
  #println("B: ", B)
  #println("M: ", M)
  #println("M:  ", M)
  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    #println("S1: ", S1)
    #println("S2: ", S2)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
      #println("S1: ", S1)
      #println("S2: ", S2)
    end
    #println("M[i,i]: ", M[i,i])
    if M[i,i].lo*M[i,i].hi > 0.0
      #println("M route 1")
      N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
      #println("Iteration 1, N[$i]: ", N[i])
    else
      #println("M route 2")
      Ntemp = N
      #println("B[i]: ", B[i])
      eD,N[i],Ntemp[i] = extProcess(N[i],X[i],M[i,i],S1,S2,B[i],rtol)
      #println("N[i]: ", N[i])
      #println("Ntemp[i]: ", Ntemp[i])
      #println("eD: ",eD)
      #println("N[i]: ",N[i])
      #println("Ntemp[i]: ",Ntemp[i])
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

  #println("exclusion:  ",exclusion)
  #println("Iflag:  ",Iflag)
  #println("Eflag:  ",Eflag)
  #println("eDflag:  ",eDflag)


  #
  while ((k<kmax) && isEqual(X,Xi,etol) == false)

    #println("Iteration $k Flag Start:")
    #println("exclusion:  ",exclusion)
    #println("Iflag:  ",Iflag)
    #println("Eflag:  ",Eflag)
    #println("eDflag:  ",eDflag)

    Xi = copy(X0)
    #println("X:  ", X)
    x_mid = mid.(X)
    #println("x_mid:  ", x_mid)
    #println("P:  ", P)
    #println("typeof(x_mid):  ", typeof(x_mid))
    #println("typeof(P):  ", typeof(P))
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
      #println("M[i,i]")
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

    println("Iteration $k Flag End:")
    println("exclusion:  ",exclusion)
    println("Iflag:  ",Iflag)
    println("Eflag:  ",Eflag)
    println("eDflag:  ",eDflag)

  end

  if exclusion
    Eflag = true
  end
  Xtemp = copy(X)
  return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
end

# Newton Gauss Siedel with Parametric Tests
function MC_KrawczykCW(X0,P,hj,h,opt,Eflag,Iflag)
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
  #println("Y: ", Y)
  #println("H: ", H)
  B = Y*H
  M = eye(nx)-Y*J
  #println("J: ", J)
  #println("Y: ", Y)
  #println("B: ", B)
  #println("M: ", M)

  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    #println("S1: ", S1)
    #println("S2: ", S2)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>=i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
      #println("S1: ", S1)
      #println("S2: ", S2)
    end
    N[i] = x_mid[i] - B[i] + S1 + S2
    #println("---Internal Krawczyk CW Iteration #1---")
    #println("N[i]: ", N[i])
    #println("X[i]: ", X[i])
    #println("Strict_XinY(N[i],X[i]): ", Strict_XinY(N[i],X[i]))
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
  #println("Inclusion: ",inclusion)

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

  #
  while ((k<kmax) && (isEqual(X,Xi,etol) == false))
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
      #println("---Internal Krawczyk CW Iteration #2 to n ---")
      #println("N[i]: ", N[i])
      #println("X[i]: ", X[i])
      #println("Strict_XinY(N[i],X[i]): ", Strict_XinY(N[i],X[i]))
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

    #println("Iflag: ", Iflag)
    #println("inclusion", inclusion)
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

  if exclusion
    #println("Eflag: ", Eflag)
    Eflag = true
  end
  return X,Eflag,Iflag,inclusionLow,inclusionHigh
end
