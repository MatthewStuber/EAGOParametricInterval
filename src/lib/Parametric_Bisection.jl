"""
--------------------------------------------------------------------------------
Function: XP_Bisection
--------------------------------------------------------------------------------
Description: Bisects box in either X or P dimension.
--------------------------------------------------------------------------------
Inputs:
h           function - h(z,p)
hj          function - Jacobian of h(z,p) w.r.t. z
CurrNode    IntervalBox - Current X node
P           IntervalBox - Current P node
NextNode    IntervalBox - New X node
NextNodeP   IntervalBox - New P node 1
Ptemp       IntervalBox - New P node 2
PIflag      Bool - Partial Inclusion Flag
Iflag       Bool - Inclusion flag
Eflag       Bool - Exclusion flag
ptol        Float64 - tolerance for P dimension size
pbisect     Bool - bisect based on relative width?
Pstart      Pstart - initial P size
--------------------------------------------------------------------------------
Returns:
bp          Int64 - Cut type flag
CurrNode    IntervalBox - Current X node
P           IntervalBox - Current P node
NextNode    IntervalBox - New X node
NextNodeP   IntervalBox - New P node 1
Ptemp       IntervalBox - New P node 2
--------------------------------------------------------------------------------
"""

function XP_Bisection(h,hj,CurrNode,P,NextNode,NextNodeP,Ptemp,PIflag,Iflag,Eflag,ptol,pbisect,Pstart)
  pPart = false
  xCut = 0.15
  pCutL = 0.0
  pCutR = 0.0
  IflagL = copy(Iflag)
  IflagR = copy(Iflag)
  EflagL = copy(Eflag)
  EflagR = copy(Eflag)
  XLeft = copy(CurrNode)
  XRight = copy(NextNode)
  PLeft = copy(P)
  PRight = copy(NextNodeP)
  Pmid = copy(Ptemp)
  j = 0
  jp = 0
  bp = 0
  j, xCut,PIflag = bisectDirX(h,CurrNode,P,xCut,PIflag)
  PIflagL = PIflag
  PIflagR = PIflag
  if (PIflag == false)
    XLeft[j] = Interval(CurrNode[j].lo,xCut)
    XRight[j] = Interval(xCut,CurrNode[j].hi)
    CurrNode = copy(XLeft)
    NextNode = copy(XRight)
    return 1,CurrNode,P,NextNode,NextNodeP,Ptemp
  else
    jp = bisectDirP(P,ptol,pbisect,Pstart)
    if (jp != -1)
      if pPart==true
        bp,jp,pCutL,pCutR = branchP(CurrNode,P,jp,pCutL,pCutR)
      else
        bp = 1
      end
      if (bp == 1)
        PLeft[jp] = Interval(P[jp].lo,mid(P[jp]))
        PRight[jp] = Interval(mid(P[jp]),NextNodeP[jp].hi)
        P = copy(PLeft)
        NextNodeP = copy(PRight)
        return 2,CurrNode,P,NextNode,NextNodeP,Ptemp
      elseif (bp == 3)
        PLeft[jp] = Interval(P[jp].lo,pCutL)
        Pmid[jp] = Interval(pCutL,pCutR)
        PRight[jp] = Interval(pCutR,P[jp].hi)
        P = copy(PLeft)
        Ptemp = copy(Pmid)
        NextNodeP = copy(PRight)
        return 3,CurrNode,P,NextNode,NextNodeP,Ptemp
      end
    end
  end
  if (false) # originally just ((jp == -1) && (bp != 1))
    xCut = xCut*CurrNode[j].lo + (1.0-xCut)*CurreNode[j].hi
    XLeft[j] = Interval(CurrNode[j].lo,xCut)
    XRight[j] = Interval(xCut,CurrNode[j].hi)
    k,rhoL = spectralR(hj,XLeft,P,rhoL)
    k,rhoR = spectralR(hj,XRight,P,rhoR)
    if (rhoL)
      XLeft,Iflag,Eflag,PIflag = partIncProc(h,hj,XLeft,P,j,Iflag,Eflag,PIflag)
    end
    if (rhoR)
      XRight,Iflag,Eflag,PIflag = partIncProc(h,hj,XRight,P,j,Iflag,Eflag,PIflag)
    end
    if ((~PIflagL) && (~PIflagR))
      PIflag = false
    end
  end

  CurrNode = copy(XLeft)
  NextNode = copy(XRight)
  return bp,CurrNode,P,NextNode,NextNodeP,Ptemp
end


"""
--------------------------------------------------------------------------------
Function: bisectDirX
--------------------------------------------------------------------------------
Description: Determine direction and position of cut in X dimension.
--------------------------------------------------------------------------------
Inputs:
h        function - h(z,p)
X        IntervalBox - for state variables
P        IntervalBox - for decision variables
x        Currently Unused
PIFlag   Bool - Flag for partial inclusion.
--------------------------------------------------------------------------------
Returns:
jmax            - The dimension to cut X in
xCut            - Cut value for X
piflagtemp      - Partial Inclusion Flag
--------------------------------------------------------------------------------
"""
function bisectDirX(h,X,P,x,PIFlag)

  piflagtemp = copy(PIFlag)
  v = 1
  iterm = 4
  jmax = 1
  k = -1
  maxL = 0.0
  disc = 20
  Xold = copy(X)
  Xnew = copy(X)
  xCut = 0.0
  maxL,jmax = findmax(diam.(X))
  if (v == 1)
    rad = diam(X[jmax])/2.0
    for l=0:(disc/2)
      xCut = X[jmax].lo+rad/(iterm) + l*(diam(X[jmax])-2.0*rad/(iterm))/(disc-1)
      Xnew[jmax] = xCut
      k = Miranda(h,Xnew,P)
      if (k>=0)
        return jmax, xCut, piflagtemp
      end
      xCut = X[jmax].hi-rad/(iterm) - l*(diam(X[jmax])-2.0*rad/(iterm))/(disc-1)
      Xnew[jmax] = xCut
      k = Miranda(h,Xnew,P)
      if (k>=0)
        return jmax, xCut, piflagtemp
      else
        Xnew[jmax] = Xold[jmax]
      end
    end
  end

  if (k == -1)
    v = 2
  end


  if (v == 2)
    iterm = iterm*4
    Xnew = copy(Xold)
    disc = 4*disc
    for i = 1:length(X)
        rad = diam(X[i])/2.0
        for l = 0:(disc/2)
          xCut = X[i].lo+rad/(iterm) + l*(diam(X[i])-2.0*rad/(iterm))/(disc-1)
          Xnew[i] = xCut
          k = Miranda(h,Xnew,P)
          if (k>=0)
            return i, xCut, piflagtemp
          end
          xCut = X[i].hi-rad/(iterm) - l*(diam(X[i])-2.0*rad/(iterm))/(disc-1)
          Xnew[i] = xCut
          k = Miranda(h,Xnew,P)
          if (k>=0)
            return i, xCut, piflagtemp
          else
            Xnew[i] = Xold[i]
          end
        end
    end
    if (k == -1)
      piflagtemp = true
    end
  end
  return jmax, xCut, piflagtemp
end

"""
--------------------------------------------------------------------------------
Function: bisectDirP
--------------------------------------------------------------------------------
Description: Determines the direction to branch P in based on absolute/relative
width and tolerance for P dimension size.
--------------------------------------------------------------------------------
Inputs:
P         IntervalBox - P Box to bisect
ptol      Float64 - tolerance for P dimension size
pbisect   Bool - bisect based on relative width?
Pstart    IntervalBox - initial P size
--------------------------------------------------------------------------------
Returns:
jmax - The dimension to branch P in.
--------------------------------------------------------------------------------
"""
function bisectDirP(P,ptol::Float64,pbisect::Bool,Pstart)
  jmax = 1
  maxL = 0.0
  if (~pbisect)
    for i=1:length(P)
      if (maxL < diam(P[i]))
        maxL = diam(P[i])
        jmax = i
      end
    end
    if (maxL < ptol)
      return -1
    end
  else
    for i=1:length(P)
      if (maxL < diam(P[i])/diam(Pstart[i]))
        maxL = diam(P[i])/diam(Pstart[i])
        jmax = i
      end
    end
    if (maxL < ptol)
      return -1
    end
  end
  return jmax
end

# Unfinish bisection based on spectral radius
function branchP(X,P,jmax,pCutL,pCutR)
  #=
  Plocal = P
  maxD = 0.0
  try
    Pf = [P[j] for j=1:length(P)]
    p = [abs(Pf[j].lo+Pf[j].hi)/2.0 for j=1:length(P)]
    Xf = [X[j] for j=1:length(X)]
    rhok = rho
    fpsecR(Xf,Pf,rho)
  catch
  end
  =#
end

"""
--------------------------------------------------------------------------------
Function:  partIncProc
--------------------------------------------------------------------------------
Description:
Partial inclusion check.
--------------------------------------------------------------------------------
Inputs:
h         function - h(z,p)
hj        function - Jacobian of h(z,p) w.r.t. z
X0        IntervalBox - state variables
P         IntervalBox - decision variables
exDim     Int64 - Dimension flag
Iflag     Bool - Inclusion flag
Eflag     Bool - Exclusion flag
PIflag    Bool - Partial Inclusion Flag
--------------------------------------------------------------------------------
Returns:
X         IntervalBox - State variables
Iflag     Bool - Inclusion flag
Eflag     Bool - Exclusion flag
PIflag    Bool - Partial Inclusion Flag
--------------------------------------------------------------------------------
"""
function partIncProc(h,hj,X0,P,exDim,Iflag,Eflag,PIflag)

  nx = length(X)
  S1 = Interval(0.0)
  S2 = Interval(0.0)
  k = 0
  incl = -1
  side = 0
  inclusion = [false for i=1:length(X)]
  exclusion = false
  Iflag = false
  Eflag = false

  X = copy(X0)
  Xpart = copy(X0)
  Xpart[exDim] = xCut
  if (X[exDim].lo == xCut)
    side = 1
  end

  x_mid = mid.(X)
  H = h(x_mid,P)
  J = hj(X,P)
  Y = Preconditioner(hj,X,P,jac="User")
  B = Y*H
  M = Y*J
  Xi = X

  for i=1:nx
    S1 = S2 = Interval(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
    if N[i] ⊆ X[i]
      inclusion[i] = true
    else
      inclusion[i] = false
    end
    if ~isdisjoint(N[i],X[i])
      if (i == exDim)
        wL = X[i].lo - N[i].lo
        wU = N[i].hi - X[i].hi
        if (side == 0)
          X[i] = Interval(X[i].lo,X[i].hi+0.25*wU)
        else
          X[i] = Interval(X[i].lo-wL*0.25,X[i].hi)
        end
      else
        X[i] = N[i] ∩ X[i]
      end
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

  Xpart = X
  if (side == 0)
    Xpart[exDim] = X[exDim].hi
    if (Xpart[exDim]<=xCut)
      PIflag = false
      return X,Iflag,Eflag,PIflag #FIX ME
    end
  else
    Xpart[exDim] = X[exDim].lo
    if (Xpart[exDim]>=xCut)
      PIflag = false
      return X,Iflag,Eflag,PIflag #FIX ME
    end
  end
  incl = Miranda(h,Xpart,P)

  return X,Iflag,Eflag,PIflag #FIX ME
end
