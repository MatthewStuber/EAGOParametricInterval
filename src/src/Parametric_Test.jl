"""
--------------------------------------------------------------------------------
Function: Miranda
--------------------------------------------------------------------------------
Description:
Runs Miranda's tests for exclusion of implicit function. Implicit function
defined by h(X,P) = 0 in (X,P) implies 0 in {h(x,p} for x in X, p in P}.
--------------------------------------------------------------------------------
Inputs:
* h:         function - Equations defining potential implicit function
* X:         Vector{Interval{Float64}} - Bounds for dependent variables
* P:         Vector{Interval{Float64}} - Bounds for independent variables
--------------------------------------------------------------------------------
Returns:
-1 if no exclusion found. Returns index at which exclusion fails otherwise.
--------------------------------------------------------------------------------
"""
function Miranda(h::Function,
                 X::Vector{Interval{Float64}},
                 P::Vector{Interval{Float64}})
  H = h(X,P)
  for i=1:length(X)
    if (H[i].lo*H[i].hi>0.0)
      return i
    end
  end
  i = -1
  return i
end

"""
--------------------------------------------------------------------------------
Function: MirandaExc
--------------------------------------------------------------------------------
Description:
Runs Miranda's tests for exclusion of implicit function. Implicit function
defined by h(X,P) = 0 in (X,P) implies 0 in {h(x,p} for x in X, p in P}.
--------------------------------------------------------------------------------
Inputs:
* h:         Function - Equations defining potential implicit function
* X:         Vector{Interval{Float64}} - Bounds for dependent variables
* P:         Vector{Interval{Float64}} - Bounds for independent variables
* Eflag_in:  Bool - Input exclusion flag.
--------------------------------------------------------------------------------
Returns:
Returns true if exclusion found. Returns input otherwise.
--------------------------------------------------------------------------------
"""
function MirandaExc(h::Function,
                    X::Vector{Interval{Float64}},
                    P::Vector{Interval{Float64}},
                    Eflag_in::Bool)
   Eflag::Bool = Eflag_in
   i::Int64 = Miranda(h,X,P)
   if (i != -1)
     Eflag = true
   end
   return Eflag
end

"""
--------------------------------------------------------------------------------
Function: partialIncTop
--------------------------------------------------------------------------------
Description:
Sets intervalbox to an the degenerate interval containing only the upper bound
for dimension i and then apply Miranda's test. Advance to next dimension and
repeat
--------------------------------------------------------------------------------
Inputs:
* h:         function - Equations defining potential implicit function
* X:         Vector{Interval{Float64}} - Bounds for dependent variables
* P:         Vector{Interval{Float64}} - Bounds for independent variables
* PIflag:    Bool - Input exclusion flag.
* incHigh:   Array{Bool} - Array that is false if potential partial inclusion
--------------------------------------------------------------------------------
Returns:
Returns tuple (k,PIflag):
* k:        Int64 - Dimension at which partial inclusion potential exists
* PIflag:   Bool - Partial inclusion flag (false if no partial inclusion)
--------------------------------------------------------------------------------
"""
function partialIncTop(h::Function,
                       X::Vector{Interval{Float64}},
                       P::Vector{Interval{Float64}},
                       PIflag::Bool,
                       incHigh::Bool)
  xU::Vector{Float64} = [X[i].hi for i=1:length(X)]
  XU::Vector{Interval{Float64}} = copy(X)
  Xprev::Vector{Interval{Float64}} = copy(X)
  k::Int64 = -1
  ku::Int64 = -1
  for i=1:length(X)
    XU[i] = Interval(xU[i])
    ku = Miranda(h,XU,P)
    if (ku == -1 && ~incHigh[i])
      k = i
      break
    else
      PIflag = false
    end
    XU[i] = Xprev[i]
  end
  return k,PIflag
end

"""
--------------------------------------------------------------------------------
Function: partialIncBot
--------------------------------------------------------------------------------
Description:
Sets intervalbox to an the degenerate interval containing only the lower bound
for dimension i and then apply Miranda's test. Advance to next dimension and
repeat
--------------------------------------------------------------------------------
Inputs:
* h:         function - Equations defining potential implicit function
* X:         Vector{Interval{Float64}} - Bounds for dependent variables
* P:         Vector{Interval{Float64}} - Bounds for independent variables
* PIflag:    Bool - Input exclusion flag.
* incLow:    Vector{Bool} - Array that is false if potential partial inclusion
--------------------------------------------------------------------------------
Returns:
Returns tuple (k,PIflag):
* k:        Dimension at which partial inclusion potential exists
* PIflag:   Partial inclusion flag (false if no partial inclusion)
--------------------------------------------------------------------------------
"""
function partialIncBot(h::Function,
                       X::Vector{Interval{Float64}},
                       P::Vector{Interval{Float64}},
                       PIflag::Bool,
                       incLow::Vector{Bool})
  xL::Vector{Float64} = [X[i].lo for i=1:length(X)]
  XL::Vector{Interval{Float64}} = copy(X)
  Xprev::Vector{Interval{Float64}} = copy(X)
  k::Int64 = -1
  kl::Int64 = -1
  for i=1:length(X)
    XL[i] = Interval(xL[i])
    kl = Miranda(h,XL,P)
    if (kl == -1 && ~incLow[i])
      k = i
      break
    else
      PIflag = false
    end
    XL[i] = Xprev[i]
  end
  return k,PIflag
end

"""
--------------------------------------------------------------------------------
Function: spectralR
--------------------------------------------------------------------------------
Description:
Computes spectral radius for novel inclusion test.
--------------------------------------------------------------------------------
Inputs:
* hj:        function - Jacobian wrt x of potential implicit function
* X:         Vector{Interval{Float64}} - Bounds for dependent variables
* P:         Vector{Interval{Float64}} - Bounds for independent variables
* rho:       Not currently used.
--------------------------------------------------------------------------------
Returns:
Returns tuple (k,eigval):
* k:        Not currently used
* eigval:   Spectral radius
--------------------------------------------------------------------------------
"""
function spectralR(hj::Function,
                    X::Vector{Interval{Float64}},
                    P::Vector{Interval{Float64}},
                    rho)
  k::Bool = true
  J::Array{Interval{Float64},2} = hj(X,P)
  Y::Array{Float64,2} = Preconditioner(hj,X,P,jac="User")
  YJ::Array{Float64,2} = abs.(Y)*diam.(J)/2
  if (length(YJ)>1)
    eigval = eigmax(YJ)
  else
    eigval = YJ[1]
  end
  return k,eigval
end

"""
--------------------------------------------------------------------------------
Function: BoundaryTest
--------------------------------------------------------------------------------
Description:
Runs novel test for inclusion of implicit function
--------------------------------------------------------------------------------
Inputs:
* h:                function - h(x,p) potential implicit function
* hj:               function - Jacobian wrt x of potential implicit function
* X0:               Vector{Interval{Float64}} - Bounds for dependent variables
* P:                Vector{Interval{Float64}} - Bounds for independent variables
* opt:              Array{Any,1} - Parameters used
* PIcert:           Bool: Partial certification flag
* PIflag:           Bool: Partial inclusion flag
* Iflag:            Bool: Inclusion flag
* Eflag:            Bool: Exclusion flag
* inclusionLow:     Vector{Bool}: Array of boolean for inclusion test
* inclusionHigh:    Vector{Bool}:  Array of boolean for inclusion test
--------------------------------------------------------------------------------
Returns:
Returns tuple (PIcert,PIflag,Iflag,Eflag):
* PIcert:           Bool: Partial certification flag
* PIflag:           Bool: Partial inclusion flag
* Iflag:            Bool: Inclusion flag
* Eflag:            Bool: Exclusion flag
--------------------------------------------------------------------------------
"""
function BoundaryTest(h::Function,
                      hj::Function,
                      X0::Vector{Interval{Float64}},
                      X::Vector{Interval{Float64}},
                      P::Vector{Interval{Float64}},
                      opt::Array{Any,1},
                      PIcert::Bool,
                      PIflag::Bool,
                      Iflag::Bool,
                      Eflag::Bool,
                      inclusionLow::Vector{Bool},
                      inclusionHigh::Vector{Bool})
  PIflagT::Bool = PIflag
  PIflagB::Bool = PIflag
  Xtemp::Vector{Interval{Float64}} = copy(X)
  Xin::Vector{Interval{Float64}} = copy(X)
  Pin::Vector{Interval{Float64}} = copy(P)
  ku::Int64,PIflagT = partialIncTop(h,Xin,Pin,PIflagT,inclusionHigh)
  kl::Int64,PIflagB = partialIncBot(h,Xin,Pin,PIflagB,inclusionLow)
  rho::Float64 = 1.01
  eDflag::Bool = false
  if (Strict_XinY(X,X0) || ((ku==-1) && (kl==-1)))
    blank,rho = spectralR(hj,Xin,Pin,rho)
    PIcert = true
    PIflag = false
    if (rho<1.0)
      Pmid::Vector{Float64} = mid.(Pin)
      if (opt[4]=="NewtonGS")
        Xnew::Vector{Interval{Float64}},Xtemp,Eflag,Iflag,eDflag = PI_NewtonGS(Xin,Pmid,hj,h,opt,Eflag,Iflag,eDflag)
      elseif (opt[4]=="KrawczykCW")
        Xnew,Eflag,Iflag,blank1,blank2 =  PI_KrawczykCW(Xin,Pmid,hj,h,opt,Eflag,Iflag)
        end
    else
      Iflag = false
    end
  end
  return PIcert,PIflag,Iflag,Eflag
end

#=
# FUTURE CODE FOR ADDITIONAL TEST
function NewBoundaryTest(Xtemp,X,P,Iflag)
  nx = length(X)
  kmax = 1000
  k = 1
  l = 0
  ztemp = 0.0
  eps = 1.0E-12
  Xf = copy(X)
  Xnew = copy(X)
  Pf = copy(P)
  l = SingleStep(Xnew,Pf)
  if (l == 0)
    k = kmax
  end
  for i = 1:nx
    r[2*i-1] = Xf[i].lo - Xnew[i].lo + eps
    r[2*i] = Xnew[i].hi - Xf[i].hi + eps
    z[2*i-1] = Xf[i].lo
    z[2*i] = Xf[i].hi
  end
  # construct Jacobian
  Jr =

end

function SingleStep()
end
=#
