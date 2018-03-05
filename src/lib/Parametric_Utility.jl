function kdel(i,j)
  i == j ? 1.0 : 0.0
end

function Expr_SubSymbols(expr,symdict)
      stack = []
      exprlist = []
      push!(stack,expr)
      cont_flag = false
      while length(stack)>0
            a = pop!(stack)
            if (isempty(exprlist))
                  push!(exprlist,a)
            end
            if ~(typeof(a)<:Number) && ~(typeof(a)==Int64)
              for i=2:length(a.args)
                  if (typeof(a.args[i])==Symbol)
                        a.args[i] = symdict[a.args[i]]
                  elseif typeof(a.args[i])==Expr
                        push!(stack,a.args[i])
                  end
              end
            end
      end
      return exprlist[1]
end

function GenerateJacobianX(ExprArr,SymXArr,SymPArr)

  # creates symbol dictionary
  symdict = Dict{Any,Any}()
  nc = length(SymXArr)
  ns = length(SymPArr)
  for i=1:nc
      symdict[SymXArr[i]] = :(z[$i])
  end
  for i=1:ns
      symdict[SymPArr[i]] = :(p[$i])
  end

  # creates calculates Jacobian symbolically and subs in symbols
  StringArr = string.(ExprArr)
  dexpr_arr = Expr(:vcat)
  for i = 1:length(StringArr)
    out = Calculus.differentiate(StringArr[i],SymXArr)
    expr_out = copy(out)
    for j = 1:length(out)
      expr_out[j] = Expr_SubSymbols(out[j],symdict)
    end
    temp_expr = Expr(:row)
    for j = 1:length(out)
      push!(temp_expr.args,out[j])
    end
    push!(dexpr_arr.args,temp_expr) # populations Jacobian storage
  end

  # puts expression object into function
  @eval Jacob(z,p) = $dexpr_arr
  return Jacob
end

function GenerateH(ExprArr, SymXArr,SymPArr)
    # creates list of symbols
    symdict = Dict{Any,Any}()
    for i=1:length(SymXArr)
        symdict[SymXArr[i]] = :(x[$i])
    end
    for i=1:length(SymPArr)
        symdict[SymPArr[i]] = :(p[$i])
    end
    # creates array of expressions
    htemp1 = deepcopy(ExprArr)
    htemp2 = []
    for i = 1:length(ExprArr)
        push!(htemp2,Expr_SubSymbols(htemp1[i],symdict))
    end
    #
    hexpr = Expr(:vect)
    for i = 1:length(htemp2)
        temp = htemp2[i]
        push!(hexpr.args,Expr(:call,:eval,temp))
    end
    @eval hnew(x,p) = $hexpr
    return hnew
end

# Unpack function
function Unpack(stack)
  cnode = pop!(stack)
  return cnode[1],cnode[2]
end

function Strict_XinY(X,Y)
  k = true
  for i=1:length(X)
    if ((X[i].lo<=Y[i].lo)||
        (X[i].hi>=Y[i].hi))
      k = false
    end
  end
  return k
end

# Tests for equality to within a tolerance atol
function isEqual(X1,X2,atol)
  out::Bool = true
  for i=1:length(X1)
    if (abs(X1[i].lo-X2[i].lo)<atol ||
        abs(X1[i].hi-X2[i].hi)<atol )
        out = false
        break
    end
  end
  return out
end

# Extended Interval Arithmetic
function extDivide(A,B,C)
  if ((A.lo == -0.0) && (A.hi == 0.0))
    B = Interval(-Inf,Inf)
    C = B
    return 0,B,C
  end
  if (A.lo == 0.0)
    B = Interval(1.0/A.hi,Inf)
    C = Interval(Inf,Inf)
    return 1,B,C
  elseif (A.hi == 0.0)
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(-Inf,-Inf)
    return 2,B,C
  else
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(1.0/A.hi,Inf)
    return 3,B,C
  end
  return -1,B,C
end

# Error & Case Handling for Extended IA in NewtonGS
function extProcess(N,X,Mii,S1,S2,B,rtol)
  v = 1
  Ntemp = copy(N)
  IML = Interval(0.0)
  IMR = Interval(0.0)
  M = (B+S1+S2)+Interval(-rtol,rtol)
  if (M.lo<=0 || M.hi>=0)
    N = Interval(-Inf,Inf)
    return 0, N, Ntemp
  end
  if (v == 1)
    k,IML,IMR = extDivide(Mii,IML,IMR)
    if (k == 1)
      N = mid(X)-M*IML
      return 0, N, Ntemp
    elseif (k == 2)
      N = mid(X)-M*IMR
      return 0, N, Ntemp
    elseif (k == 3)
      NR = mid(X)-M*IMR
      NL = mid(X)-M*IML
      if (~isdisjoint(NL,X) && isdisjoint(NR,X))
        N = NL
        return 0, N, Ntemp
      elseif (~isdisjoint(NR,X) && isdisjoint(NL,X))
        N = NR
        return 0, N, Ntemp
      elseif (~isdisjoint(NL,X) && ~isdisjoint(NR,X))
        N = NL
        Ntemp = NR
        return 1, N, Ntemp
      else
        return -1, N, Ntemp
      end
    end
  else
    N = Interval(-Inf,Inf)
    return 0, N, Ntemp
  end
end


"""
--------------------------------------------------------------------------------
Function: Preconditioner
--------------------------------------------------------------------------------
Description:
Runs Miranda's tests for exclusion of implicit function. *Implicit function
defined by h(X,P) = 0 in (X,P) implies 0 in {h(x,p} for x in X, p in P}.
--------------------------------------------------------------------------------
Inputs:
h:         function - Equations defining potential implicit function
X:         IntervalBox{N,Float64} - Bounds for dependent variables
P:         IntervalBox{N,Float64} - Bounds for independent variables
jac="User"
--------------------------------------------------------------------------------
Returns:
Preconditioning matrix.
--------------------------------------------------------------------------------
"""
function Preconditioner(h,X,P;jac="User")
  J = h(X,P)
  if (length(X)>1)
    Y = inv(mid.(J))
  else
    Y = 1.0/(mid(J[1]))
  end
  return Y
end
